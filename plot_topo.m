
clear all

% Counter for merging FFT highest values
    FFTcount = 1;

% _________INPUT DATA___________ %
        
     % Multiple files  
        start    = 600      ;
        step     = 1       ;
        last     = 600     ;

    % For topographic curves
        % Axes limit (x = X(km) ; y = topo(m))
            xmin = 2750      ;
            xmax = 3500     ;
            ymin = -12000    ;
            ymax = 5000     ;
        % Vertical exageration
            vertexg = 5     ;
            
    % For velocity curves
        % Vertical profiles (in km) 
        % Horizontal buffer for mean calculation of topography (km)
            Nmean   = 10  ;
        % Axes limit (x = Time(Myr) ; y = velocity(mm/yr)) 
            xminb   = 0    ;
            xmaxb   = 50    ;
            ymaxb   = 5     ;
            
    % LOOP FOR VERTICAL PROFILES (km)
       % for profil = 2700:1:3500 
            
% _________INPUT DATA___________ %
    

% Initial counter for gray-scale
Ngray = 0;

for num_name = start:step:last

% Total number of loops (truncated)
    TTloop  = fix(((last-start)/step)+1);
    
%File name
    nname     = 'topo'; 
    ext       = '.txt';

%read input data
    fname    =  [nname,num2str(num_name) ext];
    fid      =  fopen(fname,'r');
    A        =  textscan(fid, '%n %n %n');
    fclose(fid);
    
    timecolumn  =   A{1};
    xlength     =   A{2};
    topo        =   A{3};
% Time of file
    time = timecolumn(2,1);
    timestr = num2str(time/1e6);
% X-axe in km
    Nxlength = xlength/1000;
% Update topography to sea level
    sealevel = 10000;
    Ntopo = sealevel - topo;
% Graph setting
    x = Nxlength;
    y = Ntopo;
% Vertical exageration calculation
    vert = vertexg*1e-3;

% _____DRAWING TOPOGRAPHY 1_____ %
figure(1);
%plot(x,y,'color',[0.9,0.9,0.9]-((0.9/TTloop)*Ngray));
plot(x,y, 'k');
xlabel('Distance (km)'),ylabel('topography (m)'); title([nname ' - ' timestr ' ' 'Myr']);
axis([2700 3500 ymin ymax]); daspect([vert 1 1]); grid on; %hold on;

% Update counter for gray-scale
    Ngray = Ngray+1;
    
% VERTICAL VELOCITY CALCULATION AT THE SURFACE
% Extract mean topo value (in meters) on selected profil 1
    [c,index] = min(abs(Nxlength-profil));
    % for high-resolution area, 1 step <=> 500 m
        topo_profil = mean(Ntopo(index-((Nmean*1000)/500):index+((Nmean*1000)/500)));
% Compile topo (in mm) and time (in yr) from successive steps    
    num_nameb = num_name-(start-1);
    compil_topo(num_nameb,:) = (topo_profil.*1);%*1e3);
    compil_time(num_nameb,:) = time;
% Calculate DELTA(time) and DELTA(topo)
    Dtime = diff(compil_time);
    Dtopo = diff(compil_topo);
% Calculate vertical velocity (in mm/yr)
    VertVel = Dtopo./Dtime;
% For plotting velocity during time, remove the first value of time
    comp_time_reduc = compil_time(2:end);
   
% CREATION OF TOPOGRAPHY MATRIX FOR 3D-PLOT
% Create the initial matrix
    if Ngray < 2
        matrix_topo = zeros(length(x),TTloop);
    end
% Add topography profil for each loop (i.e. timestep)
    matrix_topo(:,Ngray) = (y);

end

%% _____DRAWING TOPO AS FUNCTION OF TIME (3D-plot) _____ %
%figure(1);
%surf((compil_time/1e6),x,matrix_topo,'EdgeColor','none');demcmap('inc',[-10000 2000],1);colorbar;
%xlabel('Time (Myr)'),ylabel('X axe (km)'),zlabel('Topography (m)');
% Scale of the plot (x = time in Myr; y = X-axis in km)
%    xlim([xminb xmaxb]);ylim([xmin xmax]);
% Orientation of 3D plot
%    view(70,40);
%   view(0,90); % View from the top


%% _____DRAWING TOPO & VELOCITY_____ %
% figure(2);
% subplot(2,1,1);
% plot(compil_time/1e6,compil_topo,'LineWidth',1.5);
% xlabel('Time (in Myr)'),ylabel('Topography (m)');
% axis([xminb xmaxb ymin ymax]);
% grid on; hold on;
% %xlim([xminb xmaxb]);

% subplot(2,1,2);
% plot((comp_time_reduc/1e6),VertVel,'LineWidth',1.5);
% xlabel('Time (in Myr)'),ylabel('Vertical velocity (mm/yr)');
% title(['Uplift/Subsidence from ' num2str(start) ' to ' num2str(last) ' steps']);
% axis([xminb xmaxb -ymaxb ymaxb]);
% grid on; hold on;


%% FAST FOURIER TRANSFORMATION
% Regular resampling of time (yr) and velocity (mm/yr)
    Reg_time = linspace(min(comp_time_reduc), max(comp_time_reduc), length(comp_time_reduc));
    Reg_veloc = interp1(comp_time_reduc,VertVel,Reg_time);
% Sampling period (depending on time step) (yr)
    Ts = mean(diff(Reg_time));
% Sampling frequency (yr-1)
    Fs = 1/Ts;                  
% Signal length
    Ls = length(Reg_time);
% Fast Fourier transformation
    FFT = fft(Reg_veloc,Ls);
% Compute 2-then-1-sided amplitude spectrum (y-axis) of FFT
    P2 = abs(FFT/Ls);
    P1 = P2(2:Ls/2+1);              % From "2" to avoid "freq = 0"
    P1(2:end-1) = 2*P1(2:end-1);
% Define frequency domain (x-axis) (in yr-1)
    Fsaxis = Fs*(1:(Ls/2))/Ls;      % From "1" to avoid "freq = 0"

% COMPILATION OF MAX-FFT CALCULATIONS
% Extract index of max value of FFT (amplitude)
    [valueFFT,indexFFT] = max(P1);
% Compile max FFT values (frequence, yr-1) in a new matrix 
    MaxFsaxis(:,FFTcount) = Fsaxis(indexFFT);
% Compile profils for calculation of FFT in a new matrix
    Compil_profil(:,FFTcount) = profil;

% Update FFT counter
    FFTcount = FFTcount+1;  

% % _____DRAWING SINGLE FFT_____ %
% figure(3);
% plot(Fsaxis,P1,'LineWidth',1.5);
% xlabel('Frequency (yr-1)'),ylabel('Amplitude');
% title(['Fast Fourier Transform from ' num2str(start) ' to ' num2str(last) ' steps']);
% grid on;
% hold on;
% xlim([0 5e-6]);

%end

% %% _____DRAWING FFT COMPILATION_____ %
% % Relevant for several time profiles
% figure(4);
% subplot(1,2,1);
% scatter(Compil_profil,((1./MaxFsaxis)/1e6),80,(1./MaxFsaxis),'filled'); grid on;
% xlabel('Location of vertical profil (km)'); ylabel('Frequency of topographic variations (Myr)'); title('Max. FFT value for different vertical profils');
% 
% subplot(1,2,2);
% histogram(((1./MaxFsaxis)/1e6),28); grid on;
% xlabel('Frequency of topographic variations (Myr)'); ylabel('Number of vertical profils'); title('Distribution of max. FFT');
