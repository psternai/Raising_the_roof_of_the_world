
% Visualising topo for active margin model
clear all; 
clf;

% Counter for pressure arrays
count = 1;

% _____INPUT DATA_____ %

    % Data filename
    % prn file
for num_name = 0:10:0

    % Box size
    xmin = 875;
    xmax = 1040;
    ymin = 5;
    ymax = 65;

    % Initial fluid velocity
    vfluid = 3e-10;
    
    % Depth of 2nd streamlines from the top of interface (km)
    str2 = 2;
    
    % Steps for streamlines
    stpstr = 10;
    
% _____INPUT DATA_____ %   

%File name
nname='clim'; 
ext='.prn';
% Filename
file = [nname,num2str(num_name)]
%read input data
fileprn=[nname,num2str(num_name) ext];
% Opening data file
fdata=fopen(fileprn,'rb'); 
% Read sizes of variables
A=fread(fdata,4,'uchar');
% Read model parameters
% Grid resolution
xnumx=fread(fdata,1,'int64');
ynumy=fread(fdata,1,'int64');
% Markers per cell
mnumx=fread(fdata,1,'int64');
mnumy=fread(fdata,1,'int64');
% Number of markers
marknum=fread(fdata,1,'int64');
% Model sizes
xsize=fread(fdata,1,'float64');
ysize=fread(fdata,1,'float64');
% Pressure value
pinit=fread(fdata,5,'float64');
% Gravity
GXKOEF=fread(fdata,1,'float64');
GYKOEF=fread(fdata,1,'float64');
% Number of rocks
rocknum=fread(fdata,1,'int32');
% Number of Boundary conditions
bondnum=fread(fdata,1,'int64');
% Stage,time
n1=fread(fdata,1,'int32');
timesum=fread(fdata,1,'float64')
% Skip rock proprties
curpos0=4+2*4+16*8+rocknum*(8*27+4);
fseek(fdata,curpos0,'bof');

% ======
%   VARIABLES
%   pr = pressure
%   vx = horizontal velocity component
%   vy = vertical velocity component
%   dro = thermal expansion??
%   exx = horizontal strain rate component
%   exy = shear strain rate
%   sxx = horizontal stress component
%   sxy = shear stress component
%   drp = ??
%   ro = density
%   nu = viscosity
%   nd = ??
%   sxxe = horizontal elastic stress component
%   sppe = ??
%   sxye = elastic shear stress component
%   exxe = horizontal elastic strain rate component
%   exye = elastic shear strain rate component
%   esp = ??
%   gg = 
%   gd = 
%   ep = 
%   et =
%   tk = temeprature
%   cp = heat capacity
%   kt = thermal conductivity
%   Ht = radioactive heat production
% ======

% Read nodes information
for i=1:1:xnumx
    for j=1:1:ynumy
        vbuf=fread(fdata,3,'float32');
        pr(j,i)=vbuf(1); 
        vx(j,i)=vbuf(2);
        vy(j,i)=vbuf(3);
        vbuf1=fread(fdata,3,'int64');
        vbuf2=fread(fdata,21,'float32');
        exx(j,i)=vbuf2(1);
        dro(j,i)=vbuf2(2);
        exy(j,i)=vbuf2(3);
        sxx(j,i)=vbuf2(4);
        drp(j,i)=vbuf2(5);
        sxy(j,i)=vbuf2(6);
        ro(j,i)=vbuf2(7);
        nu(j,i)=vbuf2(8);
        nd(j,i)=vbuf2(9);
        sxxe(j,i)=vbuf2(10);
        sppe(j,i)=vbuf2(11);
        sxye(j,i)=vbuf2(12);
        exxe(j,i)=vbuf2(13);
        exye(j,i)=vbuf2(14);
        esp(j,i)=vbuf2(15);
        mu(j,i)=vbuf2(16);
        VertProf(j,i)=vbuf2(17);
        gd(j,i)=vbuf2(18);
        ep(j,i)=vbuf2(19);
        et(j,i)=vbuf2(20);
        tk(j,i)=vbuf2(21);
        vbuf3=fread(fdata,1,'int64');
        vbuf4=fread(fdata,3,'float32');
        cp(j,i)=vbuf4(1);
        kt(j,i)=vbuf4(2);
        ht(j,i)=vbuf4(3);
    end
end

% Skip all nodes
curpos2=curpos0+(4*27+8*4)*xnumx*ynumy;
fseek(fdata,curpos2,'bof');

% Read Gridline positions
gx=fread(fdata,xnumx,'float32');
gy=fread(fdata,ynumy,'float32');
eii=ones(ynumy,xnumx)*1e-16;
sii=ones(ynumy,xnumx)*1e+4;
for i=1:1:xnumx-2
    for j=1:1:ynumy-2
        eii(j+1,i+1)=(exy(j+1,i+1)^2+((exx(j+1,i+1)+exx(j+2,i+1)+exx(j+1,i+2)+exx(j+2,i+2))/4)^2)^0.5;
        sii(j+1,i+1)=(sxy(j+1,i+1)^2+((sxx(j+1,i+1)+sxx(j+2,i+1)+sxx(j+1,i+2)+sxx(j+2,i+2))/4)^2)^0.5;
    end
end

% CALCULATION OF PRESSURE AVERAGE
    % Initial pressure
    if (count<2)
        pr_sum = zeros(size(pr));
    else
        pr_sum = pr_tt;
    end
    % Implementation of newly calculated pressure
    pr_tt = pr_sum+pr;

% CALCULATION OF Eii AVERAGE
    % Calculation of initial eii
    if (count<2)
        eii_sum = zeros(size(eii));
    else
        eii_sum = eii_tt;
    end
    % Implementation of newly calculated eii
    eii_tt = eii_sum+eii;

% CALCULATION OF Sii AVERAGE
    % Calculation of initial sii
    if (count<2)
        sii_sum = zeros(size(sii));
    else
        sii_sum = sii_tt;
    end
    % Implementation of newly calculated sii
    sii_tt = sii_sum+sii;

% Update counter
	count = count+1;
    
end

% Calculation of pressure average
pr_average = pr_tt/count;
%Calculation of Eii average
eii_average = eii_tt/count;
%Calculation of Sii average
sii_average = sii_tt/count;


%% NEW COLOR SCALE BLUE WHITE RED

bwr = [0,0,1
    0.01	,	0.01	,	1
0.02	,	0.02	,	1
0.03	,	0.03	,	1
0.04	,	0.04	,	1
0.05	,	0.05	,	1
0.06	,	0.06	,	1
0.07	,	0.07	,	1
0.08	,	0.08	,	1
0.09	,	0.09	,	1
0.1	,	0.1	,	1
0.11	,	0.11	,	1
0.12	,	0.12	,	1
0.13	,	0.13	,	1
0.14	,	0.14	,	1
0.15	,	0.15	,	1
0.16	,	0.16	,	1
0.17	,	0.17	,	1
0.18	,	0.18	,	1
0.19	,	0.19	,	1
0.2	,	0.2	,	1
0.21	,	0.21	,	1
0.22	,	0.22	,	1
0.23	,	0.23	,	1
0.24	,	0.24	,	1
0.25	,	0.25	,	1
0.26	,	0.26	,	1
0.27	,	0.27	,	1
0.28	,	0.28	,	1
0.29	,	0.29	,	1
0.3	,	0.3	,	1
0.31	,	0.31	,	1
0.32	,	0.32	,	1
0.33	,	0.33	,	1
0.34	,	0.34	,	1
0.35	,	0.35	,	1
0.36	,	0.36	,	1
0.37	,	0.37	,	1
0.38	,	0.38	,	1
0.39	,	0.39	,	1
0.4	,	0.4	,	1
0.41	,	0.41	,	1
0.42	,	0.42	,	1
0.43	,	0.43	,	1
0.44	,	0.44	,	1
0.45	,	0.45	,	1
0.46	,	0.46	,	1
0.47	,	0.47	,	1
0.48	,	0.48	,	1
0.49	,	0.49	,	1
0.5	,	0.5	,	1
0.51	,	0.51	,	1
0.52	,	0.52	,	1
0.53	,	0.53	,	1
0.54	,	0.54	,	1
0.55	,	0.55	,	1
0.56	,	0.56	,	1
0.57	,	0.57	,	1
0.58	,	0.58	,	1
0.59	,	0.59	,	1
0.6	,	0.6	,	1
0.61	,	0.61	,	1
0.62	,	0.62	,	1
0.63	,	0.63	,	1
0.64	,	0.64	,	1
0.65	,	0.65	,	1
0.66	,	0.66	,	1
0.67	,	0.67	,	1
0.68	,	0.68	,	1
0.69	,	0.69	,	1
0.7	,	0.7	,	1
0.71	,	0.71	,	1
0.72	,	0.72	,	1
0.73	,	0.73	,	1
0.74	,	0.74	,	1
0.75	,	0.75	,	1
0.76	,	0.76	,	1
0.77	,	0.77	,	1
0.78	,	0.78	,	1
0.79	,	0.79	,	1
0.8	,	0.8	,	1
0.81	,	0.81	,	1
0.82	,	0.82	,	1
0.83	,	0.83	,	1
0.84	,	0.84	,	1
0.85	,	0.85	,	1
0.86	,	0.86	,	1
0.87	,	0.87	,	1
0.88	,	0.88	,	1
0.89	,	0.89	,	1
0.9	,	0.9	,	1
0.91	,	0.91	,	1
0.92	,	0.92	,	1
0.93	,	0.93	,	1
0.94	,	0.94	,	1
0.95	,	0.95	,	1
0.96	,	0.96	,	1
0.97	,	0.97	,	1
0.98	,	0.98	,	1
0.99	,	0.99	,	1
1	,	1	,	1
1	,	0.99	,	0.99
1	,	0.98	,	0.98
1	,	0.97	,	0.97
1	,	0.96	,	0.96
1	,	0.95	,	0.95
1	,	0.94	,	0.94
1	,	0.93	,	0.93
1	,	0.92	,	0.92
1	,	0.91	,	0.91
1	,	0.9	,	0.9
1	,	0.89	,	0.89
1	,	0.88	,	0.88
1	,	0.87	,	0.87
1	,	0.86	,	0.86
1	,	0.85	,	0.85
1	,	0.84	,	0.84
1	,	0.83	,	0.83
1	,	0.82	,	0.82
1	,	0.81	,	0.81
1	,	0.8	,	0.8
1	,	0.79	,	0.79
1	,	0.78	,	0.78
1	,	0.77	,	0.77
1	,	0.76	,	0.76
1	,	0.75	,	0.75
1	,	0.74	,	0.74
1	,	0.73	,	0.73
1	,	0.72	,	0.72
1	,	0.71	,	0.71
1	,	0.7	,	0.7
1	,	0.69	,	0.69
1	,	0.68	,	0.68
1	,	0.67	,	0.67
1	,	0.66	,	0.66
1	,	0.65	,	0.65
1	,	0.64	,	0.64
1	,	0.63	,	0.63
1	,	0.62	,	0.62
1	,	0.61	,	0.61
1	,	0.6	,	0.6
1	,	0.59	,	0.59
1	,	0.58	,	0.58
1	,	0.57	,	0.57
1	,	0.56	,	0.56
1	,	0.55	,	0.55
1	,	0.54	,	0.54
1	,	0.53	,	0.53
1	,	0.52	,	0.52
1	,	0.51	,	0.51
1	,	0.5	,	0.5
1	,	0.49	,	0.49
1	,	0.48	,	0.48
1	,	0.47	,	0.47
1	,	0.46	,	0.46
1	,	0.45	,	0.45
1	,	0.44	,	0.44
1	,	0.43	,	0.43
1	,	0.42	,	0.42
1	,	0.41	,	0.41
1	,	0.40	,	0.40
1	,	0.39	,	0.39
1	,	0.38	,	0.38
1	,	0.37	,	0.37
1	,	0.36	,	0.36
1	,	0.35	,	0.35
1	,	0.34	,	0.34
1	,	0.33	,	0.33
1	,	0.32	,	0.32
1	,	0.31	,	0.31
1	,	0.30	,	0.30
1	,	0.29	,	0.29
1	,	0.28	,	0.28
1	,	0.27	,	0.27
1	,	0.26	,	0.26
1	,	0.25	,	0.25
1	,	0.24	,	0.24
1	,	0.23	,	0.23
1	,	0.22	,	0.22
1	,	0.21	,	0.21
1	,	0.20	,	0.20
1	,	0.19	,	0.19
1	,	0.18	,	0.18
1	,	0.17	,	0.17
1	,	0.16	,	0.16
1	,	0.15	,	0.15
1	,	0.14	,	0.14
1	,	0.13	,	0.13
1	,	0.12	,	0.12
1	,	0.11	,	0.11
1	,	0.10	,	0.10
1	,	0.09	,	0.09
1	,	0.08	,	0.08
1	,	0.07	,	0.07
1	,	0.06	,	0.06
1	,	0.05	,	0.05
1	,	0.04	,	0.04
1	,	0.03	,	0.03
1	,	0.02	,	0.02
1	,	0.01	,	0.01
    1,0,0];

%% SUB-SECTION FOR VISUALIZATION:

% Draw
figure(1);clf;

% CHOOSING COLOR SCALE
% flipup COMMAND = REVERSE COLOR SCALE
colormap('jet');
%colormap(bwr);

% 1ST PLOT
%subplot(211)
%fig1 = subplot(211);

%figure(kslide);
%pcolor(gx/1000,gy/1000,log10(sii_average));caxis([6 8.7]);
%pcolor(gx/1000,gy/1000,log10(sxya));caxis([6 9]);
%pcolor(gx/1000,gy/1000,log10(eii_average));caxis([-16 -12]);
%pcolor(gx/1000,gy/1000,log10(nu));caxis([18 25]);
%pcolor(gx/1000,gy/1000,ro);caxis([2500 3500]);
%pcolor(gx(200:600)/1000,gy/1000,log10(eii(:,200:600))); caxis([-18 -12]);
%pcolor(pr);%caxis([0 60])
%pcolor(gx/1000,gy/1000,tk-273.15);
%pcolor(gx/1000,gy/1000,gp); %caxis([2800 3700]);
%pcolor(gx/1000,gy/1000,vy*100*3600*24*365);caxis([-0.2 0.2]);
%pcolor(gx/1000,gy/1000,vx*100*3600*24*365);caxis([-0.05 0.05]);
pcolor(gx/1000,gy/1000,(sqrt((vx.^2)+(vy.^2)))*100*3600*24*365);caxis([0 0.1]);
shading interp; axis ij image; colorbar; title([fileprn,' ',num2str(timesum)]);
xlim([xmin xmax]);ylim([ymin ymax]); hold on;
%ISOTHERMS
isotherms = [100:100:1300]; [c,h]=contour(gx/1000,gy/1000,tk-273.15,isotherms,'w','LineWidth',0.5); hold on;
% VELOCITY FIELD
stp = 6;
quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vx(1:stp:end,1:stp:end),vy(1:stp:end,1:stp:end),6,'K');

   
% 2ND PLOT
%subplot(212)
%fig2 = subplot(212);

%figure(kslide);
%pcolor(gx/1000,gy/1000,log10(sii));caxis([6 9]);
%WITH sxy (absolute value of shear stress)
%sxya = abs(sxy);
%pcolor(gx/1000,gy/1000,log10(sxya));caxis([6 9]);
%pcolor(gx/1000,gy/1000,log10(eii));caxis([-17 -12]);
%pcolor(gx/1000,gy/1000,log10(nu));caxis([18 25]);
%pcolor(gx/1000,gy/1000,ro);caxis([2500 3500]);
%pcolor(gx(200:600)/1000,gy/1000,log10(eii(:,200:600))); caxis([-18 -12]);
%pcolor(pr);%caxis([0 60])
%pcolor(gx/1000,gy/1000,tk-273.15);
%pcolor(gx/1000,gy/1000,gp); %caxis([2800 3700]);
%pcolor(gx/1000,gy/1000,vy*100*3600*24*365);caxis([-0.1 0.1]);
%pcolor(gx/1000,gy/1000,vx*100*3600*24*365);caxis([-0.1 0.1]);
%shading interp; axis ij image; colorbar; title([fileprn,' ',num2str(timesum)]);
%xlim([xmin xmax]);ylim([ymin ymax]); hold on;
%ISOTHERMS
%isotherms = [100:100:1300];
%[c,h]=contour(gx/1000,gy/1000,tk-273.15,isotherms,'w','LineWidth',0.5);
%hold on;
% VELOCITY FIELD
%stp = 6;
%quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vx(1:stp:end,1:stp:end),vy(1:stp:end,1:stp:end),5,'K');

% =========
%   ZOOMING AREA
%xlim([730 1030])
%ylim([00 80])
% =========

% CHOOSING COLOR SCALE
% flipup COMMAND = REVERSE COLOR SCALE
%colormap(fig1,'jet');
%colormap(fig2,bwr);
%colormap(flipud(hot))

% %Image file
% filename    =  ['fig',name,num2str(nfl)];
%   print ('-dtiff', '-r300','-zbuffer ',filename);

%% Pressure gradient
    [Px Py] = gradient(pr_average);
    % Pressure gradient in Pa/m
    dpdx = (Px/500);
    dpdy = (Py/500);
    % Water velocity calculation (cf. markevps.c /* Fluid in rock */ )
        % Wt velocity coef. calculation
        vykoef = (1000*9.80665-dpdy)/(2300*9.81);
        vxkoef = ((1000*0-dpdx)/(2300*9.81));
        % Max of 2 m/s of pressure-gradient-driven water velocity
        if(vxkoef>2) 
            vxkoef=2;
        elseif(vxkoef<-2)
            vxkoef=-2;
        end
        if(vykoef>2) 
            vykoef=2;
        elseif(vykoef<-2)
            vykoef=-2;
        end
        % Final water velocity calculation (in m/yr)
        vxwtmark=vxkoef*vfluid;
        vywtmark=vykoef*vfluid;
        vwtmark=sqrt((vxwtmark.^2)+(vywtmark.^2));
        vwtmark=vwtmark*3600*24*365;

% DEFINING INITIAL LOCATION OF STREAMLINES 1
% Extract coordinates of high-strain zone (i.e. subduction channel)
[col,row] = find(eii_average>=1e-13);
ggx(:,1) = gx(row)/1e3;
ggx(:,2) = gy(col)/1e3;
datanewa = ggx(ggx(:,1)>xmin,:);
datanewa = datanewa(datanewa(:,1)<xmax,:);
datanewa = datanewa(datanewa(:,2)>ymin,:);
datanewa = datanewa(datanewa(:,2)<ymax,:);

% DEFINING INITIAL LOCATION OF STREAMLINES 2
ggx2(:,1) = gx(row)/1e3;
ggx2(:,2) = (gy(col)/1e3)+str2;
datanewb = ggx2(ggx2(:,1)>xmin,:);
datanewb = datanewb(datanewb(:,1)<xmax,:);
datanewb = datanewb(datanewb(:,2)>ymin,:);
datanewb = datanewb(datanewb(:,2)<ymax,:);

% Divide by ... the number of streamlines
datanewa = datanewa(1:3:end,:);
datanewb = datanewb(1:3:end,:);


%% Drawing Fluid flow
figure(2);
%colormap(bwr);
colormap(jet);
%colormap(flipud(hot));

%pcolor(gx/1000,gy/1000,dpdy/1e6);caxis([-0.2 0.2]);
%pcolor(gx/1000,gy/1000,log10(vwtmark));caxis([-3 -0.8]);
pcolor(gx/1000,gy/1000,log10(sii_average));caxis([6 9]);
shading interp;axis ij image;colorbar;
title([fileprn,' - ',num2str(timesum/1e6),' ','Myr','  -  Grad-Py (MPa/m)']);hold on;
%isotherms = [100:100:1300];
%[c,h]=contour(gx/1000,gy/1000,tk-273.15,isotherms,'k','LineWidth',1);hold on;
xlim([xmin xmax]);ylim([ymin ymax]);
% WATER VELOCITY FIELD
%stp=2;
%quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vxwtmark(1:stp:end,1:stp:end),vywtmark(1:stp:end,1:stp:end),250,'K');
% WATER VELOCITY STREAMLINES
    % Streamlines "base of interface"
    streamfctb = streamline(gx/1000,gy/1000,vxwtmark,vywtmark,datanewb(:,1),datanewb(:,2));
    set(streamfctb,'LineWidth',1,'Color','k');
    % Streamlines "top of interface"
    streamfcta = streamline(gx/1000,gy/1000,vxwtmark,vywtmark,datanewa(:,1),datanewa(:,2));
    set(streamfcta,'LineWidth',1,'Color','k');
    
%Image file
%    filename = ['Fig_GradPy',nname,num2str(num_name)];
%    print ('-dtiff', '-r300',filename);

% Activation if not average of pressure
%end
