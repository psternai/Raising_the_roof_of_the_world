clear all

% Counter
count = 1;

% __________ INPUT DATA __________ %

    % Unique file or multiple files
       for num_name = 10:10:610
        
    % Axis limits (km)
        Xmin = 2750 ;
        Xmax = 3500 ;
        Ymin = 0 ;
        Ymax = 300 ;
        
    % Horizontal buffer of channel properties (in km)
        buff = 5;
        
% __________ INPUT DATA __________ %

%File name
    nname     = '2Ddata'; 
    ext       = '.txt';
% File
    file     = [nname,num2str(num_name)]   
%read input data
    fname    =  [nname,num2str(num_name) ext];
    fid      =  fopen(fname,'r');
% Grid resolution
xnumx=fread(fid,1,'int64');
ynumy=fread(fid,1,'int64');
% Time
timesum=fread(fid,1,'float64');
% Read Gridline positions
gx=fread(fid,xnumx,'float32');
gy=fread(fid,ynumy,'float32');
% Read nodes information
for i=1:1:xnumx
    for j=1:1:ynumy
        nu(j,i)=fread(fid,1,'float32');
        ro(j,i)=fread(fid,1,'float32');       
        vx(j,i)=fread(fid,1,'float32');
        vy(j,i)=fread(fid,1,'float32');
        exx(j,i)=fread(fid,1,'float32');
        exy(j,i)=fread(fid,1,'float32');
        sxx(j,i)=fread(fid,1,'float32');
        sxy(j,i)=fread(fid,1,'float32');
    end
end

% Calculation of Eii and Sii
eii=ones(ynumy,xnumx)*1e-16;
sii=ones(ynumy,xnumx)*1e+4;
for i=1:1:xnumx-2
    for j=1:1:ynumy-2
        eii(j+1,i+1)=(exy(j+1,i+1)^2+((exx(j+1,i+1)+exx(j+2,i+1)+exx(j+1,i+2)+exx(j+2,i+2))/4)^2)^0.5;
        sii(j+1,i+1)=(sxy(j+1,i+1)^2+((sxx(j+1,i+1)+sxx(j+2,i+1)+sxx(j+1,i+2)+sxx(j+2,i+2))/4)^2)^0.5;
    end
end

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
    % Implementation of newly calculated eii
    sii_tt = sii_sum+sii;
    
% CALCULATION OF Sxy (SHEAR STRESS) AVERAGE
    if (count<2)
        sxy_sum = zeros(size(abs(sxy)));
    else
        sxy_sum = sxy_tt;
    end
    % Implementation of newly calculated sii
    sxy_tt = sxy_sum+abs(sxy);

% activation if integration of averaged parameter
    % Update of counter    
        count = count+1;

%Calculation of Eii average
eii_average = eii_tt/(count-1);
%Calculation of Sii average
sii_average = sii_tt/(count-1);
%Calculation of Sxy average
sxy_average = sxy_tt/(count-1);

%% NEW COLOR SCALE

% BLUE-WHITE-RED
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

%% VISUALIZATION 1

% DRAWING 1
figure(1);
%subplot(3,1,1);
%pcolor(gx/1000,gy/1000,log10(nu));caxis([17 26]);
%pcolor(gx/1000,gy/1000,ro);caxis([2200 3300]);
pcolor(gx/1000,gy/1000,log10(eii)); caxis([-18 -13]);
%pcolor(gx/1000,gy/1000,log10(eii_average));caxis([-16 -12]);
%pcolor(gx/1000,gy/1000,vy*100*3600*24*365);caxis([-0.3 0.3]);
%pcolor(gx/1000,gy/1000,log10(sii));caxis([6 8.7]);
%colormap(gca,jet);
%colormap(gca,flipud(hot));
%shading interp; axis ij image; set(gca,'fontsize',6); colorbar; 
%title([nname,num2str(num_name),' - ',num2str(timesum/1e6),' ','Myr','  -  Rho Eii Vx Vy']);hold on;
%xlim([Xmin Xmax]); ylim([Ymin Ymax]);
% VELOCITY FIELD
%stp = 5;
%quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vx(1:stp:end,1:stp:end),vy(1:stp:end,1:stp:end),7,'K');

%subplot(3,1,2);
%pcolor(gx/1000,gy/1000,log10(eii));caxis([-16 -12]);
%pcolor(gx/1000,gy/1000,log10(sii));caxis([6 8.7]);
%pcolor(gx/1000,gy/1000,log10(sii_average));caxis([6 8.7]);
%pcolor(gx/1000,gy/1000,vx*100*3600*24*365);caxis([-0.3 0.3]);
%pcolor(gx/1000,gy/1000,vy*100*3600*24*365);caxis([-0.2 0.2]);
% Color scale
%    colormap(gca,jet);
    colormap(gca,summer);
%    colormap(gca,bwr);
% Required field
    shading interp; 
    axis ij image;set(gca,'fontsize',12);hold on;xlim([Xmin Xmax]);ylim([Ymin Ymax]);% colorbar;
    c = colorbar;
    c.Label.String = 'sec.inv.str.rate (s-1)';
% VELOCITY FIELD
    stp = 5;
    quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vx(1:stp:end,1:stp:end),vy(1:stp:end,1:stp:end),2,'K');
hold off;

 % Save a .jpg figure file (yes = 1 ; no = 0)
        Exp = 1;
        filename    =  ['2Ddata_',num2str(num_name)];
% For printing file ONLY
    if Exp >= 1
        print ('-djpeg', '-r200',filename);
        num_name    =  num_name+1;
    end

    end

%subplot(3,1,3);
%pcolor(gx/1000,gy/1000,log10(sii));caxis([6 9]);
%pcolor(gx/1000,gy/1000,vy*100*3600*24*365);caxis([-0.3 0.3]);
%colormap(gca,bwr);
%shading interp; axis ij image; colorbar; set(gca,'fontsize',6); hold on;
%xlim([Xmin Xmax]); ylim([Ymin Ymax]);
% VELOCITY FIELD
%stp = 5;
%quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vx(1:stp:end,1:stp:end),vy(1:stp:end,1:stp:end),15,'K');

%subplot(4,1,4);
%pcolor(gx/1000,gy/1000,log10(eii));caxis([-17 -12]);
%pcolor(gx/1000,gy/1000,log10(sii));caxis([6 9]);
%pcolor(gx/1000,gy/1000,vx*100*3600*24*365);caxis([-0.3 0.3]);
%pcolor(gx/1000,gy/1000,vy*100*3600*24*365);caxis([-0.3 0.3]);
%colormap(gca,bwr);
%shading interp; axis ij image; colorbar; set(gca,'fontsize',6); hold on;
%xlim([Xmin Xmax]); ylim([Ymin Ymax]);
% VELOCITY FIELD
%stp = 6;
%quiver(gx(1:stp:end)/1000,gy(1:stp:end)/1000,vx(1:stp:end,1:stp:end),vy(1:stp:end,1:stp:end),8,'K');

% PRINT DATA IN OUTPUT FILE
%filename = ['Fig_',nname,num2str(num_name)];
%    print ('-dtiff', '-r300',filename);


%     
% % "end" activation if not integration of averaged parameter
% %       end
% 
% %% PROPERTIES OF SUBDUCTION INTERFACE
% 
% % Extract coordinates of high-strain zone (i.e. subduction channel)
%     [col,row] = find(eii_average>=1e-13);
%     ggx(:,1) = gx(row)/1e3;
%     ggx(:,2) = gy(col)/1e3;
%     coordChannel = ggx(ggx(:,1)>Xmin,:);
%     coordChannel = coordChannel(coordChannel(:,1)<Xmax,:);
%     coordChannel = coordChannel(coordChannel(:,2)>Ymin,:);
%     coordChannel = coordChannel(coordChannel(:,2)<Ymax,:);
% % Index of X-Y coordinates of high-strain subduction channel
%     [~,idxX] = ismember(coordChannel(:,1),gx/1e3);
%     [~,idxY] = ismember(coordChannel(:,2),gy/1e3);
% % Shear stress in high-strain subduction channel 
%     channelSxy = diag(sxy_average(idxY,idxX));
% % Vertical average of shear stress in subduction channel
%     [ud,ix,iy]=unique(coordChannel(:,1));  
%     VertmeanSxychannel = accumarray(iy,channelSxy,[],@mean);
% % For high-resolution area (node = 500 m), number of nodes for buffering the Sxy 
%     buffnode = ((buff*1000)/500);
% % Sliding window average to buffer Sxy values
%     for i = 1:length(VertmeanSxychannel)
%         meanSxychannel(i) = mean(VertmeanSxychannel(max([1 i-buffnode]):(min([i+buffnode length(VertmeanSxychannel)]))));
%     end
%         
% % Plot Sii in subduction channel
% figure(2);
% plot(ud,meanSxychannel/1e6);
% xlabel('X axis (in km)'),ylabel('Shear Stress (in MPa)'); xlim([Xmin Xmax]); grid on;