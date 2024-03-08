%Material visualization

% C_map - color map used for this visualization
%  0 - air (0)                   white           [1  1  1]
%  1 - water (1)                 azure           [0.72941 0.85882 0.94902]
%  2 - sediments  (2-sedim.)     light orange    [0.97255 0.7098 0.47451]
%  3 - sediments2 (3-prism)      brown           [0.68235 0.34118 0]
%  4 - sediments3 (4-pelagic)    orange          [1 0.50196 0]
%  5 - continental crust2 (5)    light grey      [0.75294 0.75294 0.75294]
%  6 - continental crust2 (6) -  dark grey       [0.50196 0.50196 0.50196]
%  7 - basalt (7)                dark green      [0 0.50196 0]
%  8 - gabbro (8)                light green     [0 0.84314 0]
%  9 - dry mantel1 (9)           dark blue       [0 0 0.71765]
% 10 - dry mantel2 (10)          navy blue       [0.32157 0.2 0.67451]
% 11 - hydrated mantel1 (11)     blue            [0.54118 0.72157 0.99216]
% 12 - hydrated mantel2 (12)     bright blue     [0 0.50196 1]
% 13 - serpentinized mantle (13) dark navy blue  [0 0 0.4902]
% 14 - newer appear (solid part of partially molten peridotite)                                                                 ?????
% 15 -sediments1 (15)            dark red        [0.86275 0.78431 0.23922]
% 16 -sediments1 (16)            dark red         [0.83922 0.035294
% 0.011765] ( 0.8549     0.59608     0.36078) for molten - if ever)
% 17 -   wet upper crust (17)  [0.75294 0.75294 0.75294]
% 18 -   wet lower crust (18)  [0.50196 0.50196 0.50196]
% 19 -
% 20 -
% 21 -
% 22 -
% 23 - molten sediments1 (23)     light yellow   [1 1 0.31765]
% 24 - molten sediments2 (24)     yellow         [1 0.90196 0.18824]
% 25 - molten crust1 (25)        dirty olives    [0.46667     0.46667 0.23529]
% 26 - mmolten crust2 (26)       olivs             [0.50196     0.50196 0]
% 27 - molten basalt (27)        dark violet    [0.72549 0.015686 0.78431]
% 28 - molten gabbro (28)        light violet   [0.92549 0.43922 0.99608]
% 29 -
% 30 -
% 31 - molten mantle1 (31)       red1            [0.99216 0.38824 0.30196]
% 32 - molten mantle2  (32)      red2            [0.84706 0.078431 0.15294]
% 33 -
% 34 - molten peridotite (34)    red             [1 0 0]
% 35 - molten sediments1 (36)     light yellow   [0.8549 0.59608 0.36078
% 36 - molten sediments2 (37)     yellow         [1 0.90196 0.18824]
%%
clear all

% __________ INPUT DATA __________ %

    % For loop
    for nn1 = 30:10:610

    %defining begining and end of zoomed area
        x_beg  = 0;
        x_end  = 4000;
        z_beg  = 0;
        z_end  = 300;
    
    % Steps of isotherms (°C)
        temp_step = 200;
    
    % Save a .tif figure file (yes = 1 ; no = 0)
        Exp = 1;
        
% __________ INPUT DATA __________ %

nname     =  'clim';
num_name =  nn1;

%File for composition
name      =  [nname 'c'];
ext      =  '.txt';

%open figure of temperature
name_t      =  [nname 't'];
num_name_t  = num_name;
ext_t       =  '.txt';
% Defining size of the model
x_size = 4000;
z_size = 300;
%Defining zoom
zoom = 10;

%recalculating to % of the image
if x_beg ==0
    x_beg_perc = 1;
else
    x_beg_perc = x_beg/x_size;
end
x_end_perc = x_end/x_size;
if z_beg ==0
    z_beg_perc = 1;
else
    z_beg_perc = z_beg/z_size;
end
z_end_perc = z_end/z_size;

% main loop - insert numbers of images
for i=1:1

    fname    =  [name,num2str(num_name) ext];
    fname_t     =  [name_t,num2str(num_name_t) ext_t];

    %read input data
    fid      =   fopen(fname,'r');
    A        =   textscan(fid, '%n');
    fclose(fid);
    Data     =   A{1};
    Time     =   Data(1);
    Coord_x  =   Data(2);
    Coord_z  =   Data(3);
    Data_vec =   Data(4:end);

ColorGrid = ones(Coord_z,Coord_x);
    C_map=[1 1 1;...                    % (0) air (white)
        0.72941 0.85882 0.94902;...     % (1) water (light blue)
        0.97255 0.7098 0.47451;...      % (2)
        0.97255 0.7098 0.47451;...      % (3) terrigenous sedim (light orange)
        1 0.50196 0;...                 % (4) pelagic sedim (dark orange)
        0.75294 0.75294 0.75294;...     % (5) felsic crust (grey)
        0.50196 0.50196 0.50196;...     % (6) mafic crust (dark grey)
        0 0.50196 0;...                 % (7) basalt (green)
        0 0.84314 0;...                 % (8) gabbro (light green)
        0 0 0.71765;...                 % (9) lithospheric mantle (dark blue)           
        0.32157 0.2 0.67451;...         % (10) asthenosphere (navy blue)
        0.54118 0.72157 0.99216;...     % (11) hydrated mantle (blue)
        0 0.50196 1;...                 % (12) Weak zone (bright blue)
        0 0 0.4902;...                  % (13) Serpentinite (very-dark blue)
        0 0 0;...                       % (14)
        0 0 0;...                       % (15)
        0 0 0;...                       % (16)
        0.65294 0.65294 0.65294;...     % (17) wet felsic crust (grey2)
        0.40196 0.40196 0.40196;...     % (18) wet mafic crust (dark grey2) 
        0 0 0;...                       % (19)
        0 0 0;...                       % (20)
        0 0 0;...                       % (21)
        0 0 0;...                       % (22)
        1 1 0;...                       % (23) molten terrigenous(?) sedim (yellow)
        1 1 0;...                       % (24) molten pelagic sedim (yellow)
        1 0.7 0.2;...                   % (25) molten felsic crust (light orange)
        1 0.6 0.1;...                   % (26) molten mafic crust (orange)
        0.72549 0.015686 0.78431;...    % (27) molten basalt (dark violet)
        0.92549 0.43922 0.99608;...     % (28) molten gabbro (violet)
        1 0.33 0.29;...                 % (29) molten dry peridotite (light red)
        1 0.33 0.29;...                 % (30) molten dry peridotite (light red)
        0 0 0;...                       % (31)
        1 0.33 0.29;...                 % (32) molten dry peridotite (light red)
        0 0 0;...                       % (33)
        1 0 0                           % (34) Wet peridotite melt (red)
 % (35 = NaN; 36 = molten basalt; 37 = molten felsic crust BUT I CAN'T ADD IT...) 
        1 0.6 0.1];                     % (38) Molten mafic crust (orange) 
    
    %Fill the colorgrid with the data of opened file
    num       = 1;
    ind       = 1;
    while num<length(Data_vec)
        value = Data_vec(num);

        if value==-2
            % Compressed: the next ?? indices are given the color material
            num_colors  =   Data_vec(num+1);
            material    =   Data_vec(num+2);
            ind_vec     =   ind:ind+num_colors-1;

            ind         =   ind+num_colors;
            num         =   num+3;
        else
            if value==-1
                material    =   NaN;
            else
                material    =   value;
            end
            ind_vec = ind;
            ind         =   ind+1;
            num         =   num+1;
        end

        ColorGrid(ind_vec)  =   material;

    end

    if zoom >=0
        if (z_beg ==0 & x_beg ~= 0)
            ColorGrid = ColorGrid(1:round(z_end_perc*size(ColorGrid,1)),round(x_beg_perc*size(ColorGrid,2)):round(x_end_perc*size(ColorGrid,2)));
        elseif (x_beg ==0 & z_beg ~= 0)
            ColorGrid = ColorGrid(round(z_beg_perc*size(ColorGrid,1)):round(z_end_perc*size(ColorGrid,1)),1:round(x_end_perc*size(ColorGrid,2)));
        elseif x_beg ==0 & z_beg ==0
            ColorGrid = ColorGrid(1:round(z_end_perc*size(ColorGrid,1)),1:round(x_end_perc*size(ColorGrid,2)));
        else
            ColorGrid = ColorGrid(round(z_beg_perc*size(ColorGrid,1)):round(z_end_perc*size(ColorGrid,1)),round(x_beg_perc*size(ColorGrid,2)):round(x_end_perc*size(ColorGrid,2)));
        end
    else
        ColorGrid = ColorGrid;
    end
    %save memory
    clear Data Data_vec

    %read input data of temperature
    format long g
    fid     =   fopen (fname_t,'r');
    B       =   textscan (fid, '%n')';
    fclose(fid);

    Data    =   B{1};

    Coord_x_t =   Data(2);
    Coord_z_t =   Data(3);
    Data_vec_t=   Data(4:end);
    % Convertion of the tempertaute form Kelvin to Celcius
    Data_vec_t= Data_vec_t-273;
    ColorGrid_t = ones(Coord_z_t,Coord_x_t);

    ColorGrid_t = reshape(Data_vec_t,Coord_z_t,Coord_x_t);
    
    if zoom >=0
        if (z_beg ==0 & x_beg ~= 0)
            ColorGrid_t = ColorGrid_t(1:round(z_end_perc*size(ColorGrid_t,1)),round(x_beg_perc*size(ColorGrid_t,2)):round(x_end_perc*size(ColorGrid_t,2)));
        elseif (x_beg ==0 & z_beg ~= 0)
            ColorGrid_t = ColorGrid_t(round(z_beg_perc*size(ColorGrid_t,1)):round(z_end_perc*size(ColorGrid_t,1)),1:round(x_end_perc*size(ColorGrid_t,2)));
        elseif x_beg ==0 & z_beg ==0
            ColorGrid_t = ColorGrid_t(1:round(z_end_perc*size(ColorGrid_t,1)),1:round(x_end_perc*size(ColorGrid_t,2)));
        else
            ColorGrid_t = ColorGrid_t(round(z_beg_perc*size(ColorGrid_t,1)):round(z_end_perc*size(ColorGrid_t,1)),round(x_beg_perc*size(ColorGrid_t,2)):round(x_end_perc*size(ColorGrid_t,2)));
        end
    else
        ColorGrid_t = ColorGrid_t;
    end

 % Create x and z grid - for the scale control
    if zoom >= 0
        H   =       z_end - z_beg;
        W   =       x_end - x_beg;
    else
        H       =       z_size;
        W       =       x_size;
    end
    Nx      =       size(ColorGrid,2);
    Nz      =       size(ColorGrid,1);
    dW      =       W/(Nx-1);
    dH      =       H/(Nz-1);
    x_vec   =       [x_beg:dW:x_end];
    z_vec   =       [z_beg:dH:z_end];


    % Create x and z grid - for the scale control of temperature
   if zoom >= 0
        H   =       z_end - z_beg;
        W   =       x_end - x_beg;
    else
        H       =       z_size;
        W       =       x_size;
    end
    Nx      =       size(ColorGrid_t,2);
    Nz      =       size(ColorGrid_t,1);
    dW      =       W/(Nx-1);
    dH      =       H/(Nz-1);
    x_vec_t   =       [x_beg:dW:x_end];
    z_vec_t   =       [z_beg:dH:z_end];


    % plotting
    figure(1);clf
    step=1;
    step2=1;
    colormap(C_map);
%      [C,h] = contour(x_vec_t(1:step2:end),z_vec_t(1:step2:end),ColorGrid_t(1:step2:end,1:step2:end),...
%          [100:200:1500],'r','LineWidth',0.5);
     % contour(x_vec_t(1:step2:end),z_vec_t(1:step2:end),ColorGrid_t(1:step2:end,1:step2:end),[100:200:1500],'r','LineWidth',0.5);
%     h = gcbo;
    pcolor(x_vec,z_vec,ColorGrid);
    shading flat, colormap, caxis ([0 35]), axis ij; 
    set(gca,'fontsize',5); %image
    hold on;


    [c,h]=contour(x_vec_t(1:step2:end),z_vec_t(1:step2:end),ColorGrid_t(1:step2:end,1:step2:end)/1000,...
        [100:temp_step:1500]/1000,'w','LineWidth',0.2);

    axis ij image;

%     % geting coordinates for isoterme values
%     val   =   c(1,:);
%     [i,j] =   find (val == 26);
%     t     =   ones (1,length(j));
%     val   =   c(2,:);
%     for i =   1:length(t)
%         t(i) = val(j(i));
%     end
%     t    =  ceil (t);
%     t    =  t-2;
%     % value the isoterms
%     T=100;
%     for i = 1:length(t)
%         text (15,t(i),num2str(T),'Color','w','FontSize',8,'FontAngle','italic');
%         T=T+200;
%     end

    hold off

    xlabel('Width (km)','FontSize',10);
    ylabel('Depth (km)','FontSize',10);
    title(['Time = ', num2str(round(Time/1e6,1)),' Myr'],'FontSize',20);
    
    if zoom >= 0
        filename    =  ['Fig_',name,num2str(num_name)];
    else
        filename    =  ['fig_',name,num2str(num_name)];
    end

    
% For printing file ONLY
    if Exp >= 1
        print ('-djpeg', '-r200',filename);
        num_name    =  num_name+1;
        num_name_t    =  num_name_t+1;
    end

        %save("climc340_tmp.txt","ColorGrid","-ascii")

end
    end

