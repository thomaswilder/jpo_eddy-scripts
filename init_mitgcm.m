clear all; close all

% Initial conditions for MITgcm baroclinic eddy

%% Setting up grid

% vertical grid
H = 4000; % ocean depth in m
z_levels  = ncread('ocean_vertical_grid.nc','v_grid');
z_levels(length(z_levels)+1) = H;

% vertical grid spacing
for i = 1:length(z_levels)-1 
    delta_z(i) = z_levels(i+1) - z_levels(i);
end
nz = length(delta_z); % number of vertical grid levels

% horizontal grid
nx=200; ny=200; % number of horizontal points
hx=10; hy=10; % grid spacing in km
[XG, YG] = ndgrid(-0.5*nx*hx:hx:0.5*nx*hx-hx,-0.5*ny*hy:hy:0.5*ny*hy-hy); % corner of grid cell
% [XC, YC] = ndgrid(-0.5*nx*hx+0.5*hx:hx:0.5*nx*hx-0.5*hx,...
%     -0.5*ny*hy+0.5*hy:hy:0.5*ny*hy-0.5*hy); % center of grid cell

%% Model parameters

A = 0.25; % eddy amplitude
R = 100; % eddy radius
f = 2*7.27*1e-5*sin((pi/18)*4);

g = 9.81;

D = 800; % modifies decay of temperature profile
H1 = 800; % depth of first layer
T0 = 18; % constant reference temperature
T_anom = 2.5; % temperature anomaly
alpha = 2e-4; % thermal constant

rho_1 = 1026; rho_2 = 1029; % density of each layer
N = 1e-5; % constant buoyancy frequency

% Vertical density structure
gam_e = 1;
A_bg = 3;

%% 3D temperature field

% temperature is found by inverting density using linear equation of state

% density
rho_z_b = zeros(1,nz);
for i=1:nz+1
   rho_z_b(i) = rho_1*(1-((N^2)*(-z_levels(i)))/g)+...
       (0.5)*abs(rho_2-rho_1)*(1-tanh(A_bg*(-z_levels(i)+H1)/H));
end

% background temperature
T_z_b = -(((rho_z_b-rho_1)/(alpha*rho_1))-T0);

% 3D temp field
T = zeros(nx,ny,nz);
for i = 1:nz 
    T(:,:,i) = T_anom*exp(-(XG.^2+YG.^2)/R^2)...
        *exp(gam_e*(-z_levels(i))/(D)) + T_z_b(i);
end

%% Horizontal velocities in thermal wind balance

% sea surface height using Gaussian function found at center of grid cell
eta = A*exp(-(XG.^2+YG.^2)/R^2);

% geostrophic velocities
[u, v] = geostrophic_vels(T,eta,f,alpha,hx,hy,z_levels);

%% RBCS sponge layer for temperature

% tukeywin requires signal processing toolbox

diam = nx; % Diameter of tukey window
fall = 0.1; % Falloff of function. 1 is equal to hanning window
wc=window(@tukeywin,diam,fall);
wr=window(@tukeywin,diam,fall);
[maskr,maskc] = meshgrid(wr,wc);
tukey2D = maskr.*maskc;

for k = 1:nz
    sp_mask(:,:,k) = ones(nx,ny) - tukey2D;
end

%% Topography 

h=-H*ones(nx,ny);


%% Wind vector

% vector rotates one whole period every 64 hours
nrot = 64;
uwind = zeros(nx,ny,nrot);
vwind = zeros(nx,ny,nrot);

% initial wind vector and wind speed in m/s
uwind(:,:,1) = 7; 

for i = 2:nrot
   uwind(:,:,i) = cos(-(2*pi)/nrot)*uwind(:,:,i-1)-sin(-(2*pi)/nrot)*vwind(:,:,i-1);
   vwind(:,:,i) = sin(-(2*pi)/nrot)*uwind(:,:,i-1)+cos(-(2*pi)/nrot)*vwind(:,:,i-1);
end

%% Write data to binary files for MITgcm read

ieee='b';
accuracy='real*8';

fid=fopen('data/sponge.bin','w',ieee); fwrite(fid,sp_mask,accuracy); fclose(fid);

fid=fopen('data/topog.box','w',ieee); fwrite(fid,h,accuracy); fclose(fid);

fid=fopen('data/uwind.bin','w',ieee); fwrite(fid,uwind,accuracy); fclose(fid);
fid=fopen('data/vwind.bin','w',ieee); fwrite(fid,vwind,accuracy); fclose(fid);

fid=fopen('data/eta.bin','w',ieee); fwrite(fid,eta,accuracy); fclose(fid);

fid=fopen('data/uvel.bin','w',ieee); fwrite(fid,u,accuracy); fclose(fid);
fid=fopen('data/vvel.bin','w',ieee); fwrite(fid,v,accuracy); fclose(fid);

fid=fopen('data/Temp.bin','w',ieee); fwrite(fid,T,accuracy); fclose(fid);
fid=fopen('data/Temp-bg.bin','w',ieee); fwrite(fid,T_z_b,accuracy); fclose(fid);

