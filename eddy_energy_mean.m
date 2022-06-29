clear all;

% Calculates mean energy terms for the eddy following von Storch et al (2012)
% using geostrophic horizontal velocity components.
% Calculations are either 4D or volume integrals, depending on choice.
% Can be either used locally or in a hpc environment.
% Run time will be extended for resolution greater than 10 km.
% Uses a 16 day time-mean to avoid picking up any signals in the
% instability pathways from the rotating wind

% Variables calculated:
% Km - mean kinetic energy
% Km_bt, Km_bc - barotropic, baroclinic kinetic energy
% Pm - mean potential enegry
% PmKm - conversion of potential to kinetic energy
% Wm - total wind power input

%% Setting up grid

% vertical grid
H = 4000; % ocean depth in m
z_levels  = ncread('ocean_vertical_grid.nc','v_grid');
nz = length(z_levels); % number of vertical grid levels

z_levels(length(z_levels)+1) = H;

% horizontal grid resolution
nx = 200; ny = 200; % number of horizontal points
hx = 10; hy = 10; % grid spacing in km
ht = 400; % time step


%% Model parameters

A = 25; % eddy amplitude
R = 100; % eddy radius
f = 2*7.27*1e-5*sin((pi/18)*4); % coriolis parameter
g = 9.81;

rho_1 = 1026; % reference density
T0 = 18; % constant reference temperature

A4 = (0.25*0.125*(hx*1e+3)^4*0.025)/ht; % viscous coefficient

%% Other variables

% choose wind stress and dimensions of energy terms by modifying 'str#'
% wind stress
relative = "relative";
str1 = 'nowind';

% 4D or totals
plan = "plan"; % 4D plan view
total = "total";
str2 = 'total'; 

% days
if contains(str2,total)
    % number of days for energy totals
    n = [31:1:380]; % needs to start at day 31 to run a 16 day time-mean
    nlvs = nz; % number of levels to analyse
elseif contains(str2,plan)
     n = [35 75 100 125 150 175 200 225]; % days for 4D terms
%    n = 200;
    nlvs = nz; % number of levels to analyse
end

% vertical grid spacing
for i = 1:nlvs
    delta_z(i) = z_levels(i+1) - z_levels(i);
end

% number of days used to calculate mean
pt_mean = 16;

% data paths
fn1 = 'glue_data/ACE_eta_%s_%dkm_A%d.nc';
fn2 = 'glue_data/ACE_temp_%s_%dkm_A%d.nc';
fn3 = 'glue_data/ACE_wvel_%s_%dkm_A%d.nc';
fn4 = 'glue_data/ACE_tau_%s_%dkm_A%d.nc'; 

% load background density, or rho_z_b from init_mitgcm.m
load('rho_ref'); % needs to be nz+1 levels to calculate n0(k)

% vertical derivative of background density
n0 = zeros(1,nlvs);
for k = 1:nlvs
    n0(k) = (rho_z_b(k)-rho_z_b(k+1))/delta_z(k); 
end

% start at day 11 of model iteration, first 10 days are used for adjustment
iters = ncread(sprintf(fn1,str1,hx,A),'iter');
id = find(iters == (86400/ht)*11);



%% Setup netcdf file. Attributes currently omitted.

% open netCDF file.
filename = 'data/ACE_energy_%s_%s_%dkm_A%d_mean_beta.nc';
ncid = netcdf.create(sprintf(filename,str1,str2,hx,A),'NETCDF4'); % NC_NOWRITE

if contains(str2, total)
    
    % define the dimensions of the variable.
    z_dimID = netcdf.defDim(ncid,'Z',nlvs);
    t_dimID = netcdf.defDim(ncid,'T',netcdf.getConstant('NC_UNLIMITED')); % 
    data_IDs = ([z_dimID t_dimID]);

    % define a new variable in the file.
    zID = netcdf.defVar(ncid,'Z','double',z_dimID);
    tID = netcdf.defVar(ncid,'T','double',t_dimID);

    % define output variables 
    sum_KmID = netcdf.defVar(ncid,'sum_Km','NC_FLOAT',[z_dimID t_dimID]);
    tot_KmID = netcdf.defVar(ncid,'tot_Km','NC_FLOAT',t_dimID);
    
    sum_Km_btID = netcdf.defVar(ncid,'sum_Km_bt','NC_FLOAT',[z_dimID t_dimID]);
    tot_Km_btID = netcdf.defVar(ncid,'tot_Km_bt','NC_FLOAT',t_dimID);
    
    sum_Km_bcID = netcdf.defVar(ncid,'sum_Km_bc','NC_FLOAT',[z_dimID t_dimID]);
    tot_Km_bcID = netcdf.defVar(ncid,'tot_Km_bc','NC_FLOAT',t_dimID);

    sum_PmID = netcdf.defVar(ncid,'sum_Pm','NC_FLOAT',[z_dimID t_dimID]);
    tot_PmID = netcdf.defVar(ncid,'tot_Pm','NC_FLOAT',t_dimID);

    sum_PmKmID = netcdf.defVar(ncid,'sum_PmKm','NC_FLOAT',[z_dimID t_dimID]);
    tot_PmKmID = netcdf.defVar(ncid,'tot_PmKm','NC_FLOAT',t_dimID);

    if contains(str1, relative)

     	tot_WmID = netcdf.defVar(ncid,'tot_Wm','NC_FLOAT',t_dimID);
    end
    
elseif contains(str2, plan)
    
    % define the dimensions of the variable.
    x_dimID = netcdf.defDim(ncid,'X',nx);
    y_dimID = netcdf.defDim(ncid,'Y',ny);
    z_dimID = netcdf.defDim(ncid,'Z',nlvs);
    t_dimID = netcdf.defDim(ncid,'T',netcdf.getConstant('NC_UNLIMITED')); % 
    data_IDs = ([x_dimID y_dimID z_dimID t_dimID]);

    % define a new variable in the file.
    xID = netcdf.defVar(ncid,'X','double',x_dimID);
    yID = netcdf.defVar(ncid,'Y','double',y_dimID);
    zID = netcdf.defVar(ncid,'Z','double',z_dimID);
    tID = netcdf.defVar(ncid,'T','double',t_dimID);

    % define output variables 
    KmID = netcdf.defVar(ncid,'Km','NC_FLOAT',data_IDs);
    
    Km_btID = netcdf.defVar(ncid,'Km_bt','NC_FLOAT',data_IDs);
    
    Km_bcID = netcdf.defVar(ncid,'Km_bc','NC_FLOAT',data_IDs);

    PmID = netcdf.defVar(ncid,'Pm','NC_FLOAT',data_IDs);

    PmKmID = netcdf.defVar(ncid,'PmKm','NC_FLOAT',data_IDs);

    if contains(str1, relative)

        WmID = netcdf.defVar(ncid,'Wm','NC_FLOAT',[x_dimID y_dimID t_dimID]);

    end

end

DaysID = netcdf.defVar(ncid,'Day','NC_FLOAT',t_dimID);

% leave define mode and enter data mode to write data.
netcdf.endDef(ncid);

% close the file
netcdf.close(ncid);

% end

%% Main part of script

tic

l_iter = 1;

for l = 1:length(n)
    
    dayID = id + n(l_iter) - 1;
    valID = dayID - pt_mean-1;
        
    % import data
    eta = ncread(sprintf(fn1,str1,hx,A),'ETAN',[1 1 1 valID],...
        [Inf Inf 1 2*pt_mean-1],[1 1 1 1]);
    T = ncread(sprintf(fn2,str1,hx,A),'THETA',[1 1 1 valID],...
        [Inf Inf nlvs 2*pt_mean-1],[1 1 1 1]);
    w = ncread(sprintf(fn3,str1,hx,A),'WVEL',[1 1 1 valID],...
        [Inf Inf nlvs 2*pt_mean-1],[1 1 1 1]);
    
    if contains(str1, relative)

        % Import TAU
        taux = ncread(sprintf(fn4,str1,hx,A),'EXFtaux',[1 1 1 valID],...
            [Inf Inf 1 2*pt_mean-1],[1 1 1 1]);
        tauy = ncread(sprintf(fn4,str1,hx,A),'EXFtauy',[1 1 1 valID],...
            [Inf Inf 1 2*pt_mean-1],[1 1 1 1]);
            
    end
    
    % find density at each daily output
    rho = zeros(nx,ny,size(T,3),size(T,4));
    rho_dev = zeros(nx,ny,size(T,3),size(T,4));
    for i = 1:size(T,4)
        for k = 1:size(T,3)
            % density
            rho(:,:,k,i) = rho_1*(1-2e-4*(T(:,:,k,i)-T0));
            % density deviation
            rho_dev(:,:,k,i) = rho(:,:,k,i) - rho_z_b(k);
        end
    end

    % geostrophic velocities
    [u, v] = geostrophic_uv(T,eta,f,hx,hy,z_levels);

    clear T eta
    
        
    u_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    v_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    w_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    rho_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    rho_dev_mid_mean = zeros(nx,ny,nlvs,pt_mean);


    for s = 1:pt_mean
        u_mid_mean(:,:,:,s) = mean(u(:,:,:,s:15+s),4);
        v_mid_mean(:,:,:,s) = mean(v(:,:,:,s:15+s),4);
        w_mid_mean(:,:,:,s) = mean(w(:,:,:,s:15+s),4);
        rho_mid_mean(:,:,:,s) = mean(rho(:,:,:,s:15+s),4);
        rho_dev_mid_mean(:,:,:,s) = mean(rho_dev(:,:,:,s:15+s),4);
    end

    u_mean = mean(u_mid_mean,4);
    v_mean = mean(v_mid_mean,4);
    w_mean = mean(w_mid_mean,4);
    rho_mean = mean(rho_mid_mean,4);
    rho_dev_mean = mean(rho_dev_mid_mean,4);

    % put vel means and fluctuations at centre of grid cell
    val_u_mean = zeros(nx,ny);
    val_v_mean = zeros(nx,ny);
    for k = 1:size(u,3)
        for i = 1:nx-1
            for j = 1:ny-1
                val_u_mean(i,j) = 0.5*(u_mean(i,j,k)+u_mean(i+1,j,k));
                val_v_mean(i,j) = 0.5*(v_mean(i,j,k)+v_mean(i,j+1,k));
            end
        end
        u_mean(:,:,k) = val_u_mean;
        v_mean(:,:,k) = val_v_mean;
    end
    
    % Find barotropic and baroclinic velocities
    val1 = zeros(nx,ny,nlvs); val2 = zeros(nx,ny,nlvs);
    for k = 1:size(u,3)
        val1(:,:,k) = u_mean(:,:,k).*delta_z(k);
        val2(:,:,k) = v_mean(:,:,k).*delta_z(k);
    end
    u_da = sum(val1,3)/H;
    v_da = sum(val2,3)/H;
    
    u_bt = zeros(nx,ny,nlvs); v_bt = zeros(nx,ny,nlvs);
    for k = 1:size(u,3)
        u_bt(:,:,k) = u_da;
        v_bt(:,:,k) = v_da;
    end
    
    u_bc = zeros(nx,ny,nlvs); v_bc = zeros(nx,ny,nlvs);
    for k = 1:size(u,3)
        u_bc(:,:,k) = u_mean(:,:,k) - u_bt(:,:,k);
        v_bc(:,:,k) = v_mean(:,:,k) - v_bt(:,:,k);
    end


    % -------------- calculate mean terms ---------------
    Km = zeros(nx,ny,size(u,3));
    tot_Km = zeros(1); sum_Km = zeros(1,size(u,3));
    Km_bt = zeros(nx,ny,size(u,3));
    tot_Km_bt = zeros(1); sum_Km_bt = zeros(1,size(u,3));
    Km_bc = zeros(nx,ny,size(u,3));
    tot_Km_bc = zeros(1); sum_Km_bc = zeros(1,size(u,3));
    Pm = zeros(nx,ny,size(u,3));
    tot_Pm = zeros(1); sum_Pm = zeros(1,size(u,3));
    PmKm = zeros(nx,ny,size(u,3)); sum_PmKm = zeros(1,size(u,3));
    tot_PmKm = zeros(1);
    for k = 1:size(u,3)
        % Km
        Km(:,:,k) = 0.5*(u_mean(:,:,k).^2+v_mean(:,:,k).^2);
        sum_Km(1,k) = sum(sum(Km(:,:,k)))*hx*hy*1e+6;
        % Km_bt
        Km_bt(:,:,k) = 0.5*(u_bt(:,:,k).^2+v_bt(:,:,k).^2);
        sum_Km_bt(1,k) = sum(sum(Km_bt(:,:,k)))*hx*hy*1e+6;
        % Km
        Km_bc(:,:,k) = 0.5*(u_bc(:,:,k).^2+v_bc(:,:,k).^2);
        sum_Km_bc(1,k) = sum(sum(Km_bc(:,:,k)))*hx*hy*1e+6;
        % Pm
        Pm(:,:,k) = -(g/(2*n0(k)))*(rho_dev_mean(:,:,k).^2);
        sum_Pm(1,k) = sum(sum(Pm(:,:,k)))*hx*hy*1e+6;
        % C(Pm,Km)
        PmKm(:,:,k) = -g*rho_mean(:,:,k).*w_mean(:,:,k);
        sum_PmKm(1,k) = sum(sum(PmKm(:,:,k)))*hx*hy*1e+6;
    end
    tot_Km(1) = rho_1*sum(sum_Km(1,:).*delta_z);
    tot_Km_bt(1) = rho_1*sum(sum_Km_bt(1,:).*delta_z);
    tot_Km_bc(1) = rho_1*sum(sum_Km_bc(1,:).*delta_z);
    tot_Pm(1) = sum(sum_Pm(1,:).*delta_z);
    tot_PmKm(1) = sum(sum_PmKm(1,:).*delta_z);

    if contains(str1, relative)

        taux_mid_mean = zeros(nx,ny,1,pt_mean);
        tauy_mid_mean = zeros(nx,ny,1,pt_mean);

        for s = 1:pt_mean
            taux_mid_mean(:,:,:,s) = mean(taux(:,:,:,s:15+s),4);
            tauy_mid_mean(:,:,:,s) = mean(tauy(:,:,:,s:15+s),4);
        end

        taux_mean = mean(taux_mid_mean,4);
        tauy_mean = mean(tauy_mid_mean,4);

        % work done by mean wind stress
        Wm = taux_mean.*u_mean(:,:,1)+tauy_mean.*v_mean(:,:,1);
        tot_Wm = sum(sum(Wm))*hx*hy*1e+6;

    end
            
        
        
    
    % write variables to netcdf files
    ncid = netcdf.open(sprintf(filename,str1,str2,hx,A),'NC_WRITE');
    
    if contains(str2, total)
    
        % write data to netcdf
        netcdf.putVar(ncid,sum_KmID,[0 l-1],[nlvs 1],[1 1],sum_Km);
        netcdf.putVar(ncid,tot_KmID,l-1,1,1,tot_Km);
        
        netcdf.putVar(ncid,sum_Km_btID,[0 l-1],[nlvs 1],[1 1],sum_Km_bt);
        netcdf.putVar(ncid,tot_Km_btID,l-1,1,1,tot_Km_bt);
        
        netcdf.putVar(ncid,sum_Km_bcID,[0 l-1],[nlvs 1],[1 1],sum_Km_bc);
        netcdf.putVar(ncid,tot_Km_bcID,l-1,1,1,tot_Km_bc);

        netcdf.putVar(ncid,sum_PmID,[0 l-1],[nlvs 1],[1 1],sum_Pm);
        netcdf.putVar(ncid,tot_PmID,l-1,1,1,tot_Pm);

        netcdf.putVar(ncid,sum_PmKmID,[0 l-1],[nlvs 1],[1 1],sum_PmKm);
        netcdf.putVar(ncid,tot_PmKmID,l-1,1,1,tot_PmKm);

        if contains(str1, relative)

            netcdf.putVar(ncid,tot_WmID,l-1,1,1,tot_Wm);

        end
        
    elseif contains(str2, plan)
        
        % write data to netcdf
        netcdf.putVar(ncid,KmID,[0 0 0 l-1],[nx ny nlvs 1],[1 1 1 1],Km); 
        
        netcdf.putVar(ncid,Km_btID,[0 0 0 l-1],[nx ny nlvs 1],[1 1 1 1],Km_bt); 
        
        netcdf.putVar(ncid,Km_bcID,[0 0 0 l-1],[nx ny nlvs 1],[1 1 1 1],Km_bc); 

        netcdf.putVar(ncid,PmID,[0 0 0 l-1],[nx ny nlvs 1],[1 1 1 1],Pm);  

        netcdf.putVar(ncid,PmKmID,[0 0 0 l-1],[nx ny nlvs 1],[1 1 1 1],PmKm); 

        if contains(str1, relative)

            netcdf.putVar(ncid,WmID,[0 0 l-1],[nx ny 1],[1 1 1],Wm);

        end
    
    end
    
    % write in timestamp
    netcdf.putVar(ncid,DaysID,l-1,1,1,n(l_iter));
    
    % close the file
    netcdf.close(ncid);
    
    fprintf('DAY %d\n',l);
    
    l_iter = l_iter + 1;
    
    
    
    
end
