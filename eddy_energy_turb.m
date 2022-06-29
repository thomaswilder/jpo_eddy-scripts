clear all;

% Calculates turbulent energy terms for the eddy following von Storch et al (2012)
% using geostrophic horizontal velocity components.
% Calculations are either 4D or volume integrals, depending on choice.
% Can be either used locally or in a hpc environment.
% Run time will be extended for resolution greater than 10 km.
% Uses a 16 day time-mean to avoid picking up any signals in the
% instability pathways from the rotating wind

% Variables calculated:
% Kt - turbulent kinetic energy
% Pt - turbulent potential energy
% PtPm - turbulent potential to mean potential energy
% PtKt - baroclinic pathway
% KtKm - barotropic pathway. _h and _v are horizontal and vertical components.

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
    n = [31:1:380];                  
    nlvs = nz; % number of levels to analyse
elseif contains(str2,plan)
    n = [35 75 100 125 150 175 200 225]; % days for 4D terms
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



%% Setup netcdf file

% open netCDF file.
filename = 'data/ACE_energy_%s_%s_%dkm_A%d_beta.nc';

% enables pickup of data file from last time step, used when node time limit likely to be reached on a hpc system.
% check if file exists
if isfile(sprintf(filename,str1,str2,hx,A))
    
    % get dimensions of T
    ncid = netcdf.open(sprintf(filename,str1,str2,hx,A),'NC_WRITE');
    varid = netcdf.inqVarID(ncid,'T');
    [dimname, time] = netcdf.inqDim(ncid,varid);
    netcdf.close(ncid);
    
    % then, proceed to the calculations 
    
else
    
    time = 0;

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
        sum_KtID = netcdf.defVar(ncid,'sum_Kt','NC_FLOAT',[z_dimID t_dimID]);
        tot_KtID = netcdf.defVar(ncid,'tot_Kt','NC_FLOAT',t_dimID);

        sum_PtID = netcdf.defVar(ncid,'sum_Pt','NC_FLOAT',[z_dimID t_dimID]);
        tot_PtID = netcdf.defVar(ncid,'tot_Pt','NC_FLOAT',t_dimID);

        sum_PtPmID = netcdf.defVar(ncid,'sum_PtPm','NC_FLOAT',[z_dimID t_dimID]);
        tot_PtPmID = netcdf.defVar(ncid,'tot_PtPm','NC_FLOAT',t_dimID);

        sum_PtKtID = netcdf.defVar(ncid,'sum_PtKt','NC_FLOAT',[z_dimID t_dimID]);
        tot_PtKtID = netcdf.defVar(ncid,'tot_PtKt','NC_FLOAT',t_dimID);

        sum_KtKmID = netcdf.defVar(ncid,'sum_KtKm','NC_FLOAT',[z_dimID t_dimID]);
        tot_KtKmID = netcdf.defVar(ncid,'tot_KtKm','NC_FLOAT',t_dimID);

        sum_KtKm_hID = netcdf.defVar(ncid,'sum_KtKm_h','NC_FLOAT',[z_dimID t_dimID]);
        tot_KtKm_hID = netcdf.defVar(ncid,'tot_KtKm_h','NC_FLOAT',t_dimID);

        sum_KtKm_vID = netcdf.defVar(ncid,'sum_KtKm_v','NC_FLOAT',[z_dimID t_dimID]);
        tot_KtKm_vID = netcdf.defVar(ncid,'tot_KtKm_v','NC_FLOAT',t_dimID);

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

        KtID = netcdf.defVar(ncid,'Kt','NC_FLOAT',data_IDs);

        PtID = netcdf.defVar(ncid,'Pt','NC_FLOAT',data_IDs);

        PtPmID = netcdf.defVar(ncid,'PtPm','NC_FLOAT',data_IDs);

        PtKtID = netcdf.defVar(ncid,'PtKt','NC_FLOAT',data_IDs);

        KtKmID = netcdf.defVar(ncid,'KtKm','NC_FLOAT',data_IDs);

        KtKm_hID = netcdf.defVar(ncid,'KtKm_h','NC_FLOAT',data_IDs);

        KtKm_vID = netcdf.defVar(ncid,'KtKm_v','NC_FLOAT',data_IDs);

    end

    DaysID = netcdf.defVar(ncid,'Day','NC_FLOAT',t_dimID);

    % leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);

    % close the file
    netcdf.close(ncid);

end

%% Main part of script

tic

l_iter = 1;

for l = 1:length(n)
    
    time = time + 1;
    
    dayID = id + n(l_iter) - 1;
    valID = dayID - (2*pt_mean-2);
        
    % import data
    eta = ncread(sprintf(fn1,str1,hx,A),'ETAN',[1 1 1 valID],...
        [Inf Inf 1 4*pt_mean-1],[1 1 1 1]);
    T = ncread(sprintf(fn2,str1,hx,A),'THETA',[1 1 1 valID],...
        [Inf Inf nlvs 4*pt_mean-1],[1 1 1 1]);
    w = ncread(sprintf(fn3,str1,hx,A),'WVEL',[1 1 1 valID],...
        [Inf Inf nlvs 4*pt_mean-1],[1 1 1 1]);   
    
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
    
    u_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    v_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    w_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    rho_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    rho_dev_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    
    % initialise fluctuation products
    Kt_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    urho_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    vrho_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    wrho_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    uu_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    uv_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    uw_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    vv_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    vw_anom = zeros(nx,ny,nlvs,2*pt_mean+1);
    rho_dev_anom_sq = zeros(nx,ny,nlvs,2*pt_mean+1);
    
    for r = 1:2*pt_mean+1
        
        u_mid_mean = zeros(nx,ny,nlvs,pt_mean);
        v_mid_mean = zeros(nx,ny,nlvs,pt_mean);
        w_mid_mean = zeros(nx,ny,nlvs,pt_mean);
        rho_mid_mean = zeros(nx,ny,nlvs,pt_mean);
        rho_dev_mid_mean = zeros(nx,ny,nlvs,pt_mean);
        
        
        for s = 1:pt_mean
            u_mid_mean(:,:,:,s) = mean(u(:,:,:,s+r-1:15+s+r-1),4);
            v_mid_mean(:,:,:,s) = mean(v(:,:,:,s+r-1:15+s+r-1),4);
            w_mid_mean(:,:,:,s) = mean(w(:,:,:,s+r-1:15+s+r-1),4);
            rho_mid_mean(:,:,:,s) = mean(rho(:,:,:,s+r-1:15+s+r-1),4);
            rho_dev_mid_mean(:,:,:,s) = mean(rho_dev(:,:,:,s+r-1:15+s+r-1),4);
        end
        
        u_mean = mean(u_mid_mean,4);
        v_mean = mean(v_mid_mean,4);
        w_mean = mean(w_mid_mean,4);
        rho_mean = mean(rho_mid_mean,4);
        rho_dev_mean = mean(rho_dev_mid_mean,4);
        
        u_anom(:,:,:,r) = u(:,:,:,pt_mean-1+r) - u_mean;
        v_anom(:,:,:,r) = v(:,:,:,pt_mean-1+r) - v_mean;
        w_anom(:,:,:,r) = w(:,:,:,pt_mean-1+r) - w_mean;
        rho_anom(:,:,:,r) = rho(:,:,:,pt_mean-1+r) - rho_mean;
        rho_dev_anom(:,:,:,r) = rho_dev(:,:,:,pt_mean-1+r) - rho_dev_mean;


        % put vel means and fluctuations at centre of grid cell
        val_u_anom = zeros(nx,ny);
        val_v_anom = zeros(nx,ny);
        val_u_mean = zeros(nx,ny);
        val_v_mean = zeros(nx,ny);
        for k = 1:size(u,3)
            for i = 1:nx-1
                for j = 1:ny-1
                    val_u_anom(i,j) = 0.5*(u_anom(i,j,k,r)+u_anom(i+1,j,k,r));
                    val_v_anom(i,j) = 0.5*(v_anom(i,j,k,r)+v_anom(i,j+1,k,r));
                    val_u_mean(i,j) = 0.5*(u_mean(i,j,k)+u_mean(i+1,j,k));
                    val_v_mean(i,j) = 0.5*(v_mean(i,j,k)+v_mean(i,j+1,k));
                end
            end
            u_anom(:,:,k,r) = val_u_anom;
            v_anom(:,:,k,r) = val_v_anom;
            u_mean(:,:,k) = val_u_mean;
            v_mean(:,:,k) = val_v_mean;
        end
        
        if r == pt_mean

            drmdx = zeros(nx,ny,size(u,3)); drmdy = zeros(nx,ny,size(u,3));
            dudx = zeros(nx,ny,size(u,3)); dudy = zeros(nx,ny,size(u,3));
            dvdx = zeros(nx,ny,size(u,3)); dvdy = zeros(nx,ny,size(u,3));
            for k = 1:size(u,3)
                % density
                drmdx(:,:,k) = dvald(rho_mean(:,:,k), hx, hy, 0, 'x');    
                drmdy(:,:,k) = dvald(rho_mean(:,:,k), hx, hy, 0, 'y');
                % horizontal gradients
                dudx(:,:,k) = dvald(u_mean(:,:,k), hx, hy, 0, 'x');
                dudy(:,:,k) = dvald(u_mean(:,:,k), hx, hy, 0, 'y');
                dvdx(:,:,k) = dvald(v_mean(:,:,k), hx, hy, 0, 'x');
                dvdy(:,:,k) = dvald(v_mean(:,:,k), hx, hy, 0, 'y');
            end
            % vertical gradients
            dudz = dvald(u_mean, hx, hy, 0, 'z', delta_z);
            dvdz = dvald(v_mean, hx, hy, 0, 'z', delta_z);

        end
        
        % find fluctuation products
        for k = 1:nlvs
            Kt_anom(:,:,k,r) = 0.5*(u_anom(:,:,k,r).^2+v_anom(:,:,k,r).^2);
            urho_anom(:,:,k,r) = rho_anom(:,:,k,r).*u_anom(:,:,k,r);
            vrho_anom(:,:,k,r) = rho_anom(:,:,k,r).*v_anom(:,:,k,r);
            wrho_anom(:,:,k,r) = rho_anom(:,:,k,r).*w_anom(:,:,k,r);
            uu_anom(:,:,k,r) = u_anom(:,:,k,r).^2;
            uv_anom(:,:,k,r) = u_anom(:,:,k,r).*v_anom(:,:,k,r);
            uw_anom(:,:,k,r) = u_anom(:,:,k,r).*w_anom(:,:,k,r);
            vv_anom(:,:,k,r) = v_anom(:,:,k,r).^2;
            vw_anom(:,:,k,r) = v_anom(:,:,k,r).*w_anom(:,:,k,r);
            rho_dev_anom_sq(:,:,k,r) = rho_dev_anom(:,:,k,r).^2;
        end
        
    
    end
    
    % calculate fluctuation terms
    Kt_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    urho_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    vrho_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    wrho_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    uu_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    uv_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    uw_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    vv_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    vw_anom_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    rho_dev_anom_sq_mid_mean = zeros(nx,ny,nlvs,pt_mean);
    
    for s = 1:pt_mean
        Kt_mid_mean(:,:,:,s) = mean(Kt_anom(:,:,:,s:15+s),4);
        urho_anom_mid_mean(:,:,:,s) = mean(urho_anom(:,:,:,s:15+s),4);
        vrho_anom_mid_mean(:,:,:,s) = mean(vrho_anom(:,:,:,s:15+s),4);
        wrho_anom_mid_mean(:,:,:,s) = mean(wrho_anom(:,:,:,s:15+s),4);
        uu_anom_mid_mean(:,:,:,s) = mean(uu_anom(:,:,:,s:15+s),4);
        uv_anom_mid_mean(:,:,:,s) = mean(uv_anom(:,:,:,s:15+s),4);
        uw_anom_mid_mean(:,:,:,s) = mean(uw_anom(:,:,:,s:15+s),4);
        vv_anom_mid_mean(:,:,:,s) = mean(vv_anom(:,:,:,s:15+s),4);
        vw_anom_mid_mean(:,:,:,s) = mean(vw_anom(:,:,:,s:15+s),4);
        rho_dev_anom_sq_mid_mean(:,:,:,s) = mean(rho_dev_anom_sq(:,:,:,s:15+s),4);
    end
    
    Kt = mean(Kt_mid_mean,4);
    urho_anom_mean = mean(urho_anom_mid_mean,4);
    vrho_anom_mean = mean(vrho_anom_mid_mean,4);
    wrho_anom_mean = mean(wrho_anom_mid_mean,4);
    uu_anom_mean = mean(uu_anom_mid_mean,4);
    uv_anom_mean = mean(uv_anom_mid_mean,4);
    uw_anom_mean = mean(uw_anom_mid_mean,4);
    vv_anom_mean = mean(vv_anom_mid_mean,4);
    vw_anom_mean = mean(vw_anom_mid_mean,4);
    rho_dev_anom_sq_mean = mean(rho_dev_anom_sq_mid_mean,4);
    
    % find conversion terms and sum over whole domain
    PtPm = zeros(nx,ny,size(u,3)); 
    sum_PtPm = zeros(1,size(u,3)); tot_PtPm = zeros(1);
    PtKt = zeros(nx,ny,size(u,3)); 
    sum_PtKt = zeros(1,size(u,3)); tot_PtKt = zeros(1);
    KtKm_h = zeros(nx,ny,size(u,3)); KtKm_v = zeros(nx,ny,size(u,3));
    KtKm = zeros(nx,ny,size(u,3));
    sum_KtKm = zeros(1,size(u,3)); tot_KtKm = zeros(1);
    sum_KtKm_h = zeros(1,size(u,3)); tot_KtKm_h = zeros(1);
    sum_KtKm_v = zeros(1,size(u,3)); tot_KtKm_v = zeros(1);
    tot_Kt = zeros(1); sum_Kt = zeros(1,size(u,3));
    Pt = zeros(nx,ny,size(u,3));
    tot_Pt = zeros(1); sum_Pt = zeros(1,size(u,3));
    for k = 1:size(u,3)
        % C(Pt,Pm)
        PtPm(:,:,k) = -(g/n0(k))*(urho_anom_mean(:,:,k).*drmdx(:,:,k)+...
            vrho_anom_mean(:,:,k).*drmdy(:,:,k));
        sum_PtPm(1,k) = sum(sum(PtPm(:,:,k)))*hx*hy*1e+6;
        % C(Pt,Kt)
        PtKt(:,:,k) = -g*wrho_anom_mean(:,:,k);
        sum_PtKt(1,k) = sum(sum(PtKt(:,:,k)))*hx*hy*1e+6;
        % C(Kt,Km)
        % decomposing into horizontal and vertical Reynolds shear
        KtKm_h(:,:,k) = rho_1*( uu_anom_mean(:,:,k).*dudx(:,:,k)...
            +uv_anom_mean(:,:,k).*dudy(:,:,k)+...
            uv_anom_mean(:,:,k).*dvdx(:,:,k)+...
            vv_anom_mean(:,:,k).*dvdy(:,:,k));
        KtKm_v(:,:,k) = rho_1*(uw_anom_mean(:,:,k).*dudz(:,:,k)+...
            vw_anom_mean(:,:,k).*dvdz(:,:,k));
        KtKm(:,:,k) = KtKm_h(:,:,k) + KtKm_v(:,:,k);
        sum_KtKm_h(1,k) = sum(sum(KtKm_h(:,:,k)))*hx*hy*1e+6;
        sum_KtKm_v(1,k) = sum(sum(KtKm_v(:,:,k)))*hx*hy*1e+6;
        sum_KtKm(1,k) = sum_KtKm_h(k) + sum_KtKm_v(k);
        % Kt
        sum_Kt(1,k) = sum(sum(Kt(:,:,k)))*hx*hy*1e+6;
        % Pt
        Pt(:,:,k) = -(g/(2*n0(k)))*rho_dev_anom_sq_mean(:,:,k);
        sum_Pt(1,k) = sum(sum(Pt(:,:,k)))*hx*hy*1e+6;
    end
    tot_PtPm(1) = sum(sum_PtPm(1,:).*delta_z);
    tot_PtKt(1) = sum(sum_PtKt(1,:).*delta_z);
    tot_KtKm_h(1) = sum(sum_KtKm_h(1,:).*delta_z);
    tot_KtKm_v(1) = sum(sum_KtKm_v(1,:).*delta_z);
    tot_KtKm(1) = tot_KtKm_h(1) + tot_KtKm_v(1);
    tot_Kt(1) = rho_1*sum(sum_Kt(1,:).*delta_z);
    tot_Pt(1) = sum(sum_Pt(1,:).*delta_z);
    

    
    
    % write variables to netcdf files
    ncid = netcdf.open(sprintf(filename,str1,str2,hx,A),'NC_WRITE');
    
    if contains(str2, total)
    
        % write data to netcdf
        sum_KtID = netcdf.inqVarID(ncid,'sum_Kt');
        netcdf.putVar(ncid,sum_KtID,[0 time-1],[nlvs 1],[1 1],sum_Kt);
        
        tot_KtID = netcdf.inqVarID(ncid,'tot_Kt');
        netcdf.putVar(ncid,tot_KtID,time-1,1,1,tot_Kt);

        sum_PtID = netcdf.inqVarID(ncid,'sum_Pt');
        netcdf.putVar(ncid,sum_PtID,[0 time-1],[nlvs 1],[1 1],sum_Pt);
        
        tot_PtID = netcdf.inqVarID(ncid,'tot_Pt');
        netcdf.putVar(ncid,tot_PtID,time-1,1,1,tot_Pt);

        sum_PtKtID = netcdf.inqVarID(ncid,'sum_PtKt');
        netcdf.putVar(ncid,sum_PtKtID,[0 time-1],[nlvs 1],[1 1],sum_PtKt);
        
        tot_PtKtID = netcdf.inqVarID(ncid,'tot_PtKt');
        netcdf.putVar(ncid,tot_PtKtID,time-1,1,1,tot_PtKt);

        sum_PtPmID = netcdf.inqVarID(ncid,'sum_PtPm');
        netcdf.putVar(ncid,sum_PtPmID,[0 time-1],[nlvs 1],[1 1],sum_PtPm);
        
        tot_PtPmID = netcdf.inqVarID(ncid,'tot_PtPm');
        netcdf.putVar(ncid,tot_PtPmID,time-1,1,1,tot_PtPm);

        sum_KtKmID = netcdf.inqVarID(ncid,'sum_KtKm');
        netcdf.putVar(ncid,sum_KtKmID,[0 time-1],[nlvs 1],[1 1],sum_KtKm);
        
        tot_KtKmID = netcdf.inqVarID(ncid,'tot_KtKm');
        netcdf.putVar(ncid,tot_KtKmID,time-1,1,1,tot_KtKm);

        sum_KtKm_hID = netcdf.inqVarID(ncid,'sum_KtKm_h');
        netcdf.putVar(ncid,sum_KtKm_hID,[0 time-1],[nlvs 1],[1 1],sum_KtKm_h);
        
        tot_KtKm_hID = netcdf.inqVarID(ncid,'tot_KtKm_h');
        netcdf.putVar(ncid,tot_KtKm_hID,time-1,1,1,tot_KtKm_h);

        sum_KtKm_vID = netcdf.inqVarID(ncid,'sum_KtKm_v');
        netcdf.putVar(ncid,sum_KtKm_vID,[0 time-1],[nlvs 1],[1 1],sum_KtKm_v);
        
        tot_KtKm_vID = netcdf.inqVarID(ncid,'tot_KtKm_v');
        netcdf.putVar(ncid,tot_KtKm_vID,time-1,1,1,tot_KtKm_v);
        
    elseif contains(str2, plan)
        
        % write data to netcdf 
        KtID = netcdf.inqVarID(ncid,'Kt');
        netcdf.putVar(ncid,KtID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],Kt); 

        PtID = netcdf.inqVarID(ncid,'Pt');
        netcdf.putVar(ncid,PtID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],Pt); 

        PtKtID = netcdf.inqVarID(ncid,'PtKt');
        netcdf.putVar(ncid,PtKtID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],PtKt); 

        PtPmID = netcdf.inqVarID(ncid,'PtPm');
        netcdf.putVar(ncid,PtPmID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],PtPm);

        KtKmID = netcdf.inqVarID(ncid,'KtKm');
        netcdf.putVar(ncid,KtKmID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],KtKm); 

        KtKm_hID = netcdf.inqVarID(ncid,'KtKm_h');
        netcdf.putVar(ncid,KtKm_hID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],KtKm_h); 

        KtKm_vID = netcdf.inqVarID(ncid,'KtKm_v');
        netcdf.putVar(ncid,KtKm_vID,[0 0 0 time-1],[nx ny nlvs 1],[1 1 1 1],KtKm_v);
    
    end
    
    % write in timestamp
    DaysID = netcdf.inqVarID(ncid,'Day');
    netcdf.putVar(ncid,DaysID,time-1,1,1,n(l_iter));
    
    % close the file
    netcdf.close(ncid);
    
    fprintf('DAY %d\n',l);
    
    l_iter = l_iter + 1;
    

end

toc
