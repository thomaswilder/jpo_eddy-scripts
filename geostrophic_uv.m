function [u, v] = geostrophic_uv(T,Eta,f,hx,hy,z_levels)

% Find velocities in thermal wind balance

g = 9.81;
alpha = 2e-4;

nx = size(Eta,1); ny = size(Eta,2); nz = size(T,3);

u = zeros(nx, ny, nz, size(T,4)); 
v = zeros(nx, ny, nz, size(T,4));

delta_z = zeros(1,nz);
for i = 1:length(z_levels)-1
    delta_z(i) = z_levels(i+1) - z_levels(i);
end

for l = 1:size(T,4)
    % find eta on left most corner of grid cell
    etaG = zeros(nx+1,ny+1);
    for i = 2:nx
        for j = 2:ny
            etaG(i,j) = 0.25*(Eta(i-1,j-1,1,l)+Eta(i-1,j,1,l)+...
                Eta(i,j-1,1,l)+Eta(i,j,1,l));
        end
        etaG(i,1) = 0.25*(Eta(i-1,ny,1,l)+Eta(i-1,1,1,l)+...
            Eta(i,ny,1,l)+Eta(i,1,1,l));
        etaG(i,ny+1) = 0.25*(Eta(i-1,ny,1,l)+Eta(i-1,1,1,l)+...
            Eta(i,ny,1,l)+Eta(i,1,1,l));
    end
    for j = 2:ny
        etaG(1,j) = 0.25*(Eta(nx,j-1,1,l)+Eta(nx,j,1,l)+...
            Eta(1,j-1,1,l)+Eta(1,j,1,l));
        etaG(nx+1,j) = 0.25*(Eta(nx,j-1,1,l)+Eta(nx,j,1,l)+...
            Eta(1,j-1,1,l)+Eta(1,j,1,l));
    end
    
    % surface velocities
    for i = 1:nx
        for j = 1:ny-1
            u(i,j,1,l) = -(g/f)*((etaG(i,j+1)-...
                etaG(i,j))/(hy*1e+3));
        end
        u(i,ny,1,l) = -(g/f)*((etaG(i,1)-...
                etaG(i,ny))/(hy*1e+3));
    end
    
    for j = 1:ny
        for i = 1:nx-1
            v(i,j,1,l) = (g/f)*((etaG(i+1,j)-...
                etaG(i,j))/(hx*1e+3));
        end
        v(nx,j,1,l) = (g/f)*((etaG(1,j)-...
                etaG(nx,j))/(hx*1e+3));
    end
    
    % find T on left most corner of grid cell
    TG = zeros(nx+1,ny+1,nz);
    for k = 1:nz-1
        for i = 2:nx
            for j = 2:ny
                TG(i,j,k) = 0.25*(T(i-1,j-1,k,l)+T(i-1,j,k,l)+...
                    T(i,j-1,k,l)+T(i,j,k,l));
            end
            TG(i,1,k) = 0.25*(T(i-1,ny,k,l)+T(i-1,1,k,l)+...
                T(i,ny,k,l)+T(i,1,k,l));
            TG(i,ny+1,k) = 0.25*(T(i-1,ny,k,l)+T(i-1,1,k,l)+...
                T(i,ny,k,l)+T(i,1,k,l));
        end
        for j = 2:ny
            TG(1,j,k) = 0.25*(T(nx,j-1,k,l)+T(nx,j,k,l)+...
                T(1,j-1,k,l)+T(1,j,k,l));
            TG(nx+1,j,k) = 0.25*(T(nx,j-1,k,l)+T(nx,j,k,l)+...
                T(1,j-1,k,l)+T(1,j,k,l));
        end

        TG(1,1,k) = 0.25*(T(nx,ny,k,l)+T(nx,1,k,l)+...
                    T(1,ny,k,l)+T(1,1,k,l));
        TG(1,ny+1,k) = 0.25*(T(nx,ny,k,l)+T(nx,1,k,l)+...
                    T(1,ny,k,l)+T(1,1,k,l));
        TG(nx+1,1,k) = 0.25*(T(nx,ny,k,l)+T(nx,1,k,l)+...
                    T(1,ny,k,l)+T(1,1,k,l));
        TG(nx+1,ny+1,k) = 0.25*(T(nx,ny,k,l)+T(nx,1,k,l)+...
                    T(1,ny,k,l)+T(1,1,k,l));
    end


    % temperature gradients
    dT_dx = zeros(nx,ny,nz-1);
    dT_dy = zeros(nx,ny,nz-1);
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx-1
                dT_dx(i,j,k) = ((TG(i+1,j,k)-TG(i,j,k))/(hx*1e+3));
            end
        end
        for i = 1:nx
            for j = 1:ny-1
                dT_dy(i,j,k) = ((TG(i,j+1,k)-TG(i,j,k))/(hy*1e+3));
            end
        end
    end

    % vertical shear
    dudz = zeros(nx,ny,nz);
    dvdz = zeros(nx,ny,nz);
    for k = 1:nz
        dudz(:,:,k) = -((g*alpha)/(f))*dT_dy(:,:,k)*delta_z(k);
        dvdz(:,:,k) = ((g*alpha)/(f))*dT_dx(:,:,k)*delta_z(k);
    end

    % total velocity components
    for k=2:nz-1
        u(:,:,k,l) = u(:,:,k-1,l)-dudz(:,:,k);
        v(:,:,k,l) = v(:,:,k-1,l)-dvdz(:,:,k);
    end

end


end

