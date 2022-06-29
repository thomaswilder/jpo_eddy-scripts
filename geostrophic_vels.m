function [u, v] = geostrophic_vels(T,eta,f,alpha,hx,hy,z_levels)

g = 9.81;

nx = size(T,1); ny = size(T,2); nz = size(T,3); nt = size(T,4);

u = zeros(nx,ny,nz,size(T,4)); 
v = zeros(nx,ny,nz,size(T,4));

for i = 1:length(z_levels)-1
    delta_z(i) = z_levels(i+1) - z_levels(i);
end


for l = 1:nt

    % surface geostrophic flow
    for i = 1:nx
        for j = 2:ny
            u(i,j,1,l) = -(g/f)*((eta(i,j)-...
                eta(i,j-1))/(hy*1e+3));
        end
        u(i,1,1,l) = -(g/f)*((eta(i,1)-...
                eta(i,ny))/(hy*1e+3));
    end
    
    for j = 1:ny
        for i = 2:nx
            v(i,j,1,l) = (g/f)*((eta(i,j)-...
                eta(i-1,j))/(hx*1e+3));
        end
        v(1,j,1,l) = (g/f)*((eta(1,j)-...
                eta(nx,j))/(hx*1e+3));
    end

    % temperature gradients
    dTdx = zeros(nx,ny,nz);
    dTdy = zeros(nx,ny,nz);
    for k = 1:nz
        for j = 1:ny
            for i = 2:nx
                dTdx(i,j,k) = ((T(i,j,k)-T(i-1,j,k))/(hx*1e+3));
            end
            dTdx(1,j,k) = ((T(1,j,k)-T(nx,j,k))/(hx*1e+3));
        end
        for i = 1:nx
            for j = 2:ny
                dTdy(i,j,k) = ((T(i,j,k)-T(i,j-1,k))/(hy*1e+3));
            end
            dTdy(i,1,k) = ((T(i,1,k)-T(i,ny,k))/(hy*1e+3));
        end
    end

    % thermal wind shear
    for k = 1:nz
        dudz(:,:,k) = -((g*alpha)/(f))*dTdy(:,:,k)*delta_z(k);
        dvdz(:,:,k) = ((g*alpha)/(f))*dTdx(:,:,k)*delta_z(k);
    end

    % Subtract thermal wind shear from surface velocities
    for k=2:nz
        u(:,:,k,l) = u(:,:,k-1,l)-dudz(:,:,k-1);
        v(:,:,k,l) = v(:,:,k-1,l)-dvdz(:,:,k-1);
    end

end


end

