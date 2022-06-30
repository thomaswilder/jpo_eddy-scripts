function [dval] = dvald(val, hx, hy, ht, S, delta_z)
% Finds the first derivative of 'val' in either x and y.
% Centred difference for x and y.

nx = size(val,1); ny = size(val,2);

if nargin==5
    
dval = zeros(nx,ny); 

if strcmp(S,'x') % x derivative
    for i = 2:nx-1
        for j = 1:ny
            dval(i,j) = (val(i+1,j)-val(i-1,j))/(2*hx*1e+3);
        end
    end
    
    for j = 1:ny
        dval(1,j) = (val(2,j)-val(nx,j))/(2*hx*1e+3);
        dval(nx,j) = (val(1,j)-val(nx-1,j))/(2*hx*1e+3);
    end
    
    
elseif strcmp(S,'y') % y derivative
    for i = 1:nx
        for j = 2:ny-1
            dval(i,j) = (val(i,j+1)-val(i,j-1))/(2*hy*1e+3);
        end
    end
    
    for i = 1:nx
        dval(i,1) = (val(i,2)-val(i,ny))/(2*hy*1e+3);
        dval(i,ny) = (val(i,1)-val(i,ny-1))/(2*hy*1e+3);
    end
    
elseif strcmp(S, 't')
    for i = 1:nx
        for j = 1:ny
            dval(i,j) = (val(i,j,2)-val(i,j,1))/(ht);
        end
    end
    
end

else

if strcmp(S, 'z') % z derivative
    
    nz = size(val,3);
    dval = zeros(nx,ny,nz);
    
    for k = 1:nz-1
       dval(:,:,k) = (val(:,:,k)-val(:,:,k+1))/delta_z(k);
    end
    dval(:,:,nz) = val(:,:,nz)/delta_z(nz);
    
end
end


end

