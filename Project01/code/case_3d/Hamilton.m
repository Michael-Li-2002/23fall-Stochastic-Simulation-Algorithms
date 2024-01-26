function H = Hamilton(lat,J,h)
% Hamiltonian
% 
H = 0; 
N = size(lat,1);

% for i = 1:N
%     for j = 1:N
%         mid = lat(i,j);
%         left = lat(mod(i-2,N)+1,j);
%         right = lat(mod(i,N)+1,j);
%         top = lat(i,mod(j-2,N)+1);
%         bot = lat(i,mod(j,N)+1);
%         H = H - 0.5 * J * (isSame(mid,left) + isSame(mid,right) + isSame(mid,top) + isSame(mid,bot)) - h * mid;
%     end
% end

xind = (1:N)'; yind = (1:N)'; zind = (1:N)';
right = lat(mod(xind,N)+1,yind,zind);
top = lat(xind,mod(yind,N)+1,zind);
front = lat(xind,yind,mod(zind,N)+1);
H = -J * ( sum(sum(sum(double(lat==right)))) + sum(sum(sum(double(lat==top)))) + sum(sum(sum(double(lat==front)))) ) - h*sum(sum(sum(lat)));
end

