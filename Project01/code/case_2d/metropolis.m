function [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,T)
% Metropolis for Potts model(2d case)
% 
beta = 1/kb/T;
N = size(lat,1); 
i = randi(N);
j = randi(N);
sign = 1; 

% pro = randi(q); 
pro = randi(q-1);
if pro >= lat(i,j)
    pro = pro + 1;
end
newval = pro;
mid = lat(i,j); 
left = lat(mod(i-2,N)+1,j);
right = lat(mod(i,N)+1,j);
bot = lat(i,mod(j-2,N)+1);
top = lat(i,mod(j,N)+1);
deltaH = -J * (double(pro==left) + double(pro==right) + double(pro==top) + double(pro==bot)) + J * (double(mid==left) + double(mid==right) + double(mid==top) + double(mid==bot)) - h * (pro - mid);
if rand() > min(1,exp(-beta * deltaH) )
    newval = lat(i,j);
    sign = 0; deltaH = 0;
end
position = [i,j];
end

