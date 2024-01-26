function [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,T)
% Metropolis for Potts model(2d case)
% 
beta = 1/kb/T;
N = size(lat,1); 
i = randi(N);
j = randi(N);
k = randi(N);
sign = 1; 

pro = randi(q); newval = pro;
mid = lat(i,j,k); 
left = lat(mod(i-2,N)+1,j,k);
right = lat(mod(i,N)+1,j,k);
bot = lat(i,mod(j-2,N)+1,k);
top = lat(i,mod(j,N)+1,k);
back = lat(i,j,mod(k-2,N)+1);
front = lat(i,j,mod(k,N)+1);
deltaH = -J * (double(pro==left) + double(pro==right) + double(pro==top) + double(pro==bot) + double(pro==back) + double(pro==front)) + J * (double(mid==left) + double(mid==right) + double(mid==top) + double(mid==bot) + double(mid==back) + double(mid==front)) - h * (pro - mid);
if rand() > min(1,exp(-beta * deltaH) )
    newval = lat(i,j,k);
    sign = 0; deltaH = 0;
end
position = [i,j,k];
end

