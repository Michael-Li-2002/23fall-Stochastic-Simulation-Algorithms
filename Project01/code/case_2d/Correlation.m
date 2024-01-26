function [Cormat] = Correlation(lat)
%CORRELATION 

N = size(lat,1);
Cormat = zeros(N,N,N);
ind = 1:N;
for k = 1:N
    Cormat(:,:,k) = lat.*lat(mod(ind+k-1,N)+1,:);
    Cormat(:,:,N+k) = lat.*lat(:,mod(ind+k-1,N)+1);
end
