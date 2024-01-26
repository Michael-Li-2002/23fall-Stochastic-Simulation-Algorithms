function [k,b,r2] = LeastSquare(x,y)
%LEASTSQUARE 
n = size(x,1);
A = [x,ones(n,1)];
sol = (A'*A)\(A'*y);
k = sol(1); b = sol(2);
L2err = norm( y - A * [k;b],2)^2;
SST = norm(y-mean(y),2)^2;
r2 = 1-L2err/SST;
end

