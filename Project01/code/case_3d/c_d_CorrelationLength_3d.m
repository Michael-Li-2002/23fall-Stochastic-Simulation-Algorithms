%% ASA Project 1 2100010793 ¿Óº—
% This is for calculating
% (c)(d) Correlation Length  
% ( After calculating correlation function in (c) )
% IN 3D CASE.
%% (c)(d) Correlation Length
xi = zeros(tnum,1);
for index = 1:tnum
    xdata = ([0:3])'; ydata = log(abs(Cor(index,[N,1:3])'));
    coeff = LeastSquare(xdata,ydata);
    xi(index) = -1/coeff;
end
