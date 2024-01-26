%% ASA Project 1 2100010793 ¿Óº—
% This is for calculating: 
% (c)(d) Correlation Length
% ( After calculating the Correlation Function in (c) )
%% (c)(d) Correlation Length
xi = zeros(tnum,1);     % xi(t) = correlation length at temperature t.
for index = 1:tnum
    xdata = ([0:N/4])'; ydata = log(abs(Cor(index,[N,1:N/4])'));
    coeff = LeastSquare(xdata,ydata);
    xi(index) = -1/coeff;
end