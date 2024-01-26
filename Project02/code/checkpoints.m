function points = checkpoints()
%CHECKPOINTS 计算数值解的点

rN = 3;
tN = 4;
r = ones(tN,1) * (1:rN)/(rN+1);
theta = 2*pi * (1:tN)'/(tN) * ones(1,rN);
points = [r(:).* cos(theta(:)), r(:).* sin(theta(:))];

