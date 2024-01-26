%% ASA Project 1 2100010793 Àî¼Ñ
% This is for calculating
% (b) Magnetization
% IN 3D CASE.
%% Parameters
d = 3; % dimension
q = 3; N = 20; 
J = 1; kb = 1;  T = 1; 

hstep = 0.02; hmax = 3; hnum = ceil(hmax/hstep);
hlist = hstep * (-hnum:hnum);
pretime = 10000; sampletime = 20000;
pretime2 = 100; sampletime_2 = 200;
newpretime = 40000; newsampletime = 20000;
newpretime2 = 100; newsampletime_2 = 200;
m = zeros(1,2*hnum + 1);

%% (b) Magnetization
tic;

if T > 1.85
    lat = randi(q,N,N,N);
else
    lat = ones(N,N,N);
end
M0 = sum(sum(sum(lat)));


if T > 1.65
    for hindex = 1:2*hnum +1
        tem = T; h = hlist(hindex); M = 0;
        beta = 1/kb/tem;
        if abs(hindex * hstep + hlist(1))<0.1
            thispretime = newpretime;
            thissampletime = newsampletime;
        else
            thispretime = pretime;
            thissampletime = sampletime;
        end
        for tick = 1: thispretime
            [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,tem);
            if sign == 1
                if d == 2
                    i = position(1); j = position(2);
                    M0 = M0 - lat(i,j) + newval;
                    lat(i,j) = newval;
                elseif d == 3
                    i = position(1); j = position(2); k = position(3); 
                    M0 = M0 - lat(i,j,k) + newval;
                    lat(i,j,k) = newval;
                end
            end

        end
        for tick = 1: thissampletime
            [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,tem);
            if sign == 1
                if d == 2
                    i = position(1); j = position(2);
                    M0 = M0 - lat(i,j) + newval;
                    lat(i,j) = newval;
                elseif d == 3
                    i = position(1); j = position(2); k = position(3);
                    M0 = M0 - lat(i,j,k) + newval;
                    lat(i,j,k) = newval;
                end
            end
            M = M + M0;
        end
        m(hindex) = M/thissampletime/N^d;
        disp(['finish h = ',num2str(h),'/3'])
    end
    
else
    for hindex = 1:2*hnum+1
        tem = T; h = hlist(hindex); M = 0;
        beta = 1/kb/tem;
        if abs(hindex * hstep + hlist(1))<0.1
            thispretime = newpretime2;
            thissampletime = newsampletime_2;
        else
            thispretime = pretime2;
            thissampletime = sampletime_2;
        end
        for tick = 1: thispretime
            [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
            color = lat(cluster(1,1),cluster(1,2),cluster(1,3));
            M0 = M0 + (newval-color) * clusind;
            for ind = 1:clusind
                xx = cluster(ind,1); yy = cluster(ind,2); zz = cluster(ind,3);
                lat(xx,yy,zz) = newval;
            end
        end
        for tick = 1: thissampletime
            [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
            color = lat(cluster(1,1),cluster(1,2),cluster(1,3));
            M0 = M0 + (newval-color) * clusind;
            for ind = 1:clusind
                xx = cluster(ind,1); yy = cluster(ind,2); zz = cluster(ind,3);
                lat(xx,yy,zz) = newval;
            end
            M = M + M0;
        end
        m(hindex) = M/thissampletime/N^d;
        disp(['finish h = ',num2str(h),'/3'])
    end
    
end
t = toc
%%
plot(hlist,m);