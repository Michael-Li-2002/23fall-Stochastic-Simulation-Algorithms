%% ASA Project 1 2100010793 ¿Óº—
% This is for calculating
% (a) Internal Energy and Specific Heat
% (c) Correlation Function
% IN 3D CASE.
%% Parameters
q = 3; N = 10; J = 1; kb = 1; h = 0;
% temperature list
t0 = 0.05;
tstep_1 = 0.05; t1 = 1.65; tnum_1 = ceil((t1-t0)/tstep_1);
T1 = t0 + tstep_1 * (1: tnum_1);
tstep_2 = 0.01; t2 = 1.95; tnum_2 = ceil((t2-t1)/tstep_2);
T2 = t1 + tstep_2 * (1: tnum_2);
tstep_3 = 0.05; t3 = 3; tnum_3 = ceil((t3-t2)/tstep_3);
T3 = t2 + tstep_3 * (1: tnum_3);
T = [T1,T2,T3];
tnum = tnum_1 + tnum_2 + tnum_3;
% pretime and sampletime setting
pretime3 = 8000; sampletime_3 = 280000;
pretime2 = 800; pretime1 = 800; sampletime_1 = 4000; sampletime_2 = 4000; 
u = zeros(tnum,1);
c = zeros(tnum,1); 
Cor = zeros(tnum,N);
%% Calculations: Internal Energy , Specific Heat, Correlation Function
lat = randi(q,N,N,N);         % initial lattice: uniform distribution( T = +\infty)
H0 = Hamilton(lat,J,h);       % initial Hamiltonian
Cormat = Correlation(lat);    % initial \langle \sigma_i\sigma_j \rangle
tic;
for index = tnum_1+ tnum_2 + tnum_3: -1: 1 + tnum_1 + tnum_2   % Annealing from high temperature
    tem = T(index); U = 0; C = 0;
    beta = 1/kb/tem;
    Si = zeros(N,N,N); Sij = zeros(N,N,N,3*N);
    for tick = 1: pretime3
        [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,tem);
        H1 = H0;
        if sign == 1
            H1 = H0 + deltaH;           
            i = position(1); j = position(2); k = position(3);
            lat(i,j,k) = newval;
            for kind = 1:N
                Cormat(i,j,k,kind) = lat(mod(i+kind-1,N)+1,j,k)*lat(i,j,k);
                Cormat(mod(i-kind-1,N)+1,j,k,kind) = lat(mod(i-kind-1,N)+1,j,k)*lat(i,j,k);
                Cormat(i,j,k,N+kind) = lat(i,mod(j+kind-1,N)+1,k)*lat(i,j,k);
                Cormat(i,mod(j-kind-1,N)+1,k,N+kind) = lat(i,mod(j-kind-1,N)+1,k)*lat(i,j,k);
                Cormat(i,j,k,2*N+kind) = lat(i,j,mod(k+kind-1,N)+1)*lat(i,j,k);
                Cormat(i,j,mod(k-kind-1,N)+1,2*N+kind) = lat(i,j,mod(k-kind-1,N)+1)*lat(i,j,k);
            end
        end
        H0 = H1;
    end
    for tick = 1: sampletime_3
        [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,tem);
        H1 = H0;
        if sign == 1
            H1 = H0 + deltaH;           
            i = position(1); j = position(2); k = position(3);
            lat(i,j,k) = newval;
            for kind = 1:N
                Cormat(i,j,k,kind) = lat(mod(i+kind-1,N)+1,j,k)*lat(i,j,k);
                Cormat(mod(i-kind-1,N)+1,j,k,kind) = lat(mod(i-kind-1,N)+1,j,k)*lat(i,j,k);
                Cormat(i,j,k,N+kind) = lat(i,mod(j+kind-1,N)+1,k)*lat(i,j,k);
                Cormat(i,mod(j-kind-1,N)+1,k,N+kind) = lat(i,mod(j-kind-1,N)+1,k)*lat(i,j,k);
                Cormat(i,j,k,2*N+kind) = lat(i,j,mod(k+kind-1,N)+1)*lat(i,j,k);
                Cormat(i,j,mod(k-kind-1,N)+1,2*N+kind) = lat(i,j,mod(k-kind-1,N)+1)*lat(i,j,k);
            end
        end
        H0 = H1;
        U = U + H1; C = C + H1^2;
        Si = Si + lat; Sij = Sij + Cormat;
        H0 = H1; 
    end
    u(index) = U/sampletime_3/(N^3); 
    c(index) = kb*beta^2 * (C/sampletime_3-(U/sampletime_3)^2)/(N^3);
    Nind = 1:N;
    for k = 1:N
        Cor(index,k) = 1/(2*N^3*sampletime_3)*(sum(sum(sum(Sij(:,:,:,k)))) + sum(sum(sum(Sij(:,:,:,N+k)))) + sum(sum(sum(Sij(:,:,:,2*N+k)))) ...
            - sum(sum(sum(Si.*Si(mod(Nind+k-1,N)+1,:,:))))/sampletime_3 - sum(sum(sum(Si.*Si(:,mod(Nind+k-1,N)+1,:))))/sampletime_3 - sum(sum(sum(Si.*Si(:,:,mod(Nind+k-1,N)+1))))/sampletime_3);
    end
    disp(['finish tem = ',num2str(tem),' / ',num2str(t3)])
end

for index = tnum_1 + tnum_2: -1: 1 + tnum_1
    tem = T(index); U = 0; C = 0; Si = zeros(N,N,N); Sij = zeros(N,N,N,3*N);
    beta = 1/kb/tem;
    for tick = 1: pretime2
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2); zz = cluster(ind,3);
            lat(xx,yy,zz) = newval;
        end
    end
    for tick = 1: sampletime_2
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2); zz = cluster(ind,3);
            lat(xx,yy,zz) = newval;
        end
        H1 = Hamilton(lat,J,h);
        Cormat = Correlation(lat);
        U = U + H1; C = C + H1^2; Si = Si + lat; Sij = Sij + Cormat;
        H0 = H1; 
    end
    u(index) = U/sampletime_2/(N^3); 
    c(index) = kb*beta^2 * (C/sampletime_2-(U/sampletime_2)^2)/(N^3);
    Nind = 1:N;
    for k = 1:N
        Cor(index,k) = 1/(3*N^3*sampletime_2)*(sum(sum(sum(Sij(:,:,:,k)))) + sum(sum(sum(Sij(:,:,:,N+k)))) + sum(sum(sum(Sij(:,:,:,2*N+k)))) ...
            - sum(sum(sum(Si.*Si(mod(Nind+k-1,N)+1,:,:))))/sampletime_2 - sum(sum(sum(Si.*Si(:,mod(Nind+k-1,N)+1,:))))/sampletime_2 - sum(sum(sum(Si.*Si(:,:,mod(Nind+k-1,N)+1))))/sampletime_2);
    end
    disp(['finish tem = ',num2str(tem),' / ',num2str(t3)])
end

for index = tnum_1: -1: 1
    tem = T(index); U = 0; C = 0; Si = zeros(N,N,N); Sij = zeros(N,N,N,3*N);
    beta = 1/kb/tem;
    for tick = 1: pretime1
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2); zz = cluster(ind,3);
            lat(xx,yy,zz) = newval;
        end
    end
    for tick = 1: sampletime_1
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2); zz = cluster(ind,3);
            lat(xx,yy,zz) = newval;
        end
        H1 = Hamilton(lat,J,h); Cormat = Correlation(lat);
        U = U + H1; C = C + H1^2; Si = Si + lat; Sij = Sij + Cormat;
        H0 = H1; 
    end
    u(index) = U/sampletime_1/(N^3); 
    c(index) = kb*beta^2 * (C/sampletime_1-(U/sampletime_1)^2)/(N^3);
    Nind = 1:N;
    for k = 1:N
        Cor(index,k) = 1/(3*N^3*sampletime_1)*(sum(sum(sum(Sij(:,:,:,k)))) + sum(sum(sum(Sij(:,:,:,N+k)))) + sum(sum(sum(Sij(:,:,:,2*N+k)))) ...
            - sum(sum(sum(Si.*Si(mod(Nind+k-1,N)+1,:,:))))/sampletime_1 - sum(sum(sum(Si.*Si(:,mod(Nind+k-1,N)+1,:))))/sampletime_1 - sum(sum(sum(Si.*Si(:,:,mod(Nind+k-1,N)+1))))/sampletime_1);
    end
    disp(['finish tem = ',num2str(tem),' / ',num2str(t3)])
end
