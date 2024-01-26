%% ASA Project 1 2100010793 ¿Óº—
% This is for calculating: 
% (a) Internal Energy and Specific Heat
% (c) Correlation Function
%% Parameters
q = 3; N = 20; J = 1; kb = 1; h = 0;

% Temperature list
t0 = 0.05;
tstep_1 = 0.05; t1 = 0.85; tnum_1 = ceil((t1-t0)/tstep_1);
T1 = t0 + tstep_1 * (1: tnum_1);                             % Low tempereature
tstep_2 = 0.01; t2 = 1.15; tnum_2 = ceil((t2-t1)/tstep_2);
T2 = t1 + tstep_2 * (1: tnum_2);                             % Near Critical temperature
tstep_3 = 0.05; t3 = 3; tnum_3 = ceil((t3-t2)/tstep_3);
T3 = t2 + tstep_3 * (1: tnum_3);                             % High temperature
T = [T1,T2,T3];
tnum = tnum_1 + tnum_2 + tnum_3;

% pretime: time before sampling to achieve invariant distribution
% sampletime: time for sampling and calculating
pretime3 = 8000; sampletime_3 = 150000;
pretime2 = 800; pretime1 = 800; sampletime_1 = 2500; sampletime_2 = 2500; 
u = zeros(tnum,1);   % Internal Energy      
c = zeros(tnum,1);   % Specific Heat
Cor = zeros(tnum,N); % Correlation function: Cor(t,k) = \Gamma(k) at temperature t
%% Calculations: Internal Energy , Specific Heat, Correlation function
lat = randi(q,N,N);         % initial lattice: uniform distribution( T = +\infty )
H0 = Hamilton(lat,J,h);     % Recording Hamiltonian
Cormat = Correlation(lat);  % Recording \langle \sigma_i \sigma_j \rangle
tic;
for index = tnum_1+ tnum_2 + tnum_3: -1: 1 + tnum_1 + tnum_2   % Annealing from high temperature
    tem = T(index); U = 0; C = 0;
    beta = 1/kb/tem;
    Si = zeros(N,N); Sij = zeros(N,N,2*N);
    for tick = 1: pretime3
        [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,tem);
        H1 = H0;
        if sign == 1
            H1 = H0 + deltaH;           
            lat(position(1),position(2)) = newval;
            i = position(1); j = position(2);
            for k = 1:N
                Cormat(i,j,k) = lat(mod(i+k-1,N)+1,j)*lat(i,j);
                Cormat(mod(i-k-1,N)+1,j,k) = lat(mod(i-k-1,N)+1,j)*lat(i,j);
                Cormat(i,j,N+k) = lat(i,mod(j+k-1,N)+1)*lat(i,j);
                Cormat(i,mod(j-k-1,N)+1,N+k) = lat(i,mod(j-k-1,N)+1)*lat(i,j);
            end
        end
        H0 = H1;
    end
    for tick = 1: sampletime_3
        [newval,position,sign,deltaH] = metropolis(lat,q, J, h, kb,tem);
        H1 = H0;
        if sign == 1
            H1 = H0 + deltaH;
            lat(position(1),position(2)) = newval;
            
            i = position(1); j = position(2);
            for k = 1:N
                Cormat(i,j,k) = lat(mod(i+k-1,N)+1,j)*lat(i,j);
                Cormat(mod(i-k-1,N)+1,j,k) = lat(mod(i-k-1,N)+1,j)*lat(i,j);
                Cormat(i,j,N+k) = lat(i,mod(j+k-1,N)+1)*lat(i,j);
                Cormat(i,mod(j-k-1,N)+1,N+k) = lat(i,mod(j-k-1,N)+1)*lat(i,j);
            end
        end
        U = U + H1; C = C + H1^2;
        Si = Si + lat; Sij = Sij + Cormat;
        H0 = H1; 
    end
    u(index) = U/sampletime_3/(N^2); 
    c(index) = kb*beta^2 * (C/sampletime_3-(U/sampletime_3)^2)/(N^2);
    Nind = 1:N;
    for k = 1:N
        Cor(index,k) = 1/(2*N^2*sampletime_3)*(sum(sum(Sij(:,:,k))) + sum(sum(Sij(:,:,N+k))) ...
            - sum(sum(Si.*Si(mod(Nind+k-1,N)+1,:)))/sampletime_3 - sum(sum(Si.*Si(:,mod(Nind+k-1,N)+1)))/sampletime_3);
    end
    disp(['finish tem = ',num2str(tem),' / ',num2str(t3)])
end

for index = tnum_1 + tnum_2: -1: 1 + tnum_1       % Near Critical Temperature
    tem = T(index); U = 0; C = 0; Si = zeros(N,N); Sij = zeros(N,N,2*N);
    beta = 1/kb/tem;
    for tick = 1: pretime2
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem); % Using Wolff
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2);
            lat(xx,yy) = newval;
        end
    end
    for tick = 1: sampletime_2
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2);
            lat(xx,yy) = newval;
        end
        H1 = Hamilton(lat,J,h);
        Cormat = Correlation(lat);
        U = U + H1; C = C + H1^2; Si = Si + lat; Sij = Sij + Cormat;
        H0 = H1; 
    end
    u(index) = U/sampletime_2/(N^2); 
    c(index) = kb*beta^2 * (C/sampletime_2-(U/sampletime_2)^2)/(N^2);
    Nind = 1:N;
    for k = 1:N
        Cor(index,k) = 1/(2*N^2*sampletime_2)*(sum(sum(Sij(:,:,k))) + sum(sum(Sij(:,:,N+k))) ...
            - sum(sum(Si.*Si(mod(Nind+k-1,N)+1,:)))/sampletime_2 - sum(sum(Si.*Si(:,mod(Nind+k-1,N)+1)))/sampletime_2);
        
    end
    disp(['finish tem = ',num2str(tem),' / ',num2str(t3)])
end

for index = tnum_1: -1: 1         % Low Temperature
    tem = T(index); U = 0; C = 0; Si = zeros(N,N); Sij = zeros(N,N,2*N);
    beta = 1/kb/tem;
    for tick = 1: pretime1
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2);
            lat(xx,yy) = newval;
        end
    end
    for tick = 1: sampletime_1
        [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,tem);
        for ind = 1:clusind
            xx = cluster(ind,1); yy = cluster(ind,2);
            lat(xx,yy) = newval;
        end
        H1 = Hamilton(lat,J,h); Cormat = Correlation(lat);
        U = U + H1; C = C + H1^2; Si = Si + lat; Sij = Sij + Cormat;
        H0 = H1; 
    end
    u(index) = U/sampletime_1/(N^2); 
    c(index) = kb*beta^2 * (C/sampletime_1-(U/sampletime_1)^2)/(N^2);
    Nind = 1:N;
    for k = 1:N
        Cor(index,k) = 1/(2*N^2*sampletime_1)*(sum(sum(Sij(:,:,k))) + sum(sum(Sij(:,:,N+k))) ...
         - sum(sum(Si.*Si(mod(Nind+k-1,N)+1,:)))/sampletime_1 - sum(sum(Si.*Si(:,mod(Nind+k-1,N)+1)))/sampletime_1);

    end
    disp(['finish tem = ',num2str(tem),' / ',num2str(t3)])
end
