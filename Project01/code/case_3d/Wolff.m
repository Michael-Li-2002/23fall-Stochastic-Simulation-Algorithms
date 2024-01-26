function [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,T)
%WOLFF 
beta = 1/kb/T;
N = size(lat,1); 
i = randi(N);
j = randi(N);
k = randi(N);
color = lat(i,j,k);
prob = 1 - exp(-beta*J);

not_visited = boolean(ones(N,N,N));
checklist = zeros(N^3,3);
checklist(1,1) = i; checklist(1,2) = j; checklist(1,3) = k;
cluster = zeros(N^3,3);
cluster(1,1) = i; cluster(1,2) = j; cluster(1,3) = k;

checkind = 0; checkend = 1; clusind = 1;
while checklist(checkind+1,1) ~= 0
    checkind = checkind + 1;
    i = checklist(checkind,1); j = checklist(checkind,2); k = checklist(checkind,3);
    not_visited(i,j,k) = 0;
    for dx = [-1,1]
        if not_visited(mod(i+dx-1,N)+1,j,k) && lat(mod(i+dx-1,N)+1,j,k) == color ...
                && rand() <= prob
            clusind = clusind + 1;
            cluster(clusind,1) = mod(i+dx-1,N)+1; cluster(clusind,2) = j; cluster(clusind,3) = k;
            checkend = checkend + 1;
            checklist(checkend,1) = mod(i+dx-1,N)+1; checklist(checkend,2) = j; checklist(checkend,3) = k;
            not_visited(mod(i+dx-1,N)+1,j,k) = 0;
        end
        % not_visited(mod(i+dx-1,N)+1,j) = 0;
        
        if not_visited(i,mod(j+dx-1,N)+1,k) && lat(i,mod(j+dx-1,N)+1,k) == color ...
                && rand() <= prob
            clusind = clusind + 1;
            cluster(clusind,1) = i; cluster(clusind,2) = mod(j+dx-1,N)+1; cluster(clusind,3) = k;
            checkend = checkend + 1;
            checklist(checkend,1) = i; checklist(checkend,2) = mod(j+dx-1,N)+1; checklist(checkend,3) = k;
            not_visited(i,mod(j+dx-1,N)+1,k) = 0;
        end
        % not_visited(i,mod(j+dx-1,N)+1) = 0;
        
        if not_visited(i,j,mod(k+dx-1,N)+1) && lat(i,j,mod(k+dx-1,N)+1) == color ...
                && rand() <= prob
            clusind = clusind + 1;
            cluster(clusind,1) = i; cluster(clusind,2) = j; cluster(clusind,3) = mod(k+dx-1,N)+1;
            checkend = checkend + 1;
            checklist(checkend,1) = i; checklist(checkend,2) = j; checklist(checkend,3) = mod(k+dx-1,N)+1;
            not_visited(i,j,mod(k+dx-1,N)+1) = 0;
        end
    end
    if checkind >= N^3
        break
    end
end
if h == 0
    newval = max(ceil((q-1)*rand()),1);
    if newval >= color
        newval = newval + 1;
    end
else
   % distribution = zeros(q,1);
    if h > 0
        distribution = exp(-h*beta*clusind*(q-1:-1:0));
        newval = q;
    elseif h < 0
        distribution = exp( h*beta*clusind*(1:q));
        newval = 1;
    end
    distribution = distribution./sum(distribution);
    sump = 0;
    r = rand();
    
    for i=1:q
        sump = sump + distribution(i);
        if r <= sump
            newval = i;
            break
        end
    end
% if h == 0
%     newval = max(ceil((q-1)*rand()),1);
%     if newval >= color
%         newval = newval + 1;
%     end
% else
%     distribution = zeros(q,1);
%     for i=2:color
%         distribution(i) = distribution(i-1) + exp(h*beta*clusind*(i-1-color)/2);
%     end
%     for i=color+1:q
%         distribution(i) = distribution(i-1) + exp(h*beta*clusind*(i-color)/2);
%     end
%     distribution = distribution ./ distribution(q);
%     r = rand();
%     for i=1:q
%         if r <= distribution(i)
%             if i <= color
%                 newval = max(1,i-1);
%             else
%                 newval = i;
%             end
%         end
%     end
    
end

