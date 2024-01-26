function [newval,cluster,clusind] = Wolff(lat,q, J, h, kb,T)
%WOLFF 
beta = 1/kb/T;
N = size(lat,1); 
i = randi(N);
j = randi(N);
color = lat(i,j);
prob = 1 - exp(-beta*J);

not_visited = boolean(ones(N,N));
checklist = zeros(N^2,2);
checklist(1,1) = i; checklist(1,2) = j;
cluster = zeros(N^2,2);
cluster(1,1) = i; cluster(1,2) = j;

checkind = 0; checkend = 1; clusind = 1;
while checklist(checkind+1,1) ~= 0
    checkind = checkind + 1;
    i = checklist(checkind,1); j = checklist(checkind,2);
    not_visited(i,j) = 0;
    for dx = [-1,1]
        if not_visited(mod(i+dx-1,N)+1,j) && lat(mod(i+dx-1,N)+1,j) == color ...
                && rand() <= prob
            clusind = clusind + 1;
            cluster(clusind,1) = mod(i+dx-1,N)+1; cluster(clusind,2) = j;
            checkend = checkend + 1;
            checklist(checkend,1) = mod(i+dx-1,N)+1; checklist(checkend,2) = j;
            not_visited(mod(i+dx-1,N)+1,j) = 0;
        end
        % not_visited(mod(i+dx-1,N)+1,j) = 0;
        
        if not_visited(i,mod(j+dx-1,N)+1) && lat(i,mod(j+dx-1,N)+1) == color ...
                && rand() <= prob
            clusind = clusind + 1;
            cluster(clusind,1) = i; cluster(clusind,2) = mod(j+dx-1,N)+1;
            checkend = checkend + 1;
            checklist(checkend,1) = i; checklist(checkend,2) = mod(j+dx-1,N)+1;
            not_visited(i,mod(j+dx-1,N)+1) = 0;
        end
        % not_visited(i,mod(j+dx-1,N)+1) = 0;
    end
    if checkind >= N^2
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
%     if r>=sump
%         newval = q;
%     end

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

