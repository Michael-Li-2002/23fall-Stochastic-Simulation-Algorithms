function pde = modelpde()
%MODELPDE 
pde = struct('f',@f,'b',@b,'bdval',@bdval,'inDomain',@inDomain,'exactu',@exactu,'findBdPoint',@findBdPoint);

function rhs = f(p)
    rhs = p(:,1).^2 + p(:,2).^2 + ones(size(p,1),1);
end  

function velo = b(p)
    velo = p;
end

function bd = bdval(p)
    bd = 0.5 * ones(size(p,1),1);
end

function isIn = inDomain(p)
    isIn = (p(:,1).^2 + p(:,2).^2 < 1);
end

function u = exactu(p)
    u = (p(:,1).^2 + p(:,2).^2)/2;
end

function [p,t] = findBdPoint(p0)
    lambda = sqrt(p0(:,1).^2 + p0(:,2).^2);
    p = p0 ./ lambda;
    t = 1./lambda - ones(size(p0,1),1);
end
end