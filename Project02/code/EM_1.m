function U = EM_1(pde,point,tau)
%EM_1: Euler-Maruyama, 及简单的确定 Exit Time 的方式
x0 = point; 
x1 = x0 + pde.b(x0) * tau + randn(1,2) * sqrt(tau);
intf = 0;
while pde.inDomain(x1)
    if rand() >= exp(-2*(1-sqrt(sum(x0.^2)))*(1-sqrt(sum(x1.^2)))/tau)
        intf = intf + pde.f(x0) * tau;
        x0 = x1;
        x1 = x1 + pde.b(x1) * tau + randn(1,2) * sqrt(tau);
    else
        break
    end
end
U = 1/2 - intf;
