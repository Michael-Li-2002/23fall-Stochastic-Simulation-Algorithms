function Linf_err = Linferr(pde,point,uh)
%LINFERR 
u = pde.exactu(point)';
Linf_err = max(abs(u-uh));
