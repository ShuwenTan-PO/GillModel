function [X,fval,exitflag]=Gill_d_fmincon(a1,a2,a3,a4,a5,a6,a7,a8)
% Find minimum of constrained nonlinear multivariable function
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');%,'MaxIterations',1e+4,'Algorithm','sqp'
fun=@(x)-x(1)+0*x(2);
A=[];b=[];Aeq=[];beq=[];lb=[0,0];ub=[];
% c = [];
% ceq =@(x)a1+a2*(cosh(a3*x(2))-1)^2/((sinh(a3*x(2)))^2)-a4*(cosh(a3*x(2))-1)*x(2)/sinh(a3*x(2))+a5*(x(2)^2+a6*x(1)^2/(a7*(1/sinh(a3*x(2))-coth(a3*x(2)))+x(2))^2)+a8*x(1);
% nonlcon = @(x)deal(c,ceq(x));
c = @(x)(a1+a2*(cosh(a3*x(2))-1)^2/((sinh(a3*x(2)))^2)-a4*(cosh(a3*x(2))-1)*x(2)/sinh(a3*x(2))+a5*(x(2)^2+a6*x(1)^2/(a7*(1/sinh(a3*x(2))-coth(a3*x(2)))+x(2))^2)+a8*x(1))^2;
ceq =[];
nonlcon = @(x)deal(c(x),ceq);
x0=[2e+6,1e+4];
[X,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

