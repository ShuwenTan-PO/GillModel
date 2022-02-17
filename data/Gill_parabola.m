%% This code provide the solution for the nondimensional, parobolic cross
% -section Gill's model, applied on Samoan Passage
%
% FUNCTION:
% solve max Q1 from constraint G(a+b,Q)=0
% 
% INPUTS:
% q=-1;
% w1c=[.1:.02:2 2.05:.5:10];
% Q=[0.5];
% h1c=[.01 .05:.05:2];
%
% OUTPUTS:
% Q1
%
% CALLS:
% [X,fval,exitflag]=Gill_nd_fmincon(a1,a2,a3,a4,a5,a6,a7,a8)
%
% TO NOTICE:
% 
% S.Tan, Somerville, 2019/08/18
% exclude left wall reversal flow   -- S.Tan 2019/08/18                 
%% solve from dimensionalized 
cacac
g=4e-4;% drho=0.04 based on Voet et al, 2015 Fig2 -- g=4e-4
lat=-8;
Di=1000;
Q=10e+6;
wc=[.1:.02:5];
bc=[.01 .05:.05:2];potvort=1;
dc=0.5;
rc=dc./(wc./2).^2;
Q1=nan(length(rc),length(bc));
r=6.371*10^6;
unit=2*pi*r/360;
Omega=7.292e-005;%2*pi/(24*3600);
f0=abs(2*Omega*sin(lat/180*pi));   
Q1_1=Q1;Q1_2=Q1;
aplsb=nan(length(rc),length(bc));
amnsb=nan(length(rc),length(bc));
vl=nan(length(rc),length(bc));vr=nan(length(rc),length(bc));
dl=nan(length(rc),length(bc));dr=nan(length(rc),length(bc));
Q1_int=Q1;
% solve for Q1
        for k=1:length(rc)
            for m=1:length(bc)
                    d2=Di;d3=rc(k)*f0^2/g;d4=bc(m).*1000;d5=Q;
                    Ld=sqrt(g*d2)/f0;
                    a1=g*d4-f0*d5/d2-g*d2;
                    a2=1/2*(Ld*f0*(1+2*d3*g/f0^2))^2;
                    a3=1/Ld;
                    a4=d3*g*Ld*(1+2*d3*g/f0^2);
                    a5=d3^2*g^2/2/f0^2+d3*g/4;
                    a6=(f0/d3/g/(d2+2*d3*Ld^2))^2;
                    a7=2*Ld;
                    a8=f0/2/d2;                    
                    [X,fval,exitflag]=Gill_d_fmincon(a1,a2,a3,a4,a5,a6,a7,a8);
                    G=a1+a2.*(cosh(a3.*X(2))-1).^2./((sinh(a3.*X(2))).^2)-a4.*(cosh(a3.*X(2))-1).*X(2)./sinh(a3.*X(2))+a5.*(X(2).^2+a6.*X(1).^2./(a7.*(1./sinh(a3.*X(2))-coth(a3.*X(2)))+X(2)).^2)+a8.*X(1);

% exclude left wall reversal flow   -- S.Tan 2019/08/18                 
Q1_b=X(1);w1c_b=X(2);
% now solve the critical section layer depth
x1=sinh(a3*w1c_b);
x2=cosh(a3*w1c_b);
% a+b:
x=w1c_b;
% a-b:
y=Q1_b*sqrt(a6)/(a7*(1-x2)/x1+x);
aplsb(k,m)=x;amnsb(k,m)=y;
% a and b:
a=(x+y)/2;
b=(x-y)/2;
% -a<=x0<=b:
x0=-a:(b+a)/50:b;
% d ,h, v:
h=nan(length(x0),1);
d=nan(length(x0),1);
v=nan(length(x0),1);
for i=1:length(x0)
h(i,1)=d4+d3*x0(i).^2;
d(i,1)=(d2+2*d3*Ld^2)*(sinh((x0(i)-b)/Ld)-sinh((x0(i)+a)/Ld))/x1+d2+2*d3*Ld^2;
v(i,1)=-Ld*f0*(1+2*d3*g/f0^2)*(cosh((x0(i)-b)/Ld)-cosh((x0(i)+a)/Ld))/x1-2*d3*g*x0(i)/f0;
end
% integrate to obtain Q1 -- S.Tan, IOCAS, 2020/01/06
trans = cumsum(v.*d.*[0 diff(x0)]');
Q1_int(k,m) = trans(end)*f0/g/Di/Di;

vl(k,m)=v(1);vr(k,m)=v(end);
dl(k,m)=d(1);dr(k,m)=d(end);
if v(1)<0
    Q1_1(k,m)=X(1);
else
    Q1_2(k,m)=X(1);
end
                    if X(1)<=Q&&abs(G)<1e-6
                    Q1(k,m)=X(1);w1c(k,m)=X(2);
                    end
            end
                disp(k)
end
rc=1./rc;
save Gill_parabola.mat Q1 Q1_int Q1_1 Q1_2 w1c Q bc dc rc g Di aplsb amnsb vl vr dl dr

Q1 = Q1.*f0./g./Di./Di;
a= Q1-Q1_int;
max(max(abs(a)))%Q1 and Q1_int don't differ!

