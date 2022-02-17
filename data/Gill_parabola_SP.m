%% This code predict Gill solution for Samoan Passage geometric conditions
% follow Gill_parabola.m
% wc=[0.6, 1.6, 2.8, 4.0]; bc =[0.2, 0.6, 0.8, 1.0, 1.4];
clear all
close all
clc
wc=[0.6, 1.6, 2.8, 4.0]; bc =[0.2, 0.6, 0.8, 1.0, 1.4];
potvort=1;
dc=0.5;
rc=dc./(wc./2).^2;
g=4e-4;% drho=0.04 based on Voet et al, 2015 Fig2 -- g=4e-4
r=6.371*10^6;lat=-8;
unit=2*pi*r/360;
Omega=7.292e-005;%2*pi/(24*3600);
f0=abs(2*Omega*sin(lat/180*pi));
V=nan(51,length(rc),length(bc));
H=nan(51,length(rc),length(bc));
D=nan(51,length(rc),length(bc));
X0=nan(51,length(rc),length(bc));
Q1=nan(length(rc),length(bc));
Q1_int=nan(length(rc),length(bc));
for j=1:length(bc)
for i=1:length(wc)
d2=1000;d3=rc(i)*f0^2/g;d4=bc(j)*1000;d5=10e+6;
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
if X(1)<=d5&&abs(G)<1e-6
Q1(i,j)=X(1)*abs(f0)/g/d2^2;wc1(i)=X(2);
end

Q1_b=X(1);w1c_b=X(2);
% now solve the critical section layer depth
x1=sinh(a3*w1c_b);
x2=cosh(a3*w1c_b);
% a+b:
x=w1c_b;
% a-b:
y=Q1_b*sqrt(a6)/(a7*(1-x2)/x1+x);
% a and b:
a=(x+y)/2;
b=(x-y)/2;
% -a<=x0<=b:
x0=-a:(b+a)/50:b;
% d ,h, v:
h=nan(length(x0),1);
d=nan(length(x0),1);
v=nan(length(x0),1);
for k=1:length(x0)
h(k,1)=d4+d3*x0(k).^2;
d(k,1)=(d2+2*d3*Ld^2)*(sinh((x0(k)-b)/Ld)-sinh((x0(k)+a)/Ld))/x1+d2+2*d3*Ld^2;
v(k,1)=-Ld*f0*(1+2*d3*g/f0^2)*(cosh((x0(k)-b)/Ld)-cosh((x0(k)+a)/Ld))/x1-2*d3*g*x0(k)/f0;
end

% integrate to obtain Q1 -- S.Tan, IOCAS, 2020/01/06
trans = cumsum(v.*d.*[0 diff(x0)]');
Q1_int(i,j) = trans(end)*f0/g/d2/d2;

% non-dimensionalize
V(:,i,j)=v./sqrt(g*d2);
H(:,i,j)=h./d2;
D(:,i,j)=d./d2;
X0(:,i,j)=x0./Ld;
end
end
Q=d5*abs(f0)/g/d2^2;
save Gill_parabola_SP.mat Q1 Q1_int wc wc1 Q V H D X0
