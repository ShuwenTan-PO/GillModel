%% This code provide the analitical solution for the two dimensional, parobolic cross
% -section Gill's model, applied on Samoan Passage
%
% FUNCTION:
% solve max Q1 from constraint G(a+b,Q)=0
% 
% INPUTS:
% q=-1;
% w1c=[0.6 1.2];
% Q=[0.5];
% h1c=[0.6 0.4];
% d=0.5
%
% OUTPUTS:
% Q11; Q12
%
% CALLS:
% [X,fval,exitflag]=Gill_nd_fmincon(a1,a2,a3,a4,a5,a6,a7,a8)
%
% TO NOTICE:
% This is one step further beyond the original parabolic Gill's model
% The flow tend to be critical in the western path, and if the channel
% cannot accomandate the whole inflow, the rest will enter the eastern path
% To do so, we first compute the parabolic Gill's solution at the western
% path, then compute the Bernoulli function at the boundary, assume
% Bernoulli at boundaries are the same, then substitute 
% B_mean = -|q|*Q12/2 + Bc(b) and compute the optimized Q12
%
% S.Tan, IOCAS, 2020/01/13
%% dimensional calculation - then transform to nondim
cacac
wc_w=[.1:.05:2];
wc_e=[.1:.05:2];
bc =[0.6 0.4];
potvort=-1;
dc=0.5;
rc_w=dc./(wc_w./2).^2;
rc_e=dc./(wc_e./2).^2;
g=4e-4;% drho=0.04 based on Voet et al, 2015 Fig2 -- g=4e-4
r=6.371*10^6;lat=-8;
unit=2*pi*r/360;
Omega=7.292e-005;%2*pi/(24*3600);
f0=abs(2*Omega*sin(lat/180*pi));

V=nan(2,51,length(rc_w),length(rc_e));
H=nan(2,51,length(rc_w),length(rc_e));
D=nan(2,51,length(rc_w),length(rc_e));
X0=nan(2,51,length(rc_w),length(rc_e));
%% western channel
Q1=nan(length(wc_w),length(wc_e));
Q2=nan(length(wc_w),length(wc_e));
Q1_int=nan(length(wc_w),length(wc_e));
Q2_int=nan(length(wc_w),length(wc_e));
Ber=nan(length(wc_w),length(wc_e));

for im=1:length(rc_w)
    for jm=1:length(rc_e)
        j=1;
d2=1000;d3=rc_w(im)*f0^2/g;d4=bc(j)*1000;d5=10e+6;
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
Q1(im,jm)=X(1)*abs(f0)/g/d2^2;wc1(1)=X(2);
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
% integrate to obtain Q1
trans = cumsum(v.*d.*[0 diff(x0)]');
Q1_int(im,jm) = trans(end)*f0/g/d2/d2;
% non-dimensionalize variables, western channel
V(1,:,im,jm)=v./sqrt(g*d2);
H(1,:,im,jm)=h./d2;
D(1,:,im,jm)=d./d2;
X0(1,:,im,jm)=x0./Ld;
% Bernoulli at b
% from equations 0.5356 - nondimensional: Bb/(g*d2)=1.3390
Bb=1/2*(-Ld*f0*(1+2*d3*g/f0^2)*(1-cosh(a3*wc1(1)))/sinh(a3*wc1(1))-2*d3*g*b/f0)^2+g*d4+g*d3*b^2;
% from results 0.5356
hb=d4+d3*b.^2;
db=(d2+2*d3*Ld^2)*(sinh((b-b)/Ld)-sinh((b+a)/Ld))/x1+d2+2*d3*Ld^2;
vb=-Ld*f0*(1+2*d3*g/f0^2)*(cosh((b-b)/Ld)-cosh((b+a)/Ld))/x1-2*d3*g*b/f0;
Bb=1/2*vb^2+g*(hb+db);
Ber(im,jm)=Bb;
%% eastern channel -- replace a1 to a1=g*d4-Bb;
j=2;
d2=1000;d3=rc_e(jm)*f0^2/g;d4=bc(j)*1000;d5=10e+6;
Ld=sqrt(g*d2)/f0;
a1=g*d4-Bb;
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
Q2(im,jm)=X(1)*abs(f0)/g/d2^2;wc1(2)=X(2);
end
Q2_b=X(1);w2c_b=X(2);
% now solve the critical section layer depth
x1=sinh(a3*w2c_b);
x2=cosh(a3*w2c_b);
% a+b:
x=w2c_b;
% a-b:
y=Q2_b*sqrt(a6)/(a7*(1-x2)/x1+x);
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
% integrate to obtain Q1
trans = cumsum(v.*d.*[0 diff(x0)]');
Q2_int(im,jm) = trans(end)*f0/g/d2/d2;
% non-dimensionalize variables, western channel
V(2,:,im,jm)=v./sqrt(g*d2);
H(2,:,im,jm)=h./d2;
D(2,:,im,jm)=d./d2;
X0(2,:,im,jm)=x0./Ld;
    end
end
% _w1.mat bc1=bc2=0.6
% _w2.mat bc1=0.6 bc2=0.4
save Gill_parabolic_nd_2channel_w2.mat Q1 Q2 Q1_int Q2_int rc_w rc_e wc_w wc_e bc dc Ber V H D X0
figure
plot(wc_e,Q1(find(wc_w==0.6),:));title('w_1=0.6');xlabel('w_2');ylabel('Q1');
figure
subplot(121)
contourf(wc_w,wc_e,Q1');title('Q1');xlabel('w_1');ylabel('w_2');
subplot(122)
contourf(wc_w,wc_e,Q2');title('Q2');xlabel('w_1');ylabel('w_2');

%% dimensional calculation - then transform to nondim
cacac
wc=[0.6 1.2];
bc_w=[.02:.02:1];
bc_e=[.02:.02:1];
potvort=-1;
dc=0.5;
rc=dc./(wc./2).^2;
g=4e-4;% drho=0.04 based on Voet et al, 2015 Fig2 -- g=4e-4
r=6.371*10^6;lat=-8;
unit=2*pi*r/360;
Omega=7.292e-005;%2*pi/(24*3600);
f0=abs(2*Omega*sin(lat/180*pi));

V=nan(2,51,length(bc_w),length(bc_e));
H=nan(2,51,length(bc_w),length(bc_e));
D=nan(2,51,length(bc_w),length(bc_e));
X0=nan(2,51,length(bc_w),length(bc_e));
%% western channel
Q1=nan(length(bc_w),length(bc_e));
Q2=nan(length(bc_w),length(bc_e));
Q1_int=nan(length(bc_w),length(bc_e));
Q2_int=nan(length(bc_w),length(bc_e));
Ber=nan(length(bc_w),length(bc_e));
for im=1:length(bc_w)
    for jm=1:length(bc_e)
        i=1;
d2=1000;d3=rc(i)*f0^2/g;d4=bc_w(im)*1000;d5=10e+6;
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
Q1(im,jm)=X(1)*abs(f0)/g/d2^2;wc1(1)=X(2);
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
% integrate to obtain Q1
trans = cumsum(v.*d.*[0 diff(x0)]');
Q1_int(im,jm) = trans(end)*f0/g/d2/d2;
% non-dimensionalize variables, western channel
V(1,:,im,jm)=v./sqrt(g*d2);
H(1,:,im,jm)=h./d2;
D(1,:,im,jm)=d./d2;
X0(1,:,im,jm)=x0./Ld;
% Bernoulli at b
% from equations 0.4878
Bb=1/2*(-Ld*f0*(1+2*d3*g/f0^2)*(1-cosh(a3*wc1(1)))/sinh(a3*wc1(1))-2*d3*g*b/f0)^2+g*d4+g*d3*b^2;
% from results 0.4878
hb=d4+d3*b.^2;
db=(d2+2*d3*Ld^2)*(sinh((b-b)/Ld)-sinh((b+a)/Ld))/x1+d2+2*d3*Ld^2;
vb=-Ld*f0*(1+2*d3*g/f0^2)*(cosh((b-b)/Ld)-cosh((b+a)/Ld))/x1-2*d3*g*b/f0;
Bb=1/2*vb^2+g*(hb+db);
Ber(im,jm)=Bb;
%% eastern channel -- replace a1 to a1=g*d4-Bb;
i=2;
d2=1000;d3=rc(i)*f0^2/g;d4=bc_e(jm)*1000;d5=10e+6;
Ld=sqrt(g*d2)/f0;
a1=g*d4-Bb;
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
Q2(im,jm)=X(1)*abs(f0)/g/d2^2;wc1(2)=X(2);
end
Q2_b=X(1);w2c_b=X(2);
% now solve the critical section layer depth
x1=sinh(a3*w2c_b);
x2=cosh(a3*w2c_b);
% a+b:
x=w2c_b;
% a-b:
y=Q2_b*sqrt(a6)/(a7*(1-x2)/x1+x);
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
% integrate to obtain Q1
trans = cumsum(v.*d.*[0 diff(x0)]');
Q2_int(im,jm) = trans(end)*f0/g/d2/d2;
% non-dimensionalize variables, western channel
V(2,:,im,jm)=v./sqrt(g*d2);
H(2,:,im,jm)=h./d2;
D(2,:,im,jm)=d./d2;
X0(2,:,im,jm)=x0./Ld;

    end
end
save Gill_parabolic_nd_2channel_b2.mat Q1 Q2 Q1_int Q2_int rc wc bc_w bc_e dc Ber V H D X0
% _b1.mat wc1=wc2=0.6

figure
plot(bc_e,Q1(find(bc_w==0.6),:));title('b_1=0.6');xlabel('b_2');ylabel('Q1');
figure
subplot(121)
contourf(bc_w,bc_e,Q1');title('Q1');xlabel('b_1');ylabel('b_2');
subplot(122)
contourf(bc_w,bc_e,Q2');title('Q2');xlabel('b_1');ylabel('b_2');
