%% This code compute the Gill's model solution for the nondimensional, rectangular cross-section
% 
% INPUTS:
% q=-[.01 .1 .5 1 5];
% w1c=[.01 .05 .1:.1:10];
% Q=[.01 .05 .1:.1:10];
% h1c=[.01 .05 .1:.1:10];
%
% OUTPUTS:
% dc and Q1
%
% CALLS:
%
% TO NOTICE:
% No Bernoulli jump
% 
% S.Tan, WHOI, 2019/06/07
% add separation condition
% S.Tan, IOCAS, 2020/06/15
caca
q=-[.01 .1:.1:.5 1:5];
w1c=[.01 .05 .1:.1:10];
Q=[.5 1 5 10];
h1c=[.01 .05 .1:.1:10];
% FIND THE BEST ESTIMATE
r=6.371*10^6;
unit=2*pi*r/360;
Omega=7.292e-005;%2*pi/(24*3600);
f0=2*Omega*sin(-8/180*pi);%8S
g=4e-4;% drho=0.04 based on Voet et al, 2015 Fig2 -- g=4e-4
D=1000;D_inf=1000;Ld=sqrt(g*D)/abs(f0);
best=[-D/D_inf 50e+3/Ld 10e+6*abs(f0)/g/D^2 600/D];%best estimate D_inf=D=1
d1c=nan(length(q),length(w1c),length(Q),length(h1c));
d1f=nan(length(q),length(w1c),length(Q),length(h1c));
v1c=nan(length(q),length(w1c),length(Q),length(h1c));
v1f=nan(length(q),length(w1c),length(Q),length(h1c));
Q1=nan(length(q),length(w1c),length(Q),length(h1c));
tic;
for i=1:length(q)
    for j=1:length(w1c)
        for k=1:length(Q)
            for m=1:length(h1c)
    T1c=tanh(1/2*sqrt(abs(q(i)))*w1c(j));
    a1=T1c^2;
    a2=1.5*(abs(q(i))^(-1))*(1-T1c^2);
    a4=abs(q(i))^(-2)+Q(k)-(abs(q(i))^(-1))*h1c(m)-T1c^2*(abs(q(i))^(-2))/2;
    b1=(1-T1c^2)*(abs(q(i))^(-1));
    b2=T1c^2;
    b3=1/(4*T1c^2);
    d0=1;
    fun = @(d) a1*d^2+a2*d+sqrt((b1+b2*d)*(d^3)/b3)/2-a4;
    [d, resnorm,residual, exitflag] = lsqnonlin(fun,d0);
    d1c(i,j,k,m)=d;
    Q1(i,j,k,m)=sqrt((b1+b2*d)*(d^3)/b3);
    d1f(i,j,k,m)=-Q1(i,j,k,m)./2./d1c(i,j,k,m);
    v1c(i,j,k,m)=-sqrt(abs(q(i)))*d1f(i,j,k,m)/T1c;
    v1f(i,j,k,m)=-sqrt(abs(q(i)))*T1c*(d1c(i,j,k,m)-1/abs(q(i)));
            end
        end
    end
    disp(i)
end
toc;
save Gill_rectangular.mat q w1c Q h1c d1c Q1 best d1f v1c v1f
%% separation condition
d_bar=nan(length(q),length(w1c),length(Q),length(h1c));
for i=1:length(q)
    for j=1:length(w1c)
        for k=1:length(Q)
            for m=1:length(h1c)
                B_L=abs(q(i))*Q(k)+abs(q(i))^(-1);
%                 a1=4*abs(q(i))^(-1);
%                 a2=2*abs(q(i))^(-2)*(B_L-h1c(m));
%                 a3=-2*abs(q(i))^(-1)*(B_L-h1c(m))-4*abs(q(i))^(-2);
%                 x0=1;
%                 fun = @(x) a1*x+a2*x^(-1)+a3;
%                 [x, resnorm,residual, exitflag] lsqnonlin(fun,x0);
%                 if x>0
%                     d_bar(i,j,k,m)=x;
%                 end
                A=1;
                B=-(abs(q(i))^(-1)+(B_L-h1c(m))/2);
                C=(abs(q(i))^(-1)*(B_L-h1c(m)))/2;
                d_bar_=[(-B+sqrt(B^2-4*A*C))/2/A (-B-sqrt(B^2-4*A*C))/2/A];
                l=find(d_bar_>0&d_bar_<=abs(q(i))^(-1)/2);
                if ~isempty(l)
                    d_bar(i,j,k,m)=d_bar_(l);
                end
            end
        end
    end
end
l=find(q==-1);
d_bar_=squeeze(d_bar(l,:,1,:));
Te_=sqrt(-1./(1-abs(q(l))^(-1)./d_bar_));
Te_(find(abs(imag(Te_))>0))=nan;
we_=atanh(Te_).*2./sqrt(abs(q(l)));
we_(find(abs(imag(we_))>0))=nan;
% check the Gill's function
for j=1:length(w1c)
    for m=1:length(h1c)
    x0(j,m)=h1c(m)-abs(q(l))*Q(1)-abs(q(l))^(-1); 
    end
end
G=Te_.^(-2).*d_bar_.^2+Te_.^2.*(d_bar_-abs(q(l))^(-1)).^2+2.*d_bar_.^2+2.*abs(q(l))^(-1).*d_bar_+2.*abs(q(l))^(-1).*x0;
% compute Q1
Q1_=2.*d_bar_.^2;
save Gill_rectangular.mat q w1c Q h1c d1c Q1 best d1f v1c v1f d_bar_ we_ Q1_
%%
Q1_level=[0.05 0.1:0.1:0.4 0.45];
d_bar_level=sqrt(Q1_level./2);
l=find(q==-1);
Te_level=sqrt(-1./(1-abs(q(l))^(-1)./d_bar_level));
we_level=atanh(Te_level).*2./sqrt(abs(q(l)));
k=1;
B_L=abs(q(l))*Q(k)+abs(q(l))^(-1);
h_level=B_L-2.*(d_bar_level-abs(q(l))^(-1))./(1-abs(q(l))^(-1).*d_bar_level.^(-1));
hold on
plot(we_level,h_level,'ks','markerfacecolor','k')
hold on
for i=1:length(we_level)
    text(we_level(i),h_level(i),num2str(Q1_level(i)),'fontsize',14);
end
Q1_new=squeeze(Q1(l,:,k,:));
for i=1:length(w1c)
    for j=1:length(h1c)
        Q1_pt=Q1_new(i,j);
        d_bar_pt=sqrt(Q1_pt/2);
        Te_pt=sqrt(-1./(1-abs(q(l))^(-1)./d_bar_pt));
        we_pt=atanh(Te_pt).*2./sqrt(abs(q(l)));
        h_pt=B_L-2.*(d_bar_pt-abs(q(l))^(-1))./(1-abs(q(l))^(-1).*d_bar_pt.^(-1));
        if w1c(i)>we_pt
            Q1_new(i,j)=nan;
        end
    end
end
for j=1:length(h1c)
    l=find(~isnan(Q1_new(:,j)));
    if ~isempty(l)
    Q1_new(l(end)+1:end,j)=Q1_new(l(end),j);
    end
end
save Gill_rectangular.mat q w1c Q h1c d1c Q1 best d1f v1c v1f d_bar_ we_ Q1_ Q1_new we_level h_level

for l=1:length(q)
    A = real(squeeze(Q1(l,:,1,:))');
    A(A>Q(1))=nan;
figure
contourf(w1c,h1c,A);
hold on
contour(w1c,h1c,real(squeeze(d1c(l,:,1,:))'+squeeze(d1f(l,:,1,:))'),[0,0],'w');
end
