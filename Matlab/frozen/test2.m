clear all;clc;close all
theta=[-90:0.1:90];
N=length(theta);
%Fixed Limitations
d=0.1
%lambda=2*pi/max(w);
%d=c*pi/max(w);
c=1;
M=50;
wmax=3*c*pi/d;

w=2;%linspace(1,wmax,100);
W=length(w);
y=zeros(N,W);
z=y;

wk=2;
%s=[0.65*exp(j*wk)  exp(j*wk)] % ?????????
s=[0.65  1] % ?????????
k=2;
doa=[60 -40];

%Limitations
ws=wk*d*sin(pi*doa/180)/c;

if(abs(max(ws))>pi | abs(max(doa))>90)
    %break
    warning('Aliasing !!')
end

%Model
% ws = kdsin(theta)
a=[exp(-j*ws(1)*[0:M-1]') exp(-j*ws(2)*[0:M-1]')] % steering vector?
%x=[s(:,1).*a(:,1) s(:,2).*a(:,2)];
x=a*s' + 0.1*randn(M,1)

%Beamformer response
bf=zeros(N,1);
msc=zeros(N,1);


[trasn,R]=corrmtx(x,length(x)-1,'modified');
%R = x*(x.');
[u,s,v]=svd(R);
E=v(:,1+k:end);
%R2 = x*(x');
[eigETeta, eigLambda] = eig(R);


for i=1:N

    %Steer Beam throught theta
    b=exp(-j*(wk*d*sin(pi*theta(i)/180)/c)*[0:M-1]);
    bf(i)=b*x/M;
    b=b';
    msc(i)=abs(1/(b'*E*E'*b));

end

figure
plot(theta,log10(abs(bf)))
hold on
plot(theta,log10(msc),'r')
%line([-doa(1) -doa(1)],[0 1.2])
%line([-doa(2) -doa(2)],[0 max(msc)])