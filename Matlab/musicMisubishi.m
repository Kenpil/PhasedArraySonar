clear all;
clc;

D = 7;
c = 340;
f = 40000;
fs = 90000;
sample = 0:300;
lambda = c/f;
d = lambda/2;
arrayN = 16;
arrayK = -arrayN/2:arrayN/2-1;
sinTheta = lambda * arrayK / (d * arrayN);
cosSignal = 2.5*cos(2*pi*sample*f/fs);
sinSignal = sin(2*pi*sample*f/fs);
%dirD = [-80, -60, -45, -30, -15, 0, 10, 20, 35, 50, 60, 75]/180*pi;
dirD = [ -60, -45, -15, 0, 10, 35, 50]/180*pi;
%dirD = [-45, 0, 10]/180*pi;

vTheta = zeros(arrayN, 1);
A = zeros(arrayN, D);
Rxx = zeros(arrayN, arrayN);

for k = 1:D
    sin(dirD(1, k));
    for i = 1:arrayN
        vTheta(i, 1) = exp(-1j*(i-1)*(2*pi*d/lambda)*sin(dirD(1, k)));
    end
    A(:,k) = vTheta(:,1);
end

SSH = zeros(D,D);

sample = 4;
sigma2 = 0.01; % Noise variance
n = sqrt(sigma2)*(randn(arrayN,sample) + 1j*randn(arrayN,sample))/sqrt(2);

a = cosSignal(1,1:sample);
b = 2;
c = 5;
a1 = 1;
a2 = sin(2*pi/(2*sample)*(1:sample));
a3 = sin(2*pi/(3*sample)*(1:sample));
a4 = sin(2*pi/(4*sample)*(1:sample));
a5 = sin(2*pi/(5*sample)*(1:sample));
a6 = sin(2*pi/(6*sample)*(1:sample));
a7 = cos(2*pi/(1*sample)*(1:sample));
a8 = cos(2*pi/(1.5*sample)*(1:sample));
a9 = cos(2*pi/(2.0*sample)*(1:sample));
a10 = cos(2*pi/(2.5*sample)*(1:sample));
a11 = cos(2*pi/(3.0*sample)*(1:sample));
a12 = cos(2*pi/(3.5*sample)*(1:sample));
S1 = 5 * a1.*exp(-j*2*pi*f/fs*(1:sample));
S2 = 2 * a2.*exp(-j*2*pi*f/fs*(1+b:sample+b));
S3 = 1 * a3.*exp(-j*2*pi*f/fs*(1+c:sample+c));
S4 = 1 * a4.*exp(-j*2*pi*f/fs*(1+b+c:sample+b+c));
S5 = 1 * a5.*exp(-j*2*pi*f/fs*(1-b+c:sample-b+c));
S6 = 1 * a6.*exp(-j*2*pi*f/fs*(1+2*c:sample+2*c));
S7 = 1 * a7.*exp(-j*2*pi*f/fs*(1+2*b:sample+2*b));
S8 = 1 * a8.*exp(-j*2*pi*f/fs*(1+3*c:sample+3*c));
S9 = 1 * a9.*exp(-j*2*pi*f/fs*(1-c:sample-c));
S10 = 1 * a10.*exp(-j*2*pi*f/fs*(1-2*b:sample-2*b));
S11 = 1 * a11.*exp(-j*2*pi*f/fs*(1+b-3*c:sample+b-3*c));
S12 = 1 * a12.*exp(-j*2*pi*f/fs*(1-2*b+c:sample-2*b+c));
%S = [real(S1); real(S2); real(S3); real(S4); real(S5); real(S6); real(S7); real(S8); real(S9); real(S10); real(S11); real(S12)];
S = [real(S1); real(S2); real(S3); real(S4); real(S5); real(S6); real(S7)];
%S = [real(S1); real(S2); real(S3)];
SSH = S * (S');
rankSSH = rank(SSH);
X = A * S;% + n;


% % ~Normal~
% Rxx = Rxx + X * (X)';
% Rxx = Rxx / sample;
% [E, Dim] = eig(Rxx);
% 
% theta=-90:0.5:90; %Peak search
% Pmusic = zeros(1, length(theta));
% EUnLeft = E(:,1:arrayN-D);
% for ii=1:length(theta)
%     SS=zeros(1,arrayN);
%     for jj=0:arrayN-1
%         SS(1+jj)=exp(j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
%     end
%     PP=SS*EUnLeft*EUnLeft'*SS';
%     Pmusic(ii)=abs(SS*SS'/ PP);
% end
% Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
% plot(theta,Pmusic,'-k')



 
% ~Spatial Smoothing~
arrayNSub = round((arrayN+1)/(1+0.5/sample));
%arrayNSub = 15;
subArrayN = arrayN - arrayNSub + 1;
RxxSS = zeros(arrayNSub, arrayNSub);
XbSS = zeros(arrayNSub, sample);
RxxbSS = zeros(arrayNSub, arrayNSub);
RxxSS2 = zeros(arrayNSub, arrayNSub);

for i = 1:subArrayN
    %RxxSS2 = X * (X)'
    RxxSS = RxxSS + X(i:arrayNSub+i-1,:) * X(i:arrayNSub+i-1,:)';
    XbSS = conj(flip(X(i:arrayNSub+i-1,:),1));
    RxxbSS = RxxbSS + power(XbSS * XbSS',1);
end

RxxSS = RxxSS / subArrayN;
RxxbSS = RxxbSS / subArrayN;
[ESS, DimSS] = eig((RxxSS+RxxbSS) / 2);
RxxSSCompare = imag(RxxSS - RxxSS2);
[ESS2, DimSS2] = eig(RxxSS2);
[ESS, DimSS] = eig(RxxSS);
DimSS = flip(DimSS, 1);
DimSS = flip(DimSS, 2);
ESS = flip(ESS, 1);
ESS = flip(ESS, 2);
ESSCompare = abs(ESS) - abs(ESS2);
DimSSCompare = DimSS - DimSS2;


theta=-90:0.5:90; %Peak search
Pmusic = zeros(1, length(theta));
ESSUnLeft = ESS(:,1:arrayNSub-D);

for ii=1:length(theta)
    SS=zeros(1,arrayNSub);
    for jj=0:arrayNSub-1
        SS(1+jj)=exp(1i*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
    end
    PP=SS*ESSUnLeft*ESSUnLeft'*SS';
    Pmusic(ii)=abs(SS*SS'/ PP);
end

Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
plot(theta,Pmusic,'-k')