clear all;
clc;

D = 3;
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
dirD = [1/6*pi, 0*pi, -2/3*pi];

vTheta = zeros(arrayN, 1);
A = zeros(arrayN, D);
Rxx = zeros(arrayN, arrayN);

for k = 1:D
    sin(dirD(1, k))
        for i = 1:arrayN
            vTheta(i, 1) = exp(-1j*(i-1)*(2*pi*d/lambda)*sin(dirD(1, k)));
        end
    A(:,k) = vTheta(:,1);
end

SSigma = zeros(3,1);
SSH = zeros(D,D);

snapShot = 1;
sample = 100;
for snap = 1:snapShot
    %snap = 1;
    sigma2 = 0.01; % Noise variance
    n = sqrt(sigma2)*(randn(arrayN,sample) + 1j*randn(arrayN,sample))/sqrt(2);
    
    %S = [cosSignal(1,1:sample); cosSignal(1,1:sample); cosSignal(1,1:sample)];
    %S = [cosSignal(snap)*exp(-j*2*pi*f/fs*snap); cosSignal(snap)*exp(-j*2*pi*f/fs*snap); cosSignal(snap)*exp(-j*2*pi*f/fs*snap)];
    %S1 = randn(1,sample);
    %S = [S1; S1+10; S1-3];
    a = cosSignal(1,1:sample);
    b = 2;
    c = 4;
    %S = randn(3, sample);
    %S = [a.*exp(-j*2*pi*f/fs*(1:sample)); a.*exp(-j*2*pi*f/fs*(1+b:sample+b)); a.*exp(-j*2*pi*f/fs*(1+c:sample+c))];
    %S = [cos(2*pi*f/fs*(1:sample)); cos(2*pi*f*1.1/fs*(1+b:sample+b)); cos(2*pi*f*0.9/fs*(1+c:sample+c))];
    a1 = 1;
    a2 = sin(2*pi/(2*sample)*(1:sample));
    a3 = sin(2*pi/(5*sample)*(1:sample));
    S1 = 50 * a1.*exp(-j*2*pi*f/fs*(1:sample));
    S2 = 2 * a2.*exp(-j*2*pi*f/fs*(1+b:sample+b));
    S3 = 1 * a3.*exp(-j*2*pi*f/fs*(1+c:sample+c));
    S = [S1; S2; S3];
    SSigma = SSigma + S;
    X = A * S;% + n;
    Rxx = Rxx + X * (X)';
end
SSH = S * (S');
rankSSH = rank(SSH);
SSigma = SSigma / snapShot;
Rxx = Rxx / sample;
[E, Dim] = eig(Rxx);

% PMU = zeros(1,3);
% for k = 1:D
%     upTmp = (A(:,k)') * A(:,k);
%     unTmp = (A(:,k)') * E(:,1:arrayN-D) * (E(:,1:arrayN-D)') * A(:,k);
%     PMU(1,k) = upTmp / unTmp;
% end

theta=-90:0.5:90; %Peak search
Pmusic = zeros(1, length(theta));
EUnLeft = E(:,1:arrayN-D);
for ii=1:length(theta)
    SS=zeros(1,length(arrayN));
    for jj=0:arrayN-1
        SS(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
    end
    PP=SS*EUnLeft*EUnLeft'*SS';
    Pmusic(ii)=abs(1/ PP);
end
Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
plot(theta,Pmusic,'-k')



% close all; clear all; clc;
% 
% 
% % ======= (1) TRANSMITTED SIGNALS ======= %
% 
% % Signal source directions
% az = [35;39;127]; % Azimuths
% el = zeros(size(az)); % Simple example: assume elevations zero
% M = length(az); % Number of sources
% 
% % Transmitted signals
% L = 200; % Number of data snapshots recorded by receiver
% m = randn(M,L); % Example: normally distributed random signals
% 
% % ========= (2) RECEIVED SIGNAL ========= %
% 
% % Wavenumber vectors (in units of wavelength/2)
% k = pi*[cosd(az).*cosd(el), sind(az).*cosd(el), sind(el)].';
% 
% % Array geometry [rx,ry,rz]
% N = 10; % Number of antennas
% r = [(-(N-1)/2:(N-1)/2).',zeros(N,2)]; % Assume uniform linear array
% 
% % Matrix of array response vectors
% A = exp(-1j*r*k);
% 
% % Additive noise
% sigma2 = 0.01; % Noise variance
% n = sqrt(sigma2)*(randn(N,L) + 1j*randn(N,L))/sqrt(2);
% 
% m2 = randn(M,1);
% 
% % Received signal
% x = A*m2;% + n;
% 
% % ========= (3) MUSIC ALGORITHM ========= %
% 
% % Sample covariance matrix
% Rxx = x*x';%/L;
% 
% % Eigendecompose
% [E,D] = eig(Rxx);



% %By Honghao Tang
% clc
% clear all
% format long %The data show that as long shaping scientific
% doa=[20 60]/180*pi; %Direction of arrival
% N=200;%Snapshots
% w=[pi/4 pi/3]';%Frequency
% M=10;%Number of array elements
% P=length(w); %The number of signal
% lambda=150;%Wavelength
% d=lambda/2;%Element spacing
% snr=20;%SNA
% D=zeros(P,M); %To creat a matrix with P row and M column
% for k=1:P
%     D(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]); %Assignment matrix
% end
% D=D';
% xx=2*exp(j*(w*[1:N])); %Simulate signal
% x=D*xx;
% %x=x+awgn(x,snr);%Insert Gaussian white noise
% R=x*x'; %Data covarivance matrix
% [N,V]=eig(R); %Find the eigenvalues and eigenvectors of R
% NN=N(:,1:M-P); %Estimate noise subspace
% theta=-90:0.5:90; %Peak search
% for ii=1:length(theta)
%     SS=zeros(1,length(M));
%     for jj=0:M-1
%         SS(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
%     end
%     PP=SS*NN*NN'*SS';
%     Pmusic(ii)=abs(1/ PP);
% end
% Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
% plot(theta,Pmusic,'-k')
% xlabel('angle \theta/degree')
% ylabel('spectrum function P(\theta) /dB')
% title('DOA estimation based on MUSIC algorithm ')
% grid on