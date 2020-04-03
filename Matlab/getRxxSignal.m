clear all;
clc;

D = 3;
c = 340;
f = 40000;
fs = 103000;
sample = 0:300;
lambda = c/f;
d = lambda/2;
arrayN = 16;
arrayK = -arrayN/2:arrayN/2-1;
sinTheta = lambda * arrayK / (d * arrayN);
cosSignal = 5*cos(2*pi*sample*f/fs);
sinSignal = sin(2*pi*sample*f/fs);

dirD = [1/6*pi, 0*pi, -2/3*pi];

vTheta = zeros(arrayN, 1);
A = zeros(arrayN, D);
Rxx = zeros(arrayN, arrayN);

for k = 1:D
        for i = 1:arrayN
            vTheta(i, 1) = exp(-1j*(i-1)*(2*pi*d/lambda)*sin(dirD(1, k)));
        end
        A(:,k) = vTheta(:,1);
end

SSigma = zeros(3,1);

snapShot = 1;
sample = 100;
for snap = 1:snapShot
    %snap = 1;
    sigma2 = 0.01; % Noise variance
    n = sqrt(sigma2)*(randn(arrayN,sample) + 1j*randn(arrayN,sample))/sqrt(2);
    
    S = [cosSignal(1:sample); cosSignal(1:sample); cosSignal(1:sample)];
    %S = [cosSignal(snap)*exp(-j*2*pi*f/fs*snap); cosSignal(snap)*exp(-j*2*pi*f/fs*snap); cosSignal(snap)*exp(-j*2*pi*f/fs*snap)];
    S = randn(3,sample);
    SSigma = SSigma + S;
    X = A * S + n;
    Rxx = Rxx + X * (X)';
end
SSigma = SSigma / snapShot;
Rxx = Rxx / snapShot;
[E, Dim] = eig(Rxx);

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