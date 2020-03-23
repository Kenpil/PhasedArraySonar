clear all;
clc;

comp = 0+0i;
re = 0;
im = 0i;
f = 40000;
fs = 90000;
resultLength = 1001;
point = round(resultLength/40);
arrayNum = 16;
arraivalAngle = -1/6*pi;
arraivalAngle2 = -2/6*pi;
arraivalAngle3 = -0/6*pi;
v = 340;
lambda = v/f;
d = lambda/2;
fileID = fopen('resultList.txt','w');
expFileID = fopen('phaseList.txt', 'r');
expList = importdata('phaseList.txt');
theta=[-90:1:90];

result = zeros(round(resultLength/15), arrayNum);
for i = 1:length(result)
    for j = 1:arrayNum
        result(i, j) = -0.1+0.2*rand;
    end
end

wave = 0;
wavePhase = 0;
wavePhase2 = 0;
wavePhase3 = 0;
repeat = 36;
phaseList = zeros(repeat+1, arrayNum);
phaseList2 = zeros(repeat+1, arrayNum);
for j = 0:repeat
    for i = 1:arrayNum
        wavePhase = 2*pi*(f*j/fs + (i-1)*d*sin(arraivalAngle)/lambda);
        if wavePhase > 0 && wavePhase < 16.1 * pi
            result(point+j, i) = result(point+j, i) + sin(wavePhase);
            phaseList(j+1, i) = wavePhase/(2*pi);
        end
        wavePhase2 = 2*pi*(f*j/fs + (i-1)*d*sin(arraivalAngle2)/lambda);
        if wavePhase2 > 0 && wavePhase2 < 16.1 * pi
            result(point+j, i) = result(point+j, i) + sin(wavePhase2);
            phaseList2(j+1, i) = wavePhase2/(2*pi);
        end
        wavePhase3 = 2*pi*(f*j/fs + (i-1)*d*sin(arraivalAngle3)/lambda);
        if wavePhase3 > 0 && wavePhase3 < 16.1 * pi
            result(point+j, i) = result(point+j, i) + sin(wavePhase3);
        end
    end
end

writematrix(result, 'resultList.txt');
fclose(fileID);

Rxx = zeros(arrayNum, arrayNum, length(result));
X = zeros(arrayNum, 1);
XH = zeros(1, arrayNum);
eigLambda = zeros(arrayNum, arrayNum, length(result));
eigETheta = zeros(arrayNum, arrayNum, length(result));

N = length(theta);
msc=zeros(N,1,length(result));

for l = 1:length(result)
    %X = fftshift(fft((result(l, :)))).'
    tmpa = expList(l,1);
    tmpb = 1i*expList(l,2);
     tmpc = [tmpa, tmpb];
     X = result(l,:)*(tmpa+tmpb);
    XH = X';
    Rxx(:,:,l) = X*(XH);
    [trasn,Rxx(:,:,l)]=corrmtx(X,length(X)-1,'modified');
    
    Rxx(:,:,l);
    [eigETheta(:,:,l), eigLambda(:,:,l)] = eig(Rxx(:,:,l));
    
%     [u,s,v]=svd(Rxx(:,:,l));
%     E=v(:,10:end);
%     
%     for i=1:N
%         
%         %Steer Beam throught theta
%         b=exp(-j*(d*sin(pi*theta(i)/180)/lambda)*[0:arrayNum-1]);
%         bf(i)=b*(X.')/arrayNum;
%         b=b';
%         msc(i,1,l)=abs(1/(b'*E*E'*b));
%         
%     end
end

% figure;
% plot(theta,log10(msc(:,1,40)),'r')

%eigLambda = zeros()

% kNum = -arrayNum/2:arrayNum/2-1;
% sinthetak = lambda.*kNum/(d*arrayNum);
% D = zeros(length(sinthetak),arrayNum);
% for m = 1:arrayNum
%     for l = 1:length(sinthetak)
%         D(m,l) = exp(1i*2*pi*l*d*sinthetak(m)/lambda);
%     end
% end
% 
% S = zeros(length(result),arrayNum);
% for l = 1:length(result)
%     for j = 1:arrayNum
%         S(l,j) = result(l,j)*(expList(l,1)+1i*expList(l,2));
%     end
% end
% 
% Stmp = zeros(arrayNum, 1);
% B = zeros(length(S), arrayNum);
% pow2Num = 64;
% fftB = zeros(length(S), pow2Num);
% for j = 1:length(result)
%     Stmp = (S(j, :)).';
%     B(j, :) = D*Stmp;
%     fftB(j,:) = (fft(S(j, :),pow2Num));
% end
% %B = S * (D.');
% for i = 1:length(result)
%    %B(i,:) = fftshift(B(i,:));
% end
% absfftB = abs(fftB);
% absB = abs(B);