clear all;
clc;

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
            %result(point+j, i) = result(point+j, i) + sin(wavePhase3);
        end
    end
end

writematrix(result, 'resultList.txt');
fclose(fileID);

snapShot = 8;
Rxx = zeros(arrayNum, arrayNum, length(result));
X = zeros(arrayNum, 1);
XH = zeros(1, arrayNum);
eigLambda = zeros(arrayNum, arrayNum, length(result));
eigETheta = zeros(arrayNum, arrayNum, length(result));

N = length(theta);
D = zeros(arrayNum, arrayNum);
kNum = -arrayNum/2:arrayNum/2-1;
sinThetaK = lambda*kNum/(d*arrayNum);
msc=zeros(length(result),arrayNum);
absmsc=zeros(length(result),arrayNum);
for m = 1:length(sinThetaK)
    for i = 1:arrayNum
        D(m,i) = exp(1j*(2*pi*(i-1)*d*sinThetaK(m)/lambda));
    end
end

for l = 1:length(result) - snapShot
    X(:,1) = 0;
    for i = 1:snapShot
        tmpa = expList(l+i-1,1);
        tmpb = 1i*expList(l+i-1,2);
        tmpc = [tmpa, tmpb];
        X(:,1) = X(:,1) + (result(l+i-1,:).')*(tmpa+tmpb);
    end
    XH = X';
    Rxx(:,:,l) = X*(XH);
    %[eigETheta(:,:,l), eigLambda(:,:,l)] = eig(Rxx(:,:,l));
    
    k = 2;
    [u,s,v]=svd(Rxx(:,:,l));
    E = v(:,1+k:end);
    
    for i=1:arrayNum
        atheta = D(i,:).';
        msc(l,i) = ((atheta')*atheta)/((atheta')*E*(E')*atheta);
        absmsc(l,:) = abs(msc(l,:));
    end
end