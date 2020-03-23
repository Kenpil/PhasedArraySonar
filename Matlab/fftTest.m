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
        %result(i, j) = -0.1+0.2*rand;
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
            %result(point+j, i) = result(point+j, i) + sin(wavePhase2);
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

D = zeros(arrayNum, arrayNum);
kNum = -arrayNum/2:arrayNum/2-1;
sinThetaK = lambda*kNum/(d*arrayNum);
checkD = zeros(arrayNum, arrayNum);
for m = 1:length(sinThetaK)
    for i = 1:arrayNum
        checkD(m,i) = 2*i*d*sinThetaK(m)/lambda;
        D(m,i) = exp(1j*(2*pi*i*d*sinThetaK(m)/lambda));
    end
end

B = zeros(length(result), length(sinThetaK));
S = zeros(arrayNum, 1);
absB = zeros(length(result), length(sinThetaK));
for i = 1:length(result)
    S = (result(i,:).')*(expList(i,1) + 1j*expList(i,2));
    B(i,:) = D*S;
    absB(i,:) = abs(B(i,:));
end

% figure;
% plot(kNum,log10(msc(:,1,40)),'r')