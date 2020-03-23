clear all;
clc;

% A = [1,2,3,4; 1.5,2.5,3.5,4.5];
% B = zeros(2, 4, 2);
% for i = 1:2
%     A(i,:) = A(i,:) * 2
%     B(:,:,i) = A(:,:);
% end

clear all;
clc;

comp = 0+0i;
re = 0;
im = 0i;
f = 40000;
fs = 90000;
resultLength = 901;
point = round(resultLength/32);
arrayNum = 8;
arraivalAngle = -1/6*pi;
arraivalAngle2 = 2/6*pi;
v = 340;
lambda = v/f;
d = lambda/2;
fileID = fopen('resultList.txt','w');
expFileID = fopen('phaseList.txt', 'r');
expList = importdata('phaseList.txt');

result = zeros(round(resultLength/15), arrayNum);
for i = 1:length(result)
    for j = 1:arrayNum
        result(i, j) = -0.1+0.2*rand;
    end
end

wave = 0;
wavePhase = 0;
wavePhase2 = 0;
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
        if wavePhase2 > 0 && wavePhase2 < 8.1 * pi
            result(point+j, i) = result(point+j, i) + sin(wavePhase2);
            %phaseList2(j+1, i) = wavePhase2/(2*pi);
        end
    end
end

writematrix(result, 'resultList.txt');
fclose(fileID);

[S,wo] = pmusic(result(40, :),2);
