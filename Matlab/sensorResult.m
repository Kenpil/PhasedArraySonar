clear all;
clc;

comp = 0+0i;
re = 0;
im = 0i;
f = 40000;
fs = 90000;
resultLength = 901;
point = round(resultLength/10);
arrayNum = 8;
arraivalAngle = -1/3*pi;
v = 340;
lambda = v/f;
d = lambda/2;
fileID = fopen('resultList.txt','w');
expFileID = fopen('phaseList.txt', 'r');
expList = importdata('phaseList.txt');

result = zeros(round(resultLength/4), arrayNum);

wave = 0;
wavePhase = 0;
repeat = 36;
phaseList = zeros(repeat+1, arrayNum);
for j = 0:repeat
    for i = 1:arrayNum
        wavePhase = 2*pi*(f*j/fs + (i-1)*d*sin(arraivalAngle)/lambda);
        if wavePhase > 0 && wavePhase < 16.1 * pi
            result(point+j, i) = sin(wavePhase);
             phaseList(j+1, i) = wavePhase/(2*pi);
        end
    end
end

writematrix(result, 'resultList.txt');
fclose(fileID);

kNum = -arrayNum/2:arrayNum/2-1;
sinthetak = lambda.*kNum/(d*arrayNum);
D = zeros(length(sinthetak),arrayNum);
for m = 1:arrayNum
    for i = 1:length(sinthetak)
        D(m,i) = exp(1i*2*pi*i*d*sinthetak(m)/lambda);
    end
end

S = zeros(length(result),arrayNum);
for i = 1:length(result)
    for j = 1:arrayNum
        S(i,j) = result(i,j)*(expList(i,1)+1i*expList(i,2));
    end
end

B = zeros(length(S), arrayNum);
B = S * (D.');
abs(B)