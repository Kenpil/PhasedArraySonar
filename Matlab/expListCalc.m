clear all;
clc;

comp = 0+0i;
re = 0;
im = 0i;
f = 40000;
fs = 90000;
fileID = fopen('phaseList.txt','w');
for i = 0:1:900
    comp = exp(2*pi*f*i/fs*1i);%2(pi)ft, fs90kHz,0~t~0.01s
    re = real(comp);
    im = imag(comp);
    fprintf(fileID,'%f %1f\n',re, im);
end
fclose(fileID);