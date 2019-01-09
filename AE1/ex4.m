%% Carregando o som
clear, close, clc
[y,Fs] = audioread('num02.wav');
 
f1 = 941;
f2 = 1336;

fa = 4000;
ta = 1/fa;
fim = 2;
t= 0:ta*2:fim-1;

w1 = sin(2*pi*f1*t);
w2 = sin(2*pi*f1*t);

wtot = w1 + w2;

%% Reproduzindo o som 
%sound(y,Fs)
%% Visualizando o som no DT
time = [0:length(y)-1]'/Fs;
%figure(1),plot(t',wtot'); xlabel('segundos');
figure(1),plot(time',y'); xlabel('segundos');
xlim([0 time(end)]), ylim([-1 1]);
 
%%filtro
%o filtro

p1 = 0.8*exp((800/4000)*j*pi);
p2 = 0.9*exp((941.5/4000)*j*pi);
p22 = 0.9*exp((961.5/4000)*j*pi);
p3 = 0.8*exp((1000/4000)*j*pi);

z3 = 0.95*exp((1036.5/4000)*j*pi );
z4 = 0.95*exp((600/4000)*j*pi );
z5 = 0.95*exp((1200/4000)*j*pi );

Z = [z3 -z3 z4 z4' z5 z5']';
P = [p1 p1' p2 p2' p22 p22' p3 p3']';

[num,den] = zp2tf(Z,P,1);
[h,w] = freqz(num,den);
figure(2); plot(w,abs(h)/max(abs(h)));
figure(3); zplane(num,den);

%% Visualizando o som no DF
Nfreq = length(y);
freq = linspace(0,2*pi,Nfreq)'*Fs/pi/2;
Y = fft(y,Nfreq)/Nfreq;
figure(4),plot(freq,abs(Y)); xlabel('Hertz');
xlim([0 Fs/2]);

[pico,position] = findpeaks(abs(Y));
freq_linha = 2*(position(1)*4000)/16000;
freq_coluna = 2*(position(3)*4000)/16000;
disp(freq_coluna);
disp(freq_linha);

%filtrando uma freq
teste = filter(h/max(h),1,Y);
figure(5),plot(freq,teste);
xlim([0 2000]);
