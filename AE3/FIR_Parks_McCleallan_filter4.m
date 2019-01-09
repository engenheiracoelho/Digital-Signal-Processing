% BP - Filter 4 - Parks-McCleallan
% Jessica & Leticia

clc;close all;clear all;

f1 = 1329; f2 = 1447; f3 = 1507; f4 = 1625;
fe = [f1 f2 f3 f4];

Ap = 1;
As = 40;
Gp = 0;
fs = 4000;

ds = 10^(-As/20)+0.01105;
dp = (1 - (10^(-Ap/20)))/2;
dev = [ds dp ds];
a = [0 1 0];

[n,fo,ao,w] = firpmord(fe,a,dev,fs);
n = n+24;
b = firpm(n,fo,ao,w);
b = b*(1-dp);
atraso_de_grupo = grpdelay(b)
wvtool(b);
figure(2);stem(b);title('Resposta de filtro digital');
figure(3);zplane(b,1);title('Pólos e zeros'); %1 pois não tem denominador no Z
figure(4);freqz(b,1,8000);

% % % Resposta do filtro % % %
figure(5);
[hz,wz] = freqz(b,1,8000,fs);
plot(wz,20*log10(abs(hz))); hold on;title('Filtro digital BP');
Amin = As + 20;
ylim([-Amin Gp+10]);
plot([0 fe(1) fe(1) fe(4) fe(4) fs/2],[-As -As Gp Gp -As -As], ':r');
plot([fe(2) fe(2) fe(2) fe(3) fe(3) fs/2],[-As*100 -As*100 -Ap -Ap -As*100 -As*100], ':r');