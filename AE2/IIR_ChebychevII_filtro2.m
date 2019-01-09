% BP - Filter 2 - Chebychev II
% Jessica & Leticia

close all;clear all;

%especificacoes do projeto
fp1 = 1447; %fr. de passagem inferior
fp2 = 1507; %fr. de passagem superior
fs1 = 1329; %fr. de rejeicao inferior
fs2 = 1625; %fr. de rejeicao superior
fa = 4000; %fr. de amostragem
Ap = 1; %Atenuação na fr. de passagem
As = 40; %Atenuação na fr. de rejeicao
G0 = 0;

%normalizando as fr.
wp1 = (2*pi*fp1)/fa;
wp2 = (2*pi*fp2)/fa;
ws1 = (2*pi*fs1)/fa;
ws2 = (2*pi*fs2)/fa;

lp1 = 2*fa*tan((wp1/2));
lp2 = 2*fa*tan((wp2/2));
ls1 = 2*fa*tan((ws1/2));
ls2 = 2*fa*tan((ws2/2));

%equacoes iniciais:determinacao da ordem,polos e zeros
B = lp2 - lp1;
l0 = sqrt(lp2*lp1);

ohmp = 1;
ohms = abs((-(ls2^2)+(l0^2))/(B*ls2));  % ohmS = |(-ws² + wo²)/(B*ws)|

%calculo do chebyshev tipo 2 usando as funcoes simbolicas
[N, Wn] = cheb2ord(ohmp,ohms,Ap,As,'s');
[b,a] = cheby2(N,As,Wn,'s');

syms p;
bp = poly2sym(b,'p'); %num
ap = poly2sym(a,'p'); %den
Hp(p) = vpa(bp./ap);

%visualizando o plano z de H(p)
figure(1);
zplane(b,a);
title('Pólos e zeros do Protótipo H(p)');

%visualizando o filtro com ganho normalizado
wr_p = 0:0.01:ohms+1; 
calc_p = Hp(j*wr_p);

figure(2);
hold on;
plot(wr_p,20*log10(abs(calc_p)));
plot([ohmp ohms],[-Ap -As],'r+');
plot([ohmp ohmp 0],[-Ap*100 -Ap -Ap],'--m'); % máscara da banda de passagem
plot([0 ohms ohms 10*ohms],[0 0 -As -As],'--m'); % máscara da banda de rejeição
hold off;
axis([0 6 -85 5]);
title('Magnitude de H(p) calculada');
xlabel('ohm_p');
ylabel('Magnitude (dB)');

figure(3);
[b,a] = freqs(b,a);
subplot(211);
semilogx(a,20*log10(abs(b)));
grid on;
title('Resposta em frequência do Protótipo H(p)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
semilogx(a,angle(b));
grid on;
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%passando o prototipo para analogico
%funcao de transferencia em s
syms s;
p = (((s^2) + (l0^2))./(B*s));
Hs(s) = vpa(collect(subs(Hp,p)));
s_calc =logspace(0,log10(ls2+lp2),1000);
Hs_calc = Hs(j*s_calc);

figure(4);
plot(s_calc,20*log10(abs(Hs_calc)));
hold on;
%NEW MASK
plot([lp1 ls1 l0 ls2 lp2],[-Ap -As 0 -As -Ap],'r+');
plot([100 ls1 ls1 ls2 ls2 10e5],[-As -As -Ap+1 -Ap+1 -As -As],'--m'); % máscara da banda de passagem
plot([100 lp1 lp1 lp2 lp2 10e5],[-2*As -2*As -Ap -Ap -2*As -2*As],'--m'); % máscara da banda de rejeição
title('Magnitude do Filtro Analógico H(s)');
xlabel('Frequências do Filtro em Lamda (Hz)');
ylabel('Magnitude (dB)');
ylim ([-85 5]); xlim([5e3 4e4]);
hold off;

%visualizando o plano z de H(s)
[num_s, den_s] = numden(sym(Hs(s)));
ns = sym2poly(num_s);
ds = sym2poly(den_s);
figure(5);
zplane(ns,ds);
title('Pólos e zeros do Filtro Analógico H(s)');

figure(6);
[b2,a3] = freqs(ns,ds);
subplot(211);
plot(fa*abs(a3)/pi/2,20*log10(abs(b2)));
ylim ([-85 5]); xlim([5e3 2e7]);
grid on;
title('Resposta em frequência do Filtro Analógico H(s)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
plot(fa*abs(a3)/pi/2,unwrap(angle(b2))*180/pi);
ylim ([-240 240]); xlim([5e3 2e7]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid on;
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%passando o filtro para digital
%funcao de transferencia em z usando trans. bilinear
syms z;
num_sz = 2*fa*(z-1); 
den_sz = (z+1);
s = num_sz./den_sz;
Hz(z) = vpa(collect(subs(Hs,s)));

[num_z den_z] = numden(sym(Hz(z)));
nz = sym2poly(num_z);
dz = sym2poly(den_z); 

z_calc = 0.1:0.01:pi;
Hz_calc = Hz(exp((j*z_calc))); 
%%
figure(7);
plot(((z_calc/pi)*fa/2),20*log10(abs(Hz_calc)));
hold on;
%NEW MASK
plot([fp1 fs1 fp2-30 fs2 fp2],[-Ap -As 0 -As -Ap],'r+');
plot([100 fs1 fs1 fs2 fs2 fs2 fa/2],[-As -As -Ap+1 -Ap+1 -As -As -As],'--m'); % máscara da banda de passagem
plot([fs1/pi fp1 fp1 fp2 fp2 fs2+500],[-2*As -2*As -Ap -Ap -2*As -2*As],'--m'); % máscara da banda de rejeição
title('Magnitude de H(z) calculada');
xlabel('Frequências projetadas (kHz)');
ylabel('Magnitude (dB)');
ylim ([-85 5]);xlim([500 2000]);
hold off;
%%
figure(8);
zplane(nz,dz);
title('Pólos e zeros do Filtro Digital H(z)');

figure(9);
freqz(nz,dz);
title('Magnitude x Fase de H(z)');

[b3,a3] = freqz(sym2poly(num_z),sym2poly(den_z));
figure(9);
subplot(211);
plot(fa*abs(a3)/pi/2,20*log10(abs(b3)));
ylim ([-85 5]);
grid on;
title('Resposta em frequência do Filtro Digital H(z)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
plot(fa*abs(a3)/pi/2,unwrap(angle(b3))*180/pi);
ylim ([-240 240]);
grid on;
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%calculando o atraso de grupo
figure(10);
grpdelay(nz,dz);
title('Atraso de grupo');

%equacoes
pretty(Hp)
pretty(Hs)
pretty(Hz)