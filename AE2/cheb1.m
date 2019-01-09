clear all;clc;close all;
%HP - (f1 = 1075 Hz; f2 = 1209 Hz, As = 35 dB, Chebychev 1)

%% ESPECIFICAÇÃO
fs = 1075; %fr. de passagem
fp = 1209; %fr. de rejeicao
fa = 4000; %fr. de amostragem
Ap = 1; %Atenuação em fp
As = 35; %Atenuação em fs
G0 = 0; %Ganho em db na banda de passagem
G_db = 10^(G0/20);
%% NORMALIZAÇÃO
Wp = (2*pi*fp)/fa; %fr. angular
Ws = (2*pi*fs)/fa;
lp = 2*fa*tan((Wp/2));
ls = 2*fa*tan((Ws/2));
ohmS = lp/ls;
ohmP = 1;

%% EQUAÇÕES(ordem,polos e zeros)

n = ceil((acosh(sqrt(((10^(0.1*As)-1)/(10^(0.1*Ap)-1)))))/(acosh(ohmS))); %ordem do filtro
Eps = sqrt((10^(0.1*Ap)) -1);
Angle = Eps^((-1)/n);
Fi_2 = (1/n)*(asinh(1/Eps));

%%
%Fi_1 = (((2*k) - 1)*pi)/(2*n);

for (k = 1:n)
    T_k = ((2*k-1)*pi)/(2*n);
    p_k(k) = -sinh(Fi_2)*sin(T_k)+(cosh(Fi_2)*cos(T_k))*1i;
end

den = real(poly(p_k)); %denominador
a_p = den(end);

if(mod(n,2)==0)% Ho = 1 (ímpar) e Ho = sqrt(1/(1+(Eps^2)))
    num = sqrt(1/(1+(Eps*Eps)));
else
    num = 1;
end

%% SIMBOLICO
syms p; syms s;
numS=num;
denS=poly2sym(den/a_p,p);
%% FUNÇÃO TRANSFERÊNCIA EM P
Hp(p)=symfun(numS/denS,p)*G_db;

%visualizando o plano z de H(p)
figure(1);
zplane(num,den);
title('Polos e zeros do Prototipo H(p)');

%% FILTRO COM GANHO NORMALIZADO
wr_p=0:0.01:ohmS+1;
calc_p = Hp(j*wr_p);
figure(2);
plot([ohmP ohmS],[-Ap -As]+G0,'r*');
hold on;
plot(wr_p,20*log10(abs(calc_p)));
plot([ohmP ohmP 0],[-Ap*100 -Ap -Ap]+G0,'--m'); % band pass mask
plot([0 ohmS ohmS 10*ohmS],[0 0 -As -As]+G0,'--m'); % band stop mask
axis([0.5 1.5 -40 5]);
xlabel('Omega p normalizado');
ylabel('Magnitude (dB)');
title('Magnitude do Prototipo H(p)');
hold off;

figure(3)
[b,a] = freqs(double(num),double(den));
subplot(211);
semilogx(a,20*log10(abs(b)));
xlim([0 10]);
grid on;
title('Resposta em frequência do Protótipo H(p)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
semilogx(a,angle(b));
grid on;
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%% PROTOPIPO ANALOGICO

% Funcao de transferencia em s
syms s;
Hs(s)=Hp(lp/s); % substitui s em p
wr_s=logspace(0,log10(ls+lp),1000);
[num_s_aux den_s_aux] = numden(Hs(s));
num_s = sym2poly(num_s_aux);
den_s = sym2poly(den_s_aux);

calc_s = Hs(1i.*wr_s);
%%
tf_H_s = Hs;
figure(4);
hold on;
plot([ls lp],[-As -Ap]+G0,'r*'); %%%VERIFICAR QUAL É O VALOR CORRETO
plot([0 ls ls max(wr_s)],[-As -As Ap Ap],'--m');
plot([lp lp max(wr_s)],[-As*2 -Ap -Ap],'--m');
plot(wr_s,20*log10(abs(calc_s)));
xlabel('Frequências do Filtro(H(s))');
ylabel('Magnitude (dB)');
axis([4000 14000 -45 5]);
title('Magnitude do Filtro Analogico H(s)');
hold off;

figure(5)
[b1,a1] = freqs(double(num_s),double(den_s));
subplot(211);
semilogx(a1,20*log10(abs(b1)));
ylim([-50 5]);
%xlim([0 10000]);
grid on;
title('Resposta em frequência do Filtro Analógico H(s)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
semilogx(a1,angle(b1));
grid on;    
xlim([10e2 10e5]);
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%% VISUALIZA O PLANO Z de H(s)

figure(6);
zplane(num_s,den_s);
title('Pólos e zeros do Filtro Analógico H(s)');

%% PASSANDO O FILTRO PARA DIGITAL
% Funcao de transferencia em z usando trans. bilinear
syms z;
s = 2*fa*(z-1)/(z+1);
Hz(z) = collect(subs(Hs,s));
tf_H_z = vpa(Hz);

[Nz,Dz] = numden(Hz(z));
bz = sym2poly(Nz);
az = sym2poly(Dz);

wr_z = 0.01:0.01:pi;
calc_z = Hz(exp(1i*wr_z));
figure(7);
hold on;
plot([fp fs],[-Ap -As]+G0,'r*');
plot([0 fs fs 2*fs],[-As -As Ap Ap],'--m');
plot([fp fp 2*fp],[-As*2 -Ap -Ap],'--m');
plot(fa*wr_z/pi/2,20*log10(abs(calc_z)));
xlabel('Frequências do Filtro (Hz)');xlim([fp-500 fs+500]);
ylabel('Magnitude (dB)');ylim([-40 5]);
title('Magnitude do Filtro Digital H(z)');
hold off;

figure(8);
zplane(sym2poly(Nz),sym2poly(Dz));
title('Polos e zeros do Filtro Digital H(z)');

[b2,a2] = freqz(bz,az);
figure(9);
subplot(211);
plot(fa*abs(a2)/pi/2,20*log10(abs(b2)));
ylim([-40 5]);
xlim([1000 1600]);
grid on;
title('Resposta em frequência do Filtro Digital H(z)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
plot(fa*abs(a2)/pi/2,unwrap(angle(b2))*180/pi);
xlim([1000 1600]);
grid on;
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%calculando o atraso de grupo
figure(10);
grpdelay(bz,az)
title('Atraso de grupo H(z)');
