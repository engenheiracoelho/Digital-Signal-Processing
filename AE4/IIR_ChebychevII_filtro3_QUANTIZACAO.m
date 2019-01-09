% BP - Filter 2 - Chebychev II
% Jessica & Leticia
 
close all;clear all;

%especificacoes do projeto
fp1 = 1185; %fr. de passagem inferior
fp2 = 1233; %fr. de passagem superior
fs1 = 1088; %fr. de rejeicao inferior
fs2 = 1330; %fr. de rejeicao superior
fa = 4000; %fr. de amostragem

Ap = 1; %Atenuação na fr. de passagem
As = 40; %Atenuação na fr. de rejeicao
G0 = 0;
Amin = As + 20;
fs = fa*2;

f1 = fp1; f2 = fs1; f3 = fp2; f4 = fs2;
fe = [f1 f2 f3 f4];

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
Hp(p) = vpa((1*bp)./ap);


wr_p = 0:0.01:ohms+1; 
calc_p = Hp(j*wr_p);


%passando o prototipo para analogico
%funcao de transferencia em s
syms s;
p = (((s^2) + (l0^2))./(B*s));
Hs(s) = vpa(collect(subs(Hp,p)));
s_calc =logspace(0,log10(ls2+lp2),1000);
Hs_calc = Hs(j*s_calc);

%%
%visualizando o plano z de H(s)
[num_s, den_s] = numden(sym(Hs(s)));
ns = sym2poly(num_s);
ds = sym2poly(den_s);


%figure(6);
[b2,a3] = freqs(ns,ds);


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

dz2 = dz/dz(1);
nz2 = nz/dz(1);

%z_calc = 0.1:0.0001:pi;
%Hz_calc = Hz(exp((j*z_calc))); 
[Hz_calc,z_calc] = freqz(nz,dz,40000);

[Hz_calc2,z_calc2] = freqz(nz2,dz2,20000);
Ap2 = 1.2; %Atenuação na fr. de passagem

%%
f_mask = [0 fe(2) fe(2) fe(4) fe(4) fa fa fe(1) fe(1) fe(3) fe(3)]/(fa/2);
a_mask = [-As -As Ap Ap -As -As -Amin -Amin -Ap -Ap -Amin];

%%
% Num24bsemqtz =  1.0e+125 * [2.4960    3.1438    3.2696         0   -3.2696   -3.1438   -2.4960];
% Den24bsemqtz = 1.0e+128 * [0.6537    1.2269    2.6092    2.4663    2.4544    1.0855    0.5441];
% 
% Num24bsemsos = [0.0039    0.0050    0.0051         0   -0.0051   -0.0050   -0.0039];
% Den24bsemsos = [1.0000    1.8748    3.9825    3.7606    3.7386    1.6519    0.8271];
% 
% SOS24b = [ 1.0000  0  -1.0000  1.0000 0.6148 0.9055;1.0000 0.3247 1.0000 1.0000 0.5527 0.9551;1.0000    0.9328    1.0000    1.0000    0.7072    0.9563];
% G24b = [0.0039;1.0000;1.0000;1.0000];

%funcao sem quant sem sos
[hsq,wsq] = freqz(Num24bsemqtz,Den24bsemqtz,40000);

%funcao com quant sem sos
[hquant,wquant] = freqz(Num24bsemsos,Den24bsemsos,40000);

%funcao com quant com sos
Hd = dfilt.df2sos(SOS24b,G24b);
[num_sos den_sos] = tf(Hd);
[hsos,wsos] = freqz(num_sos,den_sos,40000);
Ap = 0;
figure(7);
plot(((z_calc/pi)*fa/2),20*log10(abs(Hz_calc)));
hold on;
%plot(fa*abs(wsq)/pi/2,20*log10(abs(hsq)));
plot(fa*abs(wquant)/pi/2,20*log10(abs(hquant)));
plot(fa*abs(wsos)/pi/2,20*log10(abs(hsos)));
plot([fs1 fp1 (fp1+fp2)/2 fp2 fs2],[-As -Ap-1 0 -Ap-1 -As],'r+');
plot([0 fp1 fp1 fp2 fp2 2e3],[-2*As -2*As -Ap-1 -Ap-1 -2*As -2*As],'--m'); % máscara da banda de passagem
plot([0 fs1 fs1 fs2 fs2 2e3],[-As -As -Ap -Ap -As -As],'--m'); % máscara da banda de rejeição
title('Realização Filtro IIR');
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
ylim ([-50 5]);xlim([800 1600]);
hold off;
legend('Filtro original','Forma direta I quantizado com 24 bits', 'Forma direta I quantizado com 24 bits SOS');
%%
%figure(8);
%zplane(nz,dz);
%title('Pólos e zeros do Filtro Digital H(z)');

% figure(9);
% freqz(nz,dz);
% title('Magnitude x Fase de H(z)');

[b3,a3] = freqz(sym2poly(num_z),sym2poly(den_z));
% figure(9);
% subplot(211);
% plot(fa*abs(a3)/pi/2,20*log10(abs(b3)));
% 
% ylim ([-85 5]);
% grid on;
% title('Resposta em frequência do Filtro Digital H(z)');
% ylabel('Atenuação (dB)');
% xlabel('Frequência Normalizada');
% subplot(212);
% plot(fa*abs(a3)/pi/2,unwrap(angle(b3))*180/pi);
% ylim ([-320 320]);
% grid on;
% ylabel('Fase (graus)');
% xlabel('Frequência Normalizada');

%%
%calculando o atraso de grupo
figure(10);
grpdelay(nz,dz);
title('Atraso de grupo');

%%
num_q = nz/nz(end);
den_q = dz/dz(1);

[SOS,g] = zp2sos(num_q,den_q,1);

%equacoes
pretty(Hp)
pretty(Hs)
pretty(Hz)

% Utiliza nz e dz como nominador e denominador.
%%
% Quantização
% Forma direta sem qtz;
Nmultiplicadores = 13;
Nsomadores = 11;
Natrasos = 12; % Atraso de grupo
Nbits = 24; % Alteração para o ponto fixo. 
AHw_1 = max((2*Nmultiplicadores + Nsomadores),Natrasos) * Nbits % Área de hardware

% Forma direta com qtz;
Nmultiplicadores = 12;
Nsomadores = 11;
Natrasos = 12; % Atraso de grupo
Nbits = 24; % Alteração para o ponto fixo.
AHw_2 = max((2*Nmultiplicadores + Nsomadores),Natrasos) * Nbits % Área de hardware

% Forma direta com qtz com sos;
Nmultiplicadores = 13;
Nsomadores = 12;
Natrasos = 12; % Atraso de grupo
Nbits = 24; % Alteração para o ponto fixo.

AHw_3 = max((2*Nmultiplicadores + Nsomadores),Natrasos) * Nbits % Área de hardware
