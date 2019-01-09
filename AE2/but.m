%Filtro 1 - Butterworth
%alunas Jessica e Leticia

clear all; close all;

%especificacoes do projeto
fp = 941; %fr. de passagem
fs = 1075; %fr. de rejeicao
fa = 4000; %fr. de amostragem
Ap = 1; %Atenuação na fr. de passagem
As = 25; %Atenuação na fr. de rejeicao
G0 = 0; %Ganho em db na banda de passagem

%normalizando as fr.
Wp = (2*pi*fp)/fa; %fr. angular
Ws = (2*pi*fs)/fa;
lp = 2*fa*tan((Wp/2)); %lamdas
ls = 2*fa*tan((Ws/2));

%l = 2*fa*atan((((2*pi*(fs+fp))/fa)/2));

ohm_s = ls/lp; %omega_s
ohm_p = lp/lp; %omega_p

%equacoes iniciais:determinacao da ordem,polos e zeros

n = log10(((10^(0.1*As)) - 1)/((10^(0.1*Ap)) -1))/(2*log10(ohm_s)); %ordem do filtro
n = ceil(n)
Eps = sqrt((10^(0.1*Ap)) -1);
Angle = Eps^(-(1/n));

p_k = zeros(1,n);
for (k = 1:n)
    p_k(k) = Angle*exp(j*pi*((2*k+(n-1))/(2*n))); %testar com n+1
end

den_p = real(poly(p_k));  %denominador com os numeros apenas
a_p = den_p(end);  %pega o maior termo para dividir todos os valores depois
syms p;
D_p(p) = poly2sym(den_p/a_p,p);
polos_p = collect(D_p)
H0 = D_p(j*0);
num_p = H0*(10^(G0/20));
zeros_p = num_p

%calculo do prototipo usando funcoes simbolicas
%funcao de transferencia em p
H_p(p) = symfun(num_p/D_p,p);
tf_H_p = H_p(p)

%visualizando o plano z de H(p)
figure(1);
zplane(num_p,den_p);
title('Pólos e zeros do Protótipo H(p)');

%visualizando o filtro com ganho normalizado
wr_p=0:0.01:ohm_s+1;
calc_p = H_p(j*wr_p);

figure(2);
hold on;
plot(wr_p,20*log10(abs(calc_p)));
plot([ohm_p ohm_s],[-Ap -As]+G0,'r*');
plot([0 ohm_p ohm_p ohm_p],[-Ap -Ap -As -As*2],'--m');
plot([0 ohm_s ohm_s ohm_s],[0 0 -As -As*2],'--m');
hold off;
axis([0 1.5 -30 5]);
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
title('Magnitude do Protótipo H(p)');

figure(3)
[b,a] = freqs(double(num_p),double(den_p));
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
H_s(s)=H_p(s/lp); %substitui s em p
%wr_s=0:0.01:ohm_s+1;
wr_s=logspace(0,log10(ls+lp),1000);
calc_s = H_s(1i*wr_s);
tf_H_s = H_s
figure(4);
hold on;
plot([lp ls],[-Ap -As]+G0,'r*');
plot([lp lp 0],[-Ap*100 -Ap -Ap]+G0,'--m'); %  band pass mask
plot([0 ls ls 10*ls],[0 0 -As -As]+G0,'--m'); % band stop mask
plot(wr_s,20*log10(abs(calc_s)));
xlabel('Frequências do Filtro em Lamda (Hz)');
ylabel('Magnitude (dB)');
axis([0 10000 -30 5]);
title('Magnitude do Filtro Analógico H(s)');
hold off;

%visualizando o plano z de H(s)
D_s1(s) = poly2sym(den_p/a_p,s); %transformando os termos de d em funcao de s
D_s(s)= D_s1(s/lp);
den_s = sym2poly(D_s);
num_s = sym2poly(H0);
figure(5);
zplane(num_s,den_s);
title('Pólos e zeros do Filtro Analógico H(s)');

figure(6);
[b1,a1] = freqs(double(num_s),double(den_s));
subplot(211);
semilogx(a1,20*log10(abs(b1)));
ylim([-30 5]);
xlim([0 10000]);
grid on;
title('Resposta em frequência do Filtro Analógico H(s)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
semilogx(a1,angle(b1));
grid on;
xlim([0 10000]);
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');


%passando o filtro para digital

%funcao de transferencia em z usando trans. bilinear
syms z;
s = 2*fa*(z-1)/(z+1);
H_z(z) = collect(subs(H_s,s));
tf_H_z = vpa(H_z)

[Nz,Dz] = numden(H_z(z));
bz = sym2poly(Nz);
az = sym2poly(Dz);

wr_z = 0:0.01:pi;
calc_z = H_z(exp(1i*wr_z));
figure(7);
hold on;
plot([fp fs],[-Ap -As]+G0,'r*');
plot([fp fp 0],[-Ap*100 -Ap -Ap]+G0,'--m'); %  band pass mask
plot([0 fs fs 10*fs],[0 0 -As -As]+G0,'--m'); % band stop mask
plot(fa*wr_z/pi/2,20*log10(abs(calc_z)));
xlabel('Frequências do Filtro (Hz)');
ylabel('Magnitude (dB)');
axis([500 1200 -30 5]);
title('Magnitude do Filtro Digital H(z)');
hold off;

figure(8);
zplane(sym2poly(Nz),sym2poly(Dz));
title('Pólos e zeros do Filtro Digital H(z)');
%%
[b2,a2] = freqz(bz,az);
figure(9);
subplot(211);
plot(fa*abs(a2)/pi/2,20*log10(abs(b2)));
ylim([-30 5]);
grid on;
title('Resposta em frequência do Filtro Digital H(z)');
ylabel('Atenuação (dB)');
xlabel('Frequência Normalizada');
subplot(212);
plot(fa*abs(a2)/pi/2,unwrap(angle(b2))*180/pi);
xlim([0 1200]);
grid on;
ylabel('Fase (graus)');
xlabel('Frequência Normalizada');

%calculando o atraso de grupo
figure(10);
grpdelay(bz,az)
title('Atraso de grupo H(z)');