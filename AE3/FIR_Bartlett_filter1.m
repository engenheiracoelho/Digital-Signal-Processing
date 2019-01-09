clear all;
close all;
clc;

%% f1 = 941 Hz; f2 = 1075 Hz,As = 25 dB  // Bartlett     As = 27.48
%PASSO 1 - Escolher o tipo de janela de acordo com a atenuação do lóbulo lateral Asl e As.
%PASSO 2 - Estimar a ordem N1 do filtro considerando os parâmetros Dw
%PASSO 3 - Calcule os coeficientes clp do filtro LP , calcule os valores da janela w e obtenha a resposta ao impulso do filtro h = clp * w.
%PASSO 4 - Verifique o valor real de Dwr = wAs-wAp, e faça a correção da ordem do filtro em função do desvio constatado. N2 = N*Dwr/Dw.
%PASSO 5 - Corrija o valor de projeto dos coeficientes Clp do filtro ideal, a janela e a resposta ao impulso. Repita o PASSO 3 até 5, até obter um filtro que atenda as especificações de Dw.
%PASSO 6 - Desloque a frequência de corte wc de modo a obter o valor correto de wp. wc2 = wp + (wp-wAp).

% ESPECIFICAÇÕES
fp = 941;
fs = 1075;
fa = 4000;

wp = fp/(fa/2);
ws = fs/(fa/2);
wc = ((wp+ws)/2)*pi;

Gp = 0;
Ap = 1;
As = 25;
Amax = 10;
Amin = -As-30;

%%
%PASSO 1 - Escolher o tipo de janela de acordo com a atenuação do lóbulo lateral Asl e As.
%foi escolhida a janela de Bartlett, com As = 27,48

%%
%PASSO 2 - Estimar a ordem N1 do filtro considerando os parâmetros Dw

fcuts = [fp fs];
mags = [1 0];
devs = [0.05 0.01];

[N,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fa);

%ordem N = 67
%%
%PASSO 3 - Calcule os coeficientes clp do filtro LP , calcule os valores da janela w e obtenha a resposta ao impulso do filtro h = clp * w.

if (mod(N,2)==1) %ímpar
    L = N+2;
    M = (N+1)/2;
else %par
    L = N+1;
    M = N/2;
end

n = -M:M;

%calculo janela
w = bartlett(L);

%coeficiente clp
Clp = sin(wc.*n)./(pi.*n);%LP
Clp (M+1) = wc./pi;

Ajuste = 1;
%Ajuste = (10^(-4.3/20));
H0 = (10^(-Gp/20)*Ajuste);
w = w';

%resposta ao impulso
b = Clp.*w.*H0;

%%
%plots ate o momento
% 
% [h,W]=freqz(b',1,2048);
% figure(1);
% plot((W/pi),20*log10(abs(h)),'DisplayName','Filtro Passa-baixa com a janela' );
% hold on;
% plot([0 wp wp], [-Ap -Ap Amin], ':r');
% plot([0 ws ws pi], -[0 0 As As], ':r');
% xlim([0 1]);ylim([Amin Amax]);
% hold off;

% figure(2);stem(Clp);title('Resposta de filtro digital');
figure(3);zplane(b,1);title('Pólos e zeros'); %1 pois não tem denominador no Z
% figure(4);freqz(Clp,1,8000); title('Resposta do filtro');
% % %Ajustando o ganho
[h,w]=freqz(Clp,1,8000);
% figure(5);
% %H0 = 10^(-Ap/2/20);
% plot(w/pi,20*log10(abs(h)));hold on;
% plot(w/pi,20*log10(abs(h*H0)));
% plot([0 wp wp], [-Ap -Ap Amin], ':r');
% plot([0 ws ws pi], -[0 0 As As], ':r'); title('Resposta do filtro');
% xlim([0 1]);ylim([Amin Amax]);hold off;

%%
% figure(6)
% freqz(b,1,2048); title('Resposta do filtro com a janela');

%%
%PASSO 4 - Verifique o valor real de Dwr = wAs-wAp, e faça a correção da ordem do filtro em função do desvio constatado. N2 = N*Dwr/Dw.

Dwr = (0.5485-0.4736);Dw=(ws-wp)/pi;
N2 = ceil(N*(Dwr/Dw));
M2 = N2/2;
k = [1:N2+1];
wc2 = (ws+wp)/2 - 0.00498;
n2 = k-M2-1;

Clp2 = sin(wc2*n2)./(pi*n2);%LP
Clp2(M2+1) = wc2/pi;

w2 = bartlett(N2+1)';
b2 = Clp2.*w2*H0;

[h2,W2]=freqz(b2,1,2048);

figure(7)
subplot(211);
plot(W2*fa/2,20*log10(abs(h2*H0)));
hold on;grid on;
plot([0 fp fp], [-Ap -Ap Amin], '--r');
plot([0 fs fs fa], -[0 0 As As], '--r');
ylabel('Atenuação (dB)');
xlabel('Frequência (Hz)');
hold off;title('Filtro com a janela de bartlet ajustada ');
xlim([0 fa/2]);ylim([-50 10]);
subplot(212);
plot(W2*fa/2, unwrap(angle(h2)));
xlim([0 fa/2]);
ylabel('Fase (graus)');
xlabel('Frequência (Hz)');

%%
%%%atraso de grupo  = 118
figure(5)
grpdelay(b2)
title('Atraso de grupo');

%%
%%%aplicando zoom para verificar se esta correto
figure(6);
subplot(121);
plot(W2*fa/2,20*log10(abs(h2*H0)));
hold on;grid on;
plot([0 fp fp], [-Ap -Ap Amin], '--r');
plot([0 fs fs fa], -[0 0 As As], '--r');
ylabel('Atenuação (dB)');
xlabel('Frequência (Hz)');
hold off;title('zoom em As');
xlim([fs-1 fs+1]);ylim([-26 -24]);

subplot(122);
plot(W2*fa/2,20*log10(abs(h2*H0)));
hold on;grid on;
plot([0 fp fp], [-Ap -Ap Amin], '--r');
plot([0 fs fs fa], -[0 0 As As], '--r');
hold off;title('zoom em ap');
xlim([fp-5 fp+5]);ylim([-2 0]);
ylabel('Atenuação (dB)');
xlabel('Frequência (Hz)');
