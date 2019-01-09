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
fp = 1209;
fs = 1075;
fa = 4000;

wp = fp/(fa/2);
ws = fs/(fa/2);
wc = ((wp+ws)/2)*pi;

Gp = 0;
Ap = 1;
As = 30;

Amax = 10;
Amin = As+30;

%%
%PASSO 1 - Escolher o tipo de janela de acordo com a atenuação do lóbulo lateral Asl e As.
%foi escolhida a janela de Bartlett-hanning, com As = 40,77

%%
%PASSO 2 - Estimar a ordem N1 do filtro considerando os parâmetros Dw

fcuts = [fs fp];
mags = [1 0];
devs = [0.05 0.01];
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fa);

%ordem N = 67 porem a ordem esta um pouco alta para a largura da janela,
%vamos diminuir a ordem e ver qual é a menor que mantenha o ganho desejável
%dentro da faixa

n = 57
%%
%PASSO 3 - Escolha da janela e calculo do filtro

 if (mod(n,2)==1) %ímpar
     L = n+2;
 else %par
     L = n+1;
 end
 
%%% Kaiser
kw = kaiser(L,beta);

[hk,Wk]=freqz(kw,1,2048);
% figure(1), plot(Wk/pi,20*log10(abs(hk))); %resposta da janela
% title('janela de kaiser' );

h_fir1 = fir1(n,Wn,'high',kw);
[Hw1,w1] = freqz(h_fir1);

H0 = -(10^(-Ap/20))-0.035;
%%
figure(3);
subplot(211);
plot(w1*fa/2/pi,20*log10(abs(Hw1*H0)));
hold on;grid on;
plot([0 fs fs fa], [-As -As 0 0], '--r');
plot([fp fp fp fa], [-Amin -Ap -Ap -Ap], '--r');
ylabel('Atenuação (dB)');
xlabel('Frequência (Hz)');
hold off;title(['Filtro com a janela de kaiser N = ' num2str(n)]);
xlim([0 fa/2]);ylim([-60 10]);
subplot(212);
plot(w1*fa/2/pi, angle(Hw1)*360/pi);
ylabel('Fase (graus)');
xlabel('Frequência (Hz)');

%%% Polos e zeros do filtro
figure(4);zplane(h_fir1,1);title('Pólos e zeros'); %1 pois não tem denominador no Z

%%%atraso de grupo
figure(5)
grpdelay(h_fir1)
title('Atraso de grupo');

%%
%%%zoom nas extremas do filtro
figure(6);
subplot(121);
plot(w1*fa/2/pi,20*log10(abs(Hw1*H0)));
hold on;grid on;
plot([0 fs fs fa], [-As -As 0 0], '--r');
plot([fp fp fp fa], [-Amin -Ap -Ap -Ap], '--r');
ylabel('Atenuação (dB)');
xlabel('Frequência (Hz)');
hold off;title('zoom em As');
xlim([fs-1 fs+1]);ylim([-31 -29]);

subplot(122);
plot(w1*fa/2/pi,20*log10(abs(Hw1*H0)));
hold on;grid on;
plot([0 fs fs fa], [-As -As 0 0], '--r');
plot([fp fp fp fa], [-Amin -Ap -Ap -Ap], '--r');
hold off;title('zoom em ap');
xlim([fp-5 fp+5]);ylim([-2 0]);
ylabel('Atenuação (dB)');
xlabel('Frequência (Hz)');