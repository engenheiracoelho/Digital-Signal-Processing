%% Varia��o do Experimento 3.2 do livro:
% DINIZ, P. S. R., DA SILVA, E. A. B., e LIMA NETTO, S. Processamento Digital de Sinais: Projeto e An�lise de Sistemas. 2. ed. Porto Alegre: Bookman, 2014. 976 p. ISBN 978-8582601235.
% FILE: Ex3_2.m
 
% An�lise de sinais no dom�nio da frequ�ncia 
fs = 15;   % frequ�ncia de amostragem
f_sinal = 10;  A_sinal = 1;   % freq��ncia e amplitude do sinal 
T = 1;      % Dura��o do sinal
k_noise = 0.2;    % Intensidade do ru�do  
snr = 0;
 
time = 0 : 1/fs : (T-1/fs);
L = length(time);
freq = time * fs/T;
 
% Sinal x(n) com amplitude A_sinal e frequencia de f_sinal (Hz) 
x = A_sinal*sin(2*pi*f_sinal.*time);
 
% Adicionando um ruido com a fun��o randn
noise = k_noise*randn(1,fs*T);
x1 = x + noise;
 
% Adicionando um ruido com a fun��o awgn
x2 = awgn(x,snr);
 
% Obtendo o sinal no dom�nio da frequencia
X = abs(fft(x))/L;
X1 = abs(fft(x1))/L;
X2 = abs(fft(x2))/L;

%obtendo o sinal em db
signaldb_x = 20*log10(abs(X));
signaldb_x1 = 20*log10(abs(X1));
signaldb_x2 = 20*log10(abs(X2));
 
% Obtendo os plots dos sinais no dominio do tempo e dom�nio da frequencia
figure(1);
subplot(311);stem(time,x, 'b');hold on;
stem(time,x1, 'g');
stem(time, x2, 'r');
hold off;
title('Sinal no Dominio do Tempo'); 
legend('x(n)', 'x(n)+rand', 'x(n)awgn', 'Location','south')
xlabel('Tempo (seg)'); ylabel('x(t)');

%[peak,peakid] = findpeaks(abs(X));
%lbl = ('40 Hz');
%lb = strcat(lbl(1:length(peak)));
%lb = num2cell(lb,2) % Convert to cell array
%lbid = 1:length(lb);
subplot(312);stem(freq, (abs(X)), 'b');hold on;
stem(freq, (abs(X1)),'g');
stem(freq,(abs(X2)),'m');
hold off;
%text(freq(peakid(1)), peak(1), lb(lbid));
title('Magnitude do Sinal no Dominio da Frequencia (Linear)');
legend('X(f)', 'X(f)+rand', 'X(f)+awgn', 'Location','north'); hold off;
xlabel('Frequencia (Hz)'); ylabel('Magnitude (linear)');
ylim([0 0.7]);

subplot(313);

stem(freq, signaldb_x, 'b');hold on;
stem(freq, signaldb_x1,'g');
stem(freq, signaldb_x2','y');
hold off;

%plot(freq, signaldb_x2, 'b'); hold on ;plot(freq, signaldb_x1,'g'); plot(freq, signaldb_x,'r'); hold off;

title('Magnitude do Sinal no Dominio da Frequencia (dB)');
legend('X(f)', 'X(f)+rand', 'X(f)+awgn', 'Location','south'); hold off;
%text(freq(peakid(1)), peak(1)-5, lb(lbid));
xlabel('Frequencia (Hz)'); ylabel('Magnitude (dB)');
ylim([-100 100]);


fsa = 10; % frequ�ncia auxiliar de amostragem usada apenas para representa��o dos sinais originais
Tsa = 1/fsa;
time_aux = 0:Tsa:(1-Tsa);
figure(2);
stem(freq,abs(X),'ob');