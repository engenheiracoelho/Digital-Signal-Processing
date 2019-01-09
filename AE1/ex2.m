%% Variação do Experimento 3.2 do livro:
% DINIZ, P. S. R., DA SILVA, E. A. B., e LIMA NETTO, S. Processamento Digital de Sinais: Projeto e Análise de Sistemas. 2. ed. Porto Alegre: Bookman, 2014. 976 p. ISBN 978-8582601235.
% FILE: Ex3_2.m
 
% Análise de sinais no domínio da frequência 
fs = 200;   % frequência de amostragem
f_sinal = 110;  A_sinal = 1;   % freqüência e amplitude do sinal 
T = 1;      % Duração do sinal
k_noise = 0.3;    % Intensidade do ruído  
snr = 40;
 
time = 0 : 1/fs : (T-1/fs);
L = length(time);
freq = time * fs/T;
 
% Sinal x(n) com amplitude A_sinal e frequencia de f_sinal (Hz) 
x = A_sinal*sin(2*pi*f_sinal.*time);
 
% Adicionando um ruido com a função randn
noise = k_noise*randn(1,fs*T);
x1 = x + noise;
 
% Adicionando um ruido com a função awgn
x2 = awgn(x,snr);
 
% Obtendo o sinal no domínio da frequencia
X = abs(fft(x))/L;
X1 = abs(fft(x1))/L;
X2 = abs(fft(x2))/L;

%obtendo o sinal em db
signaldb_x = 20*log10(abs(X));
signaldb_x1 = 20*log10(abs(X1));
signaldb_x2 = 20*log10(abs(X2));
 
% Obtendo os plots dos sinais no dominio do tempo e domínio da frequencia
figure(1);
subplot(311);plot(time,x, 'b', time,x1, 'g', time, x2, 'r');
title('Sinal no Dominio do Tempo'); 
legend('x(n)', 'x(n)+rand', 'x(n)awgn', 'Location','south')
xlabel('Tempo (seg)'); ylabel('x(t)');

[peak,peakid] = findpeaks(abs(X));
lbl = ('40 Hz');
lb = strcat(lbl(1:length(peak)));
lb = num2cell(lb,2) % Convert to cell array
lbid = 1:length(lb);
subplot(312);plot(freq, (abs(X)),maxx, 'b'); hold on ;plot(freq, (abs(X1)),'g');plot(freq,(abs(X2)),'r'); 
text(freq(peakid(1)), peak(1), lb(lbid));
title('Magnitude do Sinal no Dominio da Frequencia (Linear)');
legend('X(f)', 'X(f)+rand', 'X(f)+awgn', 'Location','north'); hold off;
xlabel('Frequencia (Hz)'); ylabel('Magnitude (linear)');
ylim([0 0.7]);

subplot(313);plot(freq, signaldb_x2, 'b'); hold on ;plot(freq, signaldb_x1,'g'); plot(freq, signaldb_x,'r'); hold off;
title('Magnitude do Sinal no Dominio da Frequencia (dB)');
legend('X(f)', 'X(f)+rand', 'X(f)+awgn', 'Location','south'); hold off;
text(freq(peakid(1)), peak(1)-5, lb(lbid));
xlabel('Frequencia (Hz)'); ylabel('Magnitude (dB)');
ylim([-100 0]);