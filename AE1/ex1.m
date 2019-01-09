% Exemplificando as possiveis formas de realizar a filtragem de um sinal x(n)

clc; clear all; close all;
%% Definindo valores iniciais
%Nh = 10; Nx = 20;
Nh = 400; Nx = 10000;
% Sinal x(n) = delta de Kronecker
x = ones(1,Nx);
% A resposta ao inpulso de um sistema h(n) 
% no filtro FIR aos coeficientes b(n) = h(n) 
h = [1:Nh]; 
b = h;
%% Filtrando o sinal e medindo tempos
 
% Filtragem utilizando a convolução
% NOTE: length(y) = length(x) + length(h) -1
tic;  % iniciar a contagem do tempo
y1 = conv(x,h); 
t(1) = toc; % terminar acontagem e mostrar tempo no console
 
% filtragem utilizando a equação recursiva
% NOTE: length(y) = length(x)
tic;
y2 = filter(b,1,x);
t(2) = toc;
 
% filtragem utilizando a equação recursiva
% aumentando o tamanho de x para que length(y3) = length(y1)
x3 = [x zeros(1,length(h)-1)];
tic;
y3 = filter(h,1,x3); 
t(3) = toc;
 
length_y = length(x) + length(h) - 1;
 
% filtragem utilizando a FFT
% a y = IFFT(FFT(x)*FFT(h))
tic;
X = fft(x,length_y);
H = fft(h,length_y);
Y4 = X.*H;
y4 = ifft(Y4);
t(4) = toc;
 
% filtragem utilizando a função fftfilt
% a y = IFFT(FFT(x)*FFT(h))
 
tic
y5 = fftfilt(h,x3);
t(5) = toc;
 
disp('Comprimento do vetor de saída length(y)')
disp(['    ' num2str([length(y1) length(y2) length(y3) length(y4) length(y5)])])
disp('Tempo usado na filtragem em micro segundos')
disp(['    ' num2str(t*1e6) ' us'])
 
%%  Plotando o gráfico
subplot(411);stem(y1);
title('Filtragem utilizando a convolucao + equacao recursiva + recursiva com aumento de x');
hold on;
stem(y2,'xr');
stem(y3,'+m');
legend('y1', 'y2', 'y3')
hold off
subplot(412);stem(y1, 'ob');legend('y1')
title('Filtragem utilizando a convolucao');
subplot(413);stem(y2, 'xr'); hold on; stem(zeros(size(y1)),'.w');hold off; legend('y2')
title('Filtragem utilizando a equacao recursiva');
subplot(414);stem(y3, '+m');legend('y3')
title('Filtragem utilizando a equacao recursiva, aumentando o tamanho de x');
