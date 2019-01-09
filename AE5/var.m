clear all;

fl = 960;
fh = 1190;

banda = 1.5;
tau = 4; %ms
rej = 30; % dBm
pas = 4; %dbm

f_amostragem = 8000;

%frequencias dtmf

%%% LINHAS

%L2 -- C1, C3 
% OK:4, 6
% ERROR :5, * 
% ZERO: others


%calculo fc = ((697*1.5)/100)+2+ 697
% %fp = 697-((697*1.5)/100)-2
% f1 = 697;
% fp1 = ((f1*1.5)/100)+2+f1;
% fc1 = f1-((f1*1.5)/100)-2;

%BP 770
f2 = 770;
fp2 = ((f2*1.5)/100)+2+f2;
fc2 = f2-((f2*1.5)/100)-2;

% %BP 852
% f3 = 852;
% fp3 = ((f3*1.5)/100)+2+f3;
% fc3 = f3-((f3*1.5)/100)-2;
% 
% %BP 941
% f4 = 941;
% fp4 = ((f4*1.5)/100)+2+f4;
% fc4 = f4-((f4*1.5)/100)-2;

%%% COLUNAS
%BP 1209
f5 = 1209;
fp5 = ((f5*1.5)/100)+2+f5;
fc5 = f5-((f5*1.5)/100)-2;

% %BP 1336
% f6 = 1336;
% fp6 = ((f6*1.5)/100)+2+f6;
% fc6 = f6-((f6*1.5)/100)-2;
% 
%BP 1477
f7 = 1477;
fp7 = ((f7*1.5)/100)+2+f7;
fc7 = f7-((f7*1.5)/100)-2;

%BP 1633
% f8 = 1633;
% fp8 = ((f8*1.5)/100)+2+f8;
% fc8 = f8-((f8*1.5)/100)-2;


%%%              (OK) A frequência de amostragem f_s do sinal de entrada no sistema mostrado acima é de 8 kHz.
%%%              (OK) O circuito retificador deve se implementado pela função abs.
%%%              (OK) O circuito TC - Threshold Comparator não deve ser implementado, pois tem apenas a função de rejeitar sinais de entrada que estão acima de -3dBm ou abaixo de -40 dBm,
%%%                   Como sinais de entrada devem com frequência de amostragem de 80kHz e serão utilizados: 
%%%                  1) Os 4 Sinal de DTMF correspondentes as 4 teclas indicadas para a equipe, com duração de 1 segundo
%%%                  2) Um sinal DTMF correspondente a sequencia de teclas "123456789*0#" com duração de tom de 65ms e pausa de 65ms para cada tecla. (ver ETSI ES 201 235-2 - Specification of Dual Tone Multi-Frequency (DTMF) Transmitters and Receivers; Part 2: Transmitters)
%%%                  3) Os sinais do item (1) e (2) somados a um ruído branco, cuja relação sinal/ruído deve alterável entre 0 dB até 80 dB. 

%%%              A seleção do sinal de entrada e adição do ruido pode ser feita através de chaves manuais.
%%%              Os discriminadores das linhas e colunas não projetadas, deverão ser substituídos por chaves manuais para simular tanto a inserção dessas frequências como para gerar respostas com erro.
%%%              Caso seja inserida uma informação DTMF que ative as frequências especificadas, mas que não corresponda a um dos "Dual Tone", o discriminador deverá indicar um código de ERRO "1111" (15) na saída. Se não houver nenhuma dessas frequências, deverá ser indicado com um código de VAZIO "0000" (0). Recomenda-se que o número 0 seja representado pelo código binário "1010" (10).
%%%              Neste projeto os sinais de entrada a serem utilizados deverão ser gerados com frequência de amostragem especificada. Antes de realizar o processamento indicado no diagrama do receptor DTMF, deve ser feita uma filtragem passa-baixa com um filtro de no mínimo 2 polos do tipo Butterworth ou Chebychev 1, com fp = 4 kHz, seguido de uma subamostragem para a nova frequência de = 8 kHz.
%%%              A contante de tempo Τ = RC é equivalente a uma frequência de corte de wc = 1 / (2πΤ)
%%%              Após realizada a simulação do sistema detector DTMF, o sistema deverá ser convertido para VHDL utilizando o HDL coder.
%%%              Para comprovar o funcionamento do sistema realize a simulação final do sistema no ModelSim. 