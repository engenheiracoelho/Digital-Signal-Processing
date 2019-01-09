% BP - Filter 3 - Parks-McCleallan
% Jessica & Leticia

clc;close all;clear all;

f1 = 1088; f2 = 1185; f3 = 1233; f4 = 1330; fe = [f1 f2 f3 f4];
Ap = 1;As = 40;Gp = 0;fs = 4000;

ds = 10^(-As/20)+0.000915;dp = (1 - (10^(-Ap/20)))/2;
dev = [ds dp ds];
a = [0 1 0];
[n,fo,ao,w] = firpmord(fe,a,dev,fs);
n = n+5;
b = firpm(n,fo,ao,w);b = b*(1-dp);
atraso_de_grupo = grpdelay(b)
%wvtool(b);
%figure(2);stem(b);title('Resposta de filtro digital');
%figure(3);zplane(b,1);title('Pólos e zeros'); %1 pois não tem denominador no Z
%figure(4);freqz(b,1,8000);

%%
% % % Resposta do filtro % % %
%favor exportar no fdatool as variaveis para o plot funcionar
%funcao sem quant sem sos
[hdir,wdir] = freqz(Num_dir,Den_dir,40000);

%funcao com quant sem sos
[hsim,wsim] = freqz(Num_sim,1,40000);

figure(5);
[hz,wz] = freqz(b,1,8000,fs);
plot(wz,20*log10(abs(hz))); 
hold on;
title('Realização Filtro FIR Parks McClellan');
plot(fs*abs(wdir)/pi/2,20*log10(abs(hdir)));
plot(fs*abs(wsim)/pi/2,20*log10(abs(hsim)));
legend('Filtro original','Forma direta quantizada com 14 bits','Forma simetrica quantizada com 14 bits');
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
Amin = As + 20;
ylim([-Amin Gp+10]);
plot([0 fe(1) fe(1) fe(4) fe(4) fs/2],[-As -As Gp Gp -As -As], '--r');
plot([fe(2) fe(2) fe(2) fe(3) fe(3) fs/2],[-As*100 -As*100 -Ap -Ap -As*100 -As*100], '--r');

f_mask = [0 fe(1) fe(1) fe(4) fe(4) fs/2 fs/2 fe(2) fe(2) fe(3) fe(3)]/(fs/2);
a_mask = [-As -As Gp Gp -As -As -Amin -Amin -Ap -Ap -Amin];

%%
% Quantização
% Forma direta;
Nmultiplicadores = 65;
Nsomadores = 64;
Natrasos = 64; % Atraso de grupo
Nbits = 14; % Alteração para o ponto fixo. 
AHw_DIRETO = max((2*Nmultiplicadores + Nsomadores),Natrasos) * Nbits % Área de hardware

% Forma simétrico;
Nmultiplicadores = 33;
Nsomadores = 65;
Natrasos = 64; % Atraso de grupo
Nbits = 14; % Alteração para o ponto fixo. 

AHw_Simetrico = max((2*Nmultiplicadores + Nsomadores),Natrasos) * Nbits % Área de hardware
% OPEN THE FDATOOL
% IMPORT PROJECT
% CHANGE NUM(B) AND DENUM(1)
% CHECK IF THE ORDER AND GROUP DELAY IS OK
% IN 'edit' USE THE SYMETRIC FORM
% CREATE ALL THE MASKS
% Faz o teste com o ponto fixo para encontrar o numero de bits, o filtro tem que continuar dentro da mascara. 