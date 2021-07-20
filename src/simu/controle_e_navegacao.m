function [V,W] = controle_e_navegacao()
%% IMPORTANTE: APENAS AS VARI�VEIS DISPONIBILIZADAS EM 'VARI�VEIS DISPON�VEIS' PODEM SER UTILIZADAS


%% SOBRE O AMBIENTE SIMULADO:

% - O ROB�: Pioneer 3 DX

% - ESCALA: cada pixel na imagem do mapa corresponde a um cent�metro

% - SISTEMA DE COORDENADAS DO ROB�: posicionado no centro do rob� com o eixo X
% apontado para a frente do rob� e o eixo Y definido pela regra da m�o direita

% - AMBIENTE DE NAVEGA��O: pode ser construido livremente atrav�s de
% figuras *.bmp respeitando a escala informada acima. Pixels brancos
% representam regi�es livres de obst�culos e outras cores representam
% obst�culos no ambiente.


%% VARI�VEIS DISPON�VEIS (SOMENTE LEITURA!!! N�O ALTERAR!!!)
global Pos Pdes v_sensor s s2 angs Mapa Mapa2 i Vmax Wmax tempo tamos;
global u ind_esc;
% ATEN��O: UTILIZE AS VARI�VEIS GLOBAIS SOMENTE PARA LEITURA!!!
% ATEN��O: UTILIZE AS VARI�VEIS GLOBAIS SOMENTE PARA LEITURA!!!
% ATEN��O: UTILIZE AS VARI�VEIS GLOBAIS SOMENTE PARA LEITURA!!!

%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRI��O %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pos = vetor coluna de dimens�o 3 [x ; y ; theta] -> configura��o do rob�
% x em [cm] / y em [cm] / theta em [rad] de -pi a pi

% Pdes = vetor coluna de dimens�o 2 [x_des ; y_des] -> destino do rob�
% xdes em [cm] / ydes em [cm]

% v_sensor = vetor de dimens�o 16 com a leitura dos sensores de dist�ncia do pioneer.
% Cada elemento do vetor � a dist�ncia em cent�metros medida pelo sensor
% at� o obst�culo.

% s2 = matriz de dimens�o 2x16. Cada coluna de 's2' � um ponto (x,y)
% representando, no sistema de coordenadas do ambiente, a medi��o de um dos
% 16 sensores de dist�ncia do rob�.

% s = matriz de dimens�o 2x16. Cada coluna de 's' � um ponto (x,y)
% representando, no sistema de coordenadas do rob�, a medi��o de um dos
% 16 sensores de dist�ncia do rob�.
% ATEN��O: Os sensores s�o confi�veis para leituras a partir de 12 cm

% angs = vetor linha de dimens�o 16 cujos elementos s�o os �ngulos de cada um
% dos sensores de dist�ncia do rob�. Os �ngulos s�o medidos em radiano em
% rela��o ao eixo x positivo (frente do rob�) do sistema de coordenadas no rob�

% Mapa = matriz de dimens�o 2xN que representa o mapa do ambiente
% dispon�vel para o sistema de planejamento. O mapa � uma representa��o
% discretizada do ambiente utilizando grade regular (decomposi��o aproximada).
% 'N' � o n�mero de c�lulas ocupadas por obst�culos e as colunas da matriz 'Mapa'
% s�o as posi��es dos obst�culos (x,y) definidas no sistema de coordenadas
% do ambiente em cent�metros.
% CUIDADO! -> O MAPA PODE NAO SER PERFEITO, PODE SER DESATUALIZADO E
% PODEM HAVER OBST�CULOS M�VEIS NO AMBIENTE COMO OUTROS ROB�S, POR EXEMPLO.

% Mapa2 = matriz de mesma dimens�o que o ambiente. Representa a mesma
% informa��o que a vari�vel Mapa, mas organizada de forma diferente. O
% elemento (i,j) de Mapa2 representa uma c�lula do mapa discreizado. As
% c�lulas livre (naveg�veis) s�o representadas por 255 e as c�lulas
% ocupadas (obst�culos) s�o representadas por 0. Cada c�lula (elemento da matriz)
% representa um quadrado de 1 cent�metro de lado no ambiente
% CUIDADO! -> O MAPA PODE NAO SER PERFEITO, PODE SER DESATUALIZADO E
% PODEM HAVER OBST�CULOS M�VEIS NO AMBIENTE COMO OUTROS ROB�S, POR EXEMPLO.

% i = contador. N�mero de vezes que o algoritmo de controle e navega��o foi
% chamado. Primeira chamada tem valor i = 1.

% Vmax = velocidade linear m�xima, definida na interface gr�fica em cm/s

% Wmax = velocidade angular m�xima, definida na interface gr�fica em rad/s

%tempo = vetor de instantes de tempo. tempo(end) � o instante atual em segundos.

%tamos = per�odo de amostragem m�dio (considera as ultimas 30 itera��es)

% !!!!!! - voc� pode manipular essas vari�veis como quiser. Alguns exemplos de
% manipula��es �teis aparecem abaixo.

%% chama o planejador 

if i == 1
  dirichlet
  
end
load plano

if i == 2
  plot(vxr*1/ind_esc, vyr*1/ind_esc, 'r', 'LineWidth', 2);
  quiver(1/ind_esc*(1:1000*ind_esc), 1/ind_esc*(1:1000*ind_esc), dxn, dyn);
end


%% Dist�ncia at� o destino (d)

d = sqrt((Pdes(1)-Pos(1))^2 + (Pdes(2)-Pos(2))^2); %dist�ncia at� o destino

%% C�lculo do erro de orienta��o para o destino (theta_e)
theta_d = theta_des( round(Pos(2)*ind_esc),round(Pos(1)*ind_esc)); % �ngulo de destino de -pi a pi
theta_e = theta_d - Pos(3);

% converte theta_e para -pi a pi
if theta_e > pi, theta_e = theta_e - 2*pi; end
if theta_e < -pi, theta_e = theta_e + 2*pi; end


%% Dist�ncia do rob� at� o obst�culo mais pr�ximo (d_obs_min)
[d_obs_min, idc_min] = min(v_sensor(2:6));

%% C�lculo do erro de orienta��o para o obst�culo mais pr�ximo (theta_e_obs)
theta_obs = angs(idc_min+1); % �ngulo para o obst�culo de -pi a pi
theta_e_obs = -(theta_obs + pi);
% converte theta_e_obs para -pi a pi
if theta_e_obs > pi, theta_e_obs = theta_e_obs - 2*pi; end
if theta_e_obs < -pi, theta_e_obs = theta_e_obs + 2*pi; end


k1 = 0.5;
k2 = 1.2;
V = k1*d*cos(theta_e);
W = k2*theta_e;

plot(Pos(1), Pos(2), '.r')
 
end






    