function [V,W] = controle_e_navegacao()
%% IMPORTANTE: APENAS AS VARIÁVEIS DISPONIBILIZADAS EM 'VARIÁVEIS DISPONÍVEIS' PODEM SER UTILIZADAS


%% SOBRE O AMBIENTE SIMULADO:

% - O ROBÔ: Pioneer 3 DX

% - ESCALA: cada pixel na imagem do mapa corresponde a um centímetro

% - SISTEMA DE COORDENADAS DO ROBÔ: posicionado no centro do robô com o eixo X
% apontado para a frente do robô e o eixo Y definido pela regra da mão direita

% - AMBIENTE DE NAVEGAÇÃO: pode ser construido livremente através de
% figuras *.bmp respeitando a escala informada acima. Pixels brancos
% representam regiões livres de obstáculos e outras cores representam
% obstáculos no ambiente.


%% VARIÁVEIS DISPONÍVEIS (SOMENTE LEITURA!!! NÂO ALTERAR!!!)
global Pos Pdes v_sensor s s2 angs Mapa Mapa2 i Vmax Wmax tempo tamos;
global u ind_esc;
% ATENÇÃO: UTILIZE AS VARIÁVEIS GLOBAIS SOMENTE PARA LEITURA!!!
% ATENÇÃO: UTILIZE AS VARIÁVEIS GLOBAIS SOMENTE PARA LEITURA!!!
% ATENÇÃO: UTILIZE AS VARIÁVEIS GLOBAIS SOMENTE PARA LEITURA!!!

%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIÇÃO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pos = vetor coluna de dimensão 3 [x ; y ; theta] -> configuração do robô
% x em [cm] / y em [cm] / theta em [rad] de -pi a pi

% Pdes = vetor coluna de dimensão 2 [x_des ; y_des] -> destino do robô
% xdes em [cm] / ydes em [cm]

% v_sensor = vetor de dimensão 16 com a leitura dos sensores de distância do pioneer.
% Cada elemento do vetor é a distância em centímetros medida pelo sensor
% até o obstáculo.

% s2 = matriz de dimensão 2x16. Cada coluna de 's2' é um ponto (x,y)
% representando, no sistema de coordenadas do ambiente, a medição de um dos
% 16 sensores de distância do robô.

% s = matriz de dimensão 2x16. Cada coluna de 's' é um ponto (x,y)
% representando, no sistema de coordenadas do robô, a medição de um dos
% 16 sensores de distância do robô.
% ATENÇÃO: Os sensores são confiáveis para leituras a partir de 12 cm

% angs = vetor linha de dimensão 16 cujos elementos são os ângulos de cada um
% dos sensores de distância do robô. Os ângulos são medidos em radiano em
% relação ao eixo x positivo (frente do robô) do sistema de coordenadas no robô

% Mapa = matriz de dimensão 2xN que representa o mapa do ambiente
% disponível para o sistema de planejamento. O mapa é uma representação
% discretizada do ambiente utilizando grade regular (decomposição aproximada).
% 'N' é o número de células ocupadas por obstáculos e as colunas da matriz 'Mapa'
% são as posições dos obstáculos (x,y) definidas no sistema de coordenadas
% do ambiente em centímetros.
% CUIDADO! -> O MAPA PODE NAO SER PERFEITO, PODE SER DESATUALIZADO E
% PODEM HAVER OBSTÁCULOS MÓVEIS NO AMBIENTE COMO OUTROS ROBÔS, POR EXEMPLO.

% Mapa2 = matriz de mesma dimensão que o ambiente. Representa a mesma
% informação que a variável Mapa, mas organizada de forma diferente. O
% elemento (i,j) de Mapa2 representa uma célula do mapa discreizado. As
% células livre (navegáveis) são representadas por 255 e as células
% ocupadas (obstáculos) são representadas por 0. Cada célula (elemento da matriz)
% representa um quadrado de 1 centímetro de lado no ambiente
% CUIDADO! -> O MAPA PODE NAO SER PERFEITO, PODE SER DESATUALIZADO E
% PODEM HAVER OBSTÁCULOS MÓVEIS NO AMBIENTE COMO OUTROS ROBÔS, POR EXEMPLO.

% i = contador. Número de vezes que o algoritmo de controle e navegação foi
% chamado. Primeira chamada tem valor i = 1.

% Vmax = velocidade linear máxima, definida na interface gráfica em cm/s

% Wmax = velocidade angular máxima, definida na interface gráfica em rad/s

%tempo = vetor de instantes de tempo. tempo(end) é o instante atual em segundos.

%tamos = período de amostragem médio (considera as ultimas 30 iterações)

% !!!!!! - você pode manipular essas variáveis como quiser. Alguns exemplos de
% manipulações úteis aparecem abaixo.

%% chama o planejador 

if i == 1
  dirichlet
  
end
load plano

if i == 2
  plot(vxr*1/ind_esc, vyr*1/ind_esc, 'r', 'LineWidth', 2);
  quiver(1/ind_esc*(1:1000*ind_esc), 1/ind_esc*(1:1000*ind_esc), dxn, dyn);
end


%% Distância até o destino (d)

d = sqrt((Pdes(1)-Pos(1))^2 + (Pdes(2)-Pos(2))^2); %distância até o destino

%% Cálculo do erro de orientação para o destino (theta_e)
theta_d = theta_des( round(Pos(2)*ind_esc),round(Pos(1)*ind_esc)); % ângulo de destino de -pi a pi
theta_e = theta_d - Pos(3);

% converte theta_e para -pi a pi
if theta_e > pi, theta_e = theta_e - 2*pi; end
if theta_e < -pi, theta_e = theta_e + 2*pi; end


%% Distância do robô até o obstáculo mais próximo (d_obs_min)
[d_obs_min, idc_min] = min(v_sensor(2:6));

%% Cálculo do erro de orientação para o obstáculo mais próximo (theta_e_obs)
theta_obs = angs(idc_min+1); % ângulo para o obstáculo de -pi a pi
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






    