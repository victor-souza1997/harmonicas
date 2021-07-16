clear
close all
clc

configuracao_experimento
configuracao_robo

global u;
mapabmp = '../../contents/mapas/aula4.bmp';


tipodeplot = 3; % em terceira pessoa = 3 --- em primeira pessoa = 1  ---- sem plot = 0
tempo_max = 240; % em segundos

tempo_total = tic;
P3DX_SIMULADOR(experimento,robo,mapabmp,tipodeplot,tempo_max);
disp("tempo total")
toc(tempo_total)


tic
funcao_plotar_caminho_robo('experimento.mat')
toc
funcao_plotar_graficos('experimento.mat')
toc
