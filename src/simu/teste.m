clc
clear

mapabmp = '../../contents/mapas/aula4.bmp';
A = imread(mapabmp);
A = A(:,:,1);
A = A./(max(max(A)));
A = A;
%A = rgb2gray(A);


Mapa = ones(6,6);
Mapa(1:2,1:2) = 0;
escalonar(Mapa, 0.5);


##dim = size(Mapa)*0.5
##dx = dim(1)
##dy = dim(2)
##for j = 0:1/0.5 -1
##  for i = 0:1/0.5-1
##    Mapa(1+(j)*dy:(j+1)*dy,1+(i)*dx:(i+1)*dx)
##  end 
##end