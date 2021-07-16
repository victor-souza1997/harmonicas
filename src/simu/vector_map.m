
close all
clear

mapabmp = "mapas/map.bmp"
A = imread(mapabmp);
A = A(:,:,1);
A = A./(max(max(A)));
A = A.*255;
%A = rgb2gray(A);
Mapa2 = A;
%A2 = A(end:-1:1,:); % utilizada no plot para parecer a imagem no SC padrão

[Ay , Ax] = find(A~=255);
Mapa = [Ax,Ay]';
p_x = 40; 
p_y = 40
k = 0.01;
kr = 10;
p0 = 200;
[y,x] = find(A == 0);
for i = 1:1001
  for j = 1:1001
    M_Ua(i, j) = sqrt(k*( (i - p_y)^2 + (j - p_x)^2));
    p = min(sqrt(sum([i-x, j-y].^2,2))) + 1;

    M_Ur(i, j) = kr*(1/p - 1/p0)^2;
  endfor
endfor

surfc(M_Ur+M_Ua);