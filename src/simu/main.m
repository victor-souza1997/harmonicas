
close all
clc
clear

mapabmp = '../contents/mapa.png';
A = imread(mapabmp);
A = A(:,:,1);
A = A./(max(max(A)));
A1 = A;
A = A;
Mapa = ones(100,100);

Mapa(:,1) = 0;

Mapa(:,end) = 0;
Mapa(1,:) = 0;
Mapa(end,:) = 0;
[ny, nx] = size(Mapa);
cont = 1;
for j=1:ny
    for i =1:nx
      if Mapa(j,i) == 1
       Mapa(j,i) = cont;
       cont = cont +1;
      end
    end
end
%Mapa
A = zeros(cont - 1, cont -1);

for j=1:ny
    for i =1:nx
      if Mapa(j,i) != 0
        idc = Mapa(j,i);
        
        A(idc, idc) = -1;
   
        if Mapa(j+1, i) !=0
          A(idc, Mapa(j+1, i)) = 1/4;
        end
         if Mapa(j-1, i)!=0
          A(idc, Mapa(j-1, i)) = 1/4;
        end
        if Mapa(j,i+1)!=0
          A(idc, Mapa(j,i+1)) = 1/4;
        end
        if Mapa(j,i-1)!=0
          A(idc, Mapa(j,i-1)) = 1/4;
        end
        
      end
    end
end

%incluir destino
idc = 35;
A(idc, :) = 0;
A(idc, idc) = 1;
d = zeros(length(A),1);
d(idc) = -1;
%sistema linear de equacoes a resolver  - > A*u=d

u = A\d;
%colocar o potencial de u para cada uma das posicoes
count = 1;
for j=1:ny
    for i =1:nx
      if Mapa(j,i) !=0
       Mapa(j,i) = u(count);
       count = count +1;
      end
    end
end

%figure(1)
%[Mx, My] = meshgrid(1:length(Mapa))
%surf(Mx,My, Mapa)
%figure(2)
[dx, dy] = gradiente(Mapa);
%%parte de normalização dos vetores do gradiente
modulo = sqrt(dx.^2 + dy.^2);
dxn = - dx./modulo;
dyn = - dy./modulo;

%% fazer plot do campo vetoria
quiver(1:size(Mapa, 2), 1:size(Mapa,1), dxn, dyn, 0.5); 