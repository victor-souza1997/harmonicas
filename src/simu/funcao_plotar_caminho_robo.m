function funcao_plotar_caminho_robo(arquivo)

load(arquivo);

g = figure(2);
set(g,'name','Resultado final - trajetória executada pelo robô');
plot(P(1,:),P(2,:))
axis equal
hold on
x = P(1,end);
y = P(2,end);
theta = P(3,end);


% calculo dos pontos de interesse da plotagem
aux = linspace(0,2*pi*100,100);
vlimitex = cos(aux)';
vlimitey = sin(aux)';
xc = vlimitex*raio_robo';
yc = vlimitey*raio_robo';

xc = [xc(1:13)-7.5 ; 18-7.5 ; -8.369 ; -10.369 ;  -20.369 ; xc(41:60) ; -20.369 ; -10.369 ; -8.369 ; 18-7.5 ; xc(88:end)-7.5];
yc = [yc(1:13) ; 19 ; 19 ; 16.5 ; 16.5 ; yc(41:60) ; -16.5 ; -16.5 ; -19; -19 ; yc(88:end)];


        
pxyc = [cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[xc' ; yc'];
xc3 = pxyc(1,:)+x;
yc3 = pxyc(2,:)+y;
plot(xc3,yc3,'b')
plot([x x+18.5*cos(theta)],[y y+18.5*sin(theta)],'r');


[linhas , colunas] = size(A);
xlim([-30 colunas+30])
ylim([-30 linhas+30]) 
plot(Ax,Ay,'.','MarkerEdgeColor','k','MarkerSize',1)
plot(Pdes(1),Pdes(2),'.','MarkerEdgeColor','r','MarkerSize',20)
set(gca,'xtick',[],'ytick',[])
drawnow