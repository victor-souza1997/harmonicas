function P3DX_SIMULADOR(experimento,robo,mapabmp,tipodeplot,tempo_max)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AQUI VAI O CÓDIGO DO P3DX_SIM_CONTROL QUE VAI COMEÇAR COMO O BOTÃO %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inicialização das variáveis globais
global Pos Pdes v_sensor s s2 angs Mapa Mapa2 i Vmax Wmax tempo tamos;
global u;
colidiu = 0; % verifica se o robô colidiu para finalizar a simulação
raio_robo = robo.raio; %raio para verificar colisão e plotar o robô
tamos = experimento.tamos; % tempo de amostragem da simulação em segundos

atraso = zeros(3,round(0.25/tamos)); %vetor que modela o atraso cinemático do pioneer (~250ms)

% janela de observação para navegação em primeira pessoa
janela = robo.saturacao; %saturação do sensor

%inicialização do ruido de medição em [%]
ruido = robo.ruido;

% inicialização da posição do robô
x = experimento.rbx;%214;%416%
y = experimento.rby;%327;850%
theta = experimento.ang*pi/180; %-pi/3% 80/180*pi; inicia orientado para o eixo x do plano
Pos = [x ; y ; theta]; % é a posição atual do robô (aqui é a posição inicial).
P = Pos;  % vetor que armazena a sequência de posições do robô durante o experimento
Pvel = [ 0 ; 0]; % vetor que armazena a sequência de velocidades do robô durante o experimento
Pvel_medido = [ 0 ; 0]; % vetor que armazena a sequência de velocidades medidas do robô durante o experimento
Atual = [0; 0]; % inicialização da memória para simulação da dinâmica

%carregando o mapa do ambiente
A = imread(mapabmp);
A = A(:,:,1);
A = A./(max(max(A)));
A = A.*255;
%A = rgb2gray(A);
Mapa2 = A;
% A2 = A(end:-1:1,:); % utilizada no plot para parecer a imagem no SC padrão

[Ay , Ax] = find(A~=255);
Mapa = [Ax,Ay]';


% inicialização da posição de destino do robô
Pdes =  [800,800];%[experimento.dx ; experimento.dy ];%[800, 100]; [100,900];%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ÂNGULOS NO SISTEMA DE COORDENADAS DO ROBÔ! -> 0 graus é o eixo X (eixo de movimentação)
angs = [90 50:-20:-50 -90 -90 -130:-20:-170 170:-20:130 90]*pi/180;

% inicialização dos sensores com zero;
s = zeros(2,length(angs));
s_i = s;
s2 = s; %sensor com ruído de medição no SC do ambiente

% limite do campo de visão para primeira pessoa.
% [vlimitex,vlimitey] = scircle2(0,0,1,0); % circulo limite de visão MATLAB
aux = linspace(0,2*pi*100,100);
vlimitex = cos(aux)';
vlimitey = sin(aux)';
xc = vlimitex*raio_robo';
yc = vlimitey*raio_robo';



% calculo dos pontos de interesse do plot do robô
%[xc , yc] = scircle2(0,0,raio_robo,0); MATLAB
xc = [xc(1:13)-7.5*xc(1:13).^0 ; 18-7.5 ; -8.369 ; -10.369 ;  -20.369 ; xc(41:60) ; -20.369 ; -10.369 ; -8.369 ; 18-7.5 ; xc(88:end)-7.5*xc(88:end).^0];
yc = [yc(1:13) ; 19 ; 19 ; 16.5 ; 16.5 ; yc(41:60) ; -16.5 ; -16.5 ; -19; -19 ; yc(88:end)];


% inicialização das variáveis do loop while
Vmedido = 0;
Wmedido = 0;
dist = sqrt( (Pdes(1)-Pos(1))^2 + (Pdes(2)-Pos(2))^2 ); % distância linear para o destino
tempo = 0:tamos:tempo_max;  % controle de tempo
i = 0;  % contador

% limites de velocidade do robô
aux = robo.vmax;
if aux > 100 , aux = 100; end
Vmax = aux;         % v linear máxima em cm/s

aux = robo.wmax;
if aux > 360 , aux = 360; end
Wmax = aux*pi/180;  % v angular máxima em rad/s

% valores constantes de aceleração
aux = round(robo.acVp);
if aux > 200 , aux = 200; end
acVp = aux; % aceleração + em cm/s2

aux = round(robo.acVn);
if aux < -200 , aux = -200; end
acVn = aux; % aceleração - em cm/s2

aux = round(robo.acWp);
if aux > 300 , aux = 300; end
acWp = aux*pi/180; % aceleração + em rad/s2

aux = round(robo.acWn);
if aux < -300 , aux = -300; end
acWn = aux*pi/180; % aceleração - em rad/s2    


while  (((abs(dist) > 5) || (abs(Vmedido) > 5) || abs(Wmedido) > 0.1) && (colidiu == 0) && i*tamos<tempo_max)
      % distancia maior que 5 cm ou vlin maior q 5 cm/s ou vrot maior que 0.1 rad/s
    tic     
    % atualização das variáveis de controle de tempo
    i = i+1;          
    tempo(i+1) = i*tamos;

    %% tratamento dos sensores (antiga função inicio)
    Ps_i = s_i;
    [ymaxA , xmaxA] = size(A);
    smax = janela*(ones(1,16));
    offset = [14 18.5*ones(1,6) 14 14 18.5*ones(1,6) 14];
    smax = smax + offset;
    %s é um vetor cujos elementos representam o ponto da leitura máxima (a priori) para
    %cada um dos sensores no sistema de coordenadas do robô
    s = [smax.*cos(angs) ; smax.*sin(angs)]; % vetor no SC do robô (ainda saturado)
    s(1,1) = s(1,1) + 7;
    s(1,8) = s(1,8) + 7;
    s(1,9) = s(1,9) + -16;
    s(1,10:15) = s(1,10:15) - 7.5; 
    s(1,16) = s(1,16) + -16;
    
    Ri = [cos(Pos(3)) -sin(Pos(3)); sin(Pos(3)) cos(Pos(3))]; % matriz de rotação
    %s_i é o vetor s colocado no sistema de coordenadas do ambiente (sem ruído) (ainda saturado)
    s_i = Ri*s;
    s_i(1,:) = s_i(1,:) + Pos(1);
    s_i(2,:) = s_i(2,:) + Pos(2);
    
    pos_sensor = [7 , 18.5*cos(angs(2:7)) , 7 , -16 , 18.5*cos(angs(10:15))-7.5 , -16 ; 14 , 18.5*sin(angs(2:7)) , -14, -14 , 18.5*sin(angs(10:15)) , 14];
    pos_sensor = Ri*pos_sensor;
    pos_sensor(1,:) = pos_sensor(1,:) + Pos(1);
    pos_sensor(2,:) = pos_sensor(2,:) + Pos(2);
    
    for k=1:length(angs)
        Vsx = round(linspace(pos_sensor(1,k),s_i(1,k),janela));
        Vsy = round(linspace(pos_sensor(2,k),s_i(2,k),janela));
        Vs = [Vsx ; Vsy];
        for j=1:janela
            if (Vs(1,j)>0 && Vs(1,j)<=xmaxA && Vs(2,j)>0 && Vs(2,j)<=ymaxA) % se a leitura esta dentro do mapa
                if A(Vs(2,j),Vs(1,j)) ~= 255
                    s_i(1,k) = Vs(1,j);
                    s_i(2,k) = Vs(2,j);
                    break;
                end
            end
        end
    end
    %retorna a medição para o sistema de coordenadas do robô
    % O ROBÔ É ORIENTADO PARA O EIXO X.
    Ps_i(1,:) = s_i(1,:) - Pos(1);
    Ps_i(2,:) = s_i(2,:) - Pos(2);
    s = Ri\Ps_i;
    %adiciona ruído na medição (a resolução da medição depende da resolução do
    %mapa)
    s = s + ruido*randn(2,length(angs));
    % sensor com ruido de medição no S.C. do ambiente
    aux = Ri*s;
    s2(1,:) = aux(1,:) + Pos(1);
    s2(2,:) = aux(2,:) + Pos(2);
    
    %% tratamento dos sensores (antiga função FIM)
        
    v_sensor = sqrt( (s2(1,:)-pos_sensor(1,:)).^2 + (s2(2,:)-pos_sensor(2,:)).^2 ); 
    v_colidiu = sqrt( (s_i(1,:)-pos_sensor(1,:)).^2 + (s_i(2,:)-pos_sensor(2,:)).^2 ); % distância real entre o robô e os obstaculos
    ds0 = sort(v_colidiu);  % ds0 é a distância do sensor sem ruído com a menor leitura
    if ds0(1) <= 7 %colisão de pior caso com 
        colidiu = 1;
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%% CONTROLADOR INÍCIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%s_i = sensores no sistema de coordenadas do ambiente
%    = sem ruído e não deve ser usado pelo controlador
 
%s = sensores no sistema de coordenadas do robô
%  = com ruído adicionado. Esse pode ser utilizado pelo controlador.

%s2 = sensores no sistema de coordenadas do ambiente
%  = com ruído adicionado. Esse pode ser utilizado pelo controlador.
   [V , W] = controle_e_navegacao();    

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%% CONTROLADOR FIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Ksir = [V ; 0 ; W];
    % limitar a ação de controle desejada
    if (Ksir(1) > Vmax), Ksir(1) = Vmax; end
    if (Ksir(1) < -Vmax), Ksir(1) = -Vmax; end
    if (Ksir(3) > Wmax), Ksir(3) = Wmax; end
    if (Ksir(3) < -Wmax), Ksir(3) = -Wmax; end


    R = [cos(Pos(3)) sin(Pos(3)) 0 ; -sin(Pos(3)) cos(Pos(3)) 0 ; 0 0 1]; % matriz de rotação
    
    % modelo do atraso cinemático do pioneeer (~250 ms)
    Ksir_at = atraso(:,1);
    atraso = [atraso(:,2:end) Ksir];       
    
    % inicio da simulação da dinamica não linear do pioneer
    if ((sign(Atual(1)) == sign(Ksir_at(1))) && (abs(Atual(1)) <= abs(Ksir_at(1))))
        acV = abs(acVp)*sign(Atual(1));
        flagvn = 0;
    else
        acV = -abs(acVn)*sign(Atual(1));
        flagvn = 1;
    end
    if sign(Atual(1)) == 0
        acV = abs(acVp)*sign(Ksir_at(1));
        flagvn = 0;
    end
    
    if ((sign(Atual(2)) == sign(Ksir_at(3))) && (abs(Atual(2)) <= abs(Ksir_at(3))))
        acW = abs(acWp)*sign(Atual(2));
        flagwn = 0;
    else
        acW = -abs(acWn)*sign(Atual(2));
        flagwn = 1;
    end
    if sign(Atual(2)) == 0
        acW = abs(acWp)*sign(Ksir_at(3));
        flagwn = 0;
    end
    
    Atualn = Atual + [acV ; acW]*tamos;
    if (flagvn && (sign(Atualn(1)) ~= sign(Atual(1))) )
        Atualn(1) = 0;
    end
    if ((abs(Ksir_at(1)-Atual(1)) < abs(acV*tamos)))
        Atualn(1) = Ksir_at(1);
    end
    if (flagwn && (sign(Atualn(2)) ~= sign(Atual(2))) )
        Atualn(2) = 0;
    end
    if ((abs(Ksir_at(3)-Atual(2)) < abs(acW*tamos)))
        Atualn(2) = Ksir_at(3);
    end
    Atual = Atualn;
    % final da simulação da dinamica não linear do pioneer
    
        
    Ksi_I = R\[Atual(1) ; 0 ; Atual(2)]; % velocidade desejada no SC do ambiente
    Pos = Pos + Ksi_I*tamos; % atualização da posição do robô
    % converte theta para -pi a pi
    if Pos(3) > pi, Pos(3) = Pos(3) - 2*pi; end
    if Pos(3) < -pi, Pos(3) = Pos(3) + 2*pi; end
    Vmedido = Atual(1);
    Wmedido = Atual(2);
    

    Pvel = [Pvel [Ksir(1) ; Ksir(3)]]; % atualiza o vetor das velocidades (comandos) do robô durante o experimento.
    Pvel_medido = [Pvel_medido [Vmedido ; Wmedido]]; % atualiza o vetor das velocidades reais do robô durante o experimento.
    P = [P Pos]; % atualiza o vetor das posições do robô durante o experimento (SC do ambiente).
    
    dist = sqrt( (Pdes(1)-Pos(1))^2 + (Pdes(2)-Pos(2))^2 ); % atualiza a distância para o destino
        
   
%% PLOT DO GRÁFICO "ON LINE"
%    abs(dist)  
%    axis equal
    if tipodeplot == 1      
        if i == 1
            % plot do robô em primeira pessoa
            pxyc = [cos(pi/2) -sin(pi/2) ; sin(pi/2) cos(pi/2)]*[xc' ; yc'];
            plot(pxyc(1,:),pxyc(2,:))
            hold on
            plot([0 0],[0 18.5],'r');
            plot(vlimitex*(50),vlimitey*(50),'g')%,'g','LineWidth',1)
            plot(vlimitex*(100),vlimitey*(100),'g')%,'g','LineWidth',1)
            plot(vlimitex*(150),vlimitey*(150),'g')%,'g','LineWidth',1)
            plot(vlimitex*(200),vlimitey*(200),'g')%,'g','LineWidth',1)
            plot(vlimitex*(250),vlimitey*(250),'g')%,'g','LineWidth',1)
            axis equal;
            set(gca,'xtick',[],'ytick',[])
            xlim([-(janela+50)-10 (janela+50)+10])
            ylim([-(janela+50)-10 (janela+50)+10])
            
            % função: rotacionar 90 graus para ficar condizente com o plot do
            % robô, que é apontado para cima.
            Ps_i = [cos(pi/2) -sin(pi/2) ; sin(pi/2) cos(pi/2)]*s;
            p1 = plot(Ps_i(1,:),Ps_i(2,:),'.','MarkerEdgeColor','k','MarkerSize',7);
            % plota o destino no ambiente de navegação
            Rmapa = [cos(-Pos(3)+pi/2) sin(-Pos(3)+pi/2); -sin(-Pos(3)+pi/2) cos(-Pos(3)+pi/2)]; % matriz de rotação do mapa para o robô
            Pdesmapa = Rmapa\[Pdes(1)-Pos(1)  Pdes(2)-Pos(2)]';
            if dist > (janela+50) 
                angmapa = atan2(Pdesmapa(2),Pdesmapa(1));
                p2 = plot(cos(angmapa)*(janela+50),sin(angmapa)*(janela+50),'.','MarkerEdgeColor','r','MarkerSize',20)       
            else
                p2 = plot(Pdesmapa(1),Pdesmapa(2),'.','MarkerEdgeColor','r','MarkerSize',20)
            end
        end    
        % função: rotacionar 90 graus para ficar condizente com o plot do
        % robô, que é apontado para cima.
        Ps_i = [cos(pi/2) -sin(pi/2) ; sin(pi/2) cos(pi/2)]*s;
        set(p1,'Xdata',Ps_i(1,:),'Ydata',Ps_i(2,:));
        % plota o destino no ambiente de navegação
        Rmapa = [cos(-Pos(3)+pi/2) sin(-Pos(3)+pi/2); -sin(-Pos(3)+pi/2) cos(-Pos(3)+pi/2)]; % matriz de rotação do mapa para o robô
        Pdesmapa = Rmapa\[Pdes(1)-Pos(1)  Pdes(2)-Pos(2)]';
        if dist > (janela+50) 
            angmapa = atan2(Pdesmapa(2),Pdesmapa(1));
            set(p2,'Xdata',cos(angmapa)*(janela+50),'Ydata',sin(angmapa)*(janela+50))       
        else
            set(p2,'Xdata',Pdesmapa(1),'Ydata',Pdesmapa(2))
        end
    end
    
    if tipodeplot == 3 % plot do robô em terceira pessoa      
      if i == 1
        plot(0,0,'.')
        hold on
        image(A)
        colormap(gray(2))
        
        % plot do robô em terceira pessoa
        pxyc = [cos(Pos(3)) -sin(Pos(3)) ; sin(Pos(3)) cos(Pos(3))]*[xc' ; yc'];
        xc3 = pxyc(1,:)+Pos(1);
        yc3 = pxyc(2,:)+Pos(2);
        p1 = plot(xc3,yc3,'b')
        p2 = plot([Pos(1) Pos(1)+18.5*cos(Pos(3))],[Pos(2) Pos(2)+18.5*sin(Pos(3))],'r');
        axis equal;
        p3 = plot(s2(1,:),s2(2,:),'.','MarkerEdgeColor','k','MarkerSize',7);
        % plota o destino no ambiente de navegação
        plot(Pdes(1),Pdes(2),'.','MarkerEdgeColor','r','MarkerSize',20)
        set(gca,'xtick',[],'ytick',[])        
        xlim([0 size(A,2)])
        ylim([0 size(A,1)])        
      end      
      % plot do robô em terceira pessoa
      pxyc = [cos(Pos(3)) -sin(Pos(3)) ; sin(Pos(3)) cos(Pos(3))]*[xc' ; yc'];
      xc3 = pxyc(1,:)+Pos(1);
      yc3 = pxyc(2,:)+Pos(2);
      set(p1,'Xdata',xc3,'Ydata',yc3);
      set(p2,'Xdata',[Pos(1) Pos(1)+18.5*cos(Pos(3))],'Ydata',[Pos(2) Pos(2)+18.5*sin(Pos(3))]);
      set(p3,'Xdata',s2(1,:),'Ydata',s2(2,:));          
    end
    
    
    
drawnow
%while (toc < tamos/2) end

end
disp(tempo(end))
tempo = tempo(1:i+1);

%%%% Salvando os dados no arquivo
save -mat experimento.mat  P Pvel Pvel_medido tempo A Ax Ay Pdes raio_robo

   

end