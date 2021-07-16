function 

  %%campo potencial

  k1 = 0.09;
  k2 = 0.4;
  %V = k1*d;
  %W = k2*theta_e;
  V = k1*sqrt( (1*(Pos(2)-Pdes(2)))^2+( 1*(Pos(1)-Pdes(1)))^2); %Funcao de atracao 

  if(Pos(3) < 0) %condicional para corrigir o valor do angulo
  aux = 2*pi;
  else
  aux = 0;
  endif
  if(atan2((Pos(2) - Pdes(2)),(Pos(1) - Pdes(1))) < 0) %condicional para corrigir o valor do angulo
    aux4 = pi;
  else
    aux4 = -pi;
  endif
  if(- (Pos(3) ) + atan2((Pos(2) - Pdes(2)),(Pos(1) - Pdes(1)))+aux4 >pi)
    aux5 = -2*pi;
  else
    aux5 = 0;
    
  endif

  W = k2*( - (Pos(3) ) + atan2((Pos(2) - Pdes(2)),(Pos(1) - Pdes(1)))+aux4 + aux5); % velocidade angular apenas com a componente de atracao

  if(any(v_sensor([ 3, 6,2,7]) < 70))% caso um dos sensores dentro desse vetor tenha a distancia menor do que 70, adicionaremos uma componente de repulsao a velocidade linear e angular
    
    xy_obs = find(v_sensor == min(v_sensor([3, 6, 2,7]))); % encontrar posicao do vetor com a menor distancia para o obstaculo dentro do vetor de todos os sensores
    pos_obs = s2(:,xy_obs); %pegar posicao de medicao do sensor de menor distancia
    
    num = (Pos(2) - pos_obs(2)); %diferenca do eixo y do centro do robo para o obstaculo
    den = (Pos(1) - pos_obs(1)); %diferenca do eixo x do centro do robo para o obstaculo
    alpha_obs = atan2(num,den); %componente angular da forca de repulsao do objeto para com o robo
    if(alpha_obs < 0) %correcao do angulo
      aux3 = 2*pi;
    else
      aux3 = 0;
    endif
    th1 = (-Pos(3) + (alpha_obs + aux3));
    if(th1 > 0)
      if(th1 > pi)
        th2 = -1;
        th1 = th1 - pi;
      else
        th2 = 1;
      endif
    else
      if(th1*-1 > pi)
        th2 = -1;
        th1 = th1 + pi;
      else
        th2 = 1;
      endif 
    endif
    W =  W*k2 +0.6*th2*th1;  %
    n=500000;q0 = 10;
    V = V - sqrt(n^2*( 1/((Pos(2) - pos_obs(2))^2+(Pos(1) - pos_obs(1))^2) - 1/q0)^2*((Pos(2) - pos_obs(2))^2+(Pos(1) - pos_obs(1))^2)^-2) ;  
    
  endif


end