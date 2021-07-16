function funcao_plotar_graficos(arquivo)

load(arquivo);

g = figure(3);
set(g,'name','Evolução no tempo das velocidades linear e angular do robô');
subplot(211)
plot(tempo,Pvel(1,:))
hold on
plot(tempo,Pvel_medido(1,:),'r')
xlabel('tempo em segundos')
ylabel('V [cm/s]')
subplot(212)
plot(tempo,Pvel(2,:)*180/pi)
hold on
plot(tempo,Pvel_medido(2,:)*180/pi,'r')
xlabel('tempo em segundos')
ylabel('W [graus/s]')

g2 = figure(4);
set(g2,'name','Evolução no tempo da configuração do robô')
subplot(311)
plot(tempo,P(1,:))
hold on
xlabel('tempo em segundos')
ylabel('x [cm]')
subplot(312)
plot(tempo,P(2,:))
hold on
xlabel('tempo em segundos')
ylabel('y [cm]')
subplot(313)
plot(tempo,P(3,:)*180/pi)
hold on
xlabel('tempo em segundos')
ylabel('theta [graus]')
