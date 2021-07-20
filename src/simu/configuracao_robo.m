% Parâmetros de configuração das características do robô.

robo.ruido = 5; %desvio padrão do ruído do sensor
robo.saturacao = 10;%200; %limite do alcance sensorial em cm

robo.vmax =   100; % velocidade linear máxima cm/s (limitada a 100)
robo.acVp =   200; % aceleração linear (valor fixo positivo) positiva cm/s2 (limitada a 200)
robo.acVn =  -200; % aceleração linear (valor fixo negativo) negativa cm/s2 (limitada a -200)

robo.wmax =   360;% velocidade angular máxima graus/s (limitada a 360)
robo.acWp =   200;% aceleração angular (valor fixo positivo) positiva graus/s2 (limitada a 300)
robo.acWn =  -300;% aceleração angular (valor fixo negativo) negativa graus/s2 (limitada a -300)

%Não alterar
robo.raio = 26; %raio do robô para verificar colisão e plotar o robô [cm] - NÃO ALTERAR

