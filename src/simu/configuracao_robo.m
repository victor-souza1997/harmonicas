% Par?metros de configura??o das caracter?sticas do rob?.

robo.ruido = 5; %desvio padr?o do ru?do do sensor
robo.saturacao = 10;%200; %limite do alcance sensorial em cm

robo.vmax =   100; % velocidade linear m?xima cm/s (limitada a 100)
robo.acVp =   200; % acelera??o linear (valor fixo positivo) positiva cm/s2 (limitada a 200)
robo.acVn =  -200; % acelera??o linear (valor fixo negativo) negativa cm/s2 (limitada a -200)

robo.wmax =   360;% velocidade angular m?xima graus/s (limitada a 360)
robo.acWp =   200;% acelera??o angular (valor fixo positivo) positiva graus/s2 (limitada a 300)
robo.acWn =  -300;% acelera??o angular (valor fixo negativo) negativa graus/s2 (limitada a -300)

%N?o alterar
robo.raio = 26; %raio do rob? para verificar colis?o e plotar o rob? [cm] - N?O ALTERAR

