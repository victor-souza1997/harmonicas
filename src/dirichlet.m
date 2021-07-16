global Pdes Pos Mapa2
pkg load image
Amapa = Mapa2;

%inflar obstaculos para aumentar a seguraca e evitar coliscao

se = strel('ball', 25, 25);
Amapa = imerode(Amapa, se);

%redicionamento do mapa para escalonar
Amapa = imresize(Amapa, 0.1);

%fechar o ambiente virtualmente para o planejamento

Amapa(:,1) = 0;
Amapa(:,end) = 0;
Amapa(1,:) = 0;
Amapa(end, :) = 0;

Amapa = double(Amapa);
aux = Amapa(:);
valor = mode(aux) - abs(min(min(Amapa)));
[Ay, Ax] = find(Amapa>=valor);
[Oy, Ox] = find(Amapa<valor);
for j = 1:length(Oy)
	Amapa(Oy(j), Ox(j)) = 0; %potencial nos obstaculos
end

x_i = Pos(1)*0.1; %unidade no planejamento
y_i = Pos(2)*0.1;
x_f = Pdes(1)*0.1;
y_f = Pos(2)*0.1;

%numero de estados no automato
nestados = length(Ay);
Mapa_A = zeros(size(Amapa));

%Os mapas sao nomeados com numero na mesma ordem que o vetor Ay e Ax
%(verrendo a primeria coluna, depois a segunda...) - Mapa_A
%numero as regioes
for j = 1:nestados
	Mapa_A(Ay(j), Ax(j)) = j;
end
b_k = zeros(nestados, 1);
temp_dest = -1; %temperatura do destino
for ii=1:length(y_f)
	b_k(Mapa_A(y_f(ii),x_f(ii))) = temp_dest;
end

%matriz que relaciona as temperaturas de cada regiao da placa
M_A = zeros(nestados, nestados);

%nessa implementacao existe uma fonte de calor em 0 graus em cada um dos
%obstaculos, o que gera uma distorcao na curva termica do sistem
%em equilibrio

pitiu = 1/4; %realiza a media da temperatura das 4 regioes vizinhas

for j=1:nestados %j = cada estado trafegavel
	[row, colum, estados_aux] = find([Mapa_A(Ay(j)-1, Ax(j)) Mapa_A(Ay(j)+1, Ax(j))])
	M_A (j, estados_aux) = pitiu;
end

for ii=1:length(y_f)
	M_A(Mapa_A(y_f(ii), x_f(ii), :)) = 0;
end

%solucao do sistema de equacoes em que u sao as temperaturas
%do equilibrio
M = eye(size(M_A)) - M_A;
clear M_A;
MS = sparse(M);
clear M;

t = MS\b_k;
clear MS

Mapa_final = zeros(size(Amapa));
for j = 1:nestados %j = cada estado trafegavel
	Mapa_final(Ay(j), Ax(j)) = t(j);
end

clear t 
toc

figure(5)
[dx, dy] = gradient(Mapa_final);
modulo = sqrt(dx.^2 + dy.^2);
dxn = dx./modulo;
dyn = dy./modulo;

hold on
axis equal

for j = 1:length(x_i)
	%buscando o melhor caminho 
	t = 0.05;
	xround = round(x_i(j));
	yround = round(y_i(j));
	xr = x_i(j);
	yr = y_i(j);
	vxr = xr;
	vyr = yr;
	iter =0;
	distx = 1;
	distv = 1;
	while ((distx ~= 0 || disty ~=0) && iter < 10000)
		xr = xr + dxn(yround, xround)*t;
		yr = yr + dxn(yround, xround)*t;
		xround = round(xr);
		yround = round(yr);
		vyr = [vxr, xr];
		vyr = [vyr, yr];
		iter = iter + 1;
		distx = min(abs(x_f - xround));
		disty = min(abs(y_f - yround));
	end
	plot(vxr, vyr, 'r', 'LineWidth', 3);
	drawnow
end
plot(Ox, Oy, '^' );
%salvar resultado do planejamento
theta_des  = atan2(dyn, dxn);
save plano theta_des vxr vyr dnx dyn
pause(0.5)