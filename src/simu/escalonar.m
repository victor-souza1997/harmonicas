
function [Mapa2] = escalonar(Mapa, percentage_reduction)
  dim = size(Mapa)*percentage_reduction

  %dim = size(Mapa)*0.5
  dx = dim(1)
  dy = dim(2)
  for j = 0:1/percentage_reduction -1
    for i = 0:1/percentage_reduction-1
      Mapa(1+(j)*dy:(j+1)*dy,1+(i)*dx:(i+1)*dx)
    end 
  end
  
end

