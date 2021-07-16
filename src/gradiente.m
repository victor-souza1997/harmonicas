function [dx, dy] = gradiente(Mapa)
  [ny, nx] = size(Mapa);
  for j=1:ny
    for i=1:nx
      if j == 1
          dy(j,i) = Mapa(j+1,i) - Mapa(j, i);
      else if j == ny
          dy(j,i) = Mapa(j,i) - Mapa(j-1,i);
          else
            dy(j,i) = (Mapa(j+1,i) - Mapa(j-1,i))/2;
          end
      end
      if i == 1
        dx(j,i) = Mapa(j,i+1) - Mapa(j,i);
       else if i==nx
         dx(j,i) = Mapa(j,i) - Mapa(j,i-1);
       else
         dx(j,i) = (Mapa(j,i+1) - Mapa(j, i-1))/2;
       end
       
      end
           
    end
  end
  
  
  
endfunction