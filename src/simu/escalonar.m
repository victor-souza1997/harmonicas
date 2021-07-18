function [Mapa2] = escalonar(Mapa, percentage_reduction)
  dim = size(Mapa)*percentage_reduction%geting the size of the scan matrix
  dx = dim(1) %x dimention of scanmatrix
  dy = dim(2) %y dimetion "  "
  for j = 0:1/percentage_reduction -1
      for i = 0:1/percentage_reduction-1  
        if(find(Mapa(1+(j)*dy:(j+1)*dy,1+(i)*dx:(i+1)*dx)==0)>0), Mapa2(j+1,i+1) = 0; 
        else, Mapa2(j+1,i+1) = 1;
        end
        %Mapa(1+(j)*dy:(j+1)*dy,1+(i)*dx:(i+1)*dx); visualize the matrix being scanned
    end 
  end
end

