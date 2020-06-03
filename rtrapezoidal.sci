function [t,r]=rtrapezoidal(t0, tf, N, conIni)
    matrixSize= size(conIni,1)
    if  (matrixSize==1) then
        conIni=conIni';
    end
h=(tf-t0)/N; //Paso
t(1) = t0;
r(:,1) = conIni; //Condiciones Iniciales
for i = 1:N  //Seccion Modificable para matriz
     disp([t(i),r(1,i)',r(2,i)']) //Cambiar #Xdot
   y1 = r(:,i)+h*ecuDif(t(i), r(:,i));
   y2 = r(:,i)+(h/2)*((ecuDif(t(i), r(:,i)))+(ecuDif(t(i)+h, y1)));
   r(:,i+1) = y2;
   t(i +1) = t0 + i*h;
end
[t,r'];     
endfunction
