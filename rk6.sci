function [t,r]=rk6(t0, tf, N, conIni)
    matrixSize= size(conIni,1)
    if  (matrixSize==1) then
        conIni=conIni';
    end
h=(tf-t0)/N; //Paso
t(1) = t0;
r(:,1) = conIni; //Condiciones Iniciales
for i = 1:N  //Seccion Modificable para matriz
     disp([t(i),r(1,i)',r(2,i)']) // Cambiar #Xdot
   k1 = h*ecuDif(t(i), r(:,i));
   k2 = h*ecuDif(t(i)+h, r(:,i)+k1);
   k3 = h*ecuDif(t(i)+h/2, r(:,i)+(3*k1+k2)/8);  
   k4 = h*ecuDif(t(i)+(2*h)/3, r(:,i)+(8*k1+2*k2+8*k3)/27); 
   k5 = h*ecuDif(t(i)+((7-sqrt(21))*h)/14,r(:,i)+((3*k1*(3*sqrt(21)-7))-(8*k2*(7-sqrt(21)))+(48*k3*(7-sqrt(21)))-(3*k4*(21-sqrt(21))))/392);
   k6 = h*ecuDif(t(i)+(7+sqrt(21))*h/14, r(:,i)+(-5*k1*(231+51*sqrt(21))-40*k2*(7+sqrt(21))-320*k3*sqrt(21)+3*k4*(21+121*sqrt(21))+392*k5*(6+sqrt(21)))/1960); 
   r(:,i+1) = r(:,i)+(15*k1*(22+7*sqrt(21))+120*k2+40*k3*(7*sqrt(21)-5)-63*k4*(3*sqrt(21)-2)-14*k5*(49+9*sqrt(21))+70*k6*(7-sqrt(21)))/180;
   t(i +1) = t0 + i*h;
end
[t,r'];     
endfunction
