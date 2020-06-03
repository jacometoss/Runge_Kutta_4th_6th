function mnecudif(alg)
if(alg == 1) then
[t,r]=rk4(ti,tf,muestras,varIni);
plot(t,[r(1,:)',r(2,:)'])
legend('Xdot(1)','Xdot(2)')
xlabel('tiempo, seg')
ylabel('Y')
title('Solución Runge Kutta orden 4')
xgrid
elseif ( alg == 2) then
[t,r]=rk6(ti,tf,muestras,varIni);
plot(t,[r(1,:)',r(2,:)'])
legend('Xdot(1)','Xdot(2)')
xlabel('tiempo, seg')
ylabel('Ysol')
title('Solución Runge Kutta orden 6')
xgrid
elseif (alg == 3) then
[t,r]=rtrapezoidal(ti,tf,muestras,varIni);
plot(t,[r(1,:)',r(2,:)'])
legend('Xdot(1)','Xdot(2)')
xlabel('tiempo, seg')
ylabel('Y')
title('Solución Regla Trapezoidal Explicita')
xgrid
else
disp('Error: Seleccione un método de solución ');
end
endfunction 
