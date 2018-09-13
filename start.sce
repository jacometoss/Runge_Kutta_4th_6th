getd .;
disp('Run : runrk(4) o runrk(6) ')
function runrk(alg)
if(alg == 4)
[t,r]=rk4(ti,tf,muestras,varIni);
plot(t,[r(1,:)',r(2,:)'])
legend('Xdot(1)','Xdot(2)')
xlabel('tiempo, seg')
ylabel('Y')
title('Solución Runge Kutta orden 4')
xgrid
elseif ( alg == 6)
[t,r]=rk6(ti,tf,muestras,varIni);
plot(t,[r(1,:)',r(2,:)'])
legend('Xdot(1)','Xdot(2)')
xlabel('tiempo, seg')
ylabel('Ysol')
title('Solución Runge Kutta orden 6')
xgrid
else
disp('Error: Seleccione un método de solución ');
end
endfunction 
