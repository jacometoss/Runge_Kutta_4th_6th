# Runge Kutta 4to y 6to orden

Los métodos de Runge-Kutta (RK) son un conjunto de métodos iterativos (implícitos y explícitos) para la aproximación de soluciones de ecuaciones diferenciales ordinarias.


##  Archivos 
1. Star
    * attributes.sci
    * ecuDif.sci
    * rk4.sci (Multi ecuaciones)
    * rk6.sci (Multi ecuaciones)


##  Ejecución

De manera sencilla iniciamos el Start.sce 

** runrk(4) para ejecucar Runge-Kutta 4to Orden  o  runrk(6) para ejecutar Runge-Kutta 4to Orden **

##  Archivo : ecuDif.sci

En el archivo ecuDif es necesario especificar las ecuaciones de estado. De esta manera Xdot contiene dos ecuaciones 
```scilab
function [Xdot]=ecuDif(t,x)
  m=0.5;
  b=0.1;
  L=0.5;
  g=9.81;
  k=b/(m*L);
 Xdot=zeros(2,1); //#Xdot
 Xdot(1,1)=x(2);
 Xdot(2,1)=-(g/L)*sin(x(1))-(k/m)*x(2);
endfunction

```
