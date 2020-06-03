![VERSION](https://img.shields.io/badge/Scilab-6.0.2-lightgrey)

# Runge Kutta 4to-6to orden y Regla Trapezoidal

​	Los métodos de Runge-Kutta (RK) son un conjunto de métodos iterativos (implícitos y explícitos) para la aproximación de soluciones de ecuaciones diferenciales ordinarias. En esta documentación se muestra el de cuarto y sexto orden. Otro método de integración es la Regla Trapezoidal explicito. Para dar solución es necesario establecer las ecuaciones diferenciales y seleccionar rk4(1), rk6(2), rtrapezoidal(3). 

​	La ventaja del "scripting" es poder dar solución en forma ordenada a **n** número de ecuaciones.

## ¿Por qué utilizar este código fuente ?

Debido a su simplicidad de resolver varias ecuaciones de estado, siendo flexible para dar solución a modelos de  máquinas eléctricas.

[Motor de Inducción](https://github.com/jacometoss/Motor_de_induccion)


![GitHub Logo](https://image.ibb.co/jpL5qU/1.jpg)


##  Archivos 
1. Star
    * attributes.sci
    * ecuDif.sci
    * rk4.sci (Multi ecuaciones)
    * rk6.sci (Multi ecuaciones)
    * rtrapezoidal (Multi ecuaciones) 
    * mnecudif


##  Ejecución

De manera sencilla iniciamos el Start.sce 

**Seleccione la forma de solución : RK4(1) / RK6(2) / RTRAPEZOIDAL(3) : **

##  Archivo : ecuDif.sci

En el archivo ecuDif es necesario especificar las ecuaciones de estado, por lo que Xdot contiene dos ecuaciones de estado las cuales son impresas en los archivos rk4.sci y rk6.sci con el tiempo.
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

##  Archivo : rk6.sci

Siendo necesario cambiar la línea del archivo rk4 o rk6 que contiene  disp([t(i),r(1,i)',r(2,i)']) al incrementar las ecuaciones o únicamente al mostrar una al igual que las variables iniciames mostradas en el arhivo **attributes.sci**


```scilab
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
```

##  Archivo : rk4.sci
```scilab
function [t,r]=rk4(t0, tf, N, conIni)
    matrixSize= size(conIni,1)
    if  (matrixSize==1) then
        conIni=conIni';
    end
h=(tf-t0)/N; //Paso
t(1) = t0;
r(:,1) = conIni; //Condiciones Iniciales
for i = 1:N  //Seccion Modificable para matriz
     disp([t(i),r(1,i)',r(2,i)']) //Cambiar #Xdot
   k1 = h*ecuDif(t(i), r(:,i));
   k2 = h*ecuDif(t(i)+h/2, r(:,i)+0.5*k1);
   k3 = h*ecuDif(t(i)+h/2, r(:,i)+0.5*k2); 
   k4 = h*ecuDif(t(i)+h, r(:,i)+k3);
   r(:,i+1) = r(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
   t(i +1) = t0 + i*h;
end
[t,r'];     
endfunction
```

##  Archivo : rtrapezoidal.sci

```scilab
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
```

##  Archivo : attributes.sci

El archivo contiene los atributos de a simulación y debe establecer el número de condiciones iniciales, las cuales dependen de las variables de estado a resolver y siendo necesario modificar la impresion disp() de los arhivos rk4 y rk6
   * ti: tiempo inicial de simulación
   * tf: tiempo final de simulación
   * varIni: variables iniciales
   * muestras: muestras, 
```scilab   
global ti tf varIni muestras    
ti=0;
tf=10;
varIni=[0.01, 0.02]; //Cambiar #Xdot
muestras=10000;

```

