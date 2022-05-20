# **Métodos de integración Runge-Kutta  y Regla Trapezoidal**

El método de Runge-Kutta (RK) es un conjunto de métodos iterativos (implícitos y explícitos) para la aproximación de soluciones de ecuaciones diferenciales ordinarias, en esta documentación se muestra el de cuarto y sexto orden. Otro método agregado es la  Regla Trapezoidal que puede ser explicita, para diferenciar los métodos de solución explicito e implícito  se agregan los dos puntos siguientes:


```tex
Autor : Marco Polo Jácome Toss	
Fecha de creación :  12 de Septiembre del 2017
Licencia : GNU General Public License (GPL) 3.0
Plataforma : Scilab
Código fuente creado para dar solución numérica a múltiples ecuaciones diferenciales.
Actualizado : 04 de Abril del 2022
```

***Método de integración explicito***. En este método es posible calcular la aproximación en cada paso directamente evaluando la función f(x,y), como ejemplo de un método de integración explicito se muestra  la regla del punto medio.



<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}^{*}=y_{i}&plus;hf(x_{i},y_{i})" title="y_{i&plus;1}^{*}=y_{i}&plus;hf(x_{i},y_{i})" />

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}=y_{i}&plus;\frac{h}{2}\left&space;(&space;f(x_{i},y_{i})&plus;f(x_{i}&plus;h,y_{i&plus;1}^{*})&space;\right&space;)" title="y_{i+1}=y_{i}+\frac{h}{2}\left ( f(x_{i},y_{i})+f(x_{i}+h,y_{i+1}^{*}) \right )" />

Donde 

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}^{*}=y_{i}&plus;hf\left&space;(x_{i},y_{i}&space;\right&space;)" title="y_{i+1}^{*}=y_{i}+hf\left (x_{i},y_{i} \right )" />

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}^{}=y_{i}&plus;\frac{h}{2}\left&space;(&space;f\left&space;(x_{i},y_{i}&space;\right&space;)&plus;f\left&space;(x_{i}&plus;h,y_{i&plus;1}^{*}&space;\right&space;)&space;\right&space;),&space;0\leq&space;i\leq&space;n-1" title="y_{i+1}^{}=y_{i}+\frac{h}{2}\left ( f\left (x_{i},y_{i} \right )+f\left (x_{i}+h,y_{i+1}^{*} \right ) \right ), 0\leq i\leq n-1" />



***Método de integración implicito***. Las aproximaciones en este método vienen definidas por un sistema de ecuaciones implícito. En la siguiente ecuación se muestra la forma implícita del método de integración del punto medio.

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}=y_{i}&plus;hf\left&space;(&space;x_{i}&space;&plus;&space;\frac{h}{2},&space;\frac{y_{i}&plus;y_{i&plus;1}}{2}&space;\right&space;),&space;0\leq&space;i\leq&space;n-1" title="y_{i+1}=y_{i}+hf\left ( x_{i} + \frac{h}{2}, \frac{y_{i}+y_{i+1}}{2} \right ), 0\leq i\leq n-1" />


La ventaja del "scripting" es poder dar solución en forma ordenada a **n** número de ecuaciones diferenciales.

## Problemas de valor inicial 

Se considera el problema de valor inicial (PVI) de ecuaciones diferenciales ordinarias (EDO) :

<center><img src="https://latex.codecogs.com/gif.latex?\left.\begin{matrix}&space;y^{'}=f(x,y(x)))\\&space;y\left&space;(&space;x_{0}&space;\right&space;)=y_{0}&space;\end{matrix}\right\}&space;y,f\in&space;\mathbb{R}^{m},x\in&space;\left&space;[&space;x_{0},x_{N}&space;\right]" title="\left.\begin{matrix} y^{'}=f(x,y(x)))\\ y\left ( x_{0} \right )=y_{0} \end{matrix}\right\} y,f\in \mathbb{R}^{m},x\in \left [ x_{0},x_{N} \right]" /></center>

Lo anterior es bastante importante debido a que debemos indicar las condiciones iniciales para poder dar solución a una EDO.

## 1.1 ¿Por qué utilizar este código fuente ?

Debido a su simplicidad de resolver varias ecuaciones de estado es flexible para dar solución a modelos de  máquinas eléctricas  como para otros tipos de modelos que dependan principalmente del tiempo o en caso contrario que no involucre esta variable.

Las líneas de código siguientes pertenecen al archivo `start.sce` el cual contiene tres métodos de solución Rk4,Rk3 y Trapezoidal ya antes mencionados.

```scilab
getd .;
op=input('Seleccione la forma de solucion : RK4(1) / RK6(2) / RTRAPEZOIDAL(3) : ')
mnecudif(op)
```


![GitHub Logo](https://image.ibb.co/jpL5qU/1.jpg)

##  1.2 Lista de archivos dependientes ![VERSION](https://img.shields.io/badge/Scilab-6.0.2-lightgrey)

La siguiente lista de archivos (dependientes) muestra como se encuentra estructurado en forma general el programa y siempre debe ejecutar el  archivo `start.sce`  para poder llamar la lista restante de archivos con extensión `*.sci` .

1. Start
   * attributes.sci
   * ecuDif.sci
   * rk4.sci (Multi ecuaciones)
   * rk6.sci (Multi ecuaciones)
   * rtrapezoidal (Multi ecuaciones) 
   * mnecudif


##  1.3 Ejecución del código fuente

De manera sencilla iniciamos el *`Start.sce`* e inmediatamente se indica el tipo de método a ocupar.

```scilab
**Seleccione la forma de solución : RK4(1) / RK6(2) / RTRAPEZOIDAL(3) : **
```

##  1.4 Ecuaciones diferenciales : ecuDif.sci

En el archivo `ecuDif` es necesario especificar las ecuaciones de estado, por lo que `Xdot` contiene dos ecuaciones de estado por este motivo se tiene dos variables `Xdot(1,1)` y `Xdot(2,1)`. Estas variables se pasan a los archivos `rk4.sci` y `rk6.sci` junto con la variable de tiempo para poder dar solución a la EDO. 

En este archivo se introducen las constantes y valores necesarios de las ecuaciones de estado con sus respectivas variables de estado representados con  `Xdot`, todo es un arreglo matricial.

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

##  1.5 Método de Runge-Kutta 6to Orden : rk6.sci

Es necesario cambiar la línea del archivo rk4 o rk6 que contiene  **`disp([t(i),r(1,i)',r(2,i)'])`** al incrementar las variables y ecuaciones o únicamente  mostrar una al igual que las variables iniciales mostradas en el archivo **`attributes.sci`**.


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

##  1.6 Método de Runge Kuta 4to Orden : rk4.sci

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

##  1.7 Método Trapezoidal : rtrapezoidal.sci

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

##  1.8 Valores Iniciales de Simulación : attributes.sci

El archivo contiene los atributos de la simulación y debe establecer el número de condiciones iniciales, las cuales dependen de las variables de estado a resolver siendo necesario modificar la impresión  **`disp([t(i),r(1,i)',r(2,i)'])`**  de los archivos rk4 y rk6.

- `ti`: tiempo inicial de simulación.
- `tf`: tiempo final de simulación.
- `varIni`: variables iniciales.
- `muestras`: muestras.

Los valores iniciales mostrados en el bloque siguiente  son el tiempo inicial, tiempo final, valores iniciales en este caso para dos variables y el número de muestras.

```scilab   
global ti tf varIni muestras    
ti=0;
tf=10;
varIni=[0.01, 0.02]; //Cambiar #Xdot
muestras=10000;
```

##  1.9 Opciones de solución: mnecudif.sci

Las opciones mencionadas en la parte **1.1** pertenecen al código siguiente :

```scilab   
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
```

##  Ejemplo #1 Ecuación del péndulo

Un péndulo simple se define como una partícula de masa "m" suspendida "O" por un hilo inextensible de longitud "L" y de masa despreciable.

Aplicando la segunda Ley de Newton,  siendo la aceleración de la partícula hacia adentro de la trayectoria circular se obtiene:

<img src="https://latex.codecogs.com/gif.latex?ma_{n}=T-mg\cdot&space;\cos&space;\theta" title="ma_{n}=T-mg\cdot \cos \theta" />

Siendo la aceleración de la partícula hacia adentro dada por:

<img src="https://latex.codecogs.com/gif.latex?a_{n}=v^{2}/L" title="a_{n}=v^{2}/L" />

La aceleración tangencial de la partícula es:

<img src="https://latex.codecogs.com/gif.latex?a_{t}=dv/dt" title="a_{t}=dv/dt" />


De la segunda Ley de Newton incluyendo la fricción se obtiene:

<img src="https://latex.codecogs.com/gif.latex?ml\ddot{\theta}=-mg\sin&space;\theta&space;-kl\dot{\theta}" title="ml\ddot{\theta}=-mg\sin \theta -kl\dot{\theta}" />

El valor de k (fricción) esta dado por :

<img src="https://latex.codecogs.com/gif.latex?k=\frac{b}{mL}" title="k=\frac{b}{mL}" />


Las variables de estado son :

<img src="https://latex.codecogs.com/gif.latex?x_{1}=\theta,&space;x_{2}=\frac{d\theta}{dt}" title="x_{1}=\theta, x_{2}=\frac{d\theta}{dt}" />

Por lo tanto, las ecuaciones de estado son las siguientes:

<img src="https://latex.codecogs.com/gif.latex?\dot{x_{1}}=x_{2}" title="\dot{x_{1}}=x_{2}" />

<img src="https://latex.codecogs.com/gif.latex?\dot{x_{2}}=-\frac{g}{l}\sin&space;x_{1}&space;-&space;\frac{k}{m}x_{2}" title="\dot{x_{2}}=-\frac{g}{l}\sin x_{1} - \frac{k}{m}x_{2}" />

Estableciendo la ecuaciones en el archivo **`ecuDif.sci`** y valores correspondientes.

```scilab
function [Xdot]=ecuDif(t,x)
  m=0.5;
  b=0.1;
  L=1.5;
  g=9.81; 
  k=b/(m*L);
 Xdot=zeros(2,1); //#Xdot
 Xdot(1,1)=x(2);
 Xdot(2,1)=-(g/L)*sin(x(1))-(k/m)*x(2);
endfunction
```

La solución para`Xdot(1)` y `Xdot(2)` de las ecuaciones de estado es:

![Simulacion No.1](https://i.ibb.co/DQwzZC4/2020-06-29-00-41-49.jpg)

## Ejemplo #2 Solución analítica y numérica  

Utilizando el método de la transformada de Laplace se resolverá la ecuación diferencial siguiente:
$$ \frac{dy}{dx}+3y = 13Sin 2t \:\: \forall \:\: y\left ( 0 \right )=6 $$
Aplicando la transformada de Laplace a la Ecuación diferencial anterior se obtiene la solución para `Y(S)`.
$$ Y\left (S \right )=\frac{6s^{2}+50}{\left ( s+3 \right )\left ( s^{2}+4 \right )} $$
La solución en el dominio del tiempo se obtiene con la transformada inversa.
$$ y\left ( t \right )=8e^{-3t}-2Cos2t+3Sin2t $$
El resultado anterior se grafica con Scilab mediante el bloque siguiente:

```Scilab
t=0:0.01:1;
y=8*exp(-3*t)-2*cos(2*t)+3*sin(2*t)
plot(t,y)
xlabel("Tiempo, seg.")
ylabel('y(t)')
title('Solución de la Ecuación Diferencial')
xgrid()
```

La solución de la ecuación diferencial y gráfica se muestra a continuación:
$$ y\left ( t \right )=8e^{-3t}-2Cos2t+3Sin2t $$

![Solución ED](https://i.ibb.co/zhBRgWj/Grafica-ED-Solucion.png)

Para dar solución es necesario configurar el archivo `ecuDif` de la forma siguiente:

```scilab
function [Xdot]=ecuDif(t,x)
 Xdot=zeros(1,1); //#Xdot
 Xdot(1,1)=13*sin(2*t)-3*x(1)
endfunction
```

Modificamos las condiciones iniciales del archivo `attributes` como se muestra:

```scilab
global ti tf varIni muestras    
ti=0;
tf=1;
varIni=[6,0]; //Cambiar #Xdot
muestras=1000;
```

Se debe ajustar el archivo `mnecudif` para mostrar únicamente una solución, esta ecuación diferencial tiene una única variable.

```scilab
function mnecudif(alg)
if(alg == 1) then
[t,r]=rk4(ti,tf,muestras,varIni);
plot(t,[r(1,:)])
legend('Xdot(1)')
xlabel('tiempo, seg')
ylabel('Y')
title('Solución Runge Kutta orden 4')
xgrid
elseif ( alg == 2) then
[t,r]=rk6(ti,tf,muestras,varIni);
plot(t,[r(1,:)'])
legend('Xdot(1)')
xlabel('tiempo, seg')
ylabel('Ysol')
title('Solución Runge Kutta orden 6')
xgrid
elseif (alg == 3) then
[t,r]=rtrapezoidal(ti,tf,muestras,varIni);
plot(t,[r(1,:)'])
legend('Xdot(1)')
xlabel('tiempo, seg')
ylabel('Y')
title('Solución Regla Trapezoidal Explicita')
xgrid
else
disp('Error: Seleccione un método de solución ');
end
endfunction 
```

Finalmente al realizar estas modificaciones podrás observar la solución siguiente:

![Solución EDO](https://i.ibb.co/TvznhKN/Grafica-ED-Solucion-2.png)

## Copyright

Copyright © 2017 en adelante, Marco Polo Jácome Toss (https://jacometoss.github.io/Runge_Kutta_4th_6th/). Este programa es software libre: usted puede redistribuirlo y /o modificarlo bajo los términos de la Licencia General GNU (GNU General Public License) publicado por la Fundación para el Software Libre para la versión 3 de dicha Licencia o anterior, o cualquier versión posterior.

Este programa se distribuye con la esperanza de que sea útil pero sin ninguna garantía; incluso sin la garantía implícita de comercialización o idoneidad para  un propósito en particular.

Vea la información de Licencia de `RK4` para más detalle.
