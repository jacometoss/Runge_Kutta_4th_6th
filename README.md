# **Métodos de integración Runge-Kutta  y Regla Trapezoidal**

​	El método de Runge-Kutta (RK) es un conjunto de métodos iterativos (implícitos y explícitos) para la aproximación de soluciones de ecuaciones diferenciales ordinarias. En esta documentación se muestra el de cuarto y sexto orden. Otro método agregado es la  Regla Trapezoidal que puede ser explicito. Para diferenciar los métodos de solución explicito e implicito  se agregan los dos puntos siguientes:

```tex
Autor : Marco Polo Jácome Toss	
Fecha de creación :  12 de Septiembre del 2017
Licencia : GNU General Public License (GPL) 3.0
Plataforma : Scilab
Código fuente creado para dar solución numérica a múltiples ecuaciones diferenciales.
```

***Método de integración explicito***. En este método es posible calcular la aproximación en cada paso directamente evaluando la función f(x,y) como ejemplo de un método de integración explicito se muestra en la siguiente imagen la regla del punto medio.

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}^{*}=y_{i}&plus;hf(x_{i},y_{i})" title="y_{i+1}^{*}=y_{i}+hf(x_{i},y_{i})" />

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}=y_{i}&plus;\frac{h}{2}\left&space;(&space;f(x_{i},y_{i})&plus;f(x_{i}&plus;h,y_{i&plus;1}^{*})&space;\right&space;)" title="y_{i+1}=y_{i}+\frac{h}{2}\left ( f(x_{i},y_{i})+f(x_{i}+h,y_{i+1}^{*}) \right )" />

Donde 

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}^{*}=y_{i}&plus;hf\left&space;(x_{i},y_{i}&space;\right&space;)" title="y_{i+1}^{*}=y_{i}+hf\left (x_{i},y_{i} \right )" />

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}^{}=y_{i}&plus;\frac{h}{2}\left&space;(&space;f\left&space;(x_{i},y_{i}&space;\right&space;)&plus;f\left&space;(x_{i}&plus;h,y_{i&plus;1}^{*}&space;\right&space;)&space;\right&space;),&space;0\leq&space;i\leq&space;n-1" title="y_{i+1}^{}=y_{i}+\frac{h}{2}\left ( f\left (x_{i},y_{i} \right )+f\left (x_{i}+h,y_{i+1}^{*} \right ) \right ), 0\leq i\leq n-1" />



​	***Método de integración implicito***. Las aproximaciones en este método vienen definidas por un sistema de ecuaciones implícito. En la siguiente imagen se muestra la forma implícita del método de integración del punto medio.

<img src="https://latex.codecogs.com/gif.latex?y_{i&plus;1}=y_{i}&plus;hf\left&space;(&space;x_{i}&space;&plus;&space;\frac{h}{2},&space;\frac{y_{i}&plus;y_{i&plus;1}}{2}&space;\right&space;),&space;0\leq&space;i\leq&space;n-1" title="y_{i+1}=y_{i}+hf\left ( x_{i} + \frac{h}{2}, \frac{y_{i}+y_{i+1}}{2} \right ), 0\leq i\leq n-1" />


La ventaja del "scripting" es poder dar solución en forma ordenada a **n** número de ecuaciones diferenciales.



## **Problemas de valor inicial** 

Se considera el problema de valor inicial (PVI) de ecuaciones diferenciales ordinarias (EDO) :

<center><img src="https://latex.codecogs.com/gif.latex?\left.\begin{matrix}&space;y^{'}=f(x,y(x)))\\&space;y\left&space;(&space;x_{0}&space;\right&space;)=y_{0}&space;\end{matrix}\right\}&space;y,f\in&space;\mathbb{R}^{m},x\in&space;\left&space;[&space;x_{0},x_{N}&space;\right]" title="\left.\begin{matrix} y^{'}=f(x,y(x)))\\ y\left ( x_{0} \right )=y_{0} \end{matrix}\right\} y,f\in \mathbb{R}^{m},x\in \left [ x_{0},x_{N} \right]" /></center>

## **1.1 ¿Por qué utilizar este código fuente ?**

Debido a su simplicidad de resolver varias ecuaciones de estado es flexible para dar solución a modelos de  máquinas eléctricas  como para otros tipos de modelos que dependan principalmente del tiempo por defecto o en caso contrario que no involucre esta variable.

Las líneas de código siguientes pertenecen al archivo `start.sce` el cual contiene tres métodos de solución Rk4,Rk3 y Trapezoidal.

```scilab
getd .;
op=input('Seleccione la forma de solucion : RK4(1) / RK6(2) / RTRAPEZOIDAL(3) : ')
mnecudif(op)
```


![GitHub Logo](https://image.ibb.co/jpL5qU/1.jpg)

##  1.2 Lista de archivos dependientes ![VERSION](https://img.shields.io/badge/Scilab-6.0.2-lightgrey)

La siguiente lista de archivos dependientes muestra como se encuentra estructurado en forma general siempre se arranca toda la lista de archivos con "start.sce" por o tanto el resto de archivos es llamado para su implementación.

1. Star
    * attributes.sci
    * ecuDif.sci
    * rk4.sci (Multi ecuaciones)
    * rk6.sci (Multi ecuaciones)
    * rtrapezoidal (Multi ecuaciones) 
    * mnecudif


##  1.3 Ejecución del código fuente

De manera sencilla iniciamos el Start.sce y se indica el tipo de método a ocupar en este caso se incluye también el método trapezoidal.

```
**Seleccione la forma de solución : RK4(1) / RK6(2) / RTRAPEZOIDAL(3) : **
```



##  1.4 Ecuaciones diferenciales : ecuDif.sci

En el archivo ecuDif es necesario especificar las ecuaciones de estado, por lo que Xdot contiene dos ecuaciones de estado las cuales son impresas en los archivos rk4.sci y rk6.sci con el tiempo. En este archivo se introducen los valores necesarios de las ecuaciones de estado con sus variables de estado representados con  `Xdot`.
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

El archivo contiene los atributos de a simulación y debe establecer el número de condiciones iniciales, las cuales dependen de las variables de estado a resolver y siendo necesario modificar la impresion disp() de los arhivos rk4 y rk6
- ti: tiempo inicial de simulación
- tf: tiempo final de simulación
- varIni: variables iniciales
- muestras: muestras, 

Los valores iniciales mostrados en la parte inferior son el tiempo inicial, tiempo final, valores iniciales para dos variables y el número de muestras.

```scilab   
global ti tf varIni muestras    
ti=0;
tf=10;
varIni=[0.01, 0.02]; //Cambiar #Xdot
muestras=10000;

```

##  Ejemplo #1 Ecuación del péndulo

Un péndulo simple se define como una partícula de masa "m" suspendida de "O" por un hilo inextensible de longitud "L" y de masa despreciable

De la segunda Ley de Newton y siendo la aceleración de la partícula hacia adentro de la trayectoria circular.

<img src="https://latex.codecogs.com/gif.latex?ma_{n}=T-mg\cdot&space;\cos&space;\theta" title="ma_{n}=T-mg\cdot \cos \theta" />

Siendo la aceleración de la partícula hacia adentro dada por 

<img src="https://latex.codecogs.com/gif.latex?a_{n}=v^{2}/L" title="a_{n}=v^{2}/L" />

La aceleración tangencial de la partícula

<img src="https://latex.codecogs.com/gif.latex?a_{t}=dv/dt" title="a_{t}=dv/dt" />


De la segunda Ley de Newton incluyendo la fricción.

<img src="https://latex.codecogs.com/gif.latex?ml\ddot{\theta}=-mg\sin&space;\theta&space;-kl\dot{\theta}" title="ml\ddot{\theta}=-mg\sin \theta -kl\dot{\theta}" />

El valor de k (fricción) esta dado por :

<img src="https://latex.codecogs.com/gif.latex?k=\frac{b}{mL}" title="k=\frac{b}{mL}" />


Las variables de estado son :

<img src="https://latex.codecogs.com/gif.latex?x_{1}=\theta,&space;x_{2}=\frac{d\theta}{dt}" title="x_{1}=\theta, x_{2}=\frac{d\theta}{dt}" />

Por lo tanto las ecuaciones de estado son las siguientes :

<img src="https://latex.codecogs.com/gif.latex?\dot{x_{1}}=x_{2}" title="\dot{x_{1}}=x_{2}" />

<img src="https://latex.codecogs.com/gif.latex?\dot{x_{2}}=-\frac{g}{l}\sin&space;x_{1}&space;-&space;\frac{k}{m}x_{2}" title="\dot{x_{2}}=-\frac{g}{l}\sin x_{1} - \frac{k}{m}x_{2}" />

Estableciendo la ecuaciones en el archivo  `ecuDif.sci` y valores correspondientes.

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

La imagen de salida de la simulación con los datos anteriores se muestra con `Xdot(1)` y `Xdot(2)` 

![Simulacion No.1](https://i.ibb.co/DQwzZC4/2020-06-29-00-41-49.jpg)
