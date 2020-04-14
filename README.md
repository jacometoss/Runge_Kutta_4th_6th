# Runge Kutta 4to y 6to orden

Los métodos de Runge-Kutta (RK) son un conjunto de métodos iterativos (implícitos y explícitos) para la aproximación de soluciones de ecuaciones diferenciales ordinarias.

![GitHub Logo](https://image.ibb.co/jpL5qU/1.jpg)


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

##  Archivo : attributes.sci
El archivo contiene los atributos de a simulación y debe establecer el número de condiciones iniciales, las cuales dependen de las variables de estado a resolver y siendo necesario modificar la impresion disp() de los arhivos rk4 y rk6
   * ti: tiempo inicial de simulación
   * tf: tiempo final de simulación
   * varIni: variables iniciales
   * muestras: muestras
```scilab   
global ti tf varIni muestras    
ti=0;
tf=10;
varIni=[0.01, 0.02]; //Cambiar #Xdot
muestras=10000;

```

{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Minimal Example pandapower\n",
    "\n",
    "\n",
    "## Creating a Power System\n",
    "\n",
    "We consider the following simple 3-bus example network as a minimal example:\n",
    "\n",
    "<img src=\"pics/3bus-system.png\" width=\"50%\">\n",
    "\n",
    "The above network can be created in pandapower as follows:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandapower as pp\n",
    "\n",
    "#create empty net\n",
    "net = pp.create_empty_network()\n",
    "\n",
    "#create buses\n",
    "bus1 = pp.create_bus(net, vn_kv=20., name=\"Bus 1\")\n",
    "bus2 = pp.create_bus(net, vn_kv=0.4, name=\"Bus 2\")\n",
    "bus3 = pp.create_bus(net, vn_kv=0.4, name=\"Bus 3\")\n",
    "\n",
    "#create bus elements\n",
    "pp.create_ext_grid(net, bus=bus1, vm_pu=1.02, name=\"Grid Connection\")\n",
    "pp.create_load(net, bus=bus3, p_mw=0.100, q_mvar=0.05, name=\"Load\")\n",
    "\n",
    "#create branch elements\n",
    "trafo = pp.create_transformer(net, hv_bus=bus1, lv_bus=bus2, std_type=\"0.4 MVA 20/0.4 kV\", name=\"Trafo\")\n",
    "line = pp.create_line(net, from_bus=bus2, to_bus=bus3, length_km=0.1, std_type=\"NAYY 4x50 SE\", name=\"Line\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Data Structure\n",
    "\n",
    "Each dataframe in a pandapower net object contains the information about one pandapower element, such as line, load transformer etc."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>name</th>\n",
       "      <th>vn_kv</th>\n",
       "      <th>type</th>\n",
       "      <th>zone</th>\n",
       "      <th>in_service</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>Bus 1</td>\n",
       "      <td>20.0</td>\n",
       "      <td>b</td>\n",
       "      <td>None</td>\n",
       "      <td>True</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>Bus 2</td>\n",
       "      <td>0.4</td>\n",
       "      <td>b</td>\n",
       "      <td>None</td>\n",
       "      <td>True</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>Bus 3</td>\n",
       "      <td>0.4</td>\n",
       "      <td>b</td>\n",
       "      <td>None</td>\n",
       "      <td>True</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "    name  vn_kv type  zone  in_service\n",
       "0  Bus 1   20.0    b  None        True\n",
       "1  Bus 2    0.4    b  None        True\n",
       "2  Bus 3    0.4    b  None        True"
      ]
     },
     "execution_count": 8,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.bus"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>name</th>\n",
       "      <th>std_type</th>\n",
       "      <th>from_bus</th>\n",
       "      <th>to_bus</th>\n",
       "      <th>length_km</th>\n",
       "      <th>r_ohm_per_km</th>\n",
       "      <th>x_ohm_per_km</th>\n",
       "      <th>c_nf_per_km</th>\n",
       "      <th>g_us_per_km</th>\n",
       "      <th>max_i_ka</th>\n",
       "      <th>df</th>\n",
       "      <th>parallel</th>\n",
       "      <th>type</th>\n",
       "      <th>in_service</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>Line</td>\n",
       "      <td>NAYY 4x50 SE</td>\n",
       "      <td>1</td>\n",
       "      <td>2</td>\n",
       "      <td>0.1</td>\n",
       "      <td>0.642</td>\n",
       "      <td>0.083</td>\n",
       "      <td>210.0</td>\n",
       "      <td>0.0</td>\n",
       "      <td>0.142</td>\n",
       "      <td>1.0</td>\n",
       "      <td>1</td>\n",
       "      <td>cs</td>\n",
       "      <td>True</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   name      std_type  from_bus  to_bus  length_km  r_ohm_per_km  \\\n",
       "0  Line  NAYY 4x50 SE         1       2        0.1         0.642   \n",
       "\n",
       "   x_ohm_per_km  c_nf_per_km  g_us_per_km  max_i_ka   df  parallel type  \\\n",
       "0         0.083        210.0          0.0     0.142  1.0         1   cs   \n",
       "\n",
       "   in_service  \n",
       "0        True  "
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.line"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>name</th>\n",
       "      <th>std_type</th>\n",
       "      <th>hv_bus</th>\n",
       "      <th>lv_bus</th>\n",
       "      <th>sn_mva</th>\n",
       "      <th>vn_hv_kv</th>\n",
       "      <th>vn_lv_kv</th>\n",
       "      <th>vk_percent</th>\n",
       "      <th>vkr_percent</th>\n",
       "      <th>pfe_kw</th>\n",
       "      <th>...</th>\n",
       "      <th>tap_neutral</th>\n",
       "      <th>tap_min</th>\n",
       "      <th>tap_max</th>\n",
       "      <th>tap_step_percent</th>\n",
       "      <th>tap_step_degree</th>\n",
       "      <th>tap_pos</th>\n",
       "      <th>tap_phase_shifter</th>\n",
       "      <th>parallel</th>\n",
       "      <th>df</th>\n",
       "      <th>in_service</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>Trafo</td>\n",
       "      <td>0.4 MVA 20/0.4 kV</td>\n",
       "      <td>0</td>\n",
       "      <td>1</td>\n",
       "      <td>0.4</td>\n",
       "      <td>20.0</td>\n",
       "      <td>0.4</td>\n",
       "      <td>6.0</td>\n",
       "      <td>1.425</td>\n",
       "      <td>1.35</td>\n",
       "      <td>...</td>\n",
       "      <td>0</td>\n",
       "      <td>-2</td>\n",
       "      <td>2</td>\n",
       "      <td>2.5</td>\n",
       "      <td>0.0</td>\n",
       "      <td>0</td>\n",
       "      <td>False</td>\n",
       "      <td>1</td>\n",
       "      <td>1.0</td>\n",
       "      <td>True</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "<p>1 rows × 23 columns</p>\n",
       "</div>"
      ],
      "text/plain": [
       "    name           std_type  hv_bus  lv_bus  sn_mva  vn_hv_kv  vn_lv_kv  \\\n",
       "0  Trafo  0.4 MVA 20/0.4 kV       0       1     0.4      20.0       0.4   \n",
       "\n",
       "   vk_percent  vkr_percent  pfe_kw     ...      tap_neutral  tap_min tap_max  \\\n",
       "0         6.0        1.425    1.35     ...                0       -2       2   \n",
       "\n",
       "   tap_step_percent  tap_step_degree  tap_pos  tap_phase_shifter  parallel  \\\n",
       "0               2.5              0.0        0              False         1   \n",
       "\n",
       "    df  in_service  \n",
       "0  1.0        True  \n",
       "\n",
       "[1 rows x 23 columns]"
      ]
     },
     "execution_count": 10,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.trafo"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>name</th>\n",
       "      <th>bus</th>\n",
       "      <th>p_mw</th>\n",
       "      <th>q_mvar</th>\n",
       "      <th>const_z_percent</th>\n",
       "      <th>const_i_percent</th>\n",
       "      <th>sn_mva</th>\n",
       "      <th>scaling</th>\n",
       "      <th>in_service</th>\n",
       "      <th>type</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>Load</td>\n",
       "      <td>2</td>\n",
       "      <td>0.1</td>\n",
       "      <td>0.05</td>\n",
       "      <td>0.0</td>\n",
       "      <td>0.0</td>\n",
       "      <td>NaN</td>\n",
       "      <td>1.0</td>\n",
       "      <td>True</td>\n",
       "      <td>None</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   name  bus  p_mw  q_mvar  const_z_percent  const_i_percent  sn_mva  scaling  \\\n",
       "0  Load    2   0.1    0.05              0.0              0.0     NaN      1.0   \n",
       "\n",
       "   in_service  type  \n",
       "0        True  None  "
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.load"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Note that line and transformer are created with standard types, so thath the electric parameters of are automatically filled in from the standard type library."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Power Flow\n",
    "\n",
    "We now run a power flow:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [],
   "source": [
    "pp.runpp(net)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "And check out at the results for buses, lines an transformers:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>vm_pu</th>\n",
       "      <th>va_degree</th>\n",
       "      <th>p_mw</th>\n",
       "      <th>q_mvar</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1.020000</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>-0.107265</td>\n",
       "      <td>-0.052675</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>1.008843</td>\n",
       "      <td>-0.760126</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>0.000000</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>0.964431</td>\n",
       "      <td>0.115859</td>\n",
       "      <td>0.100000</td>\n",
       "      <td>0.050000</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "      vm_pu  va_degree      p_mw    q_mvar\n",
       "0  1.020000   0.000000 -0.107265 -0.052675\n",
       "1  1.008843  -0.760126  0.000000  0.000000\n",
       "2  0.964431   0.115859  0.100000  0.050000"
      ]
     },
     "execution_count": 13,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_bus"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>p_from_mw</th>\n",
       "      <th>q_from_mvar</th>\n",
       "      <th>p_to_mw</th>\n",
       "      <th>q_to_mvar</th>\n",
       "      <th>pl_mw</th>\n",
       "      <th>ql_mvar</th>\n",
       "      <th>i_from_ka</th>\n",
       "      <th>i_to_ka</th>\n",
       "      <th>i_ka</th>\n",
       "      <th>vm_from_pu</th>\n",
       "      <th>va_from_degree</th>\n",
       "      <th>vm_to_pu</th>\n",
       "      <th>va_to_degree</th>\n",
       "      <th>loading_percent</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>0.105392</td>\n",
       "      <td>0.050696</td>\n",
       "      <td>-0.1</td>\n",
       "      <td>-0.05</td>\n",
       "      <td>0.005392</td>\n",
       "      <td>0.000696</td>\n",
       "      <td>0.167325</td>\n",
       "      <td>0.167326</td>\n",
       "      <td>0.167326</td>\n",
       "      <td>1.008843</td>\n",
       "      <td>-0.760126</td>\n",
       "      <td>0.964431</td>\n",
       "      <td>0.115859</td>\n",
       "      <td>117.835208</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   p_from_mw  q_from_mvar  p_to_mw  q_to_mvar     pl_mw   ql_mvar  i_from_ka  \\\n",
       "0   0.105392     0.050696     -0.1      -0.05  0.005392  0.000696   0.167325   \n",
       "\n",
       "    i_to_ka      i_ka  vm_from_pu  va_from_degree  vm_to_pu  va_to_degree  \\\n",
       "0  0.167326  0.167326    1.008843       -0.760126  0.964431      0.115859   \n",
       "\n",
       "   loading_percent  \n",
       "0       117.835208  "
      ]
     },
     "execution_count": 14,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_line"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>p_hv_mw</th>\n",
       "      <th>q_hv_mvar</th>\n",
       "      <th>p_lv_mw</th>\n",
       "      <th>q_lv_mvar</th>\n",
       "      <th>pl_mw</th>\n",
       "      <th>ql_mvar</th>\n",
       "      <th>i_hv_ka</th>\n",
       "      <th>i_lv_ka</th>\n",
       "      <th>vm_hv_pu</th>\n",
       "      <th>va_hv_degree</th>\n",
       "      <th>vm_lv_pu</th>\n",
       "      <th>va_lv_degree</th>\n",
       "      <th>loading_percent</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>0.107265</td>\n",
       "      <td>0.052675</td>\n",
       "      <td>-0.105392</td>\n",
       "      <td>-0.050696</td>\n",
       "      <td>0.001873</td>\n",
       "      <td>0.001979</td>\n",
       "      <td>0.003382</td>\n",
       "      <td>0.167325</td>\n",
       "      <td>1.02</td>\n",
       "      <td>0.0</td>\n",
       "      <td>1.008843</td>\n",
       "      <td>-0.760126</td>\n",
       "      <td>29.289513</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "    p_hv_mw  q_hv_mvar   p_lv_mw  q_lv_mvar     pl_mw   ql_mvar   i_hv_ka  \\\n",
       "0  0.107265   0.052675 -0.105392  -0.050696  0.001873  0.001979  0.003382   \n",
       "\n",
       "    i_lv_ka  vm_hv_pu  va_hv_degree  vm_lv_pu  va_lv_degree  loading_percent  \n",
       "0  0.167325      1.02           0.0  1.008843     -0.760126        29.289513  "
      ]
     },
     "execution_count": 15,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_trafo"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Tap Changers\n",
    "\n",
    "We now lower the tap changer position, from position 0 to -1 and run another power flow:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {},
   "outputs": [],
   "source": [
    "net.trafo.tap_pos.at[trafo] = -1\n",
    "pp.runpp(net)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Looking at the results shows that bus voltages at the low voltage side of the transformer have increased:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>vm_pu</th>\n",
       "      <th>va_degree</th>\n",
       "      <th>p_mw</th>\n",
       "      <th>q_mvar</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1.020000</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>-0.107015</td>\n",
       "      <td>-0.052529</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>1.035301</td>\n",
       "      <td>-0.720245</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>0.000000</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>0.992135</td>\n",
       "      <td>0.109513</td>\n",
       "      <td>0.100000</td>\n",
       "      <td>0.050000</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "      vm_pu  va_degree      p_mw    q_mvar\n",
       "0  1.020000   0.000000 -0.107015 -0.052529\n",
       "1  1.035301  -0.720245  0.000000  0.000000\n",
       "2  0.992135   0.109513  0.100000  0.050000"
      ]
     },
     "execution_count": 17,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_bus"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Switches\n",
    "\n",
    "We now create an open switch at the load bus:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "0"
      ]
     },
     "execution_count": 18,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pp.create_switch(net, bus=bus3, element=line, et=\"l\", closed=False)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The open switch cuts the load bus from power supply:\n",
    "\n",
    "<img src=\"pics/3bus-system_switch.png\" width=\"8%\">"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This can be verified by running a power flow and inspecting the results. The voltage at bus 2 is given as NaN:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>vm_pu</th>\n",
       "      <th>va_degree</th>\n",
       "      <th>p_mw</th>\n",
       "      <th>q_mvar</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1.020000</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>-0.001477</td>\n",
       "      <td>0.000001</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>1.046129</td>\n",
       "      <td>-0.005637</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>0.000000</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>NaN</td>\n",
       "      <td>NaN</td>\n",
       "      <td>0.000000</td>\n",
       "      <td>0.000000</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "      vm_pu  va_degree      p_mw    q_mvar\n",
       "0  1.020000   0.000000 -0.001477  0.000001\n",
       "1  1.046129  -0.005637  0.000000  0.000000\n",
       "2       NaN        NaN  0.000000  0.000000"
      ]
     },
     "execution_count": 19,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "pp.runpp(net)\n",
    "net.res_bus"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The load does not feed in:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>p_mw</th>\n",
       "      <th>q_mvar</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>0.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "   p_mw  q_mvar\n",
       "0   0.0     0.0"
      ]
     },
     "execution_count": 20,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_load"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "And the line is in open loop operation:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>p_from_mw</th>\n",
       "      <th>q_from_mvar</th>\n",
       "      <th>p_to_mw</th>\n",
       "      <th>q_to_mvar</th>\n",
       "      <th>pl_mw</th>\n",
       "      <th>ql_mvar</th>\n",
       "      <th>i_from_ka</th>\n",
       "      <th>i_to_ka</th>\n",
       "      <th>i_ka</th>\n",
       "      <th>vm_from_pu</th>\n",
       "      <th>va_from_degree</th>\n",
       "      <th>vm_to_pu</th>\n",
       "      <th>va_to_degree</th>\n",
       "      <th>loading_percent</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1.225734e-13</td>\n",
       "      <td>-0.000001</td>\n",
       "      <td>-4.645803e-16</td>\n",
       "      <td>-5.802601e-17</td>\n",
       "      <td>1.221088e-13</td>\n",
       "      <td>-0.000001</td>\n",
       "      <td>0.000002</td>\n",
       "      <td>6.459759e-16</td>\n",
       "      <td>0.000002</td>\n",
       "      <td>1.046129</td>\n",
       "      <td>-0.005637</td>\n",
       "      <td>1.046129</td>\n",
       "      <td>-0.005649</td>\n",
       "      <td>0.001122</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "      p_from_mw  q_from_mvar       p_to_mw     q_to_mvar         pl_mw  \\\n",
       "0  1.225734e-13    -0.000001 -4.645803e-16 -5.802601e-17  1.221088e-13   \n",
       "\n",
       "    ql_mvar  i_from_ka       i_to_ka      i_ka  vm_from_pu  va_from_degree  \\\n",
       "0 -0.000001   0.000002  6.459759e-16  0.000002    1.046129       -0.005637   \n",
       "\n",
       "   vm_to_pu  va_to_degree  loading_percent  \n",
       "0  1.046129     -0.005649         0.001122  "
      ]
     },
     "execution_count": 21,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_line"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Topological Analysis\n",
    "\n",
    "The structure of the network can also be directly analyzed with the topology package. It uses an interface to the NetworkX library for graph searches. There are some predefined search algorithms, such as searching for unsupplied buses:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "{2}"
      ]
     },
     "execution_count": 22,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "import pandapower.topology as top\n",
    "top.unsupplied_buses(net)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The package correctly determines that bus 2 is cut from power supply. When we close the switch, there are no unsupplied buses anymore:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "set()"
      ]
     },
     "execution_count": 23,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.switch.closed.at[0] = True\n",
    "top.unsupplied_buses(net)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Apart from predefined search functions, it is also possible to translate the pandapower network into a NetworkX graph and run searches directly on that graph.\n",
    "\n",
    "Suppose we want to find all buses that are on the same voltage level as the load bus. We then translate the grid into a graph but excluding the transformer:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 24,
   "metadata": {},
   "outputs": [],
   "source": [
    "mg = top.create_nxgraph(net, include_trafos=False)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "And search for all buses that are connected to the load bus in that graph:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 25,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[2, 1]"
      ]
     },
     "execution_count": 25,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "list(top.connected_component(mg, 2))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The graph search finds all buses that are on the same voltage level. Searches like these can be used for feeder identification and many more applications."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Short Circuit Analysis\n",
    "\n",
    "pandapower includes a short circuit module that complies with IEC 60909. To run a short circuit analysis, we need to define short circuit parameters for the external grid:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 26,
   "metadata": {},
   "outputs": [],
   "source": [
    "net.ext_grid[\"s_sc_max_mva\"] = 100\n",
    "net.ext_grid[\"rx_max\"] = 0.1"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Now we can calculate short circuits. Here, we calculate a three phase short circuit current with a fault impedance of 2 Ohms:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandapower.shortcircuit as sc\n",
    "sc.calc_sc(net, case=\"max\", ip=True, r_fault_ohm=2.)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Initial and peak short circuit currents are given for faults at all buses:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 28,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>ikss_ka</th>\n",
       "      <th>ip_ka</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>2.534707</td>\n",
       "      <td>4.317318</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>0.126631</td>\n",
       "      <td>0.182666</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>0.122698</td>\n",
       "      <td>0.176991</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "    ikss_ka     ip_ka\n",
       "0  2.534707  4.317318\n",
       "1  0.126631  0.182666\n",
       "2  0.122698  0.176991"
      ]
     },
     "execution_count": 28,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "net.res_bus_sc"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This concludes a short walkthrough of some pandapower features. More in-depth tutorials can be found in the pandapower documentation:\n",
    "https://www.pandapower.org/start/#interactive-tutorials-"
   ]
  }
 ],
 "metadata": {
  "anaconda-cloud": {},
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
   
