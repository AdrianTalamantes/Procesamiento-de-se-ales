# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 11:04:46 2021

@author: AdrianTR
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline

#####################################################
#    Genera señal x(t) a partir de 3 funciones senoidale
#Frecuencia en Hertz (Hz) de la función x1
f1 = 15
#Frecuencia en Hertz (Hz) de la función x2
f2 = 42
#Frecuencia en Hertz (Hz) de la función x3
f3 = 120
#Ancho de banda de la señal B = fmax - fmin = 120-0=120Hz
B = f3
#Frecuencia de muestreo
fs = 100*B
#Periodo de muestreo
ts= 1/fs
#Vector de tiempos k*ts
k = np.arange(0,2000) #genera un arreglo con enteros desde 0 hasta 1999
t = k*ts
#Función senoidal x1, las amplitudes son al azar
x1=0.63*np.sin(2*np.pi*f1*t)
#Función senoidal x2
x2=0.85*np.sin(2*np.pi*f2*t)
#Función senoidal x3
x3=0.14*np.sin(2*np.pi*f3*t)
#Suma de senoides
xt=x1 +x2 +x3
#Graficos
plt.subplot(4,1,1)#4 renglones, separados por 1, en 1 columna
#titulo
plt.title('Señal Original y sus componentes')
plt.plot(t,x1)
plt.subplot(4,1,2)#4 renglones, separados por 1, en 1 columna
plt.plot(t,x2)
plt.subplot(4,1,3)#4 renglones, separados por 1, en 1 columna
plt.plot(t,x3)
plt.subplot(4,1,4)#4 renglones, separados por 1, en 1 columna
plt.plot(t,xt)
plt.xlabel('Tiempo (seg.)')

############Muestreo de la señal x(t)##########
#Tamaño del arreglo xt
L = len(xt)
#Extrae el tiempo máximo de "t"
tmax = t[L-1]
#Define nueva frecuencia de muestreo, 2B es critico, mayor es sobremuestreo, y menor es submuestreo.
fs_new = 4*B
#Número de muestras
N = int(tmax*fs_new/1.0)
#Nuevo vector de tiempos
tn = np.linspace(0, tmax, N)
#Función senoidal x1, las amplitudes son al azar
xn1=0.63*np.sin(2*np.pi*f1*tn)
#Función senoidal x2
xn2=0.85*np.sin(2*np.pi*f2*tn)
#Función senoidal x3
xn3=0.14*np.sin(2*np.pi*f3*tn)
#Suma de senoides (Señal muestreada)
xn=xn1 +xn2 +xn3
#Graficos
plt.figure(2)
plt.subplot(4,1,1)#4 renglones, separados por 1, en 1 columna
#titulo
plt.title('Señal Original y sus componentes')
plt.plot(tn,xn1)
plt.subplot(4,1,2)#4 renglones, separados por 1, en 1 columna
plt.plot(tn,xn2)
plt.subplot(4,1,3)#4 renglones, separados por 1, en 1 columna
plt.plot(tn,xn3)
plt.subplot(4,1,4)#4 renglones, separados por 1, en 1 columna
plt.stem(tn,xn)
plt.plot(tn,xn)
plt.xlabel('Tiempo (seg.)')

#### Reconstrucción de la señal, método Whittaker-Shannon######
#Nuevo periodo de muestreo (ts)
ts_new = 1/fs_new
#sen(a)/a, factor de tiempo
d = 10
#Incremento de tiempo
dt = ts_new/d
#Vector de tiempos
tr = np.arange(0,tmax, dt)
#Inicializar vector para la señal continua x(t) reconstruida
xr = np.zeros(len(tr))
#Evalúa la ecuaciónde Whittaker-Shannon
for n in range(0,N):
    xr = xr +xn[n]*np.sin(np.pi*(fs_new*tr-n))/(np.pi*(fs_new*tr-n))
    
for n in range(0,N):
    xr[n*d]=xn[n]

#Gráficos
plt.figure(3)
plt.title("Señal original vs señal reconstruida")
plt.plot(t,xt,tr,xr)
plt.xlabel("Tiempo (seg.)")


###############################
#Tarea hacer la reconstruccion de la señal con 2 metodos de interpolacion
#sugerencia:Minimos cuadrados y spline cubico natural
#Google, buscar los 2 metodos, que librerias y funciones se necesitan 
#tambien puede ser Fourier, beta spline, etc.
#Señal discreta xn, vector de tiempos discreto tn y vector de tiempos continuo tr
#Salida del metodo de interpolacion es xr.
#aprox 10 lineas de codigo (incluyendo el grafico)
#Encimar con t, xt, para checar cuanto se parece la interpolacion con la original
#Titulo de la grafica, reconstruccio  con "metodo usado"
#entregar un word, con el codigo agregado, las graficas y un comentario de lo que pasa u observemos, que metood trabajo mejor, etc. Podriamos agregar el muestreo critico, sobremuestreo y submuestreo.


######### SPLINE CUBICO
xrSc = CubicSpline(tn,xn, bc_type='natural')
#Gráficos
plt.figure(4)
plt.title("Señal original vs señal reconstruida SC")
plt.plot(t,xt,tr,xrSc(tr))
plt.xlabel("Tiempo (seg.)")

######### B-SPLINE
xrBs = BSpline(tn, xn, 4)
#Gráficos
plt.figure(5)
plt.title("Señal original vs señal reconstruida Bs")
plt.plot(t,xt,tr,xrBs(tr))
plt.xlabel("Tiempo (seg.)")

#########Akima
from scipy.interpolate import Akima1DInterpolator
xrAk = Akima1DInterpolator(tn, xn)
#Gráficos
plt.figure(6)
plt.title("Señal original vs señal reconstruida Ak")
plt.plot(t,xt,tr,xrAk(tr))
plt.xlabel("Tiempo (seg.)")
