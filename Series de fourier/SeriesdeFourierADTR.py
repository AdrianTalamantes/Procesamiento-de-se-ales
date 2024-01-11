# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 20:13:42 2021

@author: AdrianTR
"""

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import scipy.io.wavfile as rd
import sys

### Definicion de la funcion x(t)

def h(x):
    return -x+1

#    x(t)=-t+1, 0<t<1 
#Frecuencia de muestreo en HZ
fs = 100
#Periodo de muestreo
ts=1/fs
#Vector de tiempos
t = np.arange(0,1,ts)#original es  0,2,0.001
#Tamaño del vector de tiempos
L=len(t)
#Vector de magnitudes
x = np.linspace(0,1,100)
xt = h(x)
#Grafico
plt.figure(1)
plt.plot(t,xt)
plt.title('Señal -t+1')
plt.xlabel('Tiempo (seg.)')
plt.ylabel('Amplitud')


#Analisis de frecuencia

#Periodo de la señal
T=1
#Frecuencia fundamental
fo= 1/T
#coeficiente a0 (offset)
a0 = 1/2
#Numero de armonicos (componentes de frecuencia)
K=10 #La grafica nos da aproximadamente 10 maximos (ondulaciones)
#Definir los vectores para ak y bk
ak = np.zeros(2*K+1)  #Ak siempre vale 0
bk = np.zeros(2*K+1)
#Vector para las frecuencias
f = np.zeros(2*K+1)
#Definir los vectores para ck y phik
ck = np.zeros(2*K+1)
phik = np.zeros(2*K+1)
#Calculo de los coeficientes de Fourier
for k in range(-K,K+1):
    if k==0:
        ak[k+K] = a0
        ck[k+K] = a0   #c0=a0
        phik[k+K] = 0
    else:
        ak[k+K] = sys.float_info.epsilon  #ak=2.2204e-16
        bk[k+K] = (1)/(np.pi*k)
        ck[k+K] = np.sqrt(ak[k+K]**2 + bk[k+K]**2)
        
        if ak[k+K] >= 0 and bk[k+K] >= 0:
            phik[k+K] = np.arctan(np.abs(bk[k+K])/np.abs(ak[k+K]))
        if ak[k+K] < 0 and bk[k+K] >= 0:
            phik[k+K] = np.pi - np.arctan(np.abs(bk[k+K])/np.abs(ak[k+K]))
        if ak[k+K] < 0 and bk[k+K] < 0:
            phik[k+K] =  np.pi + np.arctan(np.abs(bk[k+K])/np.abs(ak[k+K]))
        if ak[k+K] >= 0 and bk[k+K] < 0:
            phik[k+K] = - np.arctan(np.abs(bk[k+K])/np.abs(ak[k+K]))

    f[k+K] = k*fo
    phik[k+K] = phik[k+K]*180/np.pi
    
#grafico
plt.figure(2)
plt.subplot(2,1,1)
plt.stem(f,ak)
plt.title('Espectro de magnitud para ak y bk')
plt.subplot(2,1,2)
plt.stem(f,bk)
plt.xlabel('Frecuencia de Hertz')

plt.figure(3)
plt.subplot(2,1,1)
plt.stem(f,ck)
plt.title('Espectro de magnitud para ck y phik')
plt.subplot(2,1,2)
plt.stem(f,phik)
plt.xlabel('Frecuencia de Hertz')

####################################################
#                 Reconstruccion de la señal
#Arreglo para recosntruir x(t)
xtt = np.zeros(L)
#Evaluamos la serie de Fourier Rectangular
for k in range(0, K+1):
    xtt = xtt + ak[k+K]*np.cos(2*np.pi*k*fo*t) + bk[k+K]*np.sin(2*np.pi*k*fo*t)
#Grafico
plt.figure(4)
plt.plot(t,xt,t,xtt)
plt.title('Señal original y reconstruida')#EFecto Gibbs
plt.xlabel('Tiempo')

#Arreglo para recosntruir x(t) POLAR
xtp = np.zeros(L)
#Matriz para almacenar las cosenoides de la reconstruccion
senoide = np.zeros((K+1,L))
#Evaluamos la serie de Fourier polar
for k in range(0, K+1):
    xtp = xtp + ck[k+K]*np.cos(2*np.pi*k*fo*t + phik[k+K]*np.pi/180)
    senoide[k,:] = ck[k+K]*np.cos(2*np.pi*k*fo*t + phik[k+K]*np.pi/180)
#Grafico
plt.figure(5)
plt.plot(t,xt,t,xtp)
plt.title('Señal original y reconstruida')#EFecto Gibbs
plt.xlabel('Tiempo')

plt.figure(6)
plt.plot(t,senoide[1,:])
plt.plot(t,senoide[2,:])
plt.plot(t,senoide[3,:])
plt.plot(t,senoide[4,:])
plt.plot(t,senoide[5,:])
plt.title('Senoide utilizadas en la reconstruccion')#EFecto Gibbs
plt.xlabel('Tiempo')

###############################################
### Transformada continua de Fourier
#Consigue valor de frecuencias mas alto
f_max = f[-1]
#Vector de frecuencias
fc = np.arange(0,f_max, 0.001)
#Evaluamos la transformada de Fourier  X(f)
Xf = np.sqrt(1-np.cos(2*np.pi*fc))/(np.sqrt(2)*np.pi*fc)
#Grafico
plt.figure(7)
plt.plot(fc,Xf,'r')
plt.stem(f[K:K+K+1],ck[K:K+K+1])
plt.title('Espectro de Magnitud')
plt.xlabel('Frecuencia de hertz')

##################################################
###### Transformada de Fourier de tiempo discreto

#Grafico de la señal x[kts]
plt.figure(8)
plt.plot(t,xt,t,xt,'.r')
plt.title('Señal -t+1 reconstruida')
plt.xlabel('Tiempo (seg.)')
plt.ylabel('Amplitud')
#Ancho de banda (Considerando caso critico en el muestreo)
B=fs/2
#Vector de frecuencias en hertz
fz = np.arange(-fs-B,fs+B,0.001)
#Frecuencia normalizada
F = fz/fs
#inicializamos vectores de soporte para la TFTD
x_real = np.zeros(len(F))
x_imag = np.zeros(len(F))
#EValuamos las TFTD
for k in range(0,L):
    x_real = x_real + xt[k]*np.cos(2*np.pi*k*F)
    x_imag = x_imag + xt[k]*np.sin(2*np.pi*k*F)
    
#Involucramos a ts en la sumatoria
x_real = x_real*ts
x_imag = x_imag*ts
#Calculamos modulo de X(f)  (|X(f)|)
Xf_mod = np.sqrt(x_real**2 + x_imag**2)
#Grafico de la señal x[kts]
plt.figure(9)
plt.plot(fz,Xf_mod)
plt.title('Espectro de magnitud continuo y periodico')
plt.xlabel('Frecuencia en Hertz')
plt.ylabel('Magnitud')

##############################################
######## Transformada de Fourier discreta

p_real = np.zeros(L)
p_imag = np.zeros(L)
Xn = np.zeros(L)
fn = np.zeros(L)
#EValuamos las TFD
for n in range(0,L):
    real = 0
    imag = 0
    for k in range(0,L):
        real = real + xt[k]*np.cos(2*np.pi*k*n/L)
        imag = imag + xt[k]*np.sin(2*np.pi*k*n/L)
    p_real[n] = real
    p_imag[n] = imag
    fn[n] = (n*fs)/L
    
#Involucramos a ts en la sumatoria
p_real = p_real*ts
p_imag = p_imag*ts
#Calculamos modulo de X(f)  (|X(f)|)
Xn_mod = np.sqrt(p_real**2 + p_imag**2)
#Grafico de la señal x[kts]
plt.figure(10)
plt.plot(fn,Xn_mod)
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia en Hertz')
plt.ylabel('Magnitud')

#Empalma resultado de la TFD con lo anterior analizado
k=0

for n in range(0,L):
    if fn[n] > 5:
        break
    else:
        k = k+1

plt.figure(7)        
plt.plot(fn[0:k],Xn_mod[0:k],'k')
