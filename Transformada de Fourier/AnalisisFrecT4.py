# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 22:36:16 2021

@author: AdrianTR
"""

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import scipy.io.wavfile as rd
import sys
from scipy import interpolate


######## Abrir la señal de ECG (Base de datos en Physionet)
file = open('datosZoo.txt',"r")
openedFile = file.read().replace('\n',' ')
data_str = openedFile.split(' ')
datos_num = list(map(float,data_str))
#Cantidad de datos en datos_num
array_length = len(datos_num)
#Arreglos para almacenar la señal original y la señal filtrada
#filtered_signal = np.zeros(int(array_length/2))
raw_signal = np.zeros(int(array_length))

#Separar señales
i=0
j=0
while(j< array_length):
    raw_signal[i] = datos_num[j]
    #filtered_signal[i] = datos_num[j+1]
    j += 1
    i += 1

#graficas
plt.figure(1)
#plt.subplot(2,1,1)
plt.plot(raw_signal)
plt.title('Señal de ECG pura de PhysioZoo de un perro')

#plt.subplot(2,1,2)
#plt.plot(filtered_signal)
###### Extraccion de un ciclo cardiaco
#Extrae ciclo cardiaco
#signal = filtered_signal[5900:6300] #para usar la otra, ponemos signal = raw_signal[5900:6300], para usar la que no es filtrada, es decir, la pura
signal = raw_signal[3:213]
# Cantidad de muestras en la señal
M = len(signal)
#Condicion de periodicidad
signal[M-1] = signal[0]
#Frecuencia de muestreo en Hertz Hz
fs = 500
#periodo de muestreo
ts = 1/fs
#vector de tiempos
t=np.arange(0,M*ts,ts)

#Tamaño del vector de tiempos
L=len(t)
#Grafico
plt.figure(2)
plt.plot(t,signal) #signal es el xt
plt.title('Ciclo cardiaco, con interpolacion')
plt.xlabel('Tiempo (min)')

###############################################################

############################################
#       Analisis en frecuencia
#Especificar Numero de muestras en la señal (siempre impar)
N = 91
#Vecto de tiempos para evaluar la interpolacion
te = np.linspace(0,t[-1],N)
#Aplica el metodo de interpolacion
inter  = interpolate.splrep(t, signal)#Argumentos, vector de tiempos y señal
yn = interpolate.splev(te,inter)
#Grafico
plt.plot(te,yn,'.r')
#Numero de trapezoides
Nt = N-1 
#Factor de las sumatorias
fact = 2/Nt
#Incremento en radianes
Inc = 2*np.pi/Nt
#Vector de angulos
phi = np.linspace(0,2*np.pi-Inc,Nt)
#Calculo del coeficiente a0
a0 = fact*np.sum(yn[0:Nt])
#Numero de armonicos
K = int(((N-1)/2)+1)#Componentes antes de que se empiece a repetir la informacion
#K=10
#Periodo de la señal
T = M*ts
#Frecuencia fundamental
fo = 1/T
#Inicializar arreglos para los coeficientes de Fourier
ak = np.zeros(2*K+1)
bk = np.zeros(2*K+1)
ck = np.zeros(2*K+1)
phik = np.zeros(2*K+1)
fk = np.zeros(2*K+1)

#Asignamos el valor de a0 a ck, c[0]=a0
ck[0] = np.abs(a0)
ak[0] = a0
#Calcula los coeficientes de Fourier
for k in range(1,K):
    #Coef rectangulares
    ak[k] = fact*np.sum(yn[0:Nt]*np.cos(k*phi))
    bk[k] = fact*np.sum(yn[0:Nt]*np.sin(k*phi))
    #Coef polares
    ck[k] = np.sqrt(ak[k]**2 + bk[k]**2)

    if ak[k] >= 0 and bk[k] >= 0:
        phik[k] = np.arctan(bk[k]/ak[k])
    if ak[k] < 0 and bk[k] >= 0:
        phik[k] = np.pi - np.arctan(bk[k]/np.abs(ak[k]))
    if ak[k] < 0 and bk[k] < 0:
        phik[k] = np.pi + np.arctan(bk[k]/np.abs(ak[k]))
    if ak[k] >= 0 and bk[k] < 0:
        phik[k] = -np.arctan(np.abs(bk[k])/np.abs(ak[k]))

    #Frecuencia de los coeficientes de Fourier
    fk[k] = k*fo
    phik[k] = phik[k]*180/np.pi
#Grafica de los espectros de la señal
plt.figure(3)
plt.subplot(2,1,1)
plt.stem(fk,ak)
plt.title('Espectro de la señal usando ak y bk')
plt.subplot(2,1,2)
plt.stem(fk,bk)
plt.xlabel('Frecuencia en Hertz')

plt.figure(4)
plt.subplot(2,1,1)
plt.stem(fk,ck)#Se utiliza para mostrar el analisis en frecuencia
plt.title('Espectro de la señal usando ck y phik')
plt.subplot(2,1,2)
plt.stem(fk,phik*180/np.pi)
plt.xlabel('Frecuencia en Hertz')

############# Reconstruccion de la señal
re_signal = np.zeros(M)
#Vector de angulos
rad_vec = np.linspace(0,2*np.pi,M)
for k in range(1,K):
    re_signal = re_signal + ak[k]*np.cos(k*rad_vec) + bk[k]*np.sin(k*rad_vec)
    
#Grafico
plt.figure(5)
plt.title('Señal original y reconstruida')
plt.plot(t,signal,t,re_signal)#para que se paresca mas, podemos aumentar el numero de muestras N
#Con N pocas muestras, es un submuestreo(traslape) es decir una perdida de informacion.


###############################################
### Transformada continua de Fourier
#Consigue valor de frecuencias mas alto
f_max = fk[K-1]
#Vector de frecuencias
fc = np.arange(0,f_max, 0.1)
#Evaluamos la transformada de Fourier  X(f)
Xf = np.sqrt(1-np.cos(2*np.pi*fc))/(np.sqrt(2)*np.pi*fc)
#Grafico
plt.figure(6)
plt.plot(fc,Xf,'r')
#plt.stem(fk[0:K+1],ck[0:K+1])
plt.title('Espectro de Magnitud')
plt.xlabel('Frecuencia de hertz')

##################################################
###### Transformada de Fourier de tiempo discreto

#Grafico de la señal x[kts]
plt.figure(7)
plt.plot(t,signal,t,signal,'.r')
plt.title('Señal Zoo discreta')
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
    x_real = x_real + signal[k]*np.cos(2*np.pi*k*F)
    x_imag = x_imag + signal[k]*np.sin(2*np.pi*k*F)
    
#Involucramos a ts en la sumatoria
x_real = x_real*ts
x_imag = x_imag*ts
#Calculamos modulo de X(f)  (|X(f)|)
Xf_mod = np.sqrt(x_real**2 + x_imag**2)
#Grafico de la señal x[kts]
plt.figure(8)
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
        real = real + signal[k]*np.cos(2*np.pi*k*n/L)
        imag = imag + signal[k]*np.sin(2*np.pi*k*n/L)
    p_real[n] = real
    p_imag[n] = imag
    fn[n] = (n*fs)/L
    
#Involucramos a ts en la sumatoria
p_real = p_real*ts
p_imag = p_imag*ts
#Calculamos modulo de X(f)  (|X(f)|)
Xn_mod = np.sqrt(p_real**2 + p_imag**2)
#Grafico de la señal x[kts]
plt.figure(9)
plt.plot(fn,Xn_mod)
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia en Hertz')
plt.ylabel('Magnitud')
