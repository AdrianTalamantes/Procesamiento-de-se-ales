# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 21:10:36 2021

@author: AdrianTR
"""

import scipy.io.wavfile as rd
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt

#Abre la señal ECG (base de datos de physionet)
file = open('ecg.txt',"r")
openedFile = file.read().replace('\n', ' ')
data_str = openedFile.split(' ')
datos_num = list(map(int,data_str))
#cantidad de datos en datos_num
array_length = len(datos_num)
#arreglos para almacenar la señal original y filtrada
filtered_signal = np.zeros(int(array_length/2))
raw_signal = np.zeros(int(array_length/2))
#separar las señales
i=0
j=0
while(j<array_length):
    raw_signal[i] = datos_num[j]
    filtered_signal[i]=datos_num[j+1]
    j+=2
    i+=1
#Graficas
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(raw_signal)
plt.title('Señal de ECG pura y filtrada')
plt.subplot(2,1,2)
plt.plot(filtered_signal)

######### Analisis de frecuencia de la señal ECG
#numero de muestras de mi señal
N1=len(raw_signal)
#Frecuencia de muestreo de la señal en Hertz
fs1 = 500 #puede cambiarse a mas
#periodo de muestre
ts1 =1/fs1
#Vector de frecuencia
k1=np.arange(0,N1)
f1 = k1*fs1/N1
#Determina TFD
Xf1 = np.abs(np.fft.fft(raw_signal))
#mitad de datos
M1 = int(np.fix(N1/2))
#Grafico
plt.figure(2)
plt.plot(f1[0:M1],Xf1[0:M1])
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')

#### Filtro pasa baja analogico
#frecuencia central del  filtro rechaza banda en Hertz
fo1=10
#frecuencia central en rad/seg
wo1=2*np.pi*fo1
#Factor de calidad
#Q=50
#ancho de banda del filtro
#B=fo/Q
#ancho de banda en rad/seg
#Br=2*np.pi*B
#Numerador de la funcion de transferencia en "s"
num1 = [wo1]
#Denominador de la funcion de transferencia en "s"
den1 = [1,wo1]
#Respuesta en frecuencia del filtro
w1,Hw1 = signal.freqs(num1,den1)
#converte w en Hertz
fhz1 = w1/(2*np.pi)
#grafico
plt.figure(3)
plt.plot(fhz1,np.abs(Hw1))
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')

###### Filtro rechazo de banda digital IIR
#Numerador de la funcion de tranferencia en "z"
numz1 = [wo1,wo1]
#denominador de la funcion de transferencia en "z"
denz1 = [(2/ts1)+wo1,(-2/ts1)+wo1]#talves le falta
#respuesta en frecuencia del filtro
Wz1,Hwz1 = signal.freqz(numz1,denz1)
#Convierte Wz en una frecuencia digital F
F1 = Wz1/(2*np.pi)
#Converte F digital en Hertz
fh1 = F1*fs1
#grafico
plt.figure(3)
plt.plot(fh1,np.abs(Hwz1))
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')

#######################################################################
#######################################################################
######### Analisis de frecuencia de la señal ECG
#numero de muestras de mi señal
N2=len(raw_signal)
#Frecuencia de muestreo de la señal en Hertz
fs2 = 500 #puede cambiarse a mas
#periodo de muestre
ts2 =1/fs2
#Vector de frecuencia
k2=np.arange(0,N2)
f2 = k2*fs2/N2
#Determina TFD
Xf3 = np.abs(np.fft.fft(raw_signal))
#mitad de datos
M2 = int(np.fix(N2/2))
#Grafico
plt.figure(2)
plt.plot(f2[0:M2],Xf3[0:M2])
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')

#### Filtro rechazo de banda analogico
#frecuencia central del  filtro rechaza banda en Hertz
fo2=52.1
#frecuencia central en rad/seg
wo2=2*np.pi*fo2
#Factor de calidad
Q2=50
#ancho de banda del filtro
B2=fo2/Q2
#ancho de banda en rad/seg
Br2=2*np.pi*B2
#Numerador de la funcion de transferencia en "s"
num2 = [1,0,wo2**2]
#Denominador de la funcion de transferencia en "s"
den2 = [1,Br2,wo2**2]
#Respuesta en frecuencia del filtro
w2,Hw2 = signal.freqs(num2,den2)
#converte w en Hertz
fhz2 = w2/(2*np.pi)
#grafico
plt.figure(3)
plt.plot(fhz2,np.abs(Hw2))
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')

###### Filtro rechazo de banda digital IIR
#Numerador de la funcion de tranferencia en "z"
numz2 = [4/ts2**2+wo2**2,-8/ts2**2+2*wo2**2,4/ts2**2+wo2**2]
#denominador de la funcion de transferencia en "z"
denz2 = [4/ts2**2+Br2*2/ts2+wo2**2,-8/ts2**2+2*wo2**2,4/ts2**2-2*Br2/ts2+wo2**2]#talves le falta
#respuesta en frecuencia del filtro
Wz2,Hwz2 = signal.freqz(numz2,denz2)
#Convierte Wz en una frecuencia digital F
F2 = Wz2/(2*np.pi)
#Converte F digital en Hertz
fh2 = F2*fs2
#grafico
plt.figure(3)
plt.plot(fh2,np.abs(Hwz2))
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')
####################################################################
####################################################################

#### Filtrado de la señal de ECG
#filtrado de la señal
filt_signal1 = signal.lfilter(numz1,denz1,raw_signal)
filt_signal = signal.lfilter(numz2,denz2,filt_signal1)
#Grafico
plt.figure(4)
plt.subplot(2,1,1)
plt.plot(raw_signal)
plt.title('Señal original vs filtrada')
plt.subplot(2,1,2)
plt.plot(filt_signal)
plt.xlabel('Muestras')

##### Analisis de frecuencia de la señal filtrada
#determina TFD
Xf2 = np.abs(np.fft.fft(filt_signal))
#Grafico
plt.figure(2)
plt.plot(f1[0:M2],Xf2[0:M2])
plt.title('Espectro de magnitud')
plt.xlabel('Frecuencia (Hz)')

######################################
#Filtro analogico pasa altas
#Filtro digital pasa altas
#Espectro original y filtrado
#señal original vs filtrada
