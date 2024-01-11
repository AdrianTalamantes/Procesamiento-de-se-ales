# -*- coding: utf-8 -*-
"""
Created on Sat Sep  18 16:40:38 2021

@author: AdrianTR
"""

##############################################################################
#                   SISTEMAS CONTINUOS Y DISCRETOS

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import scipy.io.wavfile as rd

##############################################################################
#                          MODELO CONTINUO
#Frecuencia de corte en Hertz
fc = 1200
#Proponiendo el valor de la Capacitancia 
C = 1*10**(-6)
#Determina el valor de la resistencia
R = 1/(2*np.pi*fc*C)
#Numerador de la Función de Transferencia
num = [(R*C),0]
#Denominador de la Funcion de tranferencia
den = [(R*C),1]
#Función de Transferencia
Hs = signal.lti(num,den)
#Respuesta a la función escalón unitario
t, yt = signal.step2(Hs)
#Gráfico
plt.figure(1)
plt.plot(t,yt)
plt.xlabel('Tiempo')
plt.title('Respuesta al Escalón Unitario')

##############################################################################
#                           MODELO DISCRETO
#Tamaño del arreglo de tiempos "t"
T = len(t)
#Extrae el tiempo máximo
tmax = t[T-1]
#frecuencia de muestreo propuesta
fs = 100000
#Periodo de muestreo
ts = 1/fs
#Vector de tiempos
tm = np.arange(0,tmax,ts)
#Constantes  b del filtro pasa altas digital

b = (R*C)/(R*C + ts)
#Tamaño del vector tm
K = len(tm)
#Evaluar la ecuación de diferencias  
yk = np.zeros(K)
yk_1 = 0
xk_1 = 0
xk = 1

for k in range(0,K):
    yk[k] = b*(xk -xk_1 + yk_1)
    yk_1 = yk[k]
    xk_1 = 1
    
    
#Gráfico
plt.plot(tm,yk,'.r')    

##############################################################################
#                       PRUEBA CON SEÑAL DE AUDIO

#Abre archivo de audio
fs_new, audio = rd.read('song.wav')
#Separación de canales
ch_rigth = audio[:,0]
ch_left = audio[:,1]
#Periodo de muestreo
ts_new = 1/fs_new
#Constantes a y b del filtro pasa bajas digital
#a = ts/(R*C + ts)
b = (R*C)/(R*C + ts)
#Determina el número de muestras en la señal de audio
K = len(ch_left)
#Evaluar la ecuación de diferencias  y[k] = b*(ch_left[k] -ch_left[k-1] + ykl_1)
yk_left = np.zeros(K)
yk_rigth = np.zeros(K)
ykr_1 = 0
ykl_1 = 0
xkl_1 = 0

for k in range(0,K):
    #Filtra canal izquierdo
    yk_left[k] = b*(ch_left[k] -ch_left[k-1] + ykl_1)
    ykl_1 = yk_left[k]
    xkl_1 = ch_left[k]
    
    #Filtra canal derecho
    yk_rigth[k] = b*(ch_rigth[k] -ch_rigth[k-1] + ykr_1)
    ykr_1 = yk_rigth[k]
    xkr_1 = ch_rigth[k]
    
    
#Ajusta la señal filtrada en el rango de [-1,1]
yk_left = yk_left/np.max(yk_left)
yk_rigth = yk_rigth/np.max(yk_rigth)

#Genera señal filtrada
filtered = np.zeros((K,2))
filtered[:,0] =  yk_rigth
filtered[:,1] =  yk_left

#sd.play(filtered,fs_new)

#Graficar
plt.figure(2)
plt.subplot(2, 1, 1)
plt.plot(ch_rigth)
plt.title("Canal derecho original vs filtrado")
plt.subplot(2, 1, 2)
plt.plot(yk_rigth)

plt.figure(3)
plt.subplot(2, 1, 1)
plt.plot(ch_left)
plt.title("Canal izquierdo original vs filtrado")
plt.subplot(2, 1, 2)
plt.plot(yk_left)

############# CONCLUSIONES ###############################
"""
Observamos que la grafica de respuesta al escalon unitario con el filtro pasa altas
correspode con el modelo discreto que se realizo en la tarea #4


"""
