# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 10:50:20 2021

@author: AdrianTR
"""

############# Sistemas continuos y discretos ######

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import scipy.io.wavfile as rd

############ Modelo continuo ############
#Fecuencia de corte
fc = 1200
#Proponiendo el valor de C
C = 1*10**(-6)
#Determine el valor de la resistencia
R = 1/(2*np.pi*fc*C)
#Numerador de la Funcion de tranferencia
num = [1/(R*C)]
#Denominador de la Funcion de tranferencia
den = [1,1/(R*C)]
#Funcion de transferencia
Hs = signal.lti(num,den)
#Respuesta a la funcion escalon unitario
t, yt = signal.step2(Hs) #aqui se aplica transformada inv al vector de tiempos
#Grafico
plt.figure(1)
plt.plot(t,yt)
plt.xlabel("Tiempo (seg)")
plt.title("Respuesta al Escalon unitario")

########### Modelo discreto ###############
#Tamaño del arreglo de tiempos "t"
T = len(t)
#Extrae el tiempo máximo
tmax = t[T-1]
#Frecuencia de muestreo (usando constante de tiempo thao = RC)
#fs = (R*C)/10000 <- no jalo, proponemos otra
fs = 10000
#Periodo de muestreo
ts = 1/fs
#Vector de tiempos
tm = np.arange(0,tmax,ts)
#Constantes a y b del filtro pasa bajas digital
a = ts/(R*C + ts)
b = (R*C)/(R*C + ts)
#Tamaño del vector tm
K = len(tm)
#Evaluar la ecuacion de diferencias y[k]=ax[k]+by[k-1]
yk = np.zeros(K)
yk_1 = 0
xk= 1

for k in range(0,K):
    yk[k] = a*xk + b*yk_1
    yk_1 = yk[k]
    
#Grafico
plt.plot(tm,yk,'.r')
################################################
# Prueba con señal de audio
#Abrir archivo de audio
fs_new, audio = rd.read('song.wav')
#Separacion de canales
ch_rigth = audio[:,0]
ch_left = audio[:,1]
#Periodo de muestreo (de nuevo:v), porque al abrir el archivo de audio, nos abrio una nueva fs
ts_new = 1/fs_new
#Constantes a y b del filtro pasa bajas digital
a = ts/(R*C + ts)
b = (R*C)/(R*C + ts)
#Determinar numero de muestras en la señal de audio
K = len(ch_left)
#Evaluar la ecuacion de diferencias y[k]=ax[k]+by[k-1]
yk_left = np.zeros(K)
yk_rigth = np.zeros(K)
ykr_1 = 0
ykl_1 = 0

for k in range(0,K):
    #Filtra canal izquierdo
    yk_left[k] = a*ch_left[k] + b*ykl_1
    ykl_1 = yk_left[k]
    #Filtra canal derecho
    yk_rigth[k] = a*ch_rigth[k] + b*ykr_1
    ykr_1 = yk_rigth[k]
#Ajuste la señal filtrada en el rango de [-1,1]
yk_left = yk_left/np.max(yk_left)
yk_rigth = yk_rigth/np.max(yk_rigth)
#Genera señal filtrada
filtered = np.zeros((K,2))
filtered[:,0] = yk_rigth
filtered[:,1] = yk_left
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

#Reproducir el sonido
#sd.play(filtered,fs_new)
