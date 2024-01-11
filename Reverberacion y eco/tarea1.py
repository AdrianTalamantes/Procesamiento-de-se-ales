# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 16:19:56 2021

@author: AdrianTR
"""

#######   EFECTO DE REVERBERACIÓN Y ECO

import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import scipy.io.wavfile as rd


############               PRUEBA CON SEÑAL DE AUDIO
#Abre archivo de audio
fs_new, audio = rd.read('song.wav')
#Separación de canales
ch_right = audio[:,0]
ch_left = audio[:,1]
#Periodo de muestreo
ts_new = 1/fs_new
#Determina el número de muestras en la señal de audio
K = len(ch_left)
# Factor de retraso
D = 10000
alpha = 0.8

############                ECO
#Señal con echo
yk_left_echo = np.zeros(K+D)
yk_right_echo = np.zeros(K+D)
xk_left_echo = np.concatenate((np.zeros(D),ch_left))
xk_right_echo = np.concatenate((np.zeros(D),ch_right))
#Ecuación de diferencias
for k in range(D,K+D):
    #Filtra canal izquierdo
    yk_left_echo[k] = xk_left_echo[k] + alpha*xk_left_echo[k-D]
    #Filtra canal derecho
    yk_right_echo[k] = xk_right_echo[k] + alpha*xk_right_echo[k-D]    
#Normalizar señal filtrada en el rango de [-1,1]
yk_left_echo = yk_left_echo/np.max(yk_left_echo)
yk_right_echo = yk_right_echo/np.max(yk_right_echo)
#Señal filtrada
filtered_echo = np.zeros((K+D,2))
filtered_echo[:,0] =  yk_right_echo
filtered_echo[:,1] =  yk_left_echo

#Gráfico
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(ch_right)
plt.title('Canal Derecho original vs filtrado')
plt.subplot(2,1,2)
plt.plot(yk_right_echo)
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(ch_left)
plt.title('Canal izquierdo original vs filtrado')
plt.subplot(2,1,2)
plt.plot(yk_left_echo)

#sd.play(filtered_echo,fs_new)


############                 REVERBERACIÓN

#Señal con reverberacion
yk_left_rev = np.zeros(K+D)
yk_right_rev = np.zeros(K+D)
xk_left_rev = ch_left
xk_right_rev = ch_right
#Ecuación de diferencias
for k in range(D,K+D):
    #Filtra canal izquierdo
    yk_left_rev[k] = xk_left_rev[k-D] + alpha*yk_left_rev[k-D]
    #Filtra canal derecho
    yk_right_rev[k] = xk_right_rev[k-D] + alpha*yk_right_rev[k-D]
#0 del principio
yk_left_rev = yk_left_rev[D:]
yk_right_rev = yk_right_rev[D:]
#Normalizar señal filtrada en el rango de [-1,1]
yk_left_rev = yk_left_rev/np.max(yk_left_rev)
yk_right_rev = yk_right_rev/np.max(yk_right_rev)
#Genera señal filtrada
filtered_rev = np.zeros((K,2))
filtered_rev[:,0] =  yk_right_rev
filtered_rev[:,1] =  yk_left_rev
#Gráfico
plt.figure(3)
plt.subplot(2,1,1)
plt.plot(ch_right)
plt.title('Canal Derecho original vs filtrado')
plt.subplot(2,1,2)
plt.plot(yk_right_rev)
plt.figure(4)
plt.subplot(2,1,1)
plt.plot(ch_left)
plt.title('Canal izquierdo original vs filtrado')
plt.subplot(2,1,2)
plt.plot(yk_left_rev)

sd.play(filtered_rev,fs_new)






