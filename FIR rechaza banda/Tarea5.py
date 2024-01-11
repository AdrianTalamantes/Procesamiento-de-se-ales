# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 19:56:12 2021

@author: AdrianTR
"""
from operator import sub
import scipy.io.wavfile as rd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sounddevice as sd

#################################################
#           ABRE SEñAL DE BALLENA
#Abre del archivo de la señal para analiozar
path= 'D:/UACH - Maestria/1er semestre/Procesamiento de señales/Tarea 5 parcial 3/ecg.txt'
# Abre archivo de la señal
file = open(path, "r")
openedFile = file.read().replace('\n', ' ')
data_str = openedFile.split(' ')
audio_signal = list(map(int, data_str))
fs = len(audio_signal)
# Normalizar señal
nor_signal = (audio_signal-np.mean(audio_signal))/np.std(audio_signal)
# Factor de atenuacion para dejar en el rango [-1,1] para reproducir la señal
fact = 1/6
nor_signal = nor_signal*fact
# Grafica de la señal
plt.figure(1)
plt.plot(nor_signal)
plt.title('Señal Original')
plt.xlabel('Muestras')
###################################################
#   ESPECIFICACIONES DE DISEñO DEL FILTRO FIR Rechaza banda
# Frecuencia en la banda de paso (Frecuencia de corte)
fp = 50
# Frecuencia en la banda de supresion en Hz
fbs = 500
# Frecuencia en la banda de supresion normalizada
Fp = fp/fs
# Frecuencia en la banda de supresion normalizada
Fbs = fbs/fs
# Ancho de transicion en frecuencia normalizada
dF1 = Fbs - Fp
# Ancho de transicion en frecuencia normalizada
dF2 = Fbs + Fp
# Longuitud del filtro FIR
N = int(np.fix(5.5/dF1))
# Determina si la secuencia es par o impar
even_odd = N % 2

if even_odd != 0:
    k = np.arange(-(N-1)/2,(N-1)/2+1)
else: 
    k = np.arange(-N/2-1,N/2-1)

# Frecuencia de corte normalizada
Fc1 = Fp + dF1/2
# Frecuencia de corte normalizada
Fc2 = Fp + dF2/2
# Evaluar la respuesta al impulso del filtro pasabajas
hk = (2*Fc1*np.sin(2*np.pi*k*Fc1)/(2*np.pi*k*Fc1))-(2*Fc2*np.sin(2*np.pi*k*Fc2)/(2*np.pi*k*Fc2))
if even_odd != 0:
    hk[int((N-1)/2)] = 1-(2*Fc2-Fc1)
else: 
    hk[int(N/2)+1] = 1-(2*Fc2-Fc1)

# Grafica de la respuesta al impulso
plt.figure(2)
plt.plot(hk)
plt.title('Respuesta al Impulso Truncada')
plt.xlabel('Muestras')

# Evaluar la funcion ventana
wk = 0.42 + 0.5*np.cos((2*np.pi*k)/(N-1)) + 0.08*np.cos((4*np.pi*k)/(N-1))

# Grafica a la funcion ventana
plt.figure(3)
plt.plot(wk)
plt.title('Funcion ventana tipo Blackman')
plt.xlabel('Muestras')

# Diseño del filtro FIR
hwk = hk*wk
# Grafica de la funcion ventana
plt.figure(4)
plt.plot(hwk)
plt.title('Respuesta al impulso con funcion ventana')
plt.xlabel('Muestras')

# Respuesta a la frecuencia del filtro Pasa Bajas
wz, Hz = signal.freqz(hk,1)
wz, Hwz = signal.freqz(hwk,1)

# Desnormalizar la frecuencia
F = wz/(2*np.pi)
f = F*fs

# Grafica de la respuesta al impulso
plt.figure(5)
plt.plot(f,np.abs(Hz))
plt.plot(f,np.abs(Hwz))
plt.title('Espectro de Magnitud')
plt.xlabel('Frecuencia en Hz')

# Filtrado de la señal
filt_signal = signal.lfilter(hwk,1,nor_signal)

# Grafica de la señal
plt.figure(1)
plt.plot(filt_signal)
plt.title('Señal Original')
plt.xlabel('Muestras')