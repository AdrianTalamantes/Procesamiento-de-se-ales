# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 19:56:12 2021

@author: AdrianTR
"""

import scipy.io.wavfile as rd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import sounddevice as sd

#Abrimos el archivo a analizar
fs, audio_signal = rd.read('ballena.wav')
#Normalizamos la señal
nor_signal = (audio_signal-np.mean(audio_signal))/np.std(audio_signal)
#Factor de atenuación
fact = 1/6
nor_signal = nor_signal * fact
#Grafico
plt.figure(1)
plt.plot(nor_signal)
plt.title('Señal Original')
plt.xlabel('Muestras')

#               Analisis en frecuencia de la señal
#Numero de muestras de la señal
N = len(nor_signal)
#Calcula frecuencias de las componentes de la señal
n = np.arange(0,N,dtype=float)
f = n * fs / N
#Aplicamos la transformada rapida de Fourier
Xf = np.abs(np.fft.fft(nor_signal))
#Mitad de los datos
half = int(np.fix(N/2))
#Grafica del espectro
plt.figure(2)
plt.plot(f[:half],Xf[:half])
plt.title('Espectro de la Señal')
plt.xlabel('Frecuencia (Hz)')

#       Respuesta en frecuencia del Filtro Digital Pasa Banda
#Frecuencia central
fo = np.array([31.5,63,125,250,500,1000,2000,4000,8000,16000])
#Ganancia
gain = np.array([1,10,10,100,10,10,10,10,10,10])
#Ganancia en dB
gain_dB = 10*np.log10(gain)
#Frecuencia central digital 
Fo = fo/fs
#Frecuencia de central digital en rad/s
Wo = 2*np.pi*Fo
#Factor de calidad
#Q= 1/2
Q = 10
#Ancho de banda en rad/s
Br = Fo/Q
#Constante C y Beta de la transformada de bilineal
C = np.tan(Br/2)
Beta = np.cos(Wo)
#Numerador de Hz
num = np.array([C,np.zeros(10),-C])
#Denominadr de Hz
den = np.array([C+1,-2*Beta,-C+1])
#Respuesta en frecuencia del filtro
for i in range(10):
    Wz, Hwz = signal.freqz(num[:,i],den[:,i])
    #Convierte Wz en una frecuencia digital F
    F = Wz/(2*np.pi)
    #Convierte F digital en Hertz
    fh = F*fs
    #Grafico
    plt.figure(3)
    plt.plot(fh,np.abs(Hwz))
    plt.title('Espectro de Magnitud')
    plt.xlabel('Frecuencia (Hz)')

#                   Filtrado de la señal
filt_signal = np.zeros(N)
#Filtrado de la señal
for i in range(10):
    filt_signal += signal.lfilter(num[:,i],den[:,i],nor_signal)*gain_dB[i]
#Grafico
plt.figure(4)
plt.subplot(2,1,1)
plt.plot(nor_signal)
plt.title('Señal original vs Señal filtrada')
plt.subplot(2,1,2)
plt.plot(filt_signal)
plt.xlabel('Muestras')

#           Analisis en frecuencia de la señal filtrada

#Aplicamos la transformada rapida de Fourier
Xff = np.abs(np.fft.fft(filt_signal))
#Grafica del espectro
plt.figure(5)
plt.plot(f[:half],Xf[:half])
plt.plot(f[:half],Xff[:half])
plt.title('Espectro de la Señal filtrada')
plt.xlabel('Frecuencia (Hz)')

#Reproduccion de la señal
#sd.play(nor_signal,fs)
# sd.wait()
#sd.play(filt_signal,fs)









