# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 17:36:48 2021

@author: AdrianTR
"""

import scipy.io.wavfile as rd
import numpy as np
import matplotlib.pyplot as plt

class espectrograma:
    
    def __init__(self, audio_signal, fs, frame_size, overlap, number_filters):
        self.audio = audio_signal
        self.mono = self.__Stereo()
        self.norm = self.__Norma()
        self.audio_length = len(self.norm)
        self.frame_samples = int(np.fix(frame_size * fs))
        self.overlap_samples = int(np.fix(self.frame_samples*overlap/100))
        self.n_frames = self.__Frames()
        self.half_samples = int(np.fix((frame_size * fs) / 2))
        self.Hann_window = self.__Window()
        self.freq = self.__Freq()
        self.spect = self.__Spect()
        self.BankFilter = self.__Mel_Bank(number_filters)
        self.bandFilteredSignal = self.__BandFilterSignal(number_filters)
        
    def __Stereo(self):
        n_channels = len(np.shape(self.audio))
        if n_channels > 1:
            temp1 = self.audio[:,0]
            temp2 = self.audio[:,1]
            avg = (temp1 + temp2) / 2
        else:
            avg = self.audio
        return avg

    def __Norma(self):
        mean = np.mean(self.mono)
        std = np.std(self.mono)
        norm_signal = (self.mono - mean)/std
        return norm_signal
    
    def __Frames(self):
        count = 0
        i = 0
        while(count < self.audio_length):
            count = i*self.overlap_samples + self.frame_samples
            i += 1
        return i - 1
    
    def __Window(self):
        W = np.zeros(self.frame_samples)
        mod = self.frame_samples % 2
        if mod == 0:
            end = self.half_samples
        else:
            end = self.half_samples + 1
        for n in range(-self.half_samples, end):
            W[n + self.half_samples] = 0.5 + 0.5*np.cos((2*np.pi*n)/(self.frame_samples -1))
        
        return W
    
    def __Freq(self):
        mod = self.frame_samples % 2
        if mod == 0:
            end = self.half_samples
            f = np.zeros(end)
        else:
            end = self.half_samples + 1
            f = np.zeros(end)
            
        for n in range(0, end):
            f[n] = n*fs/self.frame_samples
        return f
    
    def __Spect(self):
        frame_data = np.zeros(self.frame_samples)
        energy_matrix = np.zeros((self.half_samples,self.n_frames))

        for n in range(self.n_frames):
            frame_data = self.norm[n*self.overlap_samples : n*self.overlap_samples + self.frame_samples]
            # plt.figure(3)
            # plt.plot(frame_data)
            # plt.show()
            
            energy = np.sum(frame_data**2)
            
            if energy != 0:
                frame_data = frame_data*self.Hann_window
                #plt.figure(4)
                #plt.plot(frame_data)
                #plt.show()
                #g=1 #el breakpoint
                #Esto intenta que le de periodicidad (iniciar en 0 y terminar en 0) y darle una forma de campana, pero depende de los datos
                #Intenta eliminar las componentes de fuga
                tfd = np.fft.fft(frame_data)#Transformada rapida de Foureir en numpy
                #todas las transformadas de F, nos dan resultado complejo
                mag = np.abs(tfd)#si detecta un numero complejo, entonces aplica pitagoras, y si no, solo lo hace positivo
                #plt.figure(5)
                #plt.plot(self.freq,mag[0:self.half_samples+1])
                #plt.show()
                #g=1
                #el eje es frecuencia, un pico cerca de 2000Hz, otro en aprox 5000Hz.
                #en el segundo, se ve mas claro el 5000Hz, y otro en 2000Hz
                #En el  3ro, similar.
                #el 5to, se agrega otro pico a los 6000Hz, etc.
                #Como se analiza el ruido de una ballena, se puede decir que tiene frecuencias de los 6000Hz, 2000, etc.
                
                energy_matrix[:,n] = 10*np.log(mag[0:self.half_samples])
                #Normalizar la matriz
                # energy_matrix[:,n] = (energy_matrix[:,n] - np.mean(energy_matrix[:,n]))/np.std(energy_matrix[:,n])
                
        return energy_matrix
    
    def __Mel_Bank(self, M):
        #f_min = 0
        f_min=self.freq[0]
        #f_max=fs/2
        f_max = self.freq[-1]
        phi_min = 2595*np.log10((f_min/700)+1)
        phi_max = 2595*np.log10((f_max/700)+1)
        deltaphi = (phi_max - phi_min)/(M+1)

        phi_c = np.zeros(M+2)
        fc = np.zeros(M+2)
        for m in range(1,M+2):
            phi_c[m] = m*deltaphi
            fc[m] = 700*(10**(phi_c[m]/2595)-1)
        freqLen = len(self.freq)
        Banco_Mel = np.zeros((M, freqLen))
        for m in range(1,M+1):
            for i in range(freqLen):
                if self.freq[i] < fc[m-1]:
                    Banco_Mel[m-1,i] = 0
                if fc[m-1] <= self.freq[i] and self.freq[i] < fc[m]:
                    Banco_Mel[m-1,i] = (self.freq[i] - fc[m-1])/(fc[m]-fc[m-1])
                if fc[m] <= self.freq[i] and self.freq[i] < fc[m+1]:
                    Banco_Mel[m-1,i] = (self.freq[i] - fc[m+1])/(fc[m]-fc[m+1])
                if self.freq[i] >= fc[m+1]:
                    Banco_Mel[m-1,i] = 0
        return Banco_Mel
    
    def __BandFilterSignal(self,M):
        frame_data = np.zeros(self.frame_samples)
        filteredEnergyMatrix = np.zeros((M, self.n_frames))

        for n in range(self.n_frames):
            frame_data = self.norm[n*self.overlap_samples : n*self.overlap_samples + self.frame_samples]
            
            energy = np.sum(frame_data**2)
            
            if energy != 0:
                frame_data = frame_data*self.Hann_window
                tfd = np.fft.fft(frame_data)
                mag = np.abs(tfd)
                for m in range(M):
                    filteredEnergyMatrix[m,n] = np.sum(self.BankFilter[m] * mag[:self.half_samples])
        return filteredEnergyMatrix
            
            

#Ruta del archivo de la señal para analizar
path = 'D:/UACH - Maestria/1er semestre/Procesamiento de señales/Clase 13 Octubre/ballena.wav'
#path = 'D:/UACH - Maestria/1er semestre/Procesamiento de señales/Clase 13 Octubre/song.wav'
#Abre archivo de la señal
fs, audio_signal = rd.read(path)
#############################################################################
#       Realiza el espectrograma
#Tamaño de frame (segundos)
frame_size = 0.1
#Tamaño de traslape (porcentaje)
overlap = 60
signal = audio_signal
M = 20
#para musica podemos tomar el frame_size, mas grande como 0.05 o 0.1
#Abre archivo de la señal
#Crea objeto de la clase espectrograma
spectrog = espectrograma(signal,fs,frame_size,overlap, M)
plt.figure(1)
#plt.plot(spectrog.audio)
plt.subplot(2,1,1)
#plt.plot(spectrog.audio)
plt.plot(spectrog.mono)
plt.subplot(2,1,2)
plt.plot(spectrog.norm)

plt.figure(2)
plt.plot(spectrog.Hann_window)
plt.title('Ventana de Hann')
plt.xlabel('Muestras')
#print(spectrog.freq)#Imprimir las frecuencias
#la muestra de 11k, que sale al final, en la parte de la ballena, 
# print(spectrog.freq)
plt.figure(3)
plt.imshow(spectrog.spect[::-1], cmap=plt.get_cmap('jet'),extent=[0,spectrog.n_frames,0,fs/2],aspect="auto")
#plt.imshow(spectrog.spect[::-1], cmap = plt.get_cmap('jet'),
          # extent = [0,spectrog.n_frames,0,spectrog.freq[spectrog.half_samples-1]],
         #  aspect = 'auto')
plt.title('Espectrograma de la señal')
plt.xlabel('Número de Frames')
plt.ylabel('Frecuencia en Hertz')
#el fs/2 es porque hasta ahi llegara la frecuencia aproximadamente

#cuando graficamos el sonido de la licuadora, observamos que hay muchas frecuencias en todas aprtes, excepto en la franja de 10k-12.5k, que es donde se observa menos energia. La licuadora es un buen tipo de ruido para contaminar las señales.

#cuando graficamos la persona hablando enojada, su energia se distribuye desde 0 hasta los 12k Hz, mientras en otras partes del espectrograma, resalta mucho las partes de bajas frecuencias, debajo de los 3k Hz, en la franja azul, seguro dejo de hablar.

#con el archivo song.wav, se ve que tiene mucho contenido de bajas frecuencias, y a los 15k Hz se ve la franja que divide, por encima donde ya no se tiene frecuencias sobre los 15k Hz. Si ponemos un overlap de 10, y frame_size de 0.015 o 0.025, 0.01, se empieza a notar un poco de diferencia. 

#plt.imshow(spectrog.spect[::-1], cmap=plt.get_cmap('jet'),extent=[0,spectrog.n_frames,0,spectrog.freq[spectrog.half_samples]],aspect="auto")
#El silbido es lo que aparece en altas frecuencias (en rojo), y lo ronco en bajas frecuencias, la energia roja en 0,4000, son los picos que corresponden al silbido, luego en el tono ronco, es el costillar amarillo, con energia a baja frecuencia que se compone de mas partes de señal.
#entre mas azul, la energia tiene a cero, la amarilla es no tan grande (a medias), 
#en 200-400, vemos un costillar, de notas que se tocan a la vez
#el silbido de la ballena, es la parte de los 4000Hz, entre 0 y 200
#la parte del sonido ronco de la ballena, es la parte roja de baja frecuencia, y por eso tiene mucha energia

#Cuando cambiamos a la 04 electric (licuadora), 

plt.figure(4)
for n in range(M):
    plt.plot(spectrog.BankFilter[n])

plt.figure(5)
plt.imshow(spectrog.bandFilteredSignal[::-1], cmap = plt.get_cmap('jet'),
           extent = [0,M,0,spectrog.n_frames],
           aspect = 'auto')
plt.title('Espectrograma de la señal por bandas')
plt.xlabel('Bandas')
plt.ylabel('Frames')

