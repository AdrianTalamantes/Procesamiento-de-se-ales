# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 15:54:24 2021

@author: AdrianTR
"""

import numpy as np
import matplotlib.pyplot as plt
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
#Grafico
plt.figure(2)
plt.plot(t,signal)
plt.title('Ciclo cardiaco')
plt.xlabel('Tiempo (min)')
#Se tomo el ultimo perdiodo que es el que segun el documento, cuando la persona
#estaba en una cinta caminadora con un angulo de inclinacion
#a 1.2m/s

############################################
#       Analisis en frecuencia
#Especificar Numero de muestras en la señal (siempre impar)
N = 91
#Vecto de tiempos para evaluar la interpolacion
te = np.linspace(0,t[-1],N)
#Aplica el metodo de interpolacion
inter  = interpolate.splrep(t, signal)#Argumentos, vecto de tiempos y señal
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
#K=30
#Periodo de la señal
T = M*ts
#Frecuencia fundamental
fo = 1/T
#Inicializar arreglos para los coeficientes de Fourier
ak = np.zeros(K)
bk = np.zeros(K)
ck = np.zeros(K)
phik = np.zeros(K)
fk = np.zeros(K)
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
plt.plot(t,signal,t,re_signal)#para que se paresca mas, podemos aumentar el numero de muestras N
#Con N pocas muestras, es un submuestreo(traslape) es decir una perdida de informacion.
#puede mejorar si pudieramos mover el periodo de muestro, pero no se puede modificar, los 500Hz con que se tomo la muestra, por lo que no podremos llegar a una mejor reconstruccion, pero se acerca mucho con N=401, el desplazamiento para arriba es por a0.
#Tambien le quitamos el determinismo, y asumimos periodicidad.

#Cuando graficamos el RAW, vemos un ruido en los 50Hz. 
#La frecuencia que distribuye CFE es de 60Hz. En Europa se trabaja en 50Hz.
#El filtro que se uso para quitar ese ruido, fue un filtro pasa bajas. Por eso en la filtrada, las señales de 40Hz en adelante se hacian 0.    
