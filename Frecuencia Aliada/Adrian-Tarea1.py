# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 15:27:39 2021

@author: AdrianTR
"""

def F_aliada(rango, fs, fxn=list()):
    aliada =[]
    for fxni in fxn:
        n=0
        while fxni > rango:
           
           fxni=fxni-fs
           n+=1
        if n==0:
            aliada.append((False, fxni))
        else:
            aliada.append((True, fxni))
    return aliada

def Lista_Frec(n):
    frecuencias=[]
    print('Escriba una frecuencia y de enter hasta terminar')
    for i in range(n):
        frecuencias.append(int(input("Frecuencia " + str(i+1) + ": "))) 
    return frecuencias


# dameN = int(input("¿Cuántas frecuencias quieres probar? "))
# lista_frecuencias = Lista_Frec(dameN)
# fsoriginal= int(input("Escriba la frecuencia fs: "))

# print("Las frecuencias aliadas son: ", F_aliada(fsoriginal/2, fsoriginal, lista_frecuencias) )
