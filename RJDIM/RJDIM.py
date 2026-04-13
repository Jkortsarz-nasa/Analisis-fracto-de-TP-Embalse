# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 16:41:21 2022

@author: sabba, RM
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#plt.style.use('ggplot')
#Ingreso de datos
#file_path =""
#datos = pd.read_csv(file_path, sep = ";", header=None)
#d = datos[0].values[:] #distance from the inlet end/m
#C_exp = datos[1].values[:] #deuterium concentration/ug/g


# Fitteos de la temperatura y flujo

masa=444.7 #g de material
T=273 + 310 #K -Temperature inlet end
L=6 #m -largo del tubo
a1=7.4e-05 #inlet
b1=4.7e-10
a2=1.905e-05 #outlet
b2=3.7e-05
D=31536000*(1.17E-7)*np.exp(-8030/(1.987*T))/np.sqrt(2) #m2/años
ci=4 #ug/g -Initial hydrogen concentration
r=0.2
dx=0.01 #m
dt=r*dx*dx/D #años
nodos_x=int(L/dx)
x = np.linspace(0,6,nodos_x)

def temp(x, min=265, max=305):
    return ((max-min) / (1 + np.exp(-(2*x-6)))) + min

def flux(x):
    return (0.0212*x**4 - 0.2654*x**3 + 0.6772*x**2 + 0.9245*x - 0.0121)

x = np.linspace(0,6,nodos_x)
def diferencias(x0, A, B, C):
    masa=444.7 #g de material
    T=273 + 250 #K -Temperature inlet end
    L=6 #m -largo del tubo
    a1=7.4e-05 #inlet
    b1=4.7e-10
    a2=1.905e-05 #outlet
    b2=3.7e-05
    D=31536000*(1.17E-7)*np.exp(-8030/(1.987*T))/np.sqrt(2) #m2/años
    ci=5 #ug/g -Initial hydrogen concentration
    r=0.2
    dx=0.01 #m
    dt=r*dx*dx/D #años
    nodos_x=int(L/dx)
    x = np.linspace(0,6,nodos_x)
    def temp(x, min=265, max=305):
        return ((max-min) / (1 + np.exp(-(2*x-6)))) + min

    def flux(x):
        return (0.0212*x**4 - 0.2654*x**3 + 0.6772*x**2 + 0.9245*x - 0.0121)

    cg = ((A*temp(x)+B*flux(x) + C)*np.ones(nodos_x)).reshape(nodos_x, 1)
    #cg=np.ones(nodos_x).reshape(nodos_x, 1) #ug/g*años
    C=np.full((nodos_x,1), ci, order='F') #vector concentraacion
    b=np.zeros((nodos_x,1), order='F') #vector de términos independientes
    time=13 #years of reactor operation conditions exposure 

    #----Calculo de r en función de la temperatura a lo largo del eje del caño
    #def r_temp(x,L):
    #    return dt*(31536000*1.17E-7*np.exp(-8030/(1.987*((578-538)*x/L+538))/np.sqrt(2)))/(dx**2)

    #----Calculo de r en función de la temperatura a lo largo del eje del caño
    def r_temp(x):
        min = 250
        max = 295
        #31536000*1.17E-7
        return dt*(31536000*1.17E-6*np.exp(-8030/(1.987*(temp(x)+273)/np.sqrt(2)))/(dx**2))

    #----Armado de la matrix para el cálculo de C
    A=np.eye(nodos_x,nodos_x)
    for i in range(0,nodos_x):
        A[i,i]=(1-2*r_temp(i*dx))
    for i in range(1,nodos_x):
        A[i,i-1]=r_temp(i*dx)
    for i in range(0,nodos_x-1):
        A[i,i+1]=r_temp(i*dx)
    #Condiciones de contorno
    A[0,0]=1
    A[0,1]=0
    A[-1,-1]=1
    A[-1,-2]=0
    j=dt

    #Itero en el tiempo
    while j<=time:
        b[0]=dt*(cg[0] + (1/(masa*(a1*np.sqrt(j)+b1))))
        b[1:-1]=dt*cg[1:-1]
        b[-1]=dt*(cg[-1] + (1/(masa*(a2*np.sqrt((j))+b2))))
        C=np.dot(A,C)+0.1*b
        j=j+dt
        
    #x = np.linspace(0,6,nodos_x)
    indices = []
    index = 0
    value = x[index]
    for i in x0-1:
        while value <= i:
            index += 1
            value = x[index]
        indices.append(index)
        index = 0
        value = x[index]

    return C[indices].flatten()


y = diferencias(x,0.2,0.4,0.1)
t = temp(x)

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Distancia del Inlet [m]')
ax1.set_ylabel('Concentracion D [ug/g]', color=color)
ax1.plot(x, y, color=color)
ax1.grid(True)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Temperatura [°C]', color=color)  # we already handled the x-label with ax1
ax2.plot(x, t, color=color)
ax2.tick_params(axis='y', labelcolor=color)

ax3 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax3.spines.right.set_position(("axes", 1.15))
color = 'tab:green'
ax3.set_ylabel('Fluencia [n/m2 s]', color=color)  # we already handled the x-label with ax1
ax3.plot(x, flux(x), color=color)
ax3.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

# https://levelup.gitconnected.com/solving-2d-heat-equation-numerically-using-python-3334004aa01a
#https://towardsdatascience.com/300-times-faster-resolution-of-finite-difference-method-using-numpy-de28cdade4e1
