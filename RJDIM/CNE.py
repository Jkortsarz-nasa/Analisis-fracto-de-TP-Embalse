import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import RJDIM as RJDIM

T=273 + 250 #K -Temperature inlet end
L=6 #m -largo del tubo
D=31536000*(1.17E-7)*np.exp(-8030/(1.987*T))/np.sqrt(2) #m2/años
r=0.2
dx=0.01 #m
dt=r*dx*dx/D #años
nodos_x=int(L/dx)
x = np.linspace(0,6,nodos_x)
def temp(x, min=265, max=305):
    return ((max-min) / (1 + np.exp(-(2*x-6)))) + min

def flux(x):
    return (0.0212*x**4 - 0.2654*x**3 + 0.6772*x**2 + 0.9245*x - 0.0121)
y = RJDIM.diferencias(x,0.2,0.4,0.1)
t = temp(x)


xCNE = np.array([2.002691954,	2.454912539,	2.786002577,	3.149394218,	3.528936352,	3.876177499,3.997308046,	4.457603877,	4.998654023])
yCNE = np.array([1.485148529,	3.217821813,	4.455445587,	5.79207964,	7.277228169,	8.465345859,	8.910891173,	10.24752523,	11.485149])
plt.plot(xCNE, yCNE, 'b^', label='N1K05 Inlet')
plt.show()
p0=[0.2,0.4,0.1]

popt, pcov = curve_fit(RJDIM.diferencias, xCNE, yCNE, p0=p0,method='lm')
#x_g_data = np.linspace(0,15)
print('Resultados:')
print(popt)

#y = RJDIM.diferencias(x,0.18511348,  -0.70806712, -50.3452297)
y = RJDIM.diferencias(x,0.14599376,   1.40780342, -43.2506118)
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
