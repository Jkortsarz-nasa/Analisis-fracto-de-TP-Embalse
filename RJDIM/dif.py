import numpy as np
import fd_solver as fd
import numpy as np
import matplotlib.pyplot as plt

# Parameters
masa=444.7 #g de material
T=270 #C -Temperature inlet end
L=6 #m -largo del tubo

ar = 0.0542 # Parámetros fitteados de la tasa de ingreso de D
br = 0.0126

cg=1 #ug/g*year
ci=10 #ug/g -Initial hydrogen concentration
dt=0.01 #años
dx=0.01 #m
nodos_x=int(L/dx)

C=np.full((nodos_x,1), ci, order='F') #vector concentraacion
S=np.zeros((nodos_x,1), order='F') #vector de términos fuente
time=14

alpha =  dt*(31536000*1.17E-7*np.exp(-8030/(1.987*(T+273))))/(np.sqrt(2)*(dx**2))
print(f'Alpha: {alpha}')
K = fd.K(alpha, nodos_x)

j = 0
print('Solving model...')

Ct = []
while j <= time:
    # Actualiza el termino fuente
    S[0]=   cg + (1000/(masa*(ar*np.sqrt(j)+br)))
    S[1:-1]= cg
    S[-1]=  cg + (1000/(masa*(ar*np.sqrt(j)+br)))

    # Actualiza la concentracion en el tiempo
    C = np.dot(K, C) + dt*S
    Ct.append(C)
    j = j+dt

print('Done.')
#Gráfico
x = np.linspace(0, 6, nodos_x)
y = C

Ct = np.array(Ct)
print(np.size(Ct))
plt.plot(x, Ct[1000])
#plt.plot(x, Ct[800])
#plt.plot(x, Ct[600])
#plt.plot(x, Ct[400])
#plt.plot(x, Ct[200])
plt.xlabel('DISTANCE FROM THE INLET END[m]')
plt.ylabel('DEUTERIUM CONCENTRATION[ug/g]')
plt.grid()
plt.show()

dist = 2
print(f'Dist: {np.round(x[dist], 2)} | [D]: {np.round(y[dist][0], 2)}')