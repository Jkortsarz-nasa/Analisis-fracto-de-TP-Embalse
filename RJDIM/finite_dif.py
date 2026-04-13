import numpy as np
import matplotlib.pyplot as plt
#from tqdm import tqdm
# Diferencias finitas de ecuacion
# Ct(x, t) = D*Cxx(x, t) + S(x, t)

# h -> x -> i
# k -> t -> j

L = 6000   # Long.
T = 250 + 273.15 #Temp
#D
c =  (1.17e2*np.exp(-8030/(1.987*T)))
k =  0.1# Delta t.
h = 0.5 # Delta x.
time=10
tiempo = int(time/k)
alfa = (c*k)/(h**2)  # Converge para alfa < 1/2

# Ingreso de D
def S(x, t):
    #parametros fitteados
    a=32861.10 #inlet
    b =0.2087
    #a=7.38830313353633e-05 #inlet
    #b=4.697050219603037e-10

    r=5.373
    w=0.406
    ro=6.49*1e3
    long = 5
    scaling= 1

    if x < 1:
        return ((60*60*365*24*scaling/ ((a*np.sqrt(t) + b)*2*np.pi*r*w*long*ro)) + 1)
    else:
        return 1

#condiciones iniciales
C = 15*np.ones([int(L/h), tiempo+1])

# Condiciones de contorno constantes
#C[0, :] = 30
#C[L-1, :] = 40
print(f'alfa: {alfa}')
if alfa < 1:
    print('Solving Model...')
    for j in range(tiempo):
    #for j in tqdm(range(tiempo)):
        #Condiciones de contorno dinamicas  
        C[0, j+1] = C[0, j] + k*S(0, j)
        C[int(L/h)-1, j+1] = C[int(L/h)-1, j] + k*S(int(L/h)-1, j)

        for i in range(1,int((L/h)-1)):
            C[i, j+1] = alfa*(C[i+1, j] + C[i-1, j]) + (1 - 2*alfa)*C[i, j] + k*S(i, j)
else:
    print('Alpha value is too big!')

#plt.plot(range(int(L/h)),C)
long=50
plt.plot(range(int(long)),C[0:long,tiempo])
#plt.legend()
plt.title(f'Simulación a {tiempo*k} años')
plt.xlabel('Posicion (mm)')
plt.ylabel('C Deuterio(ug/g)')
plt.grid()
plt.show()
#print(np.round(C, 1))