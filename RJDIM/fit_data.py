import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

#funcion de fiteo
def func(t, alfa, beta):
    r=5.373
    w=0.406
    ro=6.490
    long = 5.0
    return beta * (alfa*np.sqrt(t)-np.log(1 + alfa*np.sqrt(t)))
    #return beta * (alfa*np.sqrt(t)- np.log(1 + alfa*np.sqrt(t)))
print("end")
#Bruce Inlet
xdata =[0.0,4.027,4.58,6.9,7.5,8.194,9.027,9.305,10.69,9.72,10.0,10.8,]
ydata= [0.0,59.0,55.0,80.0,72.0,79.0,90.0,59.0,67.0,105.0,86.0,97.0]
#Bruce Inlet capped
# xdata =np.array([0.0,4.027,4.58,6.9,7.5,8.194,9.027,9.72,10.0,10.8,])
# ydata= np.array([0.0,59.0,55.0,80.0,72.0,79.0,90.0,105.0,86.0,97.0])
plt.plot(xdata, ydata, 'b^', label='N1K05 Inlet')
p0=[1.,1.]

popt, pcov = curve_fit(func, xdata, ydata, p0=p0,method='lm', maxfev = 10000)
x_g_data = np.linspace(0,15)

plt.plot(x_g_data, func(x_g_data, *popt), 'r-',label='Curve Fitting: alfa=%5.3f, beta=%5.5f' % tuple(popt))
plt.xlabel('YFP')
plt.ylabel('Deuterium picked up')
plt.legend()
plt.grid()
plt.show()
print('alfa: {}'.format(popt[0]))
print('beta: {}'.format(popt[1]))
print('a: {}'.format(2/(popt[0]*popt[1])))
print('b: {}'.format(2/(popt[0]**2*popt[1])))