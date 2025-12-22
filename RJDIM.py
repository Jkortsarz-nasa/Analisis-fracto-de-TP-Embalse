import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def diferencias(x0, A, B, C, ci=5, time=13):
    # --- Constantes y malla
    masa = 444.7   # g de material
    L = 6.0        # m     #L = 6.185
    dx = 0.01      # m (target)
    # malla que incluye el extremo y respeta dx real
    nodos_x = int(round(L/dx)) + 1
    x = np.linspace(0.0, L, nodos_x)
    dx = x[1] - x[0]

    # contornos (igual que tu versión)
    a1, b1 = 7.4e-05, 4.7e-10   # inlet
    a2, b2 = 1.905e-05, 3.7e-05 # outlet

    # temperatura y flujo
    def temp(xx, min=265.0, max=305.0):
        return ((max - min) / (1.0 + np.exp(-(2.0*xx - 6.0)))) + min

    def flux(xx):
        return (0.0212*xx**4 - 0.2654*xx**3 + 0.6772*xx**2 + 0.9245*xx - 0.0121)

    # fuente cg(x) con los parámetros de ajuste (renombro C -> C0 para no pisar concentración)
    C0 = C
    cg = (A*temp(x) + B*flux(x) + C0).reshape(nodos_x, 1)  # [ug/g·año]

    # Difusividad dependiente de T(x): D(x) en m^2/año (para deuterio)
    # 31536000 convierte m^2/s -> m^2/año ; /sqrt(2) fuera del exponencial
    D_x = 31536000.0 * (1.17e-7) * np.exp(-8030.0 / (1.987 * (temp(x) + 273.0))) / np.sqrt(2.0)

    # Paso temporal a partir de un D "base" (tu criterio original)
    r = 0.2
    D_base = 31536000.0 * (1.17e-7) * np.exp(-8030.0 / (1.987 * (250.0 + 273.0))) / np.sqrt(2.0)
    dt = r * dx*dx / D_base  # años

    # Vector de concentración (estado) y término independiente
    Cconc = np.full((nodos_x, 1), float(ci))  # concentración
    b = np.zeros((nodos_x, 1))

    # r_i = D(x_i)*dt/dx^2, para llenar la matriz tridiagonal
    r_i = (D_x * dt) / (dx*dx)

    # Matriz del esquema explícito (tridiagonal con r variable por nodo)
    Acoef = np.eye(nodos_x)
    for i in range(0, nodos_x):
        Acoef[i, i] = 1.0 - 2.0 * r_i[i]
    for i in range(1, nodos_x):
        Acoef[i, i-1] = r_i[i]      # vecino izquierdo
    for i in range(0, nodos_x-1):
        Acoef[i, i+1] = r_i[i]      # vecino derecho

    # Condiciones de contorno tipo flujo (Neumann) implementadas por fila identidad
    Acoef[0, 0] = 1.0;  Acoef[0, 1] = 0.0
    Acoef[-1, -1] = 1.0; Acoef[-1, -2] = 0.0

    # Marcha temporal
    t = dt
    while t <= time:
        # términos de borde con √t coherente en ambos extremos
        b[0, 0]   = dt * ( cg[0, 0]   + 1.0 / (masa * (a1 * np.sqrt(t) + b1)) )
        b[1:-1, 0]= dt *   cg[1:-1, 0]
        b[-1, 0]  = dt * ( cg[-1, 0]  + 1.0 / (masa * (a2 * np.sqrt(t) + b2)) )

        # actualización (sin factor 0.1)
        Cconc = Acoef @ Cconc + b

        t += dt

    # Mapear x0 -> índices de malla de forma robusta
    x0 = np.asarray(x0, dtype=float)
    idx = np.searchsorted(x, x0, side='left')
    idx = np.clip(idx, 0, nodos_x - 1)

    return Cconc[idx, 0].ravel()


'''
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
'''