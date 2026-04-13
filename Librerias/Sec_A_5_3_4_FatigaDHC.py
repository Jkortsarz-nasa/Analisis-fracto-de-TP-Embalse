# ==========================
# =================================

# Sección A.2.2.1 : Fatiga y DHC

# ==========================
# =================================
# Librerias

import math
import numpy as np




def calcular_L_DHC_semi_eliptico(a, c, K_deep, K_IH):
    """
    Calcula la longitud L_DHC de iniciación de crecimiento de grieta por DHC
    para un perfil de falla semi-elíptico, según CSA N285.8-15 A.5.3.3.4.

    Parámetros:
    - a: profundidad del defecto (mm)
    - c: semilongitud del defecto (mm)
    - K_deep: factor de intensidad de tensiones en el punto más profundo (MPa√m)
    - K_IH: umbral isotérmico para iniciación de DHC (MPa√m)

    Retorna:
    - L_DHC: longitud de iniciación de DHC (mm)
    """
    if K_deep <= 0:
        return 0.0

    R = K_IH / K_deep
    rel = (a / c) ** 0.5

    # Condiciones límite según A.5.3.3.4
    if R > 1.0:
        print("Se predice que DHC no ocurrirá")
        return 0.0
    elif R <= rel:
        return 2 * c
    else:
        try:
            factor = ((a / c) * math.sqrt(1 - R**4))/((R**2) * math.sqrt(1 - ((a / c)**2)))
            return 2 * c * factor
        except ValueError:
            # En caso de dominio inválido en sqrt (por redondeos)
            return 0.0
        

def calcular_KIH(beta=None, direccion='axial'):
    """
    Calcula el umbral inferior del factor de intensidad de tensiones K_IH para DHC,
    según la sección D.6 del CSA N285.8-15.

    Parámetros:
    -----------
    beta : float (opcional)
        Ángulo de orientación de la falla en grados con respecto al eje axial (solo para dirección 'mixta')
    direccion : str
        Dirección de la falla: 'axial', 'circunferencial' o 'mixta'

    Retorna:
    --------
    K_IH : float
        Valor inferior de K_IH en MPa√m
    """
    direccion = direccion.lower()

    if direccion == 'axial':
        return 4.5  # MPa√m, según D.6-1
    elif direccion == 'circunferencial':
        return 15.0  # MPa√m, según D.6-2
    elif direccion == 'mixta':
        if beta is None:
            raise ValueError("Debe proporcionar el ángulo beta para fallas con dirección mixta.")
        if not (0 <= beta <= 90):
            raise ValueError("El ángulo beta debe estar entre 0° y 90°.")
        beta_rad = math.radians(beta)
        return 9.0 / (1.3 + 0.7 * math.cos(2 * beta_rad))  # Según ecuación D.6-3
    else:
        raise ValueError("Dirección no válida. Use 'axial', 'circunferencial' o 'mixta'.")
    

def calcular_L_DHC(a, c, K_IH, K_deep):
    """
    Calcula la longitud de crecimiento L_DHC por DHC para una falla semi-elíptica
    según la cláusula A.5.3.3.4 del CSA N285.8-15 (Ecuación A.5-29).

    Parámetros:
    -----------
    a : float
        Profundidad de la falla [mm]
    c : float
        Semi-longitud de la falla [mm]
    K_IH : float
        Valor umbral del factor de intensidad de tensiones para DHC [MPa√m]
    K_deep : float
        Valor de K_I en el punto más profundo de la falla [MPa√m]

    Retorna:
    --------
    L_DHC : float
        Longitud total de crecimiento por DHC [mm]
    """
    if K_deep <= 0:
        raise ValueError("K_deep debe ser mayor que cero para evitar división por cero.")

    ratio = K_IH / K_deep
    raiz_ac = (a / c) ** 0.5

    # Caso sin crecimiento
    if ratio > 1:
        return 0.0

    # Caso de crecimiento total
    elif ratio <= raiz_ac:
        return 2 * c

    # Caso de crecimiento parcial (Ecuación A.5-29)
    numerador = (a / c) * (1 - ratio ** 4) ** 0.5
    denominador = (ratio) ** 2 * (1 - (a / c) ** 2) ** 0.5
    L_DHC = 2 * c * (numerador / denominador)
    return L_DHC


def crecimiento_fatiga(DeltaK, DeltaN):
    """
    Calcula el crecimiento de grieta por fatiga en base a la Ecuación A.5-34 del CSA N285.8-15.

    Parámetros:
    -----------
    DeltaK : float
        Rango de factor de intensidad de tensiones [MPa√m]
    DeltaN : float
        Número de ciclos de carga
    temperatura_C : float
        Temperatura del entorno [°C]

    Retorna:
    --------
    Delta_a : float
        Incremento en la profundidad de la grieta [m]
    """

    # Selección de parámetros según temperatura (Tabla D.21)
    # if temperatura_C <= 25:
    #     C0 = 3.70e-11
    #     n0 = 3.0
    # elif temperatura_C >= 250:
    C0 = 3.438e-10
    n0 = 3.439
    # else:
    #     # Interpolación lineal si temperatura intermedia
    #     C0_25 = 3.70e-11
    #     n0_25 = 3.0
    #     C0_250 = 3.438e-10
    #     n0_250 = 3.439

    #     fraccion = (temperatura_C - 25) / (250 - 25)
    #     C0 = C0_25 + fraccion * (C0_250 - C0_25)
    #     n0 = n0_25 + fraccion * (n0_250 - n0_25)

    # Se asume DeltaK_th = 0
    Delta_a = DeltaN * C0 * (DeltaK) ** n0
    return Delta_a  # en metros



def calc_Vr_mean(T, Ti, phi, Lambda=1): 
    """
    Calcula la velocidad de crecimiento radial por DHC (Vr) media. Equation D.10-1

    Parámetros:
    ------------
    T       : float
        Temperatura ambiente en el punto de evaluación [°C]
    Ti      : float
        Temperatura de irradiación (operativa) [°C]
    phi     : float
        Flujo de neutrones equivalente (fluencia/tiempo equivalente a plena potencia) [1e16 n/(m²·s)]
    Lambda  : float
        Factor isotópico de hidrógeno [-]
        Lambda = 0 si el material no está hidrurado (supervisión),
        Lambda = 1 si el material está hidrurado (condición estándar)

    Retorna:
    ---------
    Vr : float
        Velocidad media de crecimiento radial por DHC [m/s]
    """
    exponent = (-6328.8 / (T + 273)) - 0.0218 * Ti + 0.0285 * phi + 0.382 * Lambda # revisar!!!!!!!!!!!!!!!!!
    Vr = 14.2 * np.exp(exponent)
    return Vr


def calc_Vr_bounds(T, Ti, phi, Lambda):
    """
    Calcula los límites inferior y superior de ingeniería (95%) para Vr. (Equation D.10-13 y Equation D.10-14)

    Parámetros: (idénticos a calc_Vr_mean)

    Retorna:
    ---------
    Vr_L : float
        Límite inferior de ingeniería para Vr [m/s]

    Vr_U : float
        Límite superior de ingeniería para Vr [m/s]
    """
    exponent = (-6328.8 / (T + 273)) - 0.0218 * Ti + 0.0285 * phi + 0.382 * Lambda
    Vr_L = 1.674 * np.exp(exponent)
    Vr_U = 14.20 * np.exp(exponent)
    return Vr_L
    

def calc_Va_mean(T, Ti, phi, t):
    """
    Velocidad media de crecimiento axial por DHC (Va).  Eq. D.10-9

    Parámetros
    ----------
    T   : float o array
        Temperatura local [°C]
    Ti  : float o array
        Temperatura de irradiación/operación [°C]
    phi : float o array
        Flujo equivalente (fluence/equivalent full-power time) [10^16 n/(m²·s)]
    t   : float o array
        Tiempo equivalente a plena potencia [10^3 h]  (mil horas)

    Retorna
    -------
    Va : float o array
        Velocidad axial media por DHC [m/s]
    """
    exponent = (  5951.4/(T + 273)
                - 0.018341*Ti
                + 0.01511*phi
                - 2.324/np.sqrt(t) )
    Va = 6.989 * np.exp(exponent)
    return Va


def calc_Va_bounds(T, Ti, phi, t):
    """
    Límites de ingeniería (95%) para Va.  Ecs. D.10-13 (inferior) y D.10-14 (superior)

    Parámetros: ver calc_Va_mean

    Retorna
    -------
    Va_L : float o array
        Límite inferior representativo [m/s]
    Va_U : float o array
        Límite superior representativo [m/s]
    """
    exponent = (  5951.4/(T + 273)
                - 0.018341*Ti
                + 0.01511*phi
                - 2.324/np.sqrt(t) )
    Va_L = 1.709 * np.exp(exponent)   # Eq. D.10-13
    Va_U = 6.989 * np.exp(exponent)   # Eq. D.10-14
    return Va_L, Va_U



def evaluar_DHC_cooldown_transient(
    T_list, t_list, P_list, KI_list, Ti_list, phi_list, Lambda, T_SSD, DeltaT_DHC_c, K_IH,
    a0, c0
):
    """
    Evalúa el crecimiento de grieta por DHC durante un enfriamiento.

    Parámetros:
    ------------
    T_list      : list[float] - temperaturas [°C] vs tiempo
    t_list      : list[float] - tiempos [s]
    P_list      : list[float] - presiones [MPa] vs tiempo
    KI_list     : list[float] - lista de K_I instantáneo [MPa√m] por intervalo
    Ti_list     : list[float] - temperatura de irradiación [°C] por intervalo
    phi_list    : list[float] - fluencia [1e16 n/m²s] por intervalo
    Lambda      : float       - factor de isótopo H (0 o 1)
    T_SSD       : float       - temperatura T_SSD [°C]
    DeltaT_DHC_c: float       - ΔT_DHC^c requerido [°C]
    K_IH        : float       - umbral de propagación de grieta [MPa√m]
    a0, c0      : float       - tamaño inicial de la grieta [m]

    Retorna:
    ---------
    a_final     : float - profundidad final de grieta [m]
    c_final     : float - semilongitud superficial final de grieta [m]
    """
    a = a0
    c = c0
    T_DHC_c = T_SSD - DeltaT_DHC_c

    for i in range(1, len(t_list)):
        # Determinar si se inicia crecimiento desde este punto
        if T_list[i-1] >= T_DHC_c:
            continue  # Aún no inicia enfriamiento relevante para DHC

        dt = t_list[i] - t_list[i-1]
        T_avg = 0.5 * (T_list[i-1] + T_list[i])
        Ti_avg = 0.5 * (Ti_list[i-1] + Ti_list[i])
        phi_avg = 0.5 * (phi_list[i-1] + phi_list[i])
        KI = KI_list[i]  # Se asume que KI ya fue precalculado por presión

        if KI >= K_IH:
            # Calcular velocidades DHC en este intervalo
            Vr_i = calc_Vr_mean(T_avg, Ti_avg, phi_avg, Lambda)
            Va_i = calc_Va_mean(T_avg, Ti_avg, phi_avg, t=40)  # t ficticio o estimado

            # Aplicar ecuaciones A.5-39 y A.5-40
            a += Vr_i * dt
            c += Va_i * dt

    return a, c
