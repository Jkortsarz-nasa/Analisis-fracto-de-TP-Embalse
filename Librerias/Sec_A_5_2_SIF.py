# ==========================
# =================================

# Sección A.2.2.1 : Caracterización geométrica de los defectos planares

# ==========================
# =================================
# Librerias

import pandas as pd
import math
import numpy as np
from Sec_D_3_4_2_PropMec import sigma_y_transverse_lb, sigma_y_axial_lb


from scipy.interpolate import griddata

def calcular_KI_axial(a, c, p, r_i, w, sigma_res, Fp, Fm):
    """
    Calcula el factor de intensidad de tensiones K_I para una falla axial según CSA N285.8-15 (Ecuación A.5-1).

    Parámetros:
    - a: profundidad de la falla [mm]
    - c: semi-longitud de la falla axial [mm]
    - p: presión interna [MPa]
    - r_i: radio interno del tubo [mm]
    - w: espesor de la pared [mm]
    - sigma_res: tensión residual circunferencial (rolled joint) [MPa]
    - Fp: factor de corrección geométrica para presión interna
    - Fm: factor de corrección geométrica para carga de membrana

    Retorna:
    - K_I: factor de intensidad [MPa√m]
    """
    Q = 1 + 1.464 * (a / c) ** 1.65  # A.5-2
    sigma_p = p * ((r_i / w) + 1)
    K_I = ((sigma_p * Fp) + (sigma_res * Fm)) * math.sqrt(math.pi * a / (1000.0 *Q))  # convertir a metros
    return K_I


def calcular_KI_falla_circunferencial(sigma_a, a, b, F_c):
    """
    Calcula el factor de intensidad de tensiones K_I para una falla circunferencial atravesando completamente la pared,
    según la sección A.5.2.3.2 del CSA N285.8-15.

    Parámetros:
    - sigma_a : tensión axial total aplicada (MPa)
    - a       : profundidad del defecto (mm)
    - b       : mitad de la longitud circunferencial de la falla (mm)
    - F_c     : factor geométrico para falla circunferencial (sin unidades)

    Retorna:
    - K_I : factor de intensidad de tensiones (MPa·√mm)
    """
    Q = 1 + 1.464 * (a / b) ** 1.65
    KI = sigma_a * np.sqrt(np.pi * a / (Q*1000)) * F_c
    return KI

def correccion_geometrica_superficie_libre_circunferencial(a, b, w):
    """
    Calcula el factor de corrección geométrica F_c^SP para una falla circunferencial
    en el punto de intersección con la superficie libre, según la ecuación A.5-25.

    Aplica cuando:
    - 0.20 ≤ a/b ≤ 1.0
    - 0.0 ≤ a/w ≤ 0.80

    Parámetros:
    -----------
    a : float
        Profundidad de la falla [mm]
    b : float
        Mitad de la longitud circunferencial de la falla [mm]
    w : float
        Espesor de pared del tubo [mm]

    Retorna:
    --------
    Fc_SP : float
        Factor de corrección geométrica F_c^SP
    """
    a_b = a / b
    a_w = a / w

    if not (0 <= a_b <= 1.0 and 0.0 <= a_w <= 0.80):
        raise ValueError("Parámetros fuera del rango de aplicabilidad para A.5.2.3.5")

    # Coeficientes M1, M2 y M3 (Eqs. A.5-20 a A.5-22)
    M1 = 1.13 - 0.09 * a_b
    M2 = -0.54 + 0.89 / (0.2 + a_b)
    M3 = 0.5 - (1 / (0.65 + a_b)) + 14 * (1 - a_b) ** 24

    # Ecuación A.5-25
    parte_pol = M1 + M2 * (a_w)**2 + M3 * (a_w)**4
    Fc_SP = (a_b) ** (0.5) * parte_pol * (1.1 + 0.35 * (a_w)**2)

    return Fc_SP

def calcular_KI_con_plastic_zone(tipo_falla, a, c, p, r_i, w, sigma_res, K_deep, K_surf, sigma_y):
    """
    Calcula el factor de intensidad K_I para falla axial CON corrección de zona plástica,
    según CSA N285.8-15 (Ecuaciones A.5-1 a A.5-6, Sección A.5.2.2.3).

    IMPORTANTE: Esta corrección de zona plástica SOLO debe aplicarse cuando el cálculo
    de K_I se usa para evaluar la iniciación del crecimiento de falla
    conforme a la Sección A.5.4 del código (fracture initiation evaluation).
    NO aplicar en evaluaciones de crecimiento por fatiga o DHC (A.5.3).

    Parámetros:
    - tipo_falla: "axial" o "circunferencial"
    - a: profundidad del defecto [mm]
    - c: semilongitud axial del defecto [mm]
    - p: presión interna [MPa]
    - r_i: radio interno del tubo [mm]
    - w: espesor del tubo [mm]
    - sigma_res: tensión residual [MPa]
    - K_deep: K_I en el punto más profundo de la falla [MPa√m]
    - K_surf: K_I en el punto de intersección con la superficie [MPa√m]
    - sigma_y: tensión de fluencia del material [MPa]

    Retorna:
    - K_I corregido [MPa√m] usando dimensiones efectivas a_e, c_e
    """

    # Radios de zona plástica (ecuaciones A.5-4 y A.5-6)
    r_yf = (1 / (6 * math.pi)) * (K_deep / sigma_y) ** 2
    r_yc = - (1 / (2 * math.pi)) * (K_surf / sigma_y) ** 2

    # Dimensiones efectivas con corrección plástica (A.5-3 y A.5-5)
    a_e = a + r_yf
    c_e = c + r_yc

    # Cálculo completo (A.5-1) con dimensiones corregidas
    

    if tipo_falla == "axial":
        Fp, Fm = funcion_Fp_y_Fm_DP_A_5_7_8_CORREGIDA(a_e, c_e, w, r_i, r_i + w)
        K_I_corr = calcular_KI_axial(a_e, c_e, p, r_i, w, sigma_res, Fp, Fm)
        '''
        term_presion = p * ((r_i / w) + 1) * Fp
        K_I_corr = ( term_presion + sigma_res * Fm )* math.sqrt((math.pi * a_e) / (Q_e * 1000.0))
        '''

    elif tipo_falla == "circunferencial":
        b_e=c_e
        Fc = correccion_geometrica_superficie_libre_circunferencial(a_e, b_e, w)
        K_I_corr = calcular_KI_falla_circunferencial(p, a_e, b_e, Fc)
        

    return K_I_corr


# Tabla A.2 - Coeficientes G en función de a/c y a/w
tabla_A2 = pd.DataFrame({
    "a/c":   [0.20, 0.20, 0.20, 0.40, 0.40, 0.40, 1.00, 1.00, 1.00],
    "a/w":   [0.20, 0.50, 0.80, 0.20, 0.50, 0.80, 0.20, 0.50, 0.80],
    "G0":    [1.115, 1.427, 1.872, 1.072, 1.217, 1.393, 1.015, 1.050, 1.090],
    "G1":    [0.673, 0.783, 0.960, 0.672, 0.723, 0.806, 0.715, 0.729, 0.760],
    "G2":    [0.514, 0.571, 0.671, 0.523, 0.573, 0.601, 0.588, 0.596, 0.618],
    "G3":    [0.438, 0.462, 0.529, 0.441, 0.472, 0.493, 0.512, 0.515, 0.532],
})

# Tabla A.3 estructurada como DataFrame
tabla_A3 = pd.DataFrame({
    "a/w":   [0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.500, 0.500, 0.500, 0.750, 0.750, 0.750],
    "ri/w":  [5,     10,    20,    5,     10,    20,    5,     10,    20,    5,     10,    20],
    "Fax":   [1.19,  1.20,  1.20,  1.38,  1.44,  1.45,  2.10,  2.36,  2.51,  3.30,  4.23,  5.25],
})

# Tabla A.4: Coeficientes de corrección geométrica para intersección con la superficie libre
tabla_A4 = pd.DataFrame({
    "a/c": [0.20, 0.20, 0.20, 0.40, 0.40, 0.40, 1.00, 1.00, 1.00],
    "a/w": [0.20, 0.50, 0.80, 0.20, 0.50, 0.80, 0.20, 0.50, 0.80],
    "G0":  [0.607, 0.791, 1.179, 0.777, 0.936, 1.219, 1.140, 1.219, 1.348],
    "G1":  [0.079, 0.138, 0.253, 0.125, 0.176, 0.259, 0.197, 0.221, 0.255],
    "G2":  [0.023, 0.052, 0.104, 0.043, 0.069, 0.106, 0.070, 0.085, 0.099],
    "G3":  [0.010, 0.027, 0.056, 0.021, 0.036, 0.056, 0.038, 0.044, 0.051]
})

# Tabla A.5: Geometry factor Fcr
tabla_A5 = pd.DataFrame({
    "a/w":   [0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.500, 0.500, 0.500, 0.750, 0.750, 0.750],
    "ri/w":  [5,     10,    20,    5,     10,    20,    5,     10,    20,    5,     10,    20],
    "Fcr":   [1.16,  1.19,  1.22,  1.26,  1.32,  1.36,  1.61,  1.82,  2.03,  2.15,  2.49,  2.89]
})

def interpolar_coeficientes_G(a_c_objetivo, a_w_objetivo):
    """
    Interpola los coeficientes G0, G1, G2 y G3 a partir de la Tabla A.2 del CSA N285.8-15
    según los valores dados de a/c y a/w.

    Parámetros:
    - a_c_objetivo: razón a/c del defecto
    - a_w_objetivo: razón a/w del defecto

    Retorna:
    - Diccionario con los coeficientes interpolados {G0, G1, G2, G3}
    """

    puntos = tabla_A2[["a/c", "a/w"]].values
    G0_interp = griddata(puntos, tabla_A2["G0"].values, (a_c_objetivo, a_w_objetivo), method='linear')
    G1_interp = griddata(puntos, tabla_A2["G1"].values, (a_c_objetivo, a_w_objetivo), method='linear')
    G2_interp = griddata(puntos, tabla_A2["G2"].values, (a_c_objetivo, a_w_objetivo), method='linear')
    G3_interp = griddata(puntos, tabla_A2["G3"].values, (a_c_objetivo, a_w_objetivo), method='linear')

    return {
        "G0": float(G0_interp),
        "G1": float(G1_interp),
        "G2": float(G2_interp),
        "G3": float(G3_interp),
    }


def interpolar_Fax(a_w_objetivo, ri_w_objetivo):
    """
    Interpola el valor del factor geométrico Fax a partir de la tabla A.3
    según los valores objetivo de a/w y ri/w.

    Parámetros:
    - a_w_objetivo: Relación a/w objetivo (profundidad del defecto / espesor de pared) [adimensional]
    - ri_w_objetivo: Relación ri/w objetivo (radio interno / espesor de pared) [adimensional]

    Retorna:
    - Valor interpolado de Fax [adimensional], correspondiente al punto (a/w, ri/w)
      utilizando interpolación bilineal sobre los datos de la tabla A.3.
    """
    puntos = tabla_A3[["a/w", "ri/w"]].values
    valores = tabla_A3["Fax"].values
    return float(griddata(puntos, valores, (a_w_objetivo, ri_w_objetivo), method='linear'))


def funcion_Fp_y_Fm_DP_A_5_7_8_CORREGIDA(a, c, w, ri, ro):
    """
    Calcula los factores de corrección geométrica Fp^DP y Fm^DP para una falla axial en el punto más profundo
    de penetración del defecto, según las secciones A.5.2.2.4.2, A.5.2.2.4.3 y A.5.2.2.4.4 del CSA N285.8-15.

    Aplica:
    - Ecuación A.5-7 para Fp^DP cuando 0.20 ≤ a/c ≤ 1.0 y 0.20 ≤ a/w ≤ 0.80
    - Ecuación A.5-8 para Fm^DP en el mismo dominio
    - Ecuaciones A.5-9 y A.5-10 para Fm^DP y Fp^DP cuando a/w tiende a 0
    - Ecuaciones A.5-12 y A.5-13 con interpolación de la tabla A.3 si a/c = 0

    Parámetros:
    - a : profundidad del defecto [mm]
    - c : semi-longitud axial del defecto [mm]
    - w : espesor de pared del tubo [mm]
    - ri: radio interno [mm]
    - ro: radio externo [mm]

    Retorna:
    - (Fp_DP, Fm_DP): tuple con los factores de corrección geométrica
    """

    a_c = a / c
    a_w = a / w
    ri_w = ri / w

    # Caso especial A.5.2.2.4.4: a/c = 0
    if a_c == 0:
        Fax = interpolar_Fax(a_w, ri_w)
        if Fax is None:
            return None, None
        Fp_DP = (2 * w * ro / (ro**2 - ri**2)) * Fax
        Fm_DP = 1.1261 + 0.6512 * a_w + 2.1839 * (a_w ** 2) + 3.7292 * (a_w ** 3)
        return Fp_DP, Fm_DP

    # Dominio válido para interpolación desde Tabla A.2. Caso A.5.2.2.4.2
    if 0.20 <= a_c <= 1.0 and 0.20 <= a_w <= 0.80:
        coef = interpolar_coeficientes_G(a_c, a_w)
        if None in (coef['G0'], coef['G1'], coef['G2'], coef['G3']):
            return None, None

        # Fp^DP - Ecuación A.5-7 con a/ri
        a_ri = a / ri
        parte_geom = 2 * coef['G0'] - 2 * a_ri * coef['G1'] + 3 * a_ri ** 2 * coef['G2'] - 4 * a_ri ** 3 * coef['G3']
        Fp_DP = (w * ro / (ro**2 - ri**2)) * parte_geom

        # Fm^DP - Ecuación A.5-8
        Fm_DP = coef['G0']
        return Fp_DP, Fm_DP

    # Dominio de a/w pequeño para ecuaciones A.5-9 a A.5-11
    elif 0.0 <= a_c <= 1.0 and abs(a_w) < 0.2:
        M = 1.13 - 0.07 * (a_c) ** 0.5
        Fp_DP = (2 * w * ro / (ro**2 - ri**2)) * M
        Fm_DP = M
        return Fp_DP, Fm_DP

    # Fuera del dominio válido
    return None, None


# Función para interpolar Tabla A.4
def interpolar_coeficientes_G_A4(a_c_objetivo, a_w_objetivo):
    puntos = tabla_A4[["a/c", "a/w"]].values
    G0 = float(griddata(puntos, tabla_A4["G0"].values, (a_c_objetivo, a_w_objetivo), method='linear'))
    G1 = float(griddata(puntos, tabla_A4["G1"].values, (a_c_objetivo, a_w_objetivo), method='linear'))
    G2 = float(griddata(puntos, tabla_A4["G2"].values, (a_c_objetivo, a_w_objetivo), method='linear'))
    G3 = float(griddata(puntos, tabla_A4["G3"].values, (a_c_objetivo, a_w_objetivo), method='linear'))
    return G0, G1, G2, G3

def interpolar_Fcr(a_w_objetivo, ri_w_objetivo):
    """
    Interpola el factor geométrico Fcr para una falla circunferencial de 360°
    según la Tabla A.5 del CSA N285.8-15.

    Parámetros:
    -----------
    a_w_objetivo : float
        Relación a/w de la falla circunferencial (profundidad sobre espesor).
    ri_w_objetivo : float
        Relación ri/w (radio interno sobre espesor de pared).

    Retorna:
    --------
    Fcr_interp : float
        Valor interpolado del factor geométrico Fcr.
    """
    puntos = tabla_A5[["a/w", "ri/w"]].values
    valores = tabla_A5["Fcr"].values
    Fcr_interp = griddata(puntos, valores, (a_w_objetivo, ri_w_objetivo), method='linear')

    if np.isnan(Fcr_interp):
        raise ValueError("Parámetros fuera del dominio de interpolación de la Tabla A.5")

    return float(Fcr_interp)

# Función principal para Fp^SP y Fm^SP
def funcion_Fp_y_Fm_SP_A_5_2_2_5(a, c, w, ri, ro):
    """
    Calcula los factores de corrección geométrica Fp^SP y Fm^SP para una falla axial
    en el punto de intersección con la superficie libre según CSA N285.8-15 Sección A.5.2.2.5.

    Devuelve Fp_SP, Fm_SP según:
    - A.5.2.2.5.2 si 0.20 <= a/c <= 1.0 y 0.20 <= a/w <= 0.80
    - A.5.2.2.5.3 si 0.20 <= a/c <= 1.0 y a/w tiende a 0
    """

    a_c = a / c
    a_w = a / w
    a_ri = a / ri

    # Caso 1: Dominio válido de interpolación (A.5.2.2.5.2)
    if 0.20 <= a_c <= 1.0 and 0.20 <= a_w <= 0.80:
        G0, G1, G2, G3 = interpolar_coeficientes_G_A4(a_c, a_w)
        parte_geom = 2 * G0 - 2 * a_ri * G1 + 3 * a_ri**2 * G2 - 4 * a_ri**3 * G3
        Fp_SP = (w * ro / (ro**2 - ri**2)) * parte_geom
        Fm_SP = G0
        return Fp_SP, Fm_SP

    # Caso 2: a/w → 0 (A.5.2.2.5.3)
    elif 0.0 <= a_c <= 1.0 and abs(a_w) < 0.2:          # desviacion de la norma 0.20 ≤ a/c ≤ 1.0
        M_f0 = (1.21 - 0.1 * (a_c) + 0.1 * (a_c)**4) * (a_c)**0.5
        Fp_SP = (2 * w * ro / (ro**2 - ri**2)) * M_f0
        Fm_SP = M_f0
        return Fp_SP, Fm_SP

    return None, None


def correccion_geometrica_circunferencial(a, b, w, ri, tipo="parcial"):
    """
    Calcula el factor de corrección geométrica Fc^DP para una falla circunferencial interna.

    Parámetros:
    -----------
    a : float
        Profundidad de la falla [mm]
    b : float
        Mitad de la longitud circunferencial de la falla [mm]
    w : float
        Espesor de pared del tubo [mm]
    ri : float
        Radio interno del tubo [mm]
    tipo : str
        Tipo de falla: 'parcial' o '360'

    Retorna:
    --------
    Fc_DP : float
        Factor de corrección geométrica Fc^DP
    """
    a_b = a / b
    a_w = a / w
    ri_w = ri / w

    if tipo == "360" and 0.125<=a_w<=0.750:
        # A.5.2.3.4.4 -> Fc^DP = Fcr
        return interpolar_Fcr(a_w, ri_w)

    elif 0.0<=a_b <= 0.20 and a_w < 0.2:
        # A.5.2.3.4.3 -> Fc^DP = M1
        M1 = 1.13 - 0.09 * a_b
        return M1

    elif 0.20 <= a_b <= 2.5 and 0.0 <= a_w <= 0.80:         # modificado limite superior de a_b por O11-1
        # A.5.2.3.4.2 -> Fórmula general
        M1 = 1.13 - 0.09 * a_b
        M2 = -0.54 + 0.89 / (0.2 + a_b)
        M3 = 0.5 - (1 / (0.65 + a_b)) + 14 * (1 - a_b) ** 24
        Fc_DP = M1 + M2 * (a_w) ** 2 + M3 * (a_w) ** 4
        return Fc_DP

    else:
        raise ValueError("Parámetros fuera del rango aplicable para A.5.2.3.4")


def calcular_KI_general(
    tipo_falla,              # "axial" o "circunferencial"
    punto,                   # "profundo", "superficie" o "corrección plastica"
    a,                       # profundidad del defecto [mm]
    c_o_b,                   # semi-longitud axial (c) o mitad de longitud circunferencial (b) [mm]
    p,                       # presión interna [MPa]
    r_i,                     # radio interno [mm]
    w,                       # espesor de pared [mm]
    sigma_res,              # tensión residual [MPa]
    r_o=None,                 # radio externo [mm] — requerido para fallas axiales
    tipo_circunf="parcial",   # "parcial" o "360" — requerido solo para circunferenciales
    T=None,                  # temperatura [°C] - requerido solo para corrección plástica
    sigma = 0
):
    """
    Calcula el K_I base para falla axial o circunferencial, en el punto profundo o de superficie libre.

    Parámetros:
    - tipo_falla : "axial" o "circunferencial"
    - punto      : "profundo" o "superficie"
    - a          : profundidad del defecto [mm]
    - c_o_b      : semilongitud c (axial) o b (circunferencial) [mm]
    - p          : presión interna [MPa]
    - r_i        : radio interno [mm]
    - w          : espesor de pared [mm]
    - sigma_res  : tensión residual [MPa]
    - ro         : radio externo [mm] (requerido solo para fallas axiales)
    - tipo_circunf : "parcial" o "360" (solo para fallas circunferenciales)

    Retorna:
    - K_I : factor de intensidad de tensiones [MPa·√mm]
    """

    if tipo_falla == "axial":
        c = c_o_b

        if r_o is None:
            raise ValueError("Para fallas axiales se requiere el radio externo (ro)")

        # Axial
        if punto != "corrección plastica":
            if punto == "profundo":
                Fp, Fm = funcion_Fp_y_Fm_DP_A_5_7_8_CORREGIDA(a, c, w, r_i, r_o)
            else:  # punto == "superficie"
                Fp, Fm = funcion_Fp_y_Fm_SP_A_5_2_2_5(a, c, w, r_i, r_o)

            return calcular_KI_axial(a, c, p, r_i, w, sigma_res, Fp, Fm)
        
        else:
            Fp, Fm = funcion_Fp_y_Fm_DP_A_5_7_8_CORREGIDA(a, c, w, r_i, r_o)
            K_I_Deep = calcular_KI_axial(a, c, p, r_i, w, sigma_res, Fp, Fm)
            Fp, Fm = funcion_Fp_y_Fm_SP_A_5_2_2_5(a, c, w, r_i, r_o)
            K_I_Sup = calcular_KI_axial(a, c, p, r_i, w, sigma_res, Fp, Fm)

            sigma_y = sigma_y_transverse_lb(T)

            return calcular_KI_con_plastic_zone(tipo_falla, a, c, p, r_i, w, sigma_res, K_I_Deep, K_I_Sup, sigma_y)

        # Circunferencial
    
    if tipo_falla == "circunferencial":
        b = c_o_b
        if punto != "corrección plastica":
            if punto == "profundo":
                Fc = correccion_geometrica_circunferencial(a, b, w, r_i, tipo=tipo_circunf)
            else:  # punto == "superficie"
                Fc = correccion_geometrica_superficie_libre_circunferencial(a, b, w)

            return calcular_KI_falla_circunferencial(sigma, a, b, Fc)
            
        else:
            Fc = correccion_geometrica_circunferencial(a, b, w, r_i, tipo=tipo_circunf)
            K_I_Deep = calcular_KI_falla_circunferencial(p, a, b, Fc)
            Fc = correccion_geometrica_superficie_libre_circunferencial(a, b, w)
            K_I_Sup = calcular_KI_falla_circunferencial(p, a, b, Fc)

            sigma_y = sigma_y_axial_lb(T)
            return calcular_KI_con_plastic_zone(tipo_falla, a, b, p, r_i, w, sigma_res, K_I_Deep, K_I_Sup, sigma_y)

    else:
        raise ValueError("tipo_falla debe ser 'axial' o 'circunferencial'")
    
    

def calcular_Ki_axial(T):
    """
    Calcula la tenacidad de iniciación de fractura radial-axial K_i^Axial
    para tubos de presión Zr-2.5 wt% Nb, según CSA N285.8-15 D.8.2
    (Ecuación D.8-1).

    Ecuación D.8-1:
        K_i^Axial [MPa√m] = 26.3 + 0.0227 * T(°C)

    Parámetros:
    - T: temperatura en °C


    Retorna:
    - K_i^Axial en MPa√m
    """

    return 26.3 + 0.022 * float(T)

def calcular_Ki_circ(T):
    """
    Calcula la tenacidad de iniciación de fractura radial-transversal K_i^Circ
    para tubos de presión Zr-2.5 wt% Nb, según CSA N285.8-15 D.8.3
    (Ecuación D.8-3).

    Ecuación D.8-3:
        K_i^Circ [MPa√m] = 29.9 + 0.0647 * T(°C)

    Parámetros:
    - T: temperatura en °C

    Retorna:
    - K_i^Circ en MPa√m
    """

    return 29.9 + 0.0647 * float(T)


def profile_scale(
    df_ref,                              # Tabla de referencia
    df_conv_p,                             # Tabla  con los datos para conversión de presión
    df_t,                                  # Tabla  con los datos para conversión de temperatura
):


    perfiles_por_tiempo = {}
    for j in range(len(df_conv_p)):
        Pconv_xi = []
        Tconv_xi = []
        X = []

        for i in range(len(df_ref)):
            
            # Interpolación de la presión

            Pref_xi= df_ref.loc[i,'Presión (MPa(a))']-0.1013
            Pconv_x0= df_conv_p.loc[j ,'Pressure (MPa(g))']
            Pref_x0= df_ref.loc[0,'Presión (MPa(a))']-0.1013
            Pconv_xi.append(Pref_xi*Pconv_x0/Pref_x0)
            X.append(df_ref.loc[i,'Distancia desde entrada F/C (m)'])

            # Interpolación de la temperatura

            for k in range(len(df_t)):
            
                if df_conv_p.loc[j ,'Time (sec)'] == df_t.loc[k ,'Time (sec)']:
                    Tref_xi= df_ref.loc[i,'Temperatura (°C)']
                    Tconv_x0= df_t.loc[k ,'Temperature (°C)']
                    Tref_x0= df_ref.loc[0,'Temperatura (°C)']
                    Tconv_xi.append(Tref_xi*Tconv_x0/Tref_x0)
                    break
                    

                elif df_conv_p.loc[j ,'Time (sec)'] < df_t.loc[k ,'Time (sec)']:
                    Tref_xi= df_ref.loc[i,'Temperatura (°C)']
                    Tconv_x0= df_t.loc[k-1 ,'Temperature (°C)']
                    Tref_x0= df_ref.loc[0,'Temperatura (°C)']
                    Tconv_x0= df_conv_p.loc[j ,'Time (sec)']*df_t.loc[k-1 ,'Temperature (°C)']/df_t.loc[k-1 ,'Time (sec)']
                    Tconv_xi.append(Tref_xi*Tconv_x0/Tref_x0)
                    break


        df_save = pd.DataFrame({
        "Distancia desde entrada F/C (m)": X,
        "Presión (MPa(g))": Pconv_xi,
        "Temperatura (°C)": Tconv_xi
        })

        t_key = df_conv_p.loc[j,'Time (sec)']  
        perfiles_por_tiempo[t_key] = df_save

    return perfiles_por_tiempo



