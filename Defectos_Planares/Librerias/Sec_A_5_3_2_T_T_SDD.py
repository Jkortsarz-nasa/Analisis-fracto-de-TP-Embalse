# ==========================
# =================================

# La presente librería calcula la temperatura de precipitacion de hidruros para el calculo de DHC según la sección sección A.5.3.2 de la norma CSA 285

# ==========================
# =================================
# =================================
# ========================
# Librerías

import numpy as np
import math



# =============================
# A.5.3.2.2 - Bulk Hydrogen Equivalent Concentration
# =============================

def calcular_H_eq(H_i, D):
    """
    A.5-26: Bulk hydrogen equivalent concentration at flaw location.

    H_eq = H_i + D / 2

    Parámetros:
    - H_i: concentración inicial de hidrógeno [ppm]
    - D: deuterio acumulado al final del período de evaluación [ppm]

    Retorna:
    - H_eq: concentración equivalente de hidrógeno [ppm]
    """
    return H_i + D / 2


# =============================
# A.5.3.2.3 - T_T_SSD (Temperature of terminal Solid Solubility for Dissolution)
# =============================

def calcular_T_T_SSD(H_eq, C_D=8.19e4, Q_D=34500, R=8.314):       # Constantes se obtuvieron de la sección  D.2.2
    """
    Cálculo corregido de T_SSD según la ecuación A.5-27, con coeficientes extraídos del Anexo D.2.2.

    Fórmula:
    T_SSD = (-Q_D / (R * ln(H_eq / C_D))) - 273

    Parámetros:
    - H_eq: concentración equivalente de hidrógeno [ppm]
    - C_D: constante pre-exponencial del TSSD [ppm] (Anexo D: 8.19e4)
    - Q_D: energía de activación para disolución de hidrógeno [J/mol] (Anexo D: 34500)
    - R: constante universal de los gases [J/mol·K] (8.314)

    Retorna:
    - T_SSD: temperatura en °C
    """
    if np.all(H_eq <= 0) or np.all(C_D <= 0):
        raise ValueError("H_eq y C_D deben ser mayores a 0 para cálculo logarítmico.")
    return (-Q_D / (R * math.log(H_eq / C_D))) - 273