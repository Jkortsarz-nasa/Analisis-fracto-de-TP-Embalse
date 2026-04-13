


# ==========================
# =================================

# La presente librería calcula las tensiones Hoop según la sección A.4.4 de la norma CSA 285

# ==========================
# =================================
# =================================
# ========================
# Librerías

import math
import numpy as np



# ========================
# A.4.4 - Hoop Stress
# ========================

def sigma_h_p(p, r_i, w):
    """
    A.4-16: Nominal primary hoop stress
    σ_h^p = p * r_i / w

    Parámetros:
    -----------

    - p: internal pressure [MPa]
    - r_i: inner radius of the pressure tube [mm]
    - w: wall thickness at flaw location [mm]

    
    Retorna:
    --------

    -σ_h_p : Tensión circunferencial primaria nominal [MPa]

    Returns: σ_h^p [MPa]
    """
    return p * r_i / w

def sigma_h_total(p, r_i, w, sigma_h_res):
    """
    A.4-17: Total nominal hoop stress including pressure on flaw face and residual stress
    σ_h = p * (r_i / w + 1) + σ_h^res

    Parámetros:
    -----------
    p : float
        Presión interna [MPa]
    r_i : float
        Radio interno del tubo de presión [mm]
    w : float
        Espesor de pared en la ubicación de la falla [mm]
    sigma_h_res : float
        Tensión circunferencial residual [MPa]

    Retorna:
    --------

    - σ_h : Tensión circunferencial total nominal [MPa]

    Returns: σ_h [MPa]
    """
    return p * (r_i / w + 1) + sigma_h_res

def sigma_h_fitting(p, r_i, w, chi_EF, sigma_h_res):
    """
    A.4-19: Hoop stress in rolled joint region (with fitting effect)
    σ_h = p * (1 + r_i / (w * χ_E/F)) + σ_h^res

    Parámetros:
    -----------
    p : float
        Presión interna [MPa]
    r_i : float
        Radio interno del tubo de presión [mm]
    w : float
        Espesor de pared en la ubicación de la falla [mm]
    chi_EF : float
        Factor adimensional χ_E/F (relación entre módulo de elasticidad del tubo y del fitting)
    sigma_h_res : float
        Tensión circunferencial residual [MPa]

    - chi_EF: χ_E/F ratio [dimensionless]

    Retorna:
    --------
    σ_h : float
        Tensión circunferencial total en zona de fitting [MPa]
    Returns: σ_h [MPa]
    """
    return p * (1 + r_i / (w * chi_EF)) + sigma_h_res


def momento_flector_primario(MP0: float,
                                   MP_EOL: float,
                                   t: float,
                                   t0_PT: float,
                                   t_DL: float) -> float:
    """
    A.4-13: Relajación lineal del momento flector primario (peso propio)
    M^P(t) = M^P0 - (M^P0 - M^P_EOL) * ((t - t0^PT) / t_DL)

    Parámetros:
    -----------
    MP0 : float
        Momento flector primario por peso propio al inicio de la vida de diseño del PT [kN·m]
    MP_EOL : float
        Momento flector primario por peso propio al final de la vida de diseño del PT (EOL) [kN·m]
    t : float
        Tiempo de evaluación (vida de reactor acumulada) [años EFPY u otra unidad consistente]
    t0_PT : float
        Vida del reactor al momento de instalación del PT [misma unidad de t]
    t_DL : float
        Vida de diseño del PT [misma unidad de t]. Debe ser > 0.

    Notas:
    ------
    - Modelo lineal sin recorte: si t < t0_PT, el resultado extrapola por encima de MP0;
      si t > t0_PT + t_DL, extrapola por debajo de MP_EOL.

    Retorna:
    --------
    M_P : float
        Momento flector primario por peso propio a tiempo t [kN·m]
    """
    if t_DL == 0.0:
        raise ValueError("t_DL debe ser estrictamente mayor que 0.")

    frac = (t - t0_PT) / t_DL
    M_P = MP0 - (MP0 - MP_EOL) * frac
    return M_P


def momento_flector_sismico(Mv_DBE: float, Mh_DBE: float) -> float:
    """
    A.4-14: Combinación RMS del momento flector primario por evento sísmico
    M^DBE = [ (M_v^DBE)^2 + (M_h^DBE)^2 ]^(1/2)

    Parámetros:
    -----------
    Mv_DBE : float
        Momento flector primario en el plano vertical por evento sísmico [kN·m]
    Mh_DBE : float
        Momento flector primario en el plano horizontal por evento sísmico [kN·m]

    Notas:
    ------
    - Se emplea combinación raíz-de-la-suma-de-cuadrados (RMS). El signo de cada
      componente no afecta el resultado (se eleva al cuadrado).

    Retorna:
    --------
    M_DBE : float
        Momento flector primario combinado por evento sísmico [kN·m]
    Returns: M_DBE [kN·m]
    """
    return (Mv_DBE**2 + Mh_DBE**2) ** 0.5

def momento_flector_secundario(MS0: float,
                               MS_EOL: float,
                               t: float,
                               t0_PT: float,
                               t_DL: float) -> float:
    """
    A.4-15: Relajación exponencial del momento flector secundario (desalineación)
    M^S(t) = M^S_0 * exp( ln(M^S_EOL / M^S_0) * ((t - t0^PT) / t_DL) )

    Parámetros:
    -----------
    MS0 : float
        Momento flector secundario al inicio de la vida de diseño del PT [kN·m]
    MS_EOL : float
        Momento flector secundario al final de la vida de diseño del PT (EOL) [kN·m]
    t : float
        Tiempo de evaluación (vida de reactor acumulada) [años EFPY u otra unidad consistente]
    t0_PT : float
        Vida del reactor al momento de instalación del PT [misma unidad de t]
    t_DL : float
        Vida de diseño del PT [misma unidad de t]. Debe ser > 0.

    Notas:
    ------
    - Modelo exponencial: interpola suavemente entre MS0 (en t = t0_PT) y MS_EOL (en t = t0_PT + t_DL).
    - Requiere que MS0 y MS_EOL tengan el mismo signo y que MS0 ≠ 0 para que ln(MS_EOL/MS0) sea válido.
      Si alguno fuese 0 o de signo distinto, ajustar el modelo o usar un ε pequeño según criterio de ingeniería.

    Retorna:
    --------
    M_S : float
        Momento flector secundario a tiempo t [kN·m]
    Returns: M_S [kN·m]
    """
    if t_DL == 0.0:
        raise ValueError("t_DL debe ser estrictamente mayor que 0.")
    if MS0 == 0.0:
        raise ValueError("MS0 no puede ser 0 (ln(MS_EOL/MS0) indefinido).")
    if (MS0 > 0 and MS_EOL < 0) or (MS0 < 0 and MS_EOL > 0):
        raise ValueError("MS0 y MS_EOL deben tener el mismo signo para la interpolación exponencial.")

    frac = (t - t0_PT) / t_DL
    k = math.log(MS_EOL / MS0)
    M_S = MS0 * math.exp(k * frac)
    return M_S
