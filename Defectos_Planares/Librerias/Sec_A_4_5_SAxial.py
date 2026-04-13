# ==========================
# =================================

# La presente librería calcula las tensiones Axiales según la sección A.4.5 de la norma CSA 285

# ==========================
# =================================
# =================================
# ========================
# Librerías

import math
import numpy as np





# ========================
# A.4.5 - Axial Stress
# ========================

def sigma_a_p(p, r_i, w):
    """
    A.4-20: Tensión axial primaria nominal debida a presión interna

    Fórmula:
        σ_a^p = p * r_i / (2 * w)

    Parámetros:
    -----------
    p : float
        Presión interna [MPa]
    r_i : float
        Radio interno del tubo de presión [mm]
    w : float
        Espesor de pared [mm]

    Retorna:
    --------
    σ_a^p : float
        Tensión axial primaria [MPa]
    """
    return p * r_i / (2 * w)

def sigma_a_pb(Mp, Mse, r_o, r_i):
    """
    A.4-21: Tensión axial por momento flector primario

    Fórmula:
        σ_a^pb = (Mp + Mse) * r_o / I
        I = π/4 * (r_o^4 - r_i^4)

    Parámetros:
    -----------
    Mp : float
        Momento por peso muerto [N·mm]
    Mse : float
        Momento sísmico [N·mm]
    r_o : float
        Radio externo del tubo [mm]
    r_i : float
        Radio interno del tubo [mm]

    Retorna:
    --------
    σ_a^pb : float
        Tensión axial por flexión primaria [MPa]
    """
    I = (math.pi / 4) * (r_o**4 - r_i**4)
    return (Mp + Mse) * r_o / I

def sigma_a_sb(Ms, r_o, r_i):
    """
    A.4-23: Tensión axial secundaria por flexión

    Fórmula:
        σ_a^sb = Ms * r_o / I

    Parámetros:
    -----------
    Ms : float
        Momento secundario (ej. por desplazamientos térmicos) [N·mm]
    r_o : float
        Radio externo [mm]
    r_i : float
        Radio interno [mm]

    Retorna:
    --------
    σ_a^sb : float
        Tensión axial secundaria por flexión [MPa]
    """
    I = (math.pi / 4) * (r_o**4 - r_i**4)
    return Ms * r_o / I

def sigma_a_PF(Fp, r_o, r_i):
    """
    A.4-24: Tensión axial primaria por fuerza axial

    Fórmula:
        σ_a^PF = Fp / A
        A = π * (r_o^2 - r_i^2)

    Parámetros:
    -----------
    Fp : float
        Fuerza axial primaria [N]
    r_o : float
        Radio externo [mm]
    r_i : float
        Radio interno [mm]

    Retorna:
    --------
    σ_a^PF : float
        Tensión axial primaria por fuerza axial [MPa]
    """
    A = math.pi * (r_o**2 - r_i**2)
    return Fp / A

def sigma_a_SF(Fs, r_o, r_i):
    """
    A.4-26: Tensión axial secundaria por fuerza axial

    Fórmula:
        σ_a^SF = Fs / A

    Parámetros:
    -----------
    Fs : float
        Fuerza axial secundaria (ej. por expansión térmica) [N]
    r_o : float
        Radio externo [mm]
    r_i : float
        Radio interno [mm]

    Retorna:
    --------
    σ_a^SF : float
        Tensión axial secundaria por fuerza axial [MPa]
    """
    A = math.pi * (r_o**2 - r_i**2)
    return Fs / A

def sigma_a_total(
    p,            
    r_i,                     
    w,            
    F,
    M,            
    sigma_a_res    
):
    """
    Calcula la tensión axial total (A.4-28), sumando todas las componentes:
    presión interna, momentos flectores, fuerzas axiales y esfuerzos residuales.

    Parámetros:
    -----------
    p : float
        Presión interna [MPa]
    r_i : float
        Radio interno del tubo [mm]
    r_o : float
        Radio externo del tubo [mm]
    w : float
        Espesor de pared [mm]
    Mp : float
        Momento por peso muerto [N·mm]
    Mse : float
        Momento sísmico [N·mm]
    Ms : float
        Momento secundario [N·mm]
    Fp : float
        Fuerza axial primaria [N]
    Fs : float
        Fuerza axial secundaria [N]
    sigma_a_res : float
        Tensión axial residual [MPa]

    Retorna:
    --------
    σ_a_total : float
        Tensión axial total [MPa]
    """
    r_o = r_i + w

    A = np.pi * ((r_o**2)-(r_i**2))
    I = 0.25 * np.pi * ((r_o**4)-(r_i**4))

    sigma = ((p*r_i)/(2*w)) + (F/A) + ((M*r_o)/I)
    return (
        p +
        sigma +
        sigma_a_res
    )