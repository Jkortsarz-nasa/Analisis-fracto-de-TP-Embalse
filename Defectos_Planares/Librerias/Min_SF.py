# ==========================
# =================================

# Calculo de min SF para colapso plastico e iniciación de propagación de fisuras

# ==========================
# =================================
# Librerias

import pandas as pd
import math
import numpy as np
from Sec_D_3_4_2_PropMec import sigma_y_transverse_lb, sigma_y_axial_lb, sigma_u_transverse_lb, sigma_u_axial_lb
from Sec_A_5_2_SIF import calcular_KI_general, calcular_Ki_axial
from typing import List



def sigma_hoop_colapso_plastico(
    a: float,
    w: float,
    c: float,
    r_m: float,
    sigma_y: float,
    sigma_u: float,
) -> float:
    """
    CSA N285.8-15 A.5.5.2.2 — σ_h' para colapso plástico (falla axial pasante).

    (A.5-42)  σ_h' = σ_f * [ (1 - a/w) / (1 - (1/M)*(a/w)) ]
    (A.5-43)  σ_f  = (σ_y + σ_u)/2
    (A.5-44)  M    = [ 1 + 1.255*(c²/(r_m*w)) - 0.0135*(c²/(r_m*w))² ]^(1/2)

    Unidades: a,w,c,r_m [mm]; σ_y,σ_u [MPa] -> σ_h' [MPa]
    """

    sigma_f = 0.5 * (sigma_y + sigma_u)
    aw = a / w
    z = (c * c) / (r_m * w)
    M = math.sqrt(1.0 + 1.255 * z - 0.0135 * (z ** 2))

    denom = 1.0 - (aw / M)

    sigma_h_p = sigma_f * ((1.0 - aw) / denom)
    return sigma_h_p

def calcular_Kmin_SF_axial(
    df: pd.DataFrame,
    TP,
    sigma_a_res: float,
    tabla_PT: pd.DataFrame,
    df_crecimiento: pd.DataFrame,
    corr_ax_mm: float = 2341.0,
):
    """
    Calcula K_i/K_SF para cada indicación axial, interpolando P y T en función
    de la posición axial (referida al P/T inlet con corr_ax_mm).

    Parámetros
    ----------
    df : pd.DataFrame
        Debe contener la columna "Axial" [mm].
    TP : Sequence[dict]
        Por indicación debe contener:
        - "r_final" [mm]
        - "w_final" [mm]
    sigma_a_res : float
        Tensión axial residual σ_a^res [MPa] (según unidades de tu modelo).
    tabla_PT : pd.DataFrame
        Debe contener:
        - "Distancia desde entrada F/C (m)" [m]
        - "Presión (MPa(a))" [MPa(a)]
        - "Temperatura (°C)" [°C]
    df_crecimiento : pd.DataFrame
        Debe contener:
        - "a_a" [mm]
        - "c" [mm]
    corr_ax_mm : float
        Corrección axial [mm] para convertir Axial→distancia desde P/T inlet.

    Retorna
    -------
    list[float]
        Lista con K_i / K_SF por indicación.
    """
    x = tabla_PT["Distancia desde entrada F/C (m)"].to_numpy(dtype=float)
    P = tabla_PT["Presión (MPa(a))"].to_numpy(dtype=float)
    T = tabla_PT["Temperatura (°C)"].to_numpy(dtype=float)


    min_SF: List[float] = []

    for ind in range(len(TP)):
        #last = 1000000
        # for idx in range(len(tabla_PT["Distancia desde entrada F/C (m)"])):

        xq_m = (float(df.loc[ind, "Axial"]) - corr_ax_mm) / 1000.0  # [m]

        if (df['Reactor A-Face Configuration'][ind] == 'Outlet'):
            Pi=P[::-1]
            Pq = float(np.interp(xq_m, x, Pi))
            Tq = T[0]
            #T=T[::-1]

        else:
            Pq = float(np.interp(xq_m, x, P))
            Tq = T[0]


        # Pq = P[idx]
        # Tq = T[idx]

        a_a_mm = float(df_crecimiento.loc[ind, "a_a"])
        c_mm = float(df_crecimiento.loc[ind, "c"])

        K_SF = calcular_KI_general(
            "axial",
            "corrección plastica",
            a_a_mm,
            c_mm / 2.0,
            Pq,
            float(TP[ind]["r_final"]),
            float(TP[ind]["w_final"]),
            sigma_a_res,
            float(TP[ind]["w_final"]) + float(TP[ind]["r_final"]),
            T=Tq
        )

        K_i = float(calcular_Ki_axial(Tq))#calcular_Ki_circ(T)
        new = K_i / K_SF

        # if (new<last):
        #     last = new
        #min_SF.append(last)
        min_SF.append(K_i / K_SF)

    return min_SF

def sigma_hoop_colapso_plastico_axial(
    a: float,
    w: float,
    c: float,
    r_m: float,
    sigma_y: float,
    sigma_u: float,
) -> float:
    """
    CSA N285.8-15 A.5.5.2.2 — σ_h' para colapso plástico (falla axial pasante).

    (A.5-42)  σ_h' = σ_f * [ (1 - a/w) / (1 - (1/M)*(a/w)) ]
    (A.5-43)  σ_f  = (σ_y + σ_u)/2
    (A.5-44)  M    = [ 1 + 1.255*(c²/(r_m*w)) - 0.0135*(c²/(r_m*w))² ]^(1/2)

    Unidades: a,w,c,r_m [mm]; σ_y,σ_u [MPa] -> σ_h' [MPa]
    """

    sigma_f = 0.5 * (sigma_y + sigma_u)
    aw = a / w
    z = (c * c) / (r_m * w)
    M = math.sqrt(1.0 + 1.255 * z - 0.0135 * (z ** 2))

    denom = 1.0 - (aw / M)

    sigma_h_p = sigma_f * ((1.0 - aw) / denom)
    return sigma_h_p


def calcular_Sigmamin_SF_axial(
    df: pd.DataFrame,
    TP,
    F,
    M,
    sigma_a_res: float,
    tabla_PT: pd.DataFrame,
    df_crecimiento: pd.DataFrame,
    corr_ax_mm: float = 2341.0,
    
):
    """
    Calcula sigma_h / Sigma_SF para cada indicación axial, interpolando P y T en función
    de la posición axial (referida al P/T inlet con corr_ax_mm).

    Parámetros
    ----------
    df : pd.DataFrame
        Debe contener la columna "Axial" [mm].
    TP : Sequence[dict]
        Por indicación debe contener:
        - "r_final" [mm]
        - "w_final" [mm]
    sigma_a_res : float
        Tensión axial residual σ_a^res [MPa] (según unidades de tu modelo).
    tabla_PT : pd.DataFrame
        Debe contener:
        - "Distancia desde entrada F/C (m)" [m]
        - "Presión (MPa(a))" [MPa(a)]
        - "Temperatura (°C)" [°C]
    df_crecimiento : pd.DataFrame
        Debe contener:
        - "a_a" [mm]
        - "c" [mm]
    corr_ax_mm : float
        Corrección axial [mm] para convertir Axial→distancia desde P/T inlet.

    Retorna
    -------
    list[float]
        Lista con K_i / K_SF por indicación.
    """
    x = tabla_PT["Distancia desde entrada F/C (m)"].to_numpy(dtype=float)
    P = tabla_PT["Presión (MPa(a))"].to_numpy(dtype=float)
    T = tabla_PT["Temperatura (°C)"].to_numpy(dtype=float)


    min_SF: List[float] = []

    for ind in range(len(TP)):
       
        xq_m = (float(df.loc[ind, "Axial"]) - corr_ax_mm) / 1000.0  # [m]

        if (df['Reactor A-Face Configuration'][ind] == 'Outlet'):
            Pi=P[::-1]
            Pq = float(np.interp(xq_m, x, Pi))
            Ti = T[::-1]
            Tq = float(np.interp(xq_m, x, Ti))
            

        else:
            Pq = float(np.interp(xq_m, x, P))
            Tq = float(np.interp(xq_m, x, T))


        a_a_mm = float(df_crecimiento.loc[ind, "a_a"])
        c_mm = float(df_crecimiento.loc[ind, "c"])/2
        w = TP[ind]["w_final"]
        r_i = TP[ind]["r_final"]
        r_medio = (TP[ind]["w_final"]/2) + TP[ind]["r_final"]


        sigma_y= sigma_y_transverse_lb(Tq)
        sigma_u= sigma_y_transverse_lb(Tq) 

        Sigma_SF = sigma_hoop_colapso_plastico_axial(a_a_mm, w,c_mm,r_medio,sigma_y,sigma_u)
        
        sigma_h = Pq * (r_i/w) + Pq

        min_SF.append(Sigma_SF / sigma_h)
    return min_SF

