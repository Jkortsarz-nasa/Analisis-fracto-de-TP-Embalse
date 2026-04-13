# ==========================
# =================================

# Calculo de min SF para colapso plastico e iniciación de propagación de fisuras

# ==========================
# =================================
# Librerias

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from Sec_A_5_2_SIF import calcular_KI_general
from Sec_A_5_3_4_FatigaDHC import calcular_L_DHC_semi_eliptico, calcular_KIH, calc_Vr_mean, calc_Va_mean,calc_Vr_bounds, calc_Va_bounds
from Sec_A_4_5_SAxial import sigma_a_total
from Sec_A_5_3_2_T_T_SDD import calcular_H_eq, calcular_T_T_SSD
from Sec_A_5_3_4_FatigaDHC import crecimiento_fatiga




transitorios_fatiga = {
    
    "Warmup/Cooldown": {
    "t_T_s":  [0, 900, 3000, 39000, 42000, 150000, 150001, 236400, 237000, 322800, 324000, 327600, 331800, 334800, 337800, 338700, 338760, 340800, 358800],
    "T_in_C": [65.0, 90.0, 90.0, 260.0, 260.0, 260.0, 260.5, 260.5, 260.0, 260.0, 245.0, 165.0, 165.0, 130.0, 130.0, 90.0, 90.0, 50.0, 50.0],

    # Tiempos de PRESIÓN (min -> s)
    "t_P_s":  [0, 900, 3000, 39000, 42000, 150000, 150001, 236400, 237000, 322800, 324000, 327600, 331800, 334800, 337800, 338700, 338760, 340800, 358800],
    "P_in_MPa": [9.6, 9.6, 9.6, 9.6, 11.2, 11.2, 11.2, 11.2, 11.2, 11.2, 11.2, 10.5, 10.0, 9.5, 9.5, 9.5, 7.9, 7.9, 7.9], 
        
    "SL": ["A"],
    },


    "Startup/Shutdown": {
        # Tiempos de TEMPERATURA (s)
        "t_T_s":  [0, 10, 90, 400, 401, 460, 4400],
        "T_in_C": [260, 260, 269, 269, 269, 260, 260],

        # Tiempos de PRESIÓN (s)
        "t_P_s":  [0, 180, 240, 270, 400, 401, 640, 3568, 4400],
        "P_in_MPa": [11.30, 11.60, 11.60, 11.30, 11.30, 11.30, 5.93, 11.30, 11.30],

        # Nivel de servisio
        "SL": ["A"],
    },

    "Power Manoeuvring": {
        # Presión 
        "t_P_s":    [0, 120, 1168, 1768, 1828, 1888, 1918, 2500],
        "P_in_MPa": [11.30, 9.38, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        "t_T_s":  [0, 25, 1768, 1793, 2500],
        "T_in_C": [269, 263, 269, 269, 269],

        # Nivel de servisio
        "SL": ["A"],
    },


    "Refuelling": {
    # Presión (tal cual tabla)
        "t_P_s":    [0, 120, 1168, 1768, 1828, 1888, 1918, 2500],
        "P_in_MPa": [11.30, 9.38, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura (tal cual tabla)
        "t_T_s":  [0, 206, 1768, 1974, 2500],
        "T_in_C": [269, 184, 184, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    "Loss of Offsite Load / Turbine Trip": {
        # Presión (MPa(g)) vs tiempo (s)
        "t_P_s":    [0, 10, 20, 1330, 1858, 2038, 2098, 2128, 2400],
        "P_in_MPa": [11.30, 12.10, 8.90, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura (°C) vs tiempo (s)
        "t_T_s":  [0, 10, 35, 1938, 2040, 2400],
        "T_in_C": [269, 275, 260, 269, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    "Reactor Stepback": {
    # Presión MPa(g)
        "t_P_s":    [0, 20, 1857, 2400, 2580, 2640, 2670, 3600],
        "P_in_MPa": [11.30, 7.93, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura °C
        "t_T_s":  [0, 240, 2400, 3600],
        "T_in_C": [269, 254, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    "Reactor Trip from 100% Full Power": {
    # Presión MPa(g) vs tiempo (s)
        "t_P_s":    [0, 10, 12, 57, 2044, 2400, 2580, 2640, 2670, 3600],
        "P_in_MPa": [11.30, 7.66, 7.66, 7.66, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura °C vs tiempo (s)
        "t_T_s":  [0, 30, 80, 2400, 2640, 3600],
        "T_in_C": [269, 241, 260, 269, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    # Loss of Class IV Power 
    "Loss of Class IV Power": {
    # Presión MPa(g) vs tiempo (s)
        "t_P_s":    [0, 3, 7, 42, 48, 677, 998, 2447, 3077, 3257, 3317, 3347, 3600],
        "P_in_MPa": [11.30, 12.70, 12.70, 4.73, 8.00, 8.00, 9.80, 11.30, 11.30, 11.30, 11.60, 11.30, 11.30],

        # Temperatura °C vs tiempo (s)
        "t_T_s":  [0, 2, 5, 42, 677, 998, 3077, 3157, 3317, 3600],
        "T_in_C": [269, 269, 275, 245, 260, 260, 260, 269, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    # Figura C-13 — Reactor Overpower / Loss of Regulation (RIH)
    "Reactor Overpower / Loss of Regulation": {
    # Presión MPa(g) vs tiempo (s)
        "t_P_s":    [0, 2, 12, 15, 2251, 2850, 3030, 3090, 3120, 3600],
        "P_in_MPa": [11.30, 12.70, 5.93, 7.20, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura °C vs tiempo (s)
        "t_T_s":  [0, 5, 1230, 1250, 2850, 2930, 3600],
        "T_in_C": [269, 275, 235, 260, 260, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    # Figura C-14 — Loss of Feedwater Supply from 100% FP (RIH)
    "Loss of Feedwater Supply from 100% Full Power": {
        # Presión MPa(g) vs tiempo (s)
        "t_P_s":    [0, 8, 24, 34, 1643, 2147, 2327, 2387, 2417, 2700],
        "P_in_MPa": [11.30, 12.40, 12.40, 8.35, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura °C vs tiempo (s)
        "t_T_s":  [0, 8, 32, 50, 2147, 2227, 2700],
        "T_in_C": [269, 283, 283, 260, 260, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    "Total Loss of Feedwater Supply from 100% FP, cooldown with PHT pumps": {
    # Presión MPa(g) vs tiempo (s) — puntos tal cual la tabla
        "t_P_s":    [0, 8, 24, 34, 1561, 2147, 2162, 2614, 2629, 4661, 4671, 5271, 5286, 8286, 8316, 9766, 9826, 9836, 12461, 14643, 15118, 15298, 15358, 15388, 16398],
        "P_in_MPa": [11.30, 12.40, 12.40, 8.50, 11.30, 11.30, 9.20, 9.20, 7.20, 7.20, 6.19, 6.19, 2.90, 2.90, 0.00, 0.00, 5.90, 7.30, 7.30, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

        # Temperatura °C vs tiempo (s) — puntos tal cual la tabla
        "t_T_s":  [0, 8, 32, 50, 2147, 2614, 2747, 3561, 3611, 4661, 4961, 4991, 5291, 5951, 5981, 6766, 9286, 12461, 14818, 15118, 16398],
        "T_in_C": [269, 283, 283, 260, 260, 225, 215, 177, 170, 121, 121, 85, 100, 65, 100, 27, 27, 150, 260, 269, 269],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    "Rapid Cooldown": {
    # Presión MPa(g) vs tiempo (s)
        "t_P_s":    [0, 15, 20, 60, 120, 420, 900, 1200, 1510, 4847, 4907, 4917, 7343, 9725, 10200],
        "P_in_MPa": [11.30, 11.30, 4.93, 2.89, 2.39, 1.54, 1.29, 1.29, 0.00, 0.00, 5.90, 7.30, 7.30, 11.30, 11.30],

        # Temperatura °C vs tiempo (s)
        "t_T_s":  [0, 15, 20, 60, 240, 420, 900, 1200, 1830, 2430, 4907, 7543, 9900, 10200],
        "T_in_C": [260, 260, 240, 205, 160, 139, 100, 100, 100, 27, 27, 150, 260, 260],

        # Nivel de servisio
        "SL": ["B"],
    },

    

    "Total Loss of Feedwater Supply from 100% FP, Cooldown with SDC Pump": {
    # Presión MPa(g) vs tiempo (s) — puntos de la tabla
    "t_P_s":    [0, 8, 24, 34, 2147, 2167, 3182, 4982, 4997, 5297, 5312, 8882, 8892, 11518, 13700, 14000, 14180, 14240, 15470],
    "P_in_MPa": [11.30, 12.40, 12.40,  8.50, 11.30,  7.90,  5.90,  5.90,  2.90,  2.90,  0.00,  0.00,  7.30,  7.30, 11.30, 11.30, 11.60, 11.30, 11.30],

    # Temperatura °C vs tiempo (s) — puntos de la tabla
    "t_T_s":  [0, 8, 32, 50, 2147, 2477, 2777, 3077, 4320, 8882, 11518, 13875, 14000, 15470],
    "T_in_C": [269, 283, 283, 260,  260,   95,  105,   85,   27,   27,   150,   260,   269,   269],

    # Nivel de servisio
    "SL": ["B"],
    },

   "Loss of Pressure and Inventory Control – LRVs Fail Open": {
    # Presión MPa(g) vs tiempo (s)
    "t_P_s":    [0,   95,   153,  375,  480,  670,  989, 1018,
                 1180, 1340, 1500, 1743, 2043, 2223, 2283, 2313, 3000],
    "P_in_MPa": [11.30,10.50,10.43,8.09,10.30,10.44,10.52,10.47,
                 10.70,10.82,10.88,11.30,11.30,11.60,11.60,11.30,11.30],

    # Temperatura °C vs tiempo (s)
    "t_T_s":  [0, 595, 690, 990, 1033, 1098, 1500, 2043, 2117, 3000],
    "T_in_C": [269,264,264,261, 258, 262, 261, 261, 269, 269],

    # Nivel de servisio
    "SL": ["B"],
    },

    "Fuelling Machine (F/M) D2O System Failure": {
        # Presión MPa(g) vs tiempo (s)
        "t_P_s":    [0, 0.2, 4.2, 7.2, 12.2, 36, 2243, 2253, 3347, 3362, 4433, 4463,
                    7000, 7060, 7070, 9696, 11878, 12653, 12833, 12893, 12923, 14000],
        "P_in_MPa": [11.30, 11.00, 11.20, 11.00, 11.00, 5.90, 5.90, 5.19, 5.19, 2.90, 2.90, 0.00,
                    0.00, 5.90, 7.30, 7.30, 11.30, 11.30, 11.60, 11.60, 11.30, 11.30],

    # Temperatura °C vs tiempo (s)
        "t_T_s":  [0, 4, 6, 16, 21, 36, 2243, 2543, 2573, 3173, 3473, 3533,
                4433, 4493, 5093, 7060, 9696, 12053, 12653, 12733, 14000],
        "T_in_C": [269, 271, 271, 272, 265, 253, 150, 150,  90, 100, 100,  70,
                65,  55,  27,  27, 150, 260, 260, 269, 269],

    # Nivel de servisio
    "SL": ["B"],
    }   
}

### Tablas A-6 y A-7 del Stress Report

service_conditions = {
    # ============================================================
    # SERVICE CONDITIONS LEVEL A
    # ============================================================

    "Warmup/Cooldown": {
        "service_level": "A",
        "case": 1,
        # ---- Pressure & Temperature ----
        "P_max_psia": 1654, "P_min_psia": 15, "P_range_psia": 1639,
        "T_max_F": 500, "T_min_F": 81, "T_range_F": 419,
        # ---- Axial Force (lbf) ----
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        # ---- Bending Moment (lbf·in) ----
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        # ---- Torsion (lbf·in) ----
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Startup/Shutdown": {
        "service_level": "A",
        "case": 2,
        "P_max_psia": 1697, "P_min_psia": 875, "P_range_psia": 822,
        "T_max_F": 595, "T_min_F": 500, "T_range_F": 95,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Power Manoeuvring": {
        "service_level": "A",
        "case": 3,
        "P_max_psia": 1697, "P_min_psia": 1375, "P_range_psia": 322,
        "T_max_F": 595, "T_min_F": 565, "T_range_F": 30,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Refuelling": {
        "service_level": "A",
        "case": 4,
        "P_max_psia": 1697, "P_min_psia": 1375, "P_range_psia": 322,
        "T_max_F": 595, "T_min_F": 446, "T_range_F": 149,
        "F_axial_max_lbf": 18378, "F_axial_min_lbf": -17760, "F_axial_range_lbf": 36138,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    # ============================================================
    # SERVICE CONDITIONS LEVEL B
    # ============================================================

    # Case 5 applies to several transients
    "Loss of Offsite Load / Turbine Trip": {
        "service_level": "B",
        "case": 5,
        "P_max_psia": 1857, "P_min_psia": 701, "P_range_psia": 1156,
        "T_max_F": 613, "T_min_F": 500, "T_range_F": 113,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },
    "Loss of Class IV Power": {
        "service_level": "B",
        "case": 5,
        "P_max_psia": 1857, "P_min_psia": 701, "P_range_psia": 1156,
        "T_max_F": 613, "T_min_F": 500, "T_range_F": 113,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },
    "Loss of Feedwater Supply from 100% Full Power": {
        "service_level": "B",
        "case": 5,
        "P_max_psia": 1857, "P_min_psia": 701, "P_range_psia": 1156,
        "T_max_F": 613, "T_min_F": 500, "T_range_F": 113,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },
    "Loss of Pressure and Inventory Control – LRVs Fail Open": {
        "service_level": "B",
        "case": 5,
        "P_max_psia": 1857, "P_min_psia": 701, "P_range_psia": 1156,
        "T_max_F": 613, "T_min_F": 500, "T_range_F": 113,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Reactor Stepback": {
        "service_level": "B",
        "case": 6,
        "P_max_psia": 1697, "P_min_psia": 1165, "P_range_psia": 532,
        "T_max_F": 595, "T_min_F": 540, "T_range_F": 55,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Reactor Trip from 100% Full Power": {
        "service_level": "B",
        "case": 7,
        "P_max_psia": 1857, "P_min_psia": 875, "P_range_psia": 982,
        "T_max_F": 613, "T_min_F": 469, "T_range_F": 144,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Reactor Overpower / Loss of Regulation": {
        "service_level": "B",
        "case": 7,
        "P_max_psia": 1857, "P_min_psia": 875, "P_range_psia": 982,
        "T_max_F": 613, "T_min_F": 469, "T_range_F": 144,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Rapid Cooldown": {
        "service_level": "B",
        "case": 8,
        "P_max_psia": 1654, "P_min_psia": 15, "P_range_psia": 1639,
        "T_max_F": 500, "T_min_F": 81, "T_range_F": 419,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Total Loss of Feedwater Supply from 100% FP, cooldown with PHT pumps": {
        "service_level": "B",
        "case": 9,
        "P_max_psia": 1813, "P_min_psia": 15, "P_range_psia": 1798,
        "T_max_F": 610, "T_min_F": 81, "T_range_F": 529,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Total Loss of Feedwater Supply from 100% FP, Cooldown with SDC Pump": {
        "service_level": "B",
        "case": 9,
        "P_max_psia": 1813, "P_min_psia": 15, "P_range_psia": 1798,
        "T_max_F": 610, "T_min_F": 81, "T_range_F": 529,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    "Fuelling Machine (F/M) D2O System Failure": {
        "service_level": "B",
        "case": 9,
        "P_max_psia": 1813, "P_min_psia": 15, "P_range_psia": 1798,
        "T_max_F": 610, "T_min_F": 81, "T_range_F": 529,
        "F_axial_max_lbf": 10675, "F_axial_min_lbf": -10298, "F_axial_range_lbf": 20972,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },

    # ============================================================
    # TESTING CONDITION
    # ============================================================

    "Hydrostatic, Leak Test": {
        "service_level": "Testing",
        "case": 10,
        "P_max_psia": 2500, "P_min_psia": 15, "P_range_psia": 2485,
        "T_max_F": 81, "T_min_F": 81, "T_range_F": 0,
        "F_axial_max_lbf": 2655, "F_axial_min_lbf": -2277, "F_axial_range_lbf": 4932,
        "M_bend_max_lbf_in": 8500, "M_bend_min_lbf_in": 8500, "M_bend_range_lbf_in": 0,
        "T_max_lbf_in": 8078, "T_min_lbf_in": -3742, "T_range_lbf_in": 11820,
    },
}

##### conversión de tablas A-6 y A-7 en N y N*mm

# Factores de conversión
lbf_2_N = 4.448221615
lbin_2_Nmm = 112.984829
psia_2_MPa  = 0.006894757 

# valores maximos para nivel level A y B

F = 19027 * lbf_2_N
M = 8500 * lbin_2_Nmm

# Se asume que ya definiste `service_conditions`
# con las claves ..._lbf y ..._lbf_in como antes

for cond, data in service_conditions.items():
    for key in list(data.keys()):
        # Fuerzas axiales en lbf → N
        if key.endswith("_lbf"):
            new_key = key.replace("_lbf", "_N")
            data[new_key] = data[key] * lbf_2_N

        # Momentos y torsión en lbf·in → N·mm
        elif key.endswith("_lbf_in"):
            new_key = key.replace("_lbf_in", "_Nmm")
            data[new_key] = data[key] * lbin_2_Nmm

        # Presiones en psia → MPa (absoluta)
        elif key.endswith("_psia"):
            new_key = key.replace("_psia", "_MPa")
            data[new_key] = data[key] * psia_2_MPa

        # Temperaturas en °F → °C
        elif key.endswith("_F"):
            # Ej: T_max_F → T_max_C
            new_key = key.replace("_F", "_C")
            data[new_key] = (data[key] - 32.0) * (5.0 / 9.0)



# Concentraciones de hidrógeno según INFORMACIÓN REQUERIDA PARA ISI 2024 DE CANALES COMBUSTIBLES DE CNE

df_H = pd.DataFrame({
    "LATTICE SITE": ["B13", "E06", "J01", "K12", "L10", "M03", "N20", "O07", "O11", "O13", "P17", "T16"],
    "Informe de control H (max 5 PPM)": [3, 4, 5, 3, 3, 3, 2, 3, 3, 4, 3, 4]})


import numpy as np
import pandas as pd

# Table 9 Assessment of Pressure Tube Elongation de Embalse Assessment of Pressure Tube Dimensional Data from the 2024 August Periodic Inspection

df_table9 = pd.DataFrame({
    "Lattice Site": ["B13", "E06", "J01", "K12", "L10", "M03", "N20", "O07", "O11", "O13", "P17", "T16"],
    "Length of FC (mm)": [10837, 10823, 10851, 10823, 10824, 10823, 10822, 10823, 10823, 10823, 10823, 10823],
    "Distance Between BMs L1_PT (mm)": [6182, 6172, 6193, 6174, 6175, 6172, 6171, 6175, 6174, 6182, 6174, 6171]
})

df_table2= pd.DataFrame({
    "Lattice Site": ["B13", "E06", "J01", "K12", "L10", "M03", "N20", "O07", "O11", "O13", "P17", "T16"],
    "Fuel Channel Dishing Channel Type": ["B", "C", "A", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
    "Reactor A-Face Configuration": ["Inlet", "Inlet", "Outlet", "Outlet", "Inlet", "Inlet", "Inlet", "Inlet", "Inlet", "Inlet", "Outlet", "Inlet"]
})

# -------------------------
# 0.a) Corrección axial
# -------------------------

def Correccion_axial(Table_A_1, df_table9):


    Corr_axial = {}

    for i in range(len(df_table9)):

        Corr_axial[df_table9.loc[i,"Lattice Site"]] = {
        "factor de corrección axial [mm]": float,
        }

        Corr_axial[df_table9.loc[i,"Lattice Site"]]=((df_table9.loc[i,"Length of FC (mm)"]-Table_A_1.tail(1)["Distancia desde entrada F/C (m)"]*1000)/2)

    return Corr_axial

# -------------------------
# 0.b) Inciorporar la configuración al canal
# -------------------------

def agregar_configuracion_canal(df, df_filtrado):
    """
    Agrega al DataFrame principal las columnas:
      - 'Fuel Channel Dishing Channel Type'
      - 'Reactor A-Face Configuration'
    tomando los valores de df_filtrado en base a coincidencia del número de tubo.

    Parámetros:
    -----------
    df : pd.DataFrame
        DataFrame principal que contiene una columna con el número de tubo (por defecto 'Tube No.').
    df_filtrado : pd.DataFrame
        DataFrame con columnas ['Lattice Site', 'Fuel Channel Dishing Channel Type', 'Reactor A-Face Configuration'].
    
    Devuelve:
    ----------
    pd.DataFrame : el DataFrame original con las dos columnas nuevas agregadas.
    """

    # Hacemos un merge por coincidencia de 'Lattice Site' o 'Tube No.'
    df_actualizado = df.merge(
        df_filtrado[["Lattice Site", "Fuel Channel Dishing Channel Type", "Reactor A-Face Configuration"]],
        left_on="Channel",
        right_on="Lattice Site",
        how="left"
    ).drop(columns=["Lattice Site"])

    return df_actualizado



# -------------------------
# 1) Tensiones de referencia en la ubicación del defecto
# -------------------------
def calcular_tensiones_referencia(df, Table_A_1, r_i, r_o, w_m, F, M, sigma_a_res, sigma_h_res,
                                  sigma_a_total, sigma_h_total):
    Sigma_Ax, Sigma_Hoop = [], []
    corr_ax = 2341
    save = Table_A_1
    rotado = False
    for i in range(len(w_m)):
        # modifico Table_A_1 si el canal se midió desde el outlet 

        if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet') and (rotado == False):
            primera_columna = Table_A_1.iloc[:, [0]]
            resto = Table_A_1.iloc[:, 1:][::-1].reset_index(drop=True)
            Table_A_1 = pd.concat([primera_columna, resto], axis=1)
            rotado = True
        elif (df.loc[i, 'Reactor A-Face Configuration'] == 'Inlet') and (rotado == True):
            Table_A_1 = save
            rotado = False
            

        # hallar tramo espacial más cercano al Axial (mm) -> tabla en m
        axial_mm = float(df.loc[i, "Axial"])
        j = 0
        for j in range(len(Table_A_1)):
            if Table_A_1.loc[j, "Distancia desde entrada F/C (m)"]*1000 > (axial_mm - corr_ax):
                j = max(j-1, 0)
                break


        Mse=0
        Ms=0
        Fp=5545
        Fs=0

        p_a = np.interp(axial_mm - corr_ax, Table_A_1["Distancia desde entrada F/C (m)"].to_numpy() * 1000.0, Table_A_1["Presión (MPa(a))"])
        # Tensiones de referencia (sin escalado transitorio)
        Sigma_Ax.append(
            float(sigma_a_total(p_a, r_i[i], w_m[i], F, M, sigma_a_res))
        )
        Sigma_Hoop.append(
            float(sigma_h_total(p_a, r_i[i], w_m[i], sigma_h_res))
        )
    return Sigma_Ax, Sigma_Hoop


# -------------------------
# 2) H_eq, T_T_SSD, KI_deep, L_DHC
# -------------------------
def calcular_hidrogeno_y_KI(
    df, TP, Table_A_1, df_H, t_fin, CNE,
    calcular_H_eq, calcular_T_T_SSD,
    calcular_KI_general, calcular_KIH, calcular_L_DHC_semi_eliptico,
    sigma_h_res, sigma_a_res,
    df_lattice_col="Channel",   # nombre de la columna en df que identifica el canal/sitio
    EFPH2HH = 1.03 # HH/EFPH
):
    """
    Usa df_H para obtener H0 (ppm) por canal/sitio y calcular:
    - X_by_def[i]: vector X (m) usado para ese defecto i
    - H_eq_by_def[i]: perfil H_eq(X) para el H0 de su canal
    - T_T_SSD_by_def[i]: perfil TTSSD(X) correspondiente
    - KI_deep[i], L_DHC[i], P_loc[i], T_loc[i]
    Requiere que df tenga la columna `df_lattice_col` que matchee con df_H['LATTICE SITE'].
    """

    # --- 0) Mapa H0 por lattice site, con fallback a 5 ppm ---
    H_map = (df_H.set_index("LATTICE SITE")["Informe de control H (max 5 PPM)"]
                 .fillna(5)
                 .to_dict())

    # --- 1) Caché por sitio para no recalcular Deuterium varias veces ---
    cache_por_sitio = {}   # sitio -> {"X":X, "D":D, "H_eq":H_eq, "TTSSD":TTSSD}

    # Salidas por defecto
    X_by_def, H_eq_by_def, T_T_SSD_by_def = [], [], []
    KI_deep, L_DHC_a,  L_DHC_c= [], [], []
    P_loc, T_loc = [], []


    save = Table_A_1
    rotado = False
    for i in range(len(TP)):

        if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet') and (rotado == False):
            primera_columna = Table_A_1.iloc[:, [0]]
            resto = Table_A_1.iloc[:, 1:][::-1].reset_index(drop=True)
            Table_A_1 = pd.concat([primera_columna, resto], axis=1)
            rotado = True
        elif (df.loc[i, 'Reactor A-Face Configuration'] == 'Inlet') and (rotado == True):
            Table_A_1 = save
            rotado = False


        sitio = df.loc[i, df_lattice_col]
        H0_ppm = H_map.get(sitio, 5)

        # --- Deuterio y equivalentes ---
        X_i, D_i = CNE.Deuterium(H0_ppm, t_fin*EFPH2HH/8760)
        H_eq_i = calcular_H_eq(H0_ppm, D_i)
        TTSSD_i = [float(calcular_T_T_SSD(h)) for h in H_eq_i]

        # Guardar SIEMPRE, aunque el sitio se repita
        X_by_def.append(X_i)
        H_eq_by_def.append(H_eq_i)
        T_T_SSD_by_def.append(TTSSD_i)

        # --- 3) Ubicar presión/temperatura locales del HOT en la posición del defecto ---
        axial_mm = float(df.loc[i, "Axial"])
        j = 0
        for j in range(len(Table_A_1)):
            if Table_A_1.loc[j, "Distancia desde entrada F/C (m)"]*1000 > (axial_mm - 2362.2):
                j = max(j-1, 0)
                break
        P_loc_i = float(Table_A_1.loc[j, "Presión (MPa(a))"])
        T_loc_i = float(Table_A_1.loc[j, "Temperatura (°C)"])
        P_loc.append(P_loc_i)
        T_loc.append(T_loc_i)

        # --- 4) calcular KI profundo + L_DHC ---

        # Axial
        KI_i = float(calcular_KI_general(
            "axial", "profundo",
            df.iloc[i]['Depth'], df.iloc[i]['Length']/2.0,
            P_loc_i, TP[i]['r_final'], TP[i]['w_final'],
            sigma_a_res,
            TP[i]['w_final'] + TP[i]['r_final']
        ))
        KIH_a = float(calcular_KIH('axial'))

        # Circunferencial
        KI_i = float(calcular_KI_general(
            "circunferencial", "profundo",
            df.iloc[i]['Depth'], df.iloc[i]['Width']/2.0,
            P_loc_i, TP[i]['r_final'], TP[i]['w_final'],
            sigma_h_res,
            tipo_circunf="parcial"
        ))
        KIH_c = float(calcular_KIH('circunferencial'))

        KI_deep.append(KI_i)
        L_DHC_a.append(float(calcular_L_DHC_semi_eliptico(
            df.iloc[i]['Depth'], df.iloc[i]['Length']/2.0, KI_i, KIH_a
        )))
        L_DHC_c.append(float(calcular_L_DHC_semi_eliptico(
            df.iloc[i]['Depth'], df.iloc[i]['Length']/2.0, KI_i, KIH_c
        )))

    return (
        X_by_def,            # lista de arrays X por defecto i
        H_eq_by_def,         # lista de arrays H_eq(X) por defecto i
        T_T_SSD_by_def,      # lista de arrays TTSSD(X) por defecto i
        KI_deep, 
        L_DHC_a, L_DHC_c,
        P_loc, T_loc
    )


# 2.a) Calculo de perfil de Heq por defecto

# def calcular_H_eq_defecto(
#     df, TP, df_H, t_fin, EFPH2HH,
#     df_lattice_col="Channel"   # nombre de la columna en df que identifica el canal/sitio
# ):
#     """
#     Usa df_H para obtener H0 (ppm) por canal/sitio y calcular:
#     - X_by_def[i]: vector X (m) usado para ese defecto i
#     - H_eq_by_def[i]: perfil H_eq(X) para el H0 de su canal
#     """

#     corr_ax = 2341

#     # ---  Mapa H0 por lattice site, con fallback a 5 ppm ---
#     H_map = (df_H.set_index("LATTICE SITE")["Informe de control H (max 5 PPM)"]
#                  .fillna(5)
#                  .to_dict())

  
#     H_eq_by_def = []

    

#     for i in range(len(TP)):

#         axial_mm = float(df.loc[i, "Axial"])

#         # determino los perfiles para cada tubo
#         sitio = df.loc[i, df_lattice_col]
#         H0_ppm = H_map.get(sitio, 5)

#         # --- Deuterio y equivalentes ---
#         X_i, D_i = CNE.Deuterium(H0_ppm, t_fin*EFPH2HH/8760)
#         H_eq_i = calcular_H_eq(H0_ppm, D_i)

#         if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet'):
#             H_eq_i = H_eq_i[::-1]

#         for j in range(len(X_i)):
#             if X_i[j]*1000 > (axial_mm - corr_ax):
#                 j = max(j-1, 0)
#                 print(j)
#                 break
#         H_eq_by_def.append(H_eq_i[j])

#         # Guardar SIEMPRE, aunque el sitio se repita
#         #X_by_def.append(X_i)
#         #H_eq_by_def.append(H_eq_i)
    
#     return (
#         H_eq_by_def
#         )


# 2.b) Calculo  de TTSSD por defecto

def calcular_TTSSD_defecto(H_eq_by_def):
    
    return [float(calcular_T_T_SSD(h)) for h in H_eq_by_def]

# 2.c) Calculo de L_DHC limite para caso circunferencial y axial

def calcular_L_DHC_defecto (TP, Table_A_1, df, sigma_a_res, sigma_h_res, corr_ax_mm: float = 2341.0,):

    # Ubicar presión/temperatura locales de la tabla en la posición del defecto ---
        
        P_loc = []
        T_loc = []
        L_DHC_a = []
        L_DHC_c = []
        KI_c = []
        KI_a = []
        for i in range(len(TP)):
            axial_mm = float(df.loc[i, "Axial"])

            x = Table_A_1["Distancia desde entrada F/C (m)"] * 1000

            if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet'):
                Pq = Table_A_1["Presión (MPa(a))"][::-1]
                Tq = Table_A_1["Temperatura (°C)"][::-1]
            
            else:
                Pq = Table_A_1["Presión (MPa(a))"]
                Tq = Table_A_1["Temperatura (°C)"]

            x_obj = axial_mm - corr_ax_mm

            P_loc.append(np.interp(x_obj, x, Pq))
            T_loc.append(np.interp(x_obj, x, Tq))

            P_loc_i = P_loc[0]

            #  calcular KI profundo + L_DHC 

            # Axial
            KI_i_a = float(calcular_KI_general(
                "axial", "profundo",
                df.iloc[i]['Depth'], df.iloc[i]['Length']/2.0,
                P_loc_i, TP[i]['r_final'], TP[i]['w_final'],
                sigma_a_res,
                TP[i]['w_final'] + TP[i]['r_final']
            ))
            KIH_a = float(calcular_KIH(beta=df['Angule'][i], direccion='mixta'))

            if KI_i_a>4.5:
                print("stop")


            # Circunferencial
            KI_i_c = float(calcular_KI_general(
                "circunferencial", "profundo",
                df.iloc[i]['Depth'], df.iloc[i]['Width']/2.0,
                P_loc_i, TP[i]['r_final'], TP[i]['w_final'],
                sigma_h_res,
                tipo_circunf="parcial"
            ))
            KIH_c = float(calcular_KIH('circunferencial'))

            if KI_i_c>15:
                print("stop")

            KI_c.append(KI_i_c)
            KI_a.append(KI_i_a)

            L_DHC_a.append(float(calcular_L_DHC_semi_eliptico(
                df.iloc[i]['Depth'], df.iloc[i]['Length']/2.0, KI_i_a, KIH_a
            )))
            L_DHC_c.append(float(calcular_L_DHC_semi_eliptico(
                df.iloc[i]['Depth'], df.iloc[i]['Width']/2.0, KI_i_c, KIH_c
            )))

        return (L_DHC_a, L_DHC_c, P_loc, T_loc, KI_c, KI_a)


# -------------------------
# 3) Escalado de perfiles locales (presión, temperatura) en la ubicación del defecto
# -------------------------
def escalar_perfiles_locales(Table_A_1, df_p, df_t, df):
    """
    Devuelve P_t, T_t: listas de listas [defecto][tiempo]
    Requiere df_p con 'Time (sec)','Pressure (MPa(g))'
             df_t con 'Time (sec)','Temperature (°C)'
    """
    P_t = [[] for _ in range(len(df))]
    T_t = [[] for _ in range(len(df))]


    save = Table_A_1
    rotado = False
    for i in range(len(df)):

        if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet') and (rotado == False):
            primera_columna = Table_A_1.iloc[:, [0]]
            resto = Table_A_1.iloc[:, 1:][::-1].reset_index(drop=True)
            Table_A_1 = pd.concat([primera_columna, resto], axis=1)
            rotado = True
        elif (df.loc[i, 'Reactor A-Face Configuration'] == 'Inlet') and (rotado == True):
            Table_A_1 = save
            rotado = False
        
        axial_mm = float(df.loc[i, "Axial"])
        j = 0
        for j in range(len(Table_A_1)):
            if Table_A_1.loc[j, "Distancia desde entrada F/C (m)"]*1000 > (axial_mm - 2362.2):
                j = max(j-1, 0)
                break

        P_hot_local   = float(Table_A_1.loc[j, "Presión (MPa(a))"]) - 0.1013
        P_hot_entrada = float(Table_A_1.loc[0, "Presión (MPa(a))"]) - 0.1013
        for t in range(len(df_p["Time (sec)"])):
            P_header_t = float(df_p.loc[t, "Pressure (MPa(g))"])
            P_local_t  = P_hot_local * (P_header_t / P_hot_entrada) if P_hot_entrada != 0 else 0.0
            P_t[i].append(P_local_t)

        T_hot_local   = float(Table_A_1.loc[j, "Temperatura (°C)"])
        T_hot_entrada = float(Table_A_1.loc[0, "Temperatura (°C)"])
        for t in range(len(df_t["Time (sec)"])):
            T_header_t = float(df_t.loc[t, "Temperature (°C)"])
            T_local_t  = T_hot_local * (T_header_t / T_hot_entrada) if T_hot_entrada != 0 else 0.0
            T_t[i].append(T_local_t)

    return P_t, T_t


# -------------------------
# 4) Barrido por tiempos para p_max/p_min en cada defecto
# -------------------------
def extremos_presion_perfil(perfiles_tiempo, Table_A_1, df):
    """
    perfiles_tiempo: dict {t_key: DataFrame con columnas 'Distancia...' 'Presión (MPa(g))' 'Temperatura (°C)'}
    Retorna lista de dicts 'resultados' por defecto: p_max, p_min y T asociadas.
    """
    resultados = []

    
    rotado = False
    for i in range(len(df)):
        
        '''
        if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet') and (rotado == False):
            rotado = True
        elif (df.loc[i, 'Reactor A-Face Configuration'] == 'Inlet') and (rotado == True):
            rotado = False
        '''

        p_max = -np.inf
        p_min = +np.inf
        t_max = None
        t_min = None
        axial_mm = float(df.loc[i, "Axial"])

        for tkey in perfiles_tiempo.keys():
            # ubicar fila espacial para este defecto
            j = 0
            
            perfiles_tiempo_corr = perfiles_tiempo[tkey].copy()
            if (df.loc[i, 'Reactor A-Face Configuration'] == 'Outlet'):
                primera_columna = perfiles_tiempo[tkey].iloc[:, [0]]
                resto = perfiles_tiempo[tkey].iloc[:, 1:][::-1].reset_index(drop=True)
                perfiles_tiempo_corr = pd.concat([primera_columna, resto], axis=1)
            elif (df.loc[i, 'Reactor A-Face Configuration'] == 'Inlet'):
                perfiles_tiempo_corr = perfiles_tiempo[tkey].copy()

            for j in range(len(Table_A_1)):
                if perfiles_tiempo_corr.loc[j, "Distancia desde entrada F/C (m)"]*1000 > (axial_mm - 2362.2):
                    j = max(j-1, 0)
                    break

            P_val = float(perfiles_tiempo_corr.loc[j, "Presión (MPa(g))"])
            T_val = float(perfiles_tiempo_corr.loc[j, "Temperatura (°C)"])

            if P_val > p_max:
                p_max = P_val
                t_max = T_val
            if P_val < p_min:
                p_min = P_val
                t_min = T_val

        resultados.append({
            "defecto": i,
            "canal": df.loc[i,"Channel"],
            "presion_max": p_max,
            "temperatura_de_max": t_max,
            "presion_min": p_min,
            "temperatura_de_min": t_min
        })
    return resultados


# -------------------------
# 5) Fatiga en transitorio: ΔK y Δa
# -------------------------
def crecimiento_fatiga_transitorio(df, TP, service_condition, sigma_a_res, sigma_h_res, Delta_N):
    KI_max_a, KI_min_a, Delta_KI_a, Delta_a_fatiga_a = [], [], [], []
    KI_max_c, KI_min_c, Delta_KI_c, Delta_a_fatiga_c = [], [], [], []
    for i in range(len(TP)):

        # Axial
        KI_max_i_a = float(calcular_KI_general("axial", "profundo",
                                            df.iloc[i]['a_a'], df.iloc[i]['c']/2.0,
                                            service_condition['P_max_MPa'],
                                            TP[i]['r_final'], TP[i]['w_final'],
                                            sigma_a_res,
                                            TP[i]['w_final'] + TP[i]['r_final']))
        KI_min_i_a = float(calcular_KI_general("axial", "profundo",
                                            df.iloc[i]['a_a'], df.iloc[i]['c']/2.0,
                                            service_condition['P_min_MPa'],
                                            TP[i]['r_final'], TP[i]['w_final'],
                                            sigma_a_res,
                                            TP[i]['w_final'] + TP[i]['r_final']))

        
        sigma_max = sigma_a_total(service_condition['P_max_MPa'], TP[i]['r_final'], TP[i]['w_final'], service_condition['F_axial_max_N'], service_condition['M_bend_max_Nmm'], sigma_a_res)
        sigma_min = sigma_a_total(service_condition['P_min_MPa'], TP[i]['r_final'], TP[i]['w_final'], service_condition['F_axial_min_N'], service_condition['M_bend_min_Nmm'], sigma_a_res)


        # Circunferencial
        KI_max_i_c = float(calcular_KI_general("circunferencial", "profundo",
                                            df.iloc[i]['a_c'], df.iloc[i]['b']/2.0,
                                            service_condition['P_max_MPa'],
                                            TP[i]['r_final'], TP[i]['w_final'],
                                            sigma_h_res,
                                            tipo_circunf="parcial",
                                            sigma = sigma_max))
        KI_min_i_c = float(calcular_KI_general("circunferencial", "profundo",
                                            df.iloc[i]['a_c'], df.iloc[i]['b']/2.0,
                                            service_condition['P_min_MPa'],
                                            TP[i]['r_final'], TP[i]['w_final'],
                                            sigma_h_res,
                                            tipo_circunf="parcial",
                                            sigma = sigma_min))
        

        dK_a = KI_max_i_a - KI_min_i_a
        dK_c = KI_max_i_c - KI_min_i_c
        Delta_a_a = float(crecimiento_fatiga(dK_a, Delta_N))
        Delta_a_c = float(crecimiento_fatiga(dK_c, Delta_N))

        KI_max_a.append(KI_max_i_a)
        KI_min_a.append(KI_min_i_a)
        Delta_KI_a.append(dK_a)
        Delta_a_fatiga_a.append(Delta_a_a)

        KI_max_c.append(KI_max_i_c)
        KI_min_c.append(KI_min_i_c)
        Delta_KI_c.append(dK_c)
        Delta_a_fatiga_c.append(Delta_a_c)

    return Delta_a_fatiga_a, Delta_a_fatiga_c


# -------------------------
# 6) DHC en condiciones sostenidas (A.5.3.4.3)
# -------------------------
def dhc_condiciones_sostenidas(df, TP, df_flux, Table_A_1, KI_deep_a, KI_deep_c, calcular_KIH,
                                T_i, delta_t, TTSSD):
    """
    Devuelve dict dhc_results[k] con listas Δa, Δc por defecto.

    """
    dhc_results = {}
    corr_ax = 2341
    delta_t_seconds = (delta_t) * 3600
    
    Delta_a_a_DHC,Delta_a_c_DHC, Delta_c_DHC, Delta_b_DHC, T = [], [], [], [], []

    for i in range(len(TP)):
        # ubicar posición espacial para H_eq y flujo
        axial_mm = float(df.loc[i, "Axial"])

        if (df['Reactor A-Face Configuration'][i] == 'Outlet'):
            f_pair  = df_flux["Neutron Flux (x1e16 n/m²/s)"][::-1]
            T_pair  = Table_A_1["Temperatura (°C)"][::-1]

        else:
            f_pair = df_flux["Neutron Flux (x1e16 n/m²/s)"]
            T_pair= Table_A_1["Temperatura (°C)"]

        #temperatura en el punto

        print(i)
        x_query = (axial_mm - corr_ax)/1000.0
        x_pair  = Table_A_1["Distancia desde entrada F/C (m)"].values
        T.append(float(np.interp(x_query, x_pair, T_pair)))

        # flujo neutrónico local
        x_query = (axial_mm - corr_ax)/1000.0
        x_pair  = df_flux["x (m)"].values
        f_pair  = df_flux["Neutron Flux (x1e16 n/m²/s)"].values
        flux_i  = float(np.interp(x_query, x_pair, f_pair))

        TTSSD_i = TTSSD [i]

        #Axial

        KIH_a  = float(calcular_KIH((df.loc[i, "Angule"]), direccion='mixta'))

        if (TTSSD_i < T[i]) or (KI_deep_a[i] < KIH_a):
            da_a, dc = 0.0, 0.0
        else:
            # tiempo acumulado en miles de horas para Va_mean:
            tiempo_miles_h = delta_t_seconds / (3600.0 * 1000.0)
            da_a = float(calc_Vr_bounds(T[i], T_i, flux_i, 1)) * delta_t_seconds 
            dc = float(calc_Va_bounds(T[i], T_i, flux_i, tiempo_miles_h)) * delta_t_seconds

        #Circunferencial

        KIH_c  = float(calcular_KIH(direccion='circunferencial'))

        if (TTSSD_i < T[i]) or (KI_deep_c[i] < KIH_c):
            da_c, db = 0.0, 0.0
        else:
            # tiempo acumulado en miles de horas para Va_mean:
            tiempo_miles_h = delta_t_seconds / (3600.0 * 1000.0)
            da_c = float(calc_Vr_bounds(T[i], T_i, flux_i, 1)) * delta_t_seconds
            db = float(calc_Va_bounds(T[i], T_i, flux_i, tiempo_miles_h)) * delta_t_seconds

        Delta_a_a_DHC.append(da_a)
        Delta_a_c_DHC.append(da_c)
        Delta_c_DHC.append(dc)
        Delta_b_DHC.append(db)

    dhc_results = {
        "Delta_a_a_DHC [m]": Delta_a_a_DHC,
        "Delta_a_c_DHC [m]": Delta_a_c_DHC,
        "Delta_c_DHC [m]": Delta_c_DHC,
        "Delta_b_DHC [m]": Delta_b_DHC
    }
    return dhc_results


# -------------------------
# 7) Clasificación Axial / Circunferencial
# -------------------------
def clasificar_defectos_axial_circ(df):
    df = df.copy()
    df["Circunferencial o Axial"] = np.where(df["Length"] >= df["Width"], "Axial", "Circunferencial")
    return df
 