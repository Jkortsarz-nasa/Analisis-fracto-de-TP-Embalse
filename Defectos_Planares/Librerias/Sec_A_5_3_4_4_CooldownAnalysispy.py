# ==========================
# =================================

# Sección A.5.3.4.4: Analisis de la evolución de las indicaciones durante el cooldown por DHC

# ==========================
# =================================
# Librerias

import numpy as np
import pandas as pd

from scipy.interpolate import PchipInterpolator

# Tensiones
from Sec_A_4_4_SHoop import sigma_h_total
from Sec_A_4_5_SAxial import sigma_a_total

# Solubilidad / TTSSD
from Sec_A_5_3_2_T_T_SDD import calcular_T_T_SSD

# Stress Intensity Factor (KI)
from Sec_A_5_2_SIF import calcular_KI_general

# DHC / Fatiga
from Sec_A_5_3_4_FatigaDHC import (
    calcular_KIH,
    calc_Vr_mean,
    calc_Va_mean,
)

def cooldown_en_posicion(P_in, P_out, T_in, T_out, X, L, Intlet_o_Outlet) -> dict:
    """
    Interpola linealmente el transitorio de Cool-down entre entrada y salida
    para una posición axial X en un canal de largo L.
    Parámetros
    ----------
    P_in : float
        Presión en entrada [MPa].
    P_out : float
        Presión en salida [MPa].
    T_in : float
        Temperatura en entrada [°C].
    T_out : float
        Temperatura en salida [°C].
    X : float
        Posición axial [m].
    L : float
        Largo del canal [m].
    Inlet_o_Outlet : {"Inlet", "Outlet"}
        Referencia del perfil: si X=0 corresponde a entrada ("Inlet") o a salida ("Outlet").

    Retorna
    -------
    - "P_x": Presión interpolada en X [MPa]
    - "T_x": Temperatura interpolada en X [°C]
    """


    if not (0.0 <= X <= L):
        raise ValueError("X debe estar entre 0 y L")

    frac = X / L  # 0 -> entrada, 1 -> salida

    if (Intlet_o_Outlet == "Inlet"):
        # --- Temperatura ---

        T_x = T_in + frac * (T_out - T_in)

        # --- Presión ---

        P_x = P_in + frac * (P_out - P_in)

    elif (Intlet_o_Outlet == "Outlet"):

        T_x = T_out - frac * (T_out - T_in)

        # --- Presión ---

        P_x = P_out - frac * (P_out - P_in)

    return P_x, T_x

def evaluar_cooldown_y_DHC(
    df,                                 # DataFrame con defectos: 'Axial' (mm), 'Depth', 'Width', 'Length', 'Circunferencial o Axial'
    df_crecimiento,
    X, H_eq,                            # Perfiles axiales de posición (m) y H_eq por posición
    w_m,                                # Geometría por defecto (listas/arrays, indexados por defecto)
    sigma_h_res,                        # Tensión residual hoop (MPa)
    TP,                                 # Lista/dict por defecto con 'r_final','w_final'
    df_flux,
    transitorios_DHC_cooldown,          # Tablas con la curva de cooldown  
    sigma_a_res,                        # Tensión residual axial (MPa) usada en KI_general
    circ_ax,                            # Declarar si las indicación se analizarán como axiales o circunferenciales
    F,M                                 # F = fuerza axial aplicada en el TP [N], M = momento flexor aplicada en TP [N*m]
):



    # Genero el perfil presiones a los largo del TP a en el tiempo

    Perf_P_t_cd = {}

    corr_ax = 2341 # De tabla 10 : Embalse 2024 Fuel Channel Length Measurements
    L = 6176

    for i in range(len(df)):
        clave = f"Defecto {i+1}"   # Defecto 1, Defecto 2, ...
        Perf_P_t_cd[clave] = {
            "Indic.": f"{str(df.loc[i, "Channel"])}-{str(df.loc[i, "Indic."])}",           # Número de trzabilidad de la indicación de la forma "canal"-"N° de indicación del canal"
            "Axial o circunferencial": [],                                                 # Deja registro si la indicaión fue analizada como axial o cirunferencial
            "Tiempo_P [sec]": transitorios_DHC_cooldown["Cool-down_in"]['t_P_s'].copy(),   # Eje tiempo del transitorio aplicable para la indicación
            "Presión [MPa]": [],                                                           # Presión aplicada en la zona donde se encuentra la indicación a lo largo del tiempo
            "Tiempo_T [sec]": transitorios_DHC_cooldown["Cool-down_in"]['t_T_s'].copy(),   # copia para que no compartan referencia
            "Temperature (°C)": [],                                                        # Temperatura aplicada en la zona donde se encuentra la indicación a lo largo del tiempo
            "Tensión [MPa]": [],                                                           # Tensión (Hoop o axial segun corresponda) aplicada sobre la indicaicón a lo largo del tiempo
            "Porcentaje de reducción de tensión": [],                                      # % de reducción de la tensión a lo largo del tiempo 
            "Delta T [°C/min]": [],                                                        # variación de la temperatura °C/min a lo largo del tiempo
            "H_eq": [],                                                                    # Hidrogeno equivalente en la posición de la indicación                                                      
            "T_TSDD [°C]": [],                                                             # Evolución de la temperaturo de solubilización a los largo del cooldown (sin tener en cuenta la inercia termica)
            "T_TSDD_cd [°C]":[],                                                           # Cambio de T_TSDD a lo largo de la curva de enfriamiento (teniendo en cuenta la inercia termica)
            "KI": [],                                                                      # KI(Stress intensity factor) a lo largo del cooldown
            "KIH": [],                                                                     # KIH(Limit Stress intensity factor) a lo largo del cooldown
            "Limite_alcanzado" : False,                                                    # Se alcanzó la condición KI >= KIH en en el ciclo
            "flux[x1e16 n/m²/s]": []
        }


    P_t = [[] for _ in range(len(w_m))]
    T_t = [[] for _ in range(len(w_m))]


    # Determino el perfil de P y T en la posición de la indicación     


    for i in range(len(w_m)):

        X = df.loc[i,"Axial"]-corr_ax 
                
        for t in range(len(transitorios_DHC_cooldown["Cool-down_in"]['t_T_s'])):
            P_in = transitorios_DHC_cooldown["Cool-down_in"]['P_in_MPa'][t]
            P_out = transitorios_DHC_cooldown["Cool-down_out"]['P_in_MPa'][t]
            T_in = transitorios_DHC_cooldown["Cool-down_in"]['T_in_C'][t]
            T_out = transitorios_DHC_cooldown["Cool-down_out"]['T_in_C'][t]
            Inlet_o_Outlet = df.loc[i,"Reactor A-Face Configuration"]
            P_local_t, T_local_t   = cooldown_en_posicion(P_in, P_out, T_in, T_out, X, L, Inlet_o_Outlet)

            P_t[i].append(P_local_t)
            T_t[i].append(T_local_t)


    for idx, clave in enumerate(Perf_P_t_cd):
        # Convertimos los np.float64 a float nativos para evitar problemas al guardar
        perfil_presion = [float(x) for x in P_t[idx]]
        perfil_temperatura = [float(x) for x in T_t[idx]]

        
        Perf_P_t_cd[clave]["Presión [MPa]"] = perfil_presion
        Perf_P_t_cd[clave]["Temperature (°C)"] = perfil_temperatura



    # A.5.3.4.4.3

    # Calculo la variación de tensiones en funcion de como esta clasificada el defecto
    # Defecto axial -> Tensiones hoop
    # Defectos circunferenciales -> Tensiones axiales 



    for idx, clave in enumerate(Perf_P_t_cd):

        Perf_P_t_cd[clave]["H_eq"].append(H_eq[idx])   
        
        Perf_P_t_cd[clave]["T_TSDD [°C]"].append(calcular_T_T_SSD(Perf_P_t_cd[clave]["H_eq"][0]))

        if (circ_ax =='Axial'):

            Perf_P_t_cd[clave]["Axial o circunferencial"] = "Axial"

            for idx2 in range(len(Perf_P_t_cd[clave]["Temperature (°C)"])-1):

                Perf_P_t_cd[clave]["Delta T [°C/min]"].append(60*(Perf_P_t_cd[clave]["Temperature (°C)"][idx2]-Perf_P_t_cd[clave]["Temperature (°C)"][idx2+1])/(Perf_P_t_cd[clave]["Tiempo_T [sec]"][idx2+1]-Perf_P_t_cd[clave]["Tiempo_T [sec]"][idx2]))

            for idx2 in range(len(Perf_P_t_cd[clave]["Presión [MPa]"])):

                p=Perf_P_t_cd[clave]["Presión [MPa]"][idx2]

                Perf_P_t_cd[clave]["Tensión [MPa]"].append(sigma_h_total(p, TP[idx]['r_final'], TP[idx]['w_final'], sigma_h_res))

                Perf_P_t_cd[clave]["Porcentaje de reducción de tensión"].append(((Perf_P_t_cd[clave]["Tensión [MPa]"][0]-Perf_P_t_cd[clave]["Tensión [MPa]"][idx2])/Perf_P_t_cd[clave]["Tensión [MPa]"][0])*100)      #Equation A.5-37

        

        elif(circ_ax =='Circunferencial'):

            Perf_P_t_cd[clave]["Axial o circunferencial"] = "Circunferencial"

            for idx2 in range(len(Perf_P_t_cd[clave]["Temperature (°C)"])-1):

                Perf_P_t_cd[clave]["Delta T [°C/min]"].append(60*(Perf_P_t_cd[clave]["Temperature (°C)"][idx2]-Perf_P_t_cd[clave]["Temperature (°C)"][idx2+1])/(Perf_P_t_cd[clave]["Tiempo_T [sec]"][idx2+1]-Perf_P_t_cd[clave]["Tiempo_T [sec]"][idx2]))

            for idx2 in range(len(Perf_P_t_cd[clave]["Presión [MPa]"])):

                p=Perf_P_t_cd[clave]["Presión [MPa]"][idx2]
                
                Perf_P_t_cd[clave]["Tensión [MPa]"].append(sigma_a_total(
                p,            
                TP[idx]['r_final'],                     
                TP[idx]['w_final'],             
                F,
                M,            
                sigma_a_res))    

                Perf_P_t_cd[clave]["Porcentaje de reducción de tensión"].append(((Perf_P_t_cd[clave]["Tensión [MPa]"][0]-Perf_P_t_cd[clave]["Tensión [MPa]"][idx2])/Perf_P_t_cd[clave]["Tensión [MPa]"][0])*100)      #Equation A.5-37




        '''
            For bulk hydrogen equivalent concentrations corresponding to TTSSD greater than 250 °C, or a peak
            temperature less than TTSSD , the value of ΔTC
            DHC that is used in the evaluation shall be justified. The
            peak temperature is the maximum transient or operating temperature near the flaw location. A
            value of ΔTC
            DHC equal to zero may be used.
        '''

        for idx2 in range(len(Perf_P_t_cd[clave]["Temperature (°C)"])-1):

            if (Perf_P_t_cd[clave]["T_TSDD [°C]"][0] > 250) | (max(Perf_P_t_cd[clave]["Temperature (°C)"]) < (Perf_P_t_cd[clave]["T_TSDD [°C]"][0])):

                Delta_T_DHC=0
            
            else:

                if (Perf_P_t_cd[clave]["Delta T [°C/min]"][idx2]<1):
                    Delta_T_DHC=0

                else:
                    for j in range(len(Perf_P_t_cd[clave]["Presión [MPa]"])):
                        if (Perf_P_t_cd[clave]["Tiempo_P [sec]"][j]>Perf_P_t_cd[clave]["Tiempo_T [sec]"][idx2]):
                            j=j-1

                            break
                    if (Perf_P_t_cd[clave]["Porcentaje de reducción de tensión"][j]<=20):   # DIFERE DE A.5.3.4.4.3 (b) pero no deja casos sin resolver

                        Delta_T_DHC=15


                    elif (Perf_P_t_cd[clave]["Porcentaje de reducción de tensión"][j]>=20):

                        Delta_T_DHC=20
                
            if((Perf_P_t_cd[clave]["T_TSDD [°C]"][0]-Delta_T_DHC)>Perf_P_t_cd[clave]["Temperature (°C)"][idx2]):
                idx2 = idx2 - 1
                P_1 = Perf_P_t_cd[clave]["Presión [MPa]"][idx2]
                P_2 = Perf_P_t_cd[clave]["Presión [MPa]"][idx2+1]
                T_1 = Perf_P_t_cd[clave]["Temperature (°C)"][idx2]
                T_2 = Perf_P_t_cd[clave]["Temperature (°C)"][idx2+1]
                T_TSDD = Perf_P_t_cd[clave]["T_TSDD [°C]"][0]-Delta_T_DHC
                frac = (T_TSDD  - T_1) / (T_2 - T_1)

                p_DHC = P_1 + frac * (P_2 - P_1)
                sigma_DHC = sigma_a_total(p_DHC, TP[idx]['r_final'],TP[idx]['w_final'], F, M, sigma_a_res)
                break


            
            Perf_P_t_cd[clave]["T_TSDD_cd [°C]"].append(Perf_P_t_cd[clave]["T_TSDD [°C]"][0]-Delta_T_DHC)

            

            

        if (circ_ax =="Axial"):

            Perf_P_t_cd[clave]["KI"].append(calcular_KI_general("axial", "profundo", df_crecimiento["a_a"][idx], df_crecimiento["c"][idx]/2, p_DHC, TP[idx]['r_final'],TP[idx]['w_final'], sigma_a_res, TP[idx]['w_final']+ TP[idx]['r_final']))
            Perf_P_t_cd[clave]["KIH"]=(calcular_KIH(direccion='mixta',beta = df.iloc[idx]['Angule']))

        else:

            Perf_P_t_cd[clave]["KI"].append(calcular_KI_general("circunferencial", "profundo", df_crecimiento["a_c"][idx], df_crecimiento["b"][idx]/2, p_DHC, TP[idx]['r_final'],TP[idx]['w_final'], sigma_a_res, tipo_circunf="parcial", sigma = sigma_DHC))
            Perf_P_t_cd[clave]["KIH"]=(calcular_KIH(direccion='circunferencial'))
        
        if(Perf_P_t_cd[clave]["KI"][0]>Perf_P_t_cd[clave]["KIH"]):
            Perf_P_t_cd[clave]["Limite_alcanzado"] = True
            
        
        #return Perf_P_t_cd

        #determino el flujo neutronico para el defecto

        
        Perf_P_t_cd[clave]["flux[x1e16 n/m²/s]"].append(np.interp(((df.loc[idx,"Axial"])-corr_ax)/1000, df_flux["x (m)"], df_flux["Neutron Flux (x1e16 n/m²/s)"]))



        #A.5.3.4.4.7

        if(Perf_P_t_cd[clave]["KI"][0]>Perf_P_t_cd[clave]["KIH"]):
        
            continue      # NO APLICABLE
            #A.5.3.4.4.7 (a)

            '''
            a) Divide the cool-down curve of temperature versus time into time intervals of not more than
            10 minutes starting from the temperature TC
            DHC down to the temperature at the end of the cooldown.

            Para poder cumplir con este punto se divide Perf_P_t_cd[clave]["Tiempo_T [sec]] en intervalos de 10 minutos y hago una 
            interpolación PCHIP = Piecewise Cubic Hermite Interpolating Polynomial (similar a unspline cubico)
            y adunto esotos datos a dos listas nuevas para cada defecto
            '''



            # --- Malla uniforme cada 10 min ---
            paso = 600  # 10 min en segundos
            t10 = np.arange(min(Perf_P_t_cd[clave]["Tiempo_T [sec]"]), max(Perf_P_t_cd[clave]["Tiempo_T [sec]"]) + paso, paso)

            # --- Interpolación con PCHIP (conserva forma/monotonía, evita overshoot) ---
            f_T = PchipInterpolator(Perf_P_t_cd[clave]["Tiempo_T [sec]"], Perf_P_t_cd[clave]["Temperature (°C)"])
            T10 = f_T(t10)

            # --- (Opcional) Guardar en tu diccionario ---
            try:
                Perf_P_t_cd[clave]["Tiempo_T_interp [sec]"] = t10.tolist()
                Perf_P_t_cd[clave]["Temperature_spline (°C)"] = T10.tolist()
            except NameError:
                # Si aún no tenés definido Perf_P_t_cd/clave, no falle el script:
                pass
            


            '''
            b) For time interval i, calculate the stress intensity factor KIi using the instantaneous pressure loading
                at that point in the cool-down and the instantaneous flaw depth and instantaneous flaw halflength
                at that point in the cool-down."
            '''


            for idx2 in range(len(Perf_P_t_cd[clave]["Tiempo_T_interp [sec]"])):

                if(idx2==0):
                    Perf_P_t_cd[clave]["Width"].append(df.iloc[idx]['Width'])
                    Perf_P_t_cd[clave]["Lenght"].append(df.iloc[idx]['Lenght'])



                for j in range(len(Perf_P_t_cd[clave]["Presión [MPa]"])):
                    if (Perf_P_t_cd[clave]["Tiempo_P [sec]"][j]>Perf_P_t_cd[clave]["Tiempo_T_interp [sec]"][idx2]):
                        j=j-1
                        break
            
                if (circ_ax =="Axial"):

                    Perf_P_t_cd[clave]["KI_cd"].append(calcular_KI_general("axial", "profundo", df.iloc[idx]['Depth'], df.iloc[idx]['Width']/2, Perf_P_t_cd[clave]["Presión [MPa]"][j], TP[idx]['r_final'],TP[idx]['w_final'], sigma_a_res, TP[idx]['w_final']+ TP[idx]['r_final']))
                    Perf_P_t_cd[clave]["KIH_cd"]=(calcular_KIH(direccion='axial'))

                else:

                    Perf_P_t_cd[clave]["KI_cd"].append(calcular_KI_general("circunferencial", "profundo", df.iloc[idx]['Depth'], df.iloc[idx]['Width']/2, Perf_P_t_cd[clave]["Presión [MPa]"][j], TP[idx]['r_final'],TP[idx]['w_final'], sigma_a_res, tipo_circunf="parcial"))
                    Perf_P_t_cd[clave]["KIH_cd"]=(calcular_KIH(direccion='circunferencial'))


                '''
                c) When KIi from Item b) is greater than or equal to KIH, calculate DHC crack growth during the time
                interval i in accordance with Items d) to f).
                '''

                
                Delta_a_DHC=[]
                Delta_c_DHC=[]

                if (Perf_P_t_cd[clave]["KI_cd"]>Perf_P_t_cd[clave]["KIH_cd"]):
                    
                    Delta_a_DHC.append(calc_Vr_mean(Perf_P_t_cd[clave]["Tiempo_T_interp [sec]"][idx2], T_i, Perf_P_t_cd[clave]["flux[x1e16 n/m²/s]"][idx], 1)*paso*idx2)
                    Delta_c_DHC.append(calc_Va_mean(Perf_P_t_cd[clave]["Tiempo_T_interp [sec]"][idx2], T_i, Perf_P_t_cd[clave]["flux[x1e16 n/m²/s]"][idx],(t_APP*8.76+(paso*idx2)/(60*1000))))

                else:

                    Delta_a_DHC.append(0)
                    Delta_c_DHC.append(0)

                Perf_P_t_cd[clave]["Width"].append(Perf_P_t_cd[clave]["Width"][idx2]-Delta_a_DHC[idx2])
                Perf_P_t_cd[clave]["Lenght"].append(Perf_P_t_cd[clave]['Lenght'][idx2]-Delta_a_DHC[idx2])

    return Perf_P_t_cd
