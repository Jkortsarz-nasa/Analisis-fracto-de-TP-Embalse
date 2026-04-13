# =================================
# ==========================
# A.4.2 La peresente librería calcula las variaciones del espesor y radio interno de los TP según la sección A.4.2 de la norma CSA 285
# la librería debe ejecutarse desde la función evaluar_A_4_2, este devuelve un diccionario con las datos obtenidos para cada indicación analizada
# ==========================
# =================================


# ==========================
# A.4.2.2.2 Equation A.4-1
# ==========================
def calcular_w0(w_m, C_a_SE, epsilon_ct_SE):
    """
    Ecuación A.4-1: Espesor inicial efectivo de pared del tubo de presión.
    w0 = (w_m + C_a^SE) / (1 - ε_ct^SE)

    Parámetros:
    - w_m: espesor medido al inicio del período de evaluación [mm]
    - C_a_SE: adelgazamiento por corrosión y desgaste hasta el inicio del período [mm]
    - epsilon_ct_SE: deformación de creep hasta el inicio del período [adimensional]
    """
    if epsilon_ct_SE >= 1:
        raise ValueError("ε_ct^SE debe ser menor que 1 para evitar división por cero.")
    return (w_m + C_a_SE) / (1 - epsilon_ct_SE)

# ==========================
# A.4.2.2.2 Equation A.4-2
# ==========================
def calcular_epsilon_ct_SE(epsilon_ct_total, t_SE, t_PT_instalacion, t_PT_diseño):
    """
    Ecuación A.4-2: Deformación de creep de adelgazamiento de pared hasta el inicio del período.
    ε_ct^SE = ε_ct_total * (t_SE - t_PT_instalacion) / t_PT_diseño

    Parámetros:
    - epsilon_ct_total: deformación de creep durante la vida de diseño [adim.]
    - t_SE: tiempo al inicio del período de evaluación [años]
    - t_PT_instalacion: edad del tubo al momento de instalación [años]
    - t_PT_diseño: vida de diseño del tubo [años]
    """
    return epsilon_ct_total * ((t_SE - t_PT_instalacion) / t_PT_diseño)

# ==========================
# A.4.2.2.2 Equation A.4-3
# ==========================
def calcular_Ca_SE(Ca_total, t_SE, t_PT_instalacion, t_PT_diseño):
    """
    Ecuación A.4-3: Adelgazamiento por corrosión hasta el inicio del período.
    C_a^SE = C_a_total * (t_SE - t_PT_instalacion) / t_PT_diseño

    Parámetros:
    - Ca_total: corrosión esperada durante vida de diseño [mm]
    - t_SE: tiempo al inicio del período de evaluación [años]
    - t_PT_instalacion: edad del tubo a la instalación [años]
    - t_PT_diseño: vida de diseño del tubo [años]
    """
    return Ca_total * ((t_SE - t_PT_instalacion) / t_PT_diseño)

# ==========================
# A.4.2.2.3 Equation A.4-4
# ==========================
def calcular_epsilon_ct_EV(epsilon_ct_total, t_fin, t_PT_instalacion, t_PT_diseño):
    """
    Ecuación A.4-4: Deformación de creep de adelgazamiento de pared hasta el final del período.
    ε_ct^EV = ε_ct_total * (t_fin - t_PT_instalacion) / t_PT_diseño

    Parámetros:
    - epsilon_ct_total: deformación total esperada [adim.]
    - t_fin: tiempo al final del período de evaluación [años]
    - t_PT_instalacion: edad del tubo a la instalación [años]
    - t_PT_diseño: vida de diseño del tubo [años]
    """
    return epsilon_ct_total * ((t_fin - t_PT_instalacion) / t_PT_diseño)

# ==========================
# A.4.2.2.4 Equation A.4-5
# ==========================
def calcular_Ca_EV(Ca_total, t_fin, t_PT_instalacion, t_PT_diseño):
    """
    Ecuación A.4-5: Adelgazamiento por corrosión hasta el final del período.
    C_a^EV = C_a_total * (t_fin - t_PT_instalacion) / t_PT_diseño

    Parámetros:
    - Ca_total: corrosión total estimada para vida de diseño [mm]
    - t_fin: tiempo al final del período de evaluación [años]
    - t_PT_instalacion: edad del tubo a la instalación [años]
    - t_PT_diseño: vida de diseño del tubo [años]
    """
    return Ca_total * ((t_fin - t_PT_instalacion) / t_PT_diseño)

# ==========================
# A.4.2.2.5 Equation A.4-6
# ==========================
def calcular_w_final(w0, epsilon_ct_EV, Ca_EV):
    """
    Ecuación A.4-6: Espesor final de la pared del tubo de presión.
    w = (1 - ε_ct^EV) * w0 - C_a^EV

    Parámetros:
    - w0: espesor corregido al inicio de vida útil [mm]
    - epsilon_ct_EV: deformación axial de creep hasta fin del período [adim.]
    - Ca_EV: corrosión acumulada hasta fin del período [mm]
    """
    return (1 - epsilon_ct_EV) * w0 - Ca_EV

# ==========================
# A.4.2.3.3 Equation A.4-7
# ==========================
def calcular_rm0(r_inicial, epsilon_ct_SE, epsilon_cr_SE, w0, Ca_SE):
    """
    Ecuación A.4-7: Radio medio inicial de vida útil del tubo.
    rm0 = r_meas - 0.5 * (1 - ε_ct^SE)/(1 + ε_ct^SE) * (w0 - C_a^SE)

    Parámetros:
    - r_inicial: radio interno medido al inicio del período [mm]
    - epsilon_ct_SE: deformación por creep axial hasta inicio [adim.]
    - epsilon_cr_SE: deformación por creep diametral hasta inicio [adim.]
    - w0: espesor corregido inicial [mm]
    - Ca_SE: corrosión acumulada hasta inicio del período [mm]
    """
    numerador = r_inicial + 0.5 * (1 - epsilon_ct_SE) * w0 - Ca_SE
    denominador = 1 + epsilon_cr_SE
    return numerador / denominador

# ==========================
# A.4.2.3.3 Equation A.4-8
# ==========================
def calcular_epsilon_cr_SE(epsilon_cr_total, t_SE, t_PT_instalacion, t_PT_diseño):
    """
    Ecuación A.4-8: Deformación diametral por creep hasta el inicio del período.
    ε_cr^SE = ε_cr_total * (t_SE - t_PT_instalacion) / t_PT_diseño

    Parámetros:
    - epsilon_cr_total: deformación diametral total esperada [adim.]
    - t_SE: tiempo al inicio del período [años]
    - t_PT_instalacion: tiempo de instalación del tubo [años]
    - t_PT_diseño: vida de diseño [años]
    """
    return epsilon_cr_total * ((t_SE - t_PT_instalacion) / t_PT_diseño)

# ==========================
# A.4.2.3.4 Equation A.4-9
# ==========================
def calcular_epsilon_cr_EV(epsilon_cr_total, t_fin, t_PT_instalacion, t_PT_diseño):
    """
    Ecuación A.4-9: Deformación diametral por creep hasta el final del período.
    ε_cr^EV = ε_cr_total * (t_fin - t_PT_instalacion) / t_PT_diseño

    Parámetros:
    - epsilon_cr_total: deformación diametral total [adim.]
    - t_fin: tiempo al final del período [años]
    - t_PT_instalacion: instalación del tubo [años]
    - t_PT_diseño: vida de diseño [años]
    """
    return epsilon_cr_total * ((t_fin - t_PT_instalacion) / t_PT_diseño)

# ==========================
# A.4.2.3.6 Equation A.4-10
# ==========================
def calcular_r_final(rm0, epsilon_cr_EV, epsilon_ct_EV, w0, Ca_EV):
    """
    Ecuación A.4-10 (CSA N285.8-15): radio interno final del tubo al final del período.

    r_i = (1 + ε_cr^EV) * r_m0 - 0.5 * (1 - ε_ct^EV) * w0 + C_a^EV

    Parámetros
    ----------
    rm0 : float
        Radio medio corregido al inicio [mm].
    epsilon_cr_EV : float
        Deformación diametral final (creep/creep+growth, según aplique) [adim.].
    epsilon_ct_EV : float
        Deformación axial final [adim.].
    w0 : float
        Espesor corregido inicial [mm].
    Ca_EV : float
        Corrosión/desgaste acumulado al final [mm] .

    Retorna
    -------
    float
        Radio interno final r_i [mm].
    """
    return (1 + epsilon_cr_EV) * rm0 - 0.5 * (1 - epsilon_ct_EV) * w0 + Ca_EV


def calcular_A_4_2_2_2(w_m, Ca_total, epsilon_ct_total, t_SE, t_PT_instalacion, t_PT_diseño):
    """
    A.4.2.2.2: Cálculo del espesor de pared corregido al inicio de vida útil del tubo.
    Retorna:
    - w0: espesor corregido
    - Ca_SE: corrosión acumulada al inicio
    - epsilon_ct_SE: creep acumulado al inicio
    """
    epsilon_ct_SE = calcular_epsilon_ct_SE(epsilon_ct_total, t_SE, t_PT_instalacion, t_PT_diseño)
    Ca_SE = calcular_Ca_SE(Ca_total, t_SE, t_PT_instalacion, t_PT_diseño)
    w0 = calcular_w0(w_m, Ca_SE, epsilon_ct_SE)
    return w0, Ca_SE, epsilon_ct_SE


def calcular_A_4_2_2_3(epsilon_ct_total, Ca_total, t_fin, t_PT_instalacion, t_PT_diseño):
    """
    A.4.2.2.3: Cálculo del creep y corrosión hasta el final del período.
    Retorna:
    - epsilon_ct_EV: creep acumulado al final
    - Ca_EV: corrosión acumulada al final
    """
    epsilon_ct_EV = calcular_epsilon_ct_EV(epsilon_ct_total, t_fin, t_PT_instalacion, t_PT_diseño)
    Ca_EV = calcular_Ca_EV(Ca_total, t_fin, t_PT_instalacion, t_PT_diseño)
    return epsilon_ct_EV, Ca_EV


def calcular_A_4_2_2_5(w0, epsilon_ct_EV, Ca_EV):
    """
    A.4.2.2.5: Cálculo del espesor final del tubo.
    Retorna:
    - w: espesor final
    """
    return calcular_w_final(w0, epsilon_ct_EV, Ca_EV)


def calcular_A_4_2_3_3(r_meas, epsilon_ct_SE, epsilon_cr_SE, w0, Ca_SE):
    """
    A.4.2.3.3: Cálculo del radio medio inicial del tubo.
    Retorna:
    - rm0: radio medio inicial
    """
    return calcular_rm0(r_meas, epsilon_ct_SE, epsilon_cr_SE, w0, Ca_SE)


def calcular_A_4_2_3_4(epsilon_cr_total, t_SE, t_fin, t_PT_instalacion, t_PT_diseño):
    """
    A.4.2.3.4: Cálculo del creep diametral acumulado.
    Retorna:
    - epsilon_cr_SE: hasta el inicio
    - epsilon_cr_EV: hasta el final
    """
    epsilon_cr_SE = calcular_epsilon_cr_SE(epsilon_cr_total, t_SE, t_PT_instalacion, t_PT_diseño)
    epsilon_cr_EV = calcular_epsilon_cr_EV(epsilon_cr_total, t_fin, t_PT_instalacion, t_PT_diseño)
    return epsilon_cr_SE, epsilon_cr_EV


def calcular_A_4_2_3_6(rm0, epsilon_cr_EV, epsilon_ct_EV, w0, Ca_EV):
    """
    A.4.2.3.6: Cálculo del radio interno final del tubo.
    Retorna:
    - ri: radio interno final
    """
    return calcular_r_final(rm0, epsilon_cr_EV, epsilon_ct_EV, w0, Ca_EV)

def evaluar_A_4_2(
    w_m, Ca_total, epsilon_ct_total, epsilon_cr_total,
    r_meas,
    t_SE, t_fin, t_PT_instalacion, t_PT_diseño
):
    """
    Ejecuta todo el flujo de cálculo de la sección A.4.2 del CSA N285.8-15.

    Parámetros:
    - w_m: espesor medido al inicio [mm]
    - Ca_total: corrosión total esperada en vida útil [mm]
    - epsilon_ct_total: deformación de creep axial total [adim.]
    - epsilon_cr_total: deformación de creep diametral total [adim.]
    - r_meas: radio interno medido al inicio [mm]
    - t_SE: tiempo al inicio del período de evaluación [años]
    - t_fin: tiempo al final del período de evaluación [años]
    - t_PT_instalacion: tiempo de instalación del tubo [años]
    - t_PT_diseño: vida de diseño del tubo [años]

    Retorna:
    - Diccionario con:
        - "w0": Espesor de pared corregido al inicio del período [mm]
        - "w_final": Espesor de pared proyectado al final del período [mm]
        - "Ca_SE": Corrosión acumulada al inicio del período [mm]
        - "Ca_EV": Corrosión acumulada al final del período [mm]
        - "epsilon_ct_SE": Deformación axial por creep acumulada al inicio [adim.]
        - "epsilon_ct_EV": Deformación axial por creep acumulada al final [adim.]
        - "rm0": Radio medio corregido al inicio del período [mm]
        - "epsilon_cr_SE": Deformación diametral por creep acumulada al inicio [adim.]
        - "epsilon_cr_EV": Deformación diametral por creep acumulada al final [adim.]
        - "r_final": Radio interno corregido al final del período [mm]
    """

    # A.4.2.2.2 - Espesor corregido inicial
    w0, Ca_SE, epsilon_ct_SE = calcular_A_4_2_2_2(w_m, Ca_total, epsilon_ct_total, t_SE, t_PT_instalacion, t_PT_diseño)

    # A.4.2.2.3 - Creep y corrosión al final
    epsilon_ct_EV, Ca_EV = calcular_A_4_2_2_3(epsilon_ct_total, Ca_total, t_fin, t_PT_instalacion, t_PT_diseño)

    # A.4.2.2.5 - Espesor final
    w_final = calcular_A_4_2_2_5(w0, epsilon_ct_EV, Ca_EV)

    # A.4.2.3.4 - Creep diametral
    epsilon_cr_SE, epsilon_cr_EV = calcular_A_4_2_3_4(epsilon_cr_total, t_SE, t_fin, t_PT_instalacion, t_PT_diseño)

    # A.4.2.3.3 - Radio medio inicial
    rm0 = calcular_A_4_2_3_3(r_meas, epsilon_ct_SE, epsilon_cr_SE, w0, Ca_SE)

    # A.4.2.3.6 - Radio interno final
    r_final = calcular_A_4_2_3_6(rm0, epsilon_cr_EV, epsilon_ct_EV, w0, Ca_EV)

    return {
        "w0": w0,                           # [mm] espesor corregido al inicio
        "w_final": w_final,                 # [mm] espesor final al fin del período
        "Ca_SE": Ca_SE,                     # [mm] corrosión acumulada al inicio
        "Ca_EV": Ca_EV,                     # [mm] corrosión acumulada al final
        "epsilon_ct_SE": epsilon_ct_SE,     # [-] creep axial acumulado al inicio
        "epsilon_ct_EV": epsilon_ct_EV,     # [-] creep axial acumulado al final
        "rm0": rm0,                         # [mm] radio medio corregido al inicio
        "epsilon_cr_SE": epsilon_cr_SE,     # [-] creep diametral al inicio
        "epsilon_cr_EV": epsilon_cr_EV,     # [-] creep diametral al final
        "r_final": r_final                  # [mm] radio interno final proyectado
    }
