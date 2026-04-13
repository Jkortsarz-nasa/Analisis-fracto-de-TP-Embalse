# ==========================
# =================================

# La presente librería calcula las propiedades mecanicas del material en el sentido axial y circunferencial según las secciones D.3.4.2 y D.3.4.3 de la norma CSA 285

# ==========================
# =================================



def sigma_y_axial_lb(T_C: float) -> float:
    """
    Calcula la tensión de fluencia axial límite inferior σ_y^H [MPa]
    según Eq. D.3-3 (CSA N285.8-15).

    Parámetros
    ----------
    T_C : float
        Temperatura en grados Celsius [°C].

    Retorna
    -------
    float
        σ_y^H axial límite inferior en MPa.
    """
    return 758.0 - 0.595 * T_C


def sigma_y_transverse_lb(T_C: float) -> float:
    """
    Calcula la tensión de fluencia transversal límite inferior σ_y^H [MPa]
    según Eq. D.3-5 (CSA N285.8-15).

    Parámetros
    ----------
    T_C : float
        Temperatura en grados Celsius [°C].

    Retorna
    -------
    float
        σ_y^H transversal límite inferior en MPa.
    """
    return 988.0 - 1.154 * T_C

def sigma_u_axial_lb(T_C: float) -> float:
    """
    Calcula la tensión última axial límite inferior σ_u^H [MPa]
    según Eq. D.3-4 (CSA N285.8-15).

    Parámetros
    ----------
    T_C : float
        Temperatura en grados Celsius [°C].

    Retorna
    -------
    float
        σ_u^H axial límite inferior en MPa.
    """
    return 907.0 - 0.722 * T_C


def sigma_u_transverse_lb(T_C: float) -> float:
    """
    Calcula la tensión última transversal límite inferior σ_u^H [MPa]
    según Eq. D.3-6 (CSA N285.8-15).

    Parámetros
    ----------
    T_C : float
        Temperatura en grados Celsius [°C].

    Retorna
    -------
    float
        σ_u^H transversal límite inferior en MPa.
    """
    return 1021.0 - 1.245 * T_C






