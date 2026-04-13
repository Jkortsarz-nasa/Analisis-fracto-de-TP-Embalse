# ==========================
# =================================

# Sección A.2.2.1 : Caracterización geométrica de los defectos planares

# ==========================
# =================================
# Librerias

import pandas as pd



def clasificar_origen_falla_v2(tipo_fisico, orientacion=None):
    """
    Clasifica el tipo de falla según la Tabla 1 de CSA N285.8-15 y A.2.1
    Retorna: tipo de evaluación normativa correspondiente
    """
    tipo_fisico = tipo_fisico.lower()
    if tipo_fisico in [
        "bearing pad fretting", "debris fretting", "scores", "scratches",
        "crevice corrosion", "wall thinning", "mechanical damage"
    ]:
        return "Volumetric (seguir A.2.3)"
    elif tipo_fisico in ["crack", "lamination", "stringer", "manufacturing flaw"]:
        if orientacion is not None and abs(orientacion) < 10:
            return "Laminar (seguir A.2.2.3 y IWA-3360)"
        else:
            return "Planar (seguir A.2.2.2)"
    else:
        return "Tipo no definido (evaluación de ingeniería requerida)"


def caracterizar_falla_superficial(a, long_axial, long_circ):
    """
    Caracteriza una falla planar superficial según A.2.2.2.2
    """
    return {
        "Tipo": "Falla planar superficial",
        "a (profundidad)": a,
        "2c (longitud axial)": long_axial,
        "2b (longitud circunferencial)": long_circ
    }


def evaluar_multiples_fallas(a1, a2, Sz, Sθ, c1, c2, b1, b2):
    """
    Evalúa si múltiples fallas superficiales deben combinarse (A.2.2.2.3.1)
    """
    if Sz <= 0.5 * max(a1, a2) and Sθ <= 0.5 * max(a1, a2):
        return {
            "Combinación requerida": True,
            "2c combinado": 2 * c1 + Sz + 2 * c2,
            "2b combinado": 2 * b1 + Sθ + 2 * b2,
            "a combinado": max(a1, a2)
        }
    else:
        return {
            "Combinación requerida": False,
            "Nota": "Las fallas se evalúan individualmente"
        }


def reclasificar_falla_por_proximidad(a, S):
    """
    Reclasifica una falla subsuperficial como superficial si S < 0.4a (A.2.2.3 + IWA-3320)
    """
    if S < 0.4 * a:
        return "Falla subsuperficial reclasificada como superficial (S < 0.4a)"
    else:
        return "Falla subsuperficial se mantiene como subsuperficial"


def evaluar_fallas_laminares_combinadas(W1, l1, W2, l2, S1, S2, H, a1, a2):
    """
    Evalúa si dos fallas laminares deben combinarse según IWA-3360(b),
    y devuelve geometría combinada si corresponde (IWA-3360(c)).
    """
    criterio = 0.37 * min(W1, l1, W2, l2)
    criterio_h = 0.17 * min(W1, l1, W2, l2)

    cumple_S1 = S1 <= criterio
    cumple_S2 = S2 <= criterio
    cumple_H = H <= criterio_h

    combinar = cumple_S1 and cumple_S2 and cumple_H

    resultado = {
        "Combinar fallas laminares": combinar,
        "S1 cumple": cumple_S1,
        "S2 cumple": cumple_S2,
        "H cumple": cumple_H,
        "Umbral S1/S2": round(criterio, 3),
        "Umbral H": round(criterio_h, 3)
    }

    if combinar:
        W_comb = W1 + S1 + W2
        l_comb = l1 + S2 + l2
        a_comb = max(a1, a2)
        area_combinada = 0.75 * (W_comb * l_comb)

        resultado.update({
            "Ancho combinado (W)": W_comb,
            "Largo combinado (l)": l_comb,
            "Profundidad combinada (a)": a_comb,
            "Área proyectada equivalente": round(area_combinada, 3)
        })

    return resultado

# Mostrar resumen de todas las funciones definidas
import pandas as pd
from IPython.display import display

funciones = pd.DataFrame({
    "Función": [
        "clasificar_origen_falla_v2",
        "caracterizar_falla_superficial",
        "evaluar_multiples_fallas",
        "reclasificar_falla_por_proximidad",
        "evaluar_fallas_laminares_combinadas"
    ],
    "Descripción": [
        "Clasifica el tipo de falla según su naturaleza y orientación (A.2.1, Tabla 1)",
        "Genera parámetros geométricos para falla superficial (A.2.2.2.2)",
        "Determina si múltiples fallas superficiales deben combinarse (A.2.2.2.3.1)",
        "Reclasifica una falla subsuperficial si está cerca de la superficie (A.2.2.3 + IWA-3320)",
        "Evalúa combinación de fallas laminares y devuelve geometría combinada (IWA-3360(b)(c))"
    ]
})

display(funciones)



def evaluar_interaccion_fallas_por_canal_v2(canal, lista_fallas):
    """
    Versión mejorada: evalúa combinación entre fallas planas y laminares, y siempre reporta
    dimensiones combinadas si corresponde.
    """
    combinaciones = []
    n = len(lista_fallas)

    for i in range(n):
        for j in range(i+1, n):
            falla_i = lista_fallas[i]
            falla_j = lista_fallas[j]
            tipo_i = falla_i["tipo"]
            tipo_j = falla_j["tipo"]

            if tipo_i == "Planar" and tipo_j == "Planar":
                Sz = falla_i["Sz"].get(j, None)
                Sθ = falla_i["Sθ"].get(j, None)
                if Sz is not None and Sθ is not None:
                    resultado = evaluar_multiples_fallas(
                        a1=falla_i["a"], a2=falla_j["a"],
                        Sz=Sz, Sθ=Sθ,
                        c1=falla_i["c"], c2=falla_j["c"],
                        b1=falla_i["b"], b2=falla_j["b"]
                    )
                    resultado.update({
                        "Canal": canal,
                        "Falla A": falla_i["id"],
                        "Falla B": falla_j["id"],
                        "Tipo": "Planar"
                    })
                    # renombrar dimensiones combinadas para uniformidad
                    if resultado["Combinación requerida"]:
                        resultado["Ancho combinado (2b)"] = resultado.pop("2b combinado")
                        resultado["Largo combinado (2c)"] = resultado.pop("2c combinado")
                        resultado["Profundidad combinada (a)"] = resultado.pop("a combinado")
                    combinaciones.append(resultado)

            elif tipo_i == "Laminar" and tipo_j == "Laminar":
                S1 = falla_i["Sz"].get(j, None)
                S2 = falla_i["Sθ"].get(j, None)
                H  = falla_i["H"].get(j, None)
                if S1 is not None and S2 is not None and H is not None:
                    resultado = evaluar_fallas_laminares_combinadas(
                        W1=falla_i["W"], l1=falla_i["l"],
                        W2=falla_j["W"], l2=falla_j["l"],
                        S1=S1, S2=S2, H=H,
                        a1=falla_i["a"], a2=falla_j["a"]
                    )
                    resultado.update({
                        "Canal": canal,
                        "Falla A": falla_i["id"],
                        "Falla B": falla_j["id"],
                        "Tipo": "Laminar"
                    })
                    # renombrar para compatibilidad
                    if resultado["Combinar fallas laminares"]:
                        resultado["Ancho combinado (2b)"] = resultado.pop("Ancho combinado (W)")
                        resultado["Largo combinado (2c)"] = resultado.pop("Largo combinado (l)")
                        resultado["Profundidad combinada (a)"] = resultado.pop("Profundidad combinada (a)")
                    combinaciones.append(resultado)

            elif "subsuperficial" in tipo_i.lower() or "subsuperficial" in tipo_j.lower():
                combinaciones.append({
                    "Canal": canal,
                    "Falla A": falla_i["id"],
                    "Falla B": falla_j["id"],
                    "Tipo": "Subsuperficial",
                    "Combinación requerida": False,
                    "Nota": "No se requiere combinación normativa para fallas subsuperficiales"
                })

    return combinaciones