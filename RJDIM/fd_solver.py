from matplotlib.pyplot import axis
from scipy.sparse import diags
import numpy as np


def K(alpha, n):
    # Elementos de las diagonales
    C = alpha*np.ones(n-1)
    C[-1] = 0

    D = (1-2*alpha)*np.ones(n)
    D[0] = 1
    D[-1] = 1

    E = alpha*np.ones(n-1)
    E[0] = 0

    k = [C,D,E]
    offset = [-1, 0, 1]
    A = diags(k,offset).toarray()

    return(A)