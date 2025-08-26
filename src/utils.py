import numpy as np

def feddes_alpha(h, h1=-10.0, h2=-25.0, h3=-400.0, h4=-8000.0):
    """
    Resposta de estresse radicular (Feddes). h em cm (negativo em sucção).
    Retorna α∈[0,1]. Ajuste limiares via config se quiser.
    """
    # Convenção: valores típicos (ex.: cultura comum); ajuste conforme o caso
    if h >= h1:               # saturação/anaerobiose
        return 0.0
    if h1 > h >= h2:
        # ramp up de 0 -> 1
        return (h1 - h) / (h1 - h2 + 1e-12)
    if h2 > h >= h3:
        return 1.0
    if h3 > h >= h4:
        # ramp down de 1 -> 0
        return (h - h4) / (h3 - h4 + 1e-12)
    return 0.0

def solve_tridiagonal(a, b, c, d):
    """
    Resolve o sistema tridiagonal Ax = d, onde:
      - a: subdiagonal (len n)   [a[0] é ignorado]
      - b: diagonal   (len n)
      - c: superdiag  (len n)    [c[-1] é ignorado]
      - d: lado direito (len n)
    Retorna x (len n).
    Implementação: algoritmo de Thomas (O(n)).
    """
    import numpy as np

    a = np.asarray(a, dtype=float).copy()
    b = np.asarray(b, dtype=float).copy()
    c = np.asarray(c, dtype=float).copy()
    d = np.asarray(d, dtype=float).copy()

    n = b.size
    if any(arr.size != n for arr in (a, c, d)):
        raise ValueError("a, b, c, d devem ter o mesmo tamanho")

    # forward sweep
    for i in range(1, n):
        denom = b[i-1]
        if abs(denom) < 1e-18:
            raise ZeroDivisionError("Diagonal principal com zero durante Thomas")
        w = a[i] / denom
        b[i] = b[i] - w * c[i-1]
        d[i] = d[i] - w * d[i-1]

    # back substitution
    x = np.zeros(n, dtype=float)
    if abs(b[-1]) < 1e-18:
        raise ZeroDivisionError("Diagonal principal com zero no último passo do Thomas")
    x[-1] = d[-1] / b[-1]
    for i in range(n-2, -1, -1):
        if abs(b[i]) < 1e-18:
            raise ZeroDivisionError("Diagonal principal com zero na retro-substituição")
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]

    return x