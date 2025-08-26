# solute.py
import numpy as np
from material import fk
from watflow import build_sink_profile
from utils import solve_tridiagonal

def solve_solute_transport(state, config, materials):
    """
    Transporte 1D de soluto (adveção-dispersão) com remoção por extração radicular:
      ∂(θc)/∂t = -∂(q c)/∂z + ∂/∂z(D ∂c/∂z) - S c
    onde S (1/dia) vem de build_sink_profile() ou state['sink'].
    """
    NumNP = int(config.get('NumNP', 100))
    dz    = float(config.get('dz', 1.0))
    dt    = float(config.get('dt', 0.1))

    alphaL = float(config.get('alphaL', 1.0))  # cm
    Dm     = float(config.get('Dm', 1e-3))     # cm^2/dia
    c_top  = config.get('c_top', None)
    bc_c_bot = str(config.get('bc_c_bot', 'neumann')).lower()
    c_bot  = float(config.get('c_bot', 0.0))

    # Hídrico
    h     = state['h']
    theta = state['theta']

    # Concentração
    if 'c' not in state or state['c'] is None:
        state['c'] = np.zeros(NumNP, dtype=float)
    c_old = state['c'].astype(float).copy()

    # Material
    Par = materials[0]['Par']

    # Fluxos q nas interfaces
    K = np.array([fk(0, hi, Par) for hi in h], dtype=float)
    q_int = np.zeros(NumNP-1)
    for i in range(NumNP-1):
        Ki = 0.5 * (K[i] + K[i+1])
        q_int[i] = -Ki * (h[i+1] - h[i]) / dz  # cm/d

    # Dispersão efetiva D_int
    theta_int = 0.5 * (theta[:-1] + theta[1:])
    D_int = alphaL * np.abs(q_int) / np.maximum(theta_int, 1e-12) + Dm

    # SINK (1/dia)
    S = build_sink_profile(config, state)

    # Sistema tridiagonal
    a = np.zeros(NumNP)
    b = np.zeros(NumNP)
    c_sup = np.zeros(NumNP)
    rhs = np.zeros(NumNP)

    # Difusão/Dispersão implícita + termo temporal
    for i in range(1, NumNP-1):
        Dp  = D_int[i]
        Dm_ = D_int[i-1]
        a[i]   += Dm_ / (dz*dz)
        c_sup[i] += Dp  / (dz*dz)
        b[i]   += (theta[i] / dt) + a[i] + c_sup[i]
        rhs[i]  = theta[i] * c_old[i] / dt

        # Termo de remoção por extração: + S[i] * c_new  -> soma S[i] na diagonal
        b[i] += S[i]

    # Advecção upwind semi-implícita
    for i in range(1, NumNP-1):
        qR = q_int[i]
        qL = q_int[i-1]
        if qR >= 0:
            b[i] +=  qR / dz
        else:
            c_sup[i] += -qR / dz
        if qL >= 0:
            a[i] +=  qL / dz
        else:
            b[i] += -qL / dz

    # CC topo
    if c_top is not None:
        b[0] = 1.0
        rhs[0] = float(c_top)
    else:
        b[0] = 1.0
        rhs[0] = c_old[0]

    # CC base
    if bc_c_bot == 'dirichlet':
        b[-1] = 1.0
        rhs[-1] = float(c_bot)
    else:
        b[-1] = 1.0
        rhs[-1] = c_old[-1]

    # Resolve
    c_new = solve_tridiagonal(a, b, c_sup, rhs)

    state['c'] = c_new
    # opcional manter para pós-processo
    state['sink'] = S
    return {'c': c_new.copy(), 'q_int': q_int.copy(), 'D_int': D_int.copy(), 'sink': S.copy()}
