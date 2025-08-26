# solute.py
import numpy as np
from material import fk
from watflow import solve_tridiagonal


def solve_solute_transport(state, config, materials):
    """
    Transporte 1D de soluto (adveção-dispersão) em meio poroso não saturado.
    Esquema: difusão/dispersão implícita + advecção upwind semi-implícita.

    Requer no estado:
      - state['h']      : pressão matricial (cm)
      - state['theta']  : teor de água (cm3/cm3)
      - state['c']      : concentração (inicializa se não existir)

    Parâmetros em config (com defaults):
      - NumNP (int)           : número de nós
      - dz (float)            : espaçamento (cm)
      - dt (float)            : passo de tempo (dia)
      - alphaL (float)        : dispersividade longitudinal (cm)
      - Dm (float)            : difusão molecular (cm^2/dia)
      - c_top (float|None)    : Dirichlet no topo (se None, aplica condição “natural” simples)
      - bc_c_bot (str)        : 'neumann' (padrão) ou 'dirichlet'
      - c_bot (float)         : valor de Dirichlet na base, se aplicável
    """

    NumNP = int(config.get("NumNP", 100))
    dz = float(config.get("dz", 1.0))
    dt = float(config.get("dt", 0.1))

    alphaL = float(config.get("alphaL", 1.0))  # cm
    Dm = float(config.get("Dm", 1e-3))         # cm^2/dia
    c_top = config.get("c_top", None)
    bc_c_bot = str(config.get("bc_c_bot", "neumann")).lower()
    c_bot = float(config.get("c_bot", 0.0))

    # Estado hídrico
    h = state["h"]
    theta = state["theta"]

    # Estado de concentração (inicializa se necessário)
    if "c" not in state or state["c"] is None:
        state["c"] = np.zeros(NumNP, dtype=float)
    c_old = state["c"].astype(float).copy()

    # Parâmetros de material (assume único material por enquanto)
    Par = materials[0]["Par"]

    # Condutividade hidráulica nos nós
    K = np.array([fk(0, hi, Par) for hi in h], dtype=float)

    # Fluxos de Darcy nas interfaces (i+1/2): q = -K * dh/dz
    q_int = np.zeros(NumNP - 1, dtype=float)
    for i in range(NumNP - 1):
        Ki = 0.5 * (K[i] + K[i + 1])
        q_int[i] = -Ki * (h[i + 1] - h[i]) / dz  # cm/dia (positivo para cima)

    # Dispersão efetiva nas interfaces: D = alphaL * |q| / theta + Dm
    theta_int = 0.5 * (theta[:-1] + theta[1:])
    D_int = alphaL * np.abs(q_int) / np.maximum(theta_int, 1e-12) + Dm

    # Montagem do sistema tridiagonal: A * c_new = rhs
    a = np.zeros(NumNP, dtype=float)  # subdiagonal
    b = np.zeros(NumNP, dtype=float)  # diagonal
    c_sup = np.zeros(NumNP, dtype=float)  # superdiagonal
    rhs = np.zeros(NumNP, dtype=float)

    # Difusão/Dispersão implícita
    for i in range(1, NumNP - 1):
        Dp = D_int[i]       # entre i e i+1
        Dm_ = D_int[i - 1]  # entre i-1 e i

        a[i] += Dm_ / (dz * dz)
        c_sup[i] += Dp / (dz * dz)
        b[i] += (theta[i] / dt) + a[i] + c_sup[i]

    # Advecção upwind semi-implícita (- d/dz (q c))
    for i in range(1, NumNP - 1):
        qR = q_int[i]     # interface direita (i -> i+1)
        qL = q_int[i - 1] # interface esquerda (i-1 -> i)

        # Contribuição da interface direita
        if qR >= 0:      # fluxo para cima (usa c_i)
            b[i] += qR / dz
        else:            # fluxo para baixo (usa c_{i+1})
            c_sup[i] += -qR / dz

        # Contribuição da interface esquerda
        if qL >= 0:      # fluxo para cima (usa c_{i-1})
            a[i] += qL / dz
        else:            # fluxo para baixo (usa c_i)
            b[i] += -qL / dz

        # Termo temporal
        rhs[i] = theta[i] * c_old[i] / dt

    # Condição de contorno no topo
    if c_top is not None:  # Dirichlet
        b[0] = 1.0
        rhs[0] = float(c_top)
    else:
        # Aproximação simples: manter valor anterior no topo
        b[0] = 1.0
        rhs[0] = c_old[0]

    # Condição de contorno na base
    if bc_c_bot == "dirichlet":
        b[-1] = 1.0
        rhs[-1] = float(c_bot)
    else:
        # Neumann (gradiente zero): mantém valor anterior
        b[-1] = 1.0
        rhs[-1] = c_old[-1]

    # Resolve sistema tridiagonal
    c_new = solve_tridiagonal(a, b, c_sup, rhs)

    # Atualiza estado e retorna
    state["c"] = c_new
    return {
        "c": c_new.copy(),
        "q_int": q_int.copy(),
        "D_int": D_int.copy(),
    }
