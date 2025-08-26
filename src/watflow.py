import numpy as np
from material import fk, theta_vg, c_vg
from sink import build_sink_profile
from utils import solve_tridiagonal

def solve_water_flow(state, config, materials):
    """
    Solver simplificado 1D para a equação de Richards (implícito, Newton).
    Inclui termo de extração S (1/dia): ∂θ/∂t = -∂q/∂z - S
    """
    NumNP = config.get("NumNP", 100)
    dz = config.get("dz", 1.0)
    dt = config.get("dt", 0.1)
    MaxIt = config.get("MaxIt", 20)
    TolH = config.get("TolH", 1e-6)

    Par = materials[0]["Par"]  # supõe material único
    h = state["h"].copy()
    theta_old = state["theta"].copy()
    theta = theta_vg(h, Par)

    # Condições de contorno simples
    h_top = config.get("h_top", -100.0)  # Dirichlet no topo (cm)
    q_bot = config.get("q_bot", 0.0)     # Neumann na base (cm/dia)

    # --- SINK (1/dia) ---
    S = build_sink_profile(config, state)  # tamanho NumNP

    for _ in range(MaxIt):
        # propriedades
        K = np.array([fk(0, hi, Par) for hi in h], dtype=float)
        C = c_vg(h, Par)

        # coeficientes do sistema tridiagonal
        a = np.zeros(NumNP)
        b = np.zeros(NumNP)
        c = np.zeros(NumNP)
        r = np.zeros(NumNP)

        # nós internos
        for i in range(1, NumNP - 1):
            Ki_plus = 0.5 * (K[i] + K[i + 1])
            Ki_minus = 0.5 * (K[i] + K[i - 1])

            # fluxos (cm/dia)
            flux_plus = -Ki_plus * (h[i + 1] - h[i]) / dz
            flux_minus = -Ki_minus * (h[i] - h[i - 1]) / dz

            # MASSA: ∂θ/∂t = -∂q/∂z - S
            # Resíduo desejando zero: (θ-θ_old)/dt + ∂q/∂z + S = 0
            div_q = (flux_plus - flux_minus) / dz
            r[i] = (theta[i] - theta_old[i]) / dt + div_q + S[i]

            # Coefs lineares básicos (aprox. um passo de Newton simples)
            a[i] =  Ki_minus / (dz * dz)
            c[i] =  Ki_plus  / (dz * dz)
            b[i] = C[i] / dt + a[i] + c[i]
            # Termo sink atua como +S[i] na equação: adiciona a diagonal
            b[i] += 0.0  # (aqui S entra apenas no resíduo; para fully implicit em h, precisaria dS/dh)

        # Topo: Dirichlet
        b[0] = 1.0
        r[0] = h[0] - h_top

        # Base: Neumann em q_bot
        Kb = K[-1]
        b[-1] = 1.0
        r[-1] = (h[-1] - h[-2]) / dz + q_bot / max(Kb, 1e-12)

        # Resolve
        dh = solve_tridiagonal(a, b, c, -r)

        h += dh
        theta = theta_vg(h, Par)

        if np.linalg.norm(dh, ord=np.inf) < TolH:
            break

    state["h"] = h
    state["theta"] = theta
    # Opcional: armazenar S
    state["sink"] = S
    return {"h": h.copy(), "theta": theta.copy(), "sink": S.copy()}


