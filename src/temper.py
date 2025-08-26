# temper.py
import numpy as np
from watflow import solve_tridiagonal
from material import fk  # para obter K(h) e calcular q nas interfaces


def solve_heat_transport(state, config, materials):
    """
    Transporte 1D de calor em meio poroso não saturado (condução + advecção pela água).
    Esquema: condução implícita + termo advectivo upwind semi-implícito.

    Requer no estado:
      - state['h']      : pressão matricial (cm)
      - state['theta']  : teor de água (cm3/cm3)
      - state['T']      : temperatura (°C) (inicializa se não existir)

    Parâmetros em config (com defaults):
      — malha/tempo:
        - NumNP (int)         : número de nós
        - dz (float)          : espaçamento (cm)
        - dt (float)          : passo de tempo (dia)
      — propriedades térmicas:
        - lambda_dry (W/m/K)  : condutividade térmica do solo seco
        - lambda_sat (W/m/K)  : condutividade térmica do solo saturado
        - rho_w (kg/m3)       : densidade da água
        - c_w (J/kg/K)        : calor específico da água
        - rho_s (kg/m3)       : densidade dos sólidos (solo)
        - c_s (J/kg/K)        : calor específico dos sólidos
        - alphaT (m)          : dispersividade térmica (difusividade mecânica)
      — condições de contorno:
        - T_top (float|None)  : Dirichlet no topo (°C). Se None, mantém valor anterior.
        - bc_T_bot (str)      : 'neumann' (gradiente zero) ou 'dirichlet'
        - T_bot (float)       : valor para Dirichlet na base, se aplicável

    Observações de unidade:
      - O código usa entrada hidráulica em cm e dia; propriedades térmicas em SI.
      - Convertemos para SI (m, s) internamente para a parte térmica.
    """

    # -------------------------
    # 1) Leitura de parâmetros
    # -------------------------
    NumNP = int(config.get("NumNP", 100))
    dz_cm = float(config.get("dz", 1.0))           # cm
    dt_day = float(config.get("dt", 0.1))          # dia

    # Conversões p/ SI
    dz = dz_cm / 100.0                              # m
    dt = dt_day * 24.0 * 3600.0                     # s

    # Propriedades térmicas (SI)
    lam_dry = float(config.get("lambda_dry", 0.25))  # W/m/K
    lam_sat = float(config.get("lambda_sat", 2.5))   # W/m/K
    rho_w   = float(config.get("rho_w", 1000.0))     # kg/m3
    c_w     = float(config.get("c_w", 4180.0))       # J/kg/K
    rho_s   = float(config.get("rho_s", 2650.0))     # kg/m3
    c_s     = float(config.get("c_s", 800.0))        # J/kg/K
    alphaT  = float(config.get("alphaT", 0.01))      # m (dispersividade térmica)

    # Condições de contorno
    T_top    = config.get("T_top", None)                     # °C
    bc_T_bot = str(config.get("bc_T_bot", "neumann")).lower()
    T_bot    = float(config.get("T_bot", 20.0))              # °C

    # Estado
    h = state["h"]
    theta = state["theta"]

    if "T" not in state or state["T"] is None:
        state["T"] = np.full(NumNP, 20.0, dtype=float)
    T_old = state["T"].astype(float).copy()

    # Parâmetros de material (assume único)
    Par = materials[0]["Par"]
    # Qr, Qs
    Qr, Qs = float(Par[0]), float(Par[1])

    # ---------------------------------------
    # 2) Fluxos e propriedades nas interfaces
    # ---------------------------------------
    # Condutividade hidráulica e fluxos (em cm/d); converter para m/s
    K = np.array([fk(0, hi, Par) for hi in h], dtype=float)
    q_int_cm_day = np.zeros(NumNP - 1, dtype=float)
    for i in range(NumNP - 1):
        Ki = 0.5 * (K[i] + K[i + 1])
        q_int_cm_day[i] = -Ki * (h[i + 1] - h[i]) / dz_cm

    q_int = q_int_cm_day * (1.0 / 100.0) / (24.0 * 3600.0)  # m/s

    # Condutividade térmica λ(θ) (método simples: mistura linear pela saturação efetiva)
    Se = np.clip((theta - Qr) / max(Qs - Qr, 1e-12), 0.0, 1.0)  # adimensional
    lam_node = lam_dry + (lam_sat - lam_dry) * Se               # W/m/K
    lam_int = 0.5 * (lam_node[:-1] + lam_node[1:])              # nas interfaces

    # Dispersão térmica (advecção de calor -> termo efetivo adicional de difusão)
    # λ_eff = λ + ρ_w c_w α_T |q|
    lam_eff_int = lam_int + rho_w * c_w * alphaT * np.abs(q_int)  # W/m/K

    # Capacidade térmica volumétrica (J/m3/K): C_th = ρ_s c_s (1 - porosidade) + ρ_w c_w θ
    # Assumindo porosidade ~ Qs (volume de poros); sólido ocupa (1 - Qs)
    C_th = rho_s * c_s * (1.0 - Qs) + rho_w * c_w * theta  # vetor nos nós

    # --------------------------------------------------
    # 3) Montagem do sistema tridiagonal para T^{n+1}
    # --------------------------------------------------
    a = np.zeros(NumNP, dtype=float)    # subdiagonal
    b = np.zeros(NumNP, dtype=float)    # diagonal
    c_sup = np.zeros(NumNP, dtype=float)  # superdiagonal
    rhs = np.zeros(NumNP, dtype=float)

    # Condução (implícita) + termo temporal
    for i in range(1, NumNP - 1):
        kR = lam_eff_int[i]       # entre i e i+1
        kL = lam_eff_int[i - 1]   # entre i-1 e i

        a[i] += kL / (dz * dz)
        c_sup[i] += kR / (dz * dz)
        b[i] += C_th[i] / dt + a[i] + c_sup[i]
        rhs[i] = C_th[i] * T_old[i] / dt

    # Advecção de calor: termo -ρ_w c_w * d/dz (q T)
    # upwind semi-implícito pelas interfaces (similar ao soluto)
    for i in range(1, NumNP - 1):
        qR = q_int[i]     # i -> i+1 (m/s)
        qL = q_int[i - 1] # i-1 -> i

        adv_coef = rho_w * c_w / dz  # J/m3/K/s

        if qR >= 0:  # fluxo para cima: usa T_i
            b[i] += adv_coef * qR
        else:        # fluxo para baixo: usa T_{i+1}
            c_sup[i] += -adv_coef * qR

        if qL >= 0:  # (i-1 -> i): usa T_{i-1}
            a[i] += adv_coef * qL
        else:        # (i -> i-1): usa T_i
            b[i] += -adv_coef * qL

    # Condições de contorno
    # Topo
    if T_top is not None:  # Dirichlet
        b[0] = 1.0
        rhs[0] = float(T_top)
    else:
        # Mantém valor anterior (aproximação suave)
        b[0] = 1.0
        rhs[0] = T_old[0]

    # Base
    if bc_T_bot == "dirichlet":
        b[-1] = 1.0
        rhs[-1] = float(T_bot)
    else:
        # Neumann: gradiente zero -> mantém valor anterior (alternativa simples)
        b[-1] = 1.0
        rhs[-1] = T_old[-1]

    # ---------------------------
    # 4) Resolver e atualizar T
    # ---------------------------
    T_new = solve_tridiagonal(a, b, c_sup, rhs)
    state["T"] = T_new

    return {
        "T": T_new.copy(),
        "lambda_int": lam_int.copy(),
        "q_int": q_int.copy(),
    }
