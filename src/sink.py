import numpy as np
from utils import feddes_alpha

def build_sink_profile(config, state):
    """
    Constrói S[i] em 1/dia (termo volumétrico de extração) distribuindo Tpot (cm/dia)
    ao longo do perfil de raízes (beta_root), com estresse α(h).
    - Se já existir state['sink'], usa-o diretamente (array de tamanho NumNP, 1/dia).
    - Caso contrário, usa config['Tpot'] (cm/d), config['beta_root'] (array que soma 1),
      e depth/shape implícitos. Se beta_root não existir, assume uniforme até z_root.
    """
    NumNP = int(config.get("NumNP", len(state.get("h", [])) or 100))
    dz = float(config.get("dz", 1.0))  # cm
    h = state["h"]

    # 1) Se veio pronto no estado, use
    if isinstance(state.get("sink", None), np.ndarray):
        S = state["sink"].astype(float)
        return S

    # 2) Monte a partir de config
    Tpot = float(config.get("Tpot", 0.0))  # cm/dia (positivo para extração)
    if abs(Tpot) < 1e-12:
        return np.zeros(NumNP, dtype=float)

    # perfil de raízes
    beta = np.array(config.get("beta_root", []), dtype=float)
    if beta.size != NumNP:
        # Cria beta uniforme até z_root
        z_root = float(config.get("z_root", NumNP * dz))  # cm de profundidade radicular
        n_root = max(1, min(NumNP, int(round(z_root / dz))))
        beta = np.zeros(NumNP, dtype=float)
        beta[:n_root] = 1.0 / n_root

    # estresse (Feddes) opcional
    fed = config.get("feddes_params", None)
    if fed and isinstance(fed, (list, tuple)) and len(fed) == 4:
        h1, h2, h3, h4 = map(float, fed)
        alpha = np.array([feddes_alpha(hi, h1, h2, h3, h4) for hi in h], dtype=float)
    else:
        alpha = np.ones(NumNP, dtype=float)

    # S tem unidades [1/dia]: (Tpot cm/d) * (beta adim) * (alpha adim) / (dz cm)
    S = (Tpot * beta * alpha) / max(dz, 1e-12)
    return S

