# time_control.py
from typing import Callable, Dict, Any, Optional, List
import numpy as np


def run_time_loop(
    config: Dict[str, Any],
    profile: Dict[str, Any],
    materials: List[Dict[str, Any]],
    state: Dict[str, Any],
    solver_water: Callable[[Dict[str, Any], Dict[str, Any], List[Dict[str, Any]]], Dict[str, np.ndarray]],
    solver_solute: Optional[Callable[[Dict[str, Any], Dict[str, Any], List[Dict[str, Any]]], Dict[str, np.ndarray]]] = None,
    solver_heat: Optional[Callable[[Dict[str, Any], Dict[str, Any], List[Dict[str, Any]]], Dict[str, np.ndarray]]] = None,
) -> Dict[str, Any]:
    """
    Orquestra a simulação no tempo.

    Args:
        config: dicionário de configuração (dt, t_end, NumNP, etc.).
        profile: perfil/condições iniciais (não usado diretamente aqui; mantido por compatibilidade).
        materials: lista de materiais (parametrizações hidráulicas/termo-físicas).
        state: estado mutável do modelo; deve conter ao menos 'h' e 'theta'. Pode conter 'c' e 'T'.
        solver_water: função que atualiza o estado hídrico (obrigatória).
        solver_solute: função que atualiza o transporte de soluto (opcional).
        solver_heat: função que atualiza o transporte de calor (opcional).

    Returns:
        results: dict com séries temporais e coleções nodais (para posterior escrita).
    """

    # Tempo
    t = float(state.get("t", 0.0))
    dt = float(config.get("dt", 0.1))         # dia
    t_end = float(config.get("t_end", 1.0))   # dia
    save_every = int(config.get("save_every", 1))  # salvar a cada N passos
    NumNP = int(config.get("NumNP", len(state.get("h", [])) or 100))

    # Acumuladores de saída
    time_levels = []      # registros por passo (compacto)
    nodes_series = []     # h nodal
    conc_series = []      # c nodal (se houver soluto)
    temp_series = []      # T nodal (se houver calor)
    balance = []          # espaço para balanços se quiser acumular
    check_msgs = []

    step = 0
    while t < t_end - 1e-12:
        step += 1

        # --- 1) Água ---
        try:
            res_w = solver_water(state, config, materials)
        except Exception as e:
            check_msgs.append(f"[t={t:.6g} d] Falha no solver de água: {e}")
            break

        # --- 2) Soluto (opcional) ---
        res_c = None
        if solver_solute is not None:
            try:
                res_c = solver_solute(state, config, materials)
            except Exception as e:
                check_msgs.append(f"[t={t:.6g} d] Falha no solver de soluto: {e}")
                # Decide: parar ou continuar sem soluto. Aqui vamos parar:
                break

        # --- 3) Calor (opcional) ---
        res_T = None
        if solver_heat is not None:
            try:
                res_T = solver_heat(state, config, materials)
            except Exception as e:
                check_msgs.append(f"[t={t:.6g} d] Falha no solver de calor: {e}")
                break

        # Atualiza tempo
        t += dt
        state["t"] = t

        # --- 4) Coleta de resultados ---
        # Registro compacto do passo
        h_now = res_w["h"] if "h" in res_w else state["h"]
        time_levels.append({
            "step": step,
            "t": t,
            "h_top": float(h_now[0]),
            "h_mid": float(h_now[len(h_now)//2]),
            "h_bot": float(h_now[-1]),
        })

        # Salvar séries nodais a cada N passos (reduz I/O)
        if step % save_every == 0:
            nodes_series.append({
                "t": t,
                **{f"h_{i}": float(v) for i, v in enumerate(h_now)}
            })

            if res_c is not None and "c" in res_c:
                c_now = res_c["c"]
                conc_series.append({
                    "t": t,
                    **{f"c_{i}": float(v) for i, v in enumerate(c_now)}
                })

            if res_T is not None and "T" in res_T:
                T_now = res_T["T"]
                temp_series.append({
                    "t": t,
                    **{f"T_{i}": float(v) for i, v in enumerate(T_now)}
                })

        # (Opcional) espaço para balanços/fluxos agregados por passo
        # Ex.: balance.append({...})

    # Mensagem de checagem
    if not check_msgs:
        check = "Simulação executada com sucesso."
    else:
        check = "\n".join(check_msgs)

    results = {
        "time_levels": time_levels,
        "nodes": nodes_series,
        "check": check,
    }

    if conc_series:
        results["concentrations"] = conc_series
    if temp_series:
        results["temperatures"] = temp_series
    if balance:
        results["balance"] = balance

    return results
