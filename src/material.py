# material.py
import numpy as np

EPS = 1e-12


def theta_vg(h, Par):
    """
    Curva de retenção van Genuchten.
    Par = [Qr, Qs, Alfa, n, Ks, (opcional) BPar]
    h   : pressão matricial (cm). Aceita escalar ou array.
    """
    Qr, Qs, Alfa, n = float(Par[0]), float(Par[1]), float(Par[2]), float(Par[3])
    m = 1.0 - 1.0 / n
    h_abs = np.maximum(np.abs(h), EPS)
    Se = (1.0 + (Alfa * h_abs) ** n) ** (-m)
    return Qr + (Qs - Qr) * Se


def c_vg(h, Par):
    """
    Capacitância hidráulica aproximada: C = dθ/dh.
    Usa derivada numérica central para funcionar para qualquer Par.
    """
    h = np.asarray(h, dtype=float)
    eps = 1e-6
    return (theta_vg(h + eps, Par) - theta_vg(h - eps, Par)) / (2.0 * eps)


def fk(iModel: int, h, Par):
    """
    Condutividade hidráulica não saturada K(h).
    Implementação básica para van Genuchten (iModel=0 ou 3) e fallback.
    Aceita escalar ou array para h.
    Par = [Qr, Qs, Alfa, n, Ks, (opcional) BPar]
    """
    Qr, Qs, Alfa, n = float(Par[0]), float(Par[1]), float(Par[2]), float(Par[3])
    Ks = max(float(Par[4]), 1e-37)
    BPar = float(Par[5]) if len(Par) > 5 else 0.5

    h = np.asarray(h, dtype=float)

    if iModel in (0, 3):  # van Genuchten clássico e variante com he=const
        m = 1.0 - 1.0 / n
        Se = (1.0 + (Alfa * np.maximum(np.abs(h), EPS)) ** n) ** (-m)
        # Mualem-van Genuchten:
        # K = Ks * Se^B * [1 - (1 - Se^(1/m))^m]^2
        term = 1.0 - (1.0 - np.power(Se, 1.0 / m)) ** m
        K = Ks * np.power(Se, BPar) * np.power(term, 2.0)
        return K

    # Brooks-Corey (esqueleto simples; expanda conforme necessário)
    if iModel == 2:
        # Espera Par extra: Lambda (Par[5]), hb (Par[6])
        Lambda = float(Par[5])
        hb = max(float(Par[6]), EPS)
        Se = (np.maximum(np.abs(h), EPS) / hb) ** (-Lambda)
        Se = np.clip(Se, 0.0, 1.0)
        return Ks * Se ** ((2 + 3 * Lambda) / Lambda)

    # Fallback: retorna Ks (saturado)
    return Ks * np.ones_like(h, dtype=float)


def load_materials(config: dict):
    """
    Stub de carregamento de materiais. Substitua pela leitura real do Selector/Profile.
    Retorna lista de materiais; cada material é um dict com 'model' e 'Par'.
    """
    # Exemplo: único material van Genuchten
    mats = [{
        "model": 0,
        # [Qr, Qs, Alfa (1/cm), n, Ks (cm/dia), BPar]
        "Par": [0.05, 0.45, 0.01, 1.6, 10.0, 0.5],
    }]
    return mats
