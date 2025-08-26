
import numpy as np
from material import fk, theta_vg, c_vg

def solve_water_flow(state, config, materials):
    NumNP = config.get("NumNP", 100)
    dz = config.get("dz", 1.0)
    dt = config.get("dt", 0.1)
    MaxIt = config.get("MaxIt", 20)
    TolH = config.get("TolH", 1e-6)

    Par = materials[0]["Par"]
    h = state["h"].copy()
    theta = theta_vg(h, Par)

    h_top = config.get("h_top", -100.0)
    q_bot = config.get("q_bot", 0.0)

    for newton in range(MaxIt):
        K = np.array([fk(0, hi, Par) for hi in h])
        C = c_vg(h, Par)

        a = np.zeros(NumNP)
        b = np.zeros(NumNP)
        c = np.zeros(NumNP)
        r = np.zeros(NumNP)

        for i in range(1, NumNP-1):
            Ki_plus = 0.5 * (K[i] + K[i+1])
            Ki_minus = 0.5 * (K[i] + K[i-1])

            flux_plus = -Ki_plus * (h[i+1] - h[i]) / dz
            flux_minus = -Ki_minus * (h[i] - h[i-1]) / dz

            r[i] = (theta[i] - state["theta"][i]) / dt - (flux_plus - flux_minus) / dz

            a[i] =  Ki_minus / (dz*dz)
            c[i] =  Ki_plus  / (dz*dz)
            b[i] = C[i] / dt + a[i] + c[i]

        b[0] = 1.0
        r[0] = h[0] - h_top

        Kb = K[-1]
        b[-1] = 1.0
        r[-1] = (h[-1] - h[-2]) / dz + q_bot / max(Kb, 1e-12)

        dh = solve_tridiagonal(a, b, c, -r)
        h += dh
        theta = theta_vg(h, Par)

        if np.linalg.norm(dh, ord=np.inf) < TolH:
            break

    state["h"] = h
    state["theta"] = theta
    return {"h": h.copy(), "theta": theta.copy()}

def solve_tridiagonal(a, b, c, d):
    n = len(b)
    cp = np.zeros(n)
    dp = np.zeros(n)
    x = np.zeros(n)

    cp[0] = c[0]/b[0]
    dp[0] = d[0]/b[0]
    for i in range(1, n):
        denom = b[i] - a[i]*cp[i-1]
        cp[i] = c[i]/denom if i < n-1 else 0.0
        dp[i] = (d[i] - a[i]*dp[i-1]) / denom

    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i]*x[i+1]
    return x
