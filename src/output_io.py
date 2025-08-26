import pandas as pd
from pathlib import Path

def write_outputs(output_dir: str, results: dict):
    """
    Salva resultados da simulação em arquivos CSV/logs.
    results: dicionário com diferentes resultados numéricos
    """

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # 1. Arquivo de checagem
    if "check" in results:
        with open(out_path / "check.log", "w") as f:
            f.write(results["check"])

    # 2. Informações da simulação
    if "run_info" in results:
        pd.DataFrame(results["run_info"]).to_csv(out_path / "run_info.csv", index=False)

    # 3. Saídas por nível de tempo
    if "time_levels" in results:
        pd.DataFrame(results["time_levels"]).to_csv(out_path / "time_levels.csv", index=False)

    # 4. Informações nodais
    if "nodes" in results:
        pd.DataFrame(results["nodes"]).to_csv(out_path / "nodes.csv", index=False)

    # 5. Balanço
    if "balance" in results:
        pd.DataFrame(results["balance"]).to_csv(out_path / "balance.csv", index=False)

    # 6. Observações
    if "observations" in results:
        pd.DataFrame(results["observations"]).to_csv(out_path / "observations.csv", index=False)

    # 7. Perfis finais
    if "profile" in results:
        pd.DataFrame(results["profile"]).to_csv(out_path / "profile_out.csv", index=False)

    print(f"✅ Resultados salvos em {out_path}")
