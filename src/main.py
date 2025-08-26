"""
HYDRUS 1D - Python version (esqueleto inicial)

Modelo numérico de fluxo de água, transporte de calor e solutos em meio poroso não saturado.
Baseado no código Fortran original (Simunek, Sejna, van Genuchten).
"""

import sys
import signal
import numpy as np

from input_io import read_inputs
from output_io import write_outputs
from material import load_materials
from watflow import solve_water_flow
from solute import solve_solute_transport
from temper import solve_heat_transport
from time_control import run_time_loop


def signal_handler(sig, frame):
    print("⚠️ Execução interrompida pelo usuário (CTRL+C).")
    sys.exit(0)


def main():
    # Captura interrupção (CTRL+C)
    signal.signal(signal.SIGINT, signal_handler)

    print("=== HYDRUS 1D (Python) ===")

    # 1. Ler arquivos de entrada
    config, profile, options = read_inputs("data/input/")

    # 2. Carregar propriedades de materiais
    materials = load_materials(config)

    # 3. Inicializar estado do modelo
    state = {
        "water": None,
        "solute": None,
        "heat": None,
    }

    # 4. Executar simulação no tempo
    results = run_time_loop(
        config=config,
        profile=profile,
        materials=materials,
        state=state,
        solver_water=solve_water_flow,
        solver_solute=solve_solute_transport,
        solver_heat=solve_heat_transport,
    )

    # 5. Escrever saídas
    write_outputs("data/output/", results)

    print("✅ Simulação concluída.")


if __name__ == "__main__":
    main()
