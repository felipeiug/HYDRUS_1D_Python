from pathlib import Path

def read_inputs(input_dir: str):
    """
    Lê arquivos de entrada do HYDRUS (versão Python).
    Retorna: config (dict), profile (dict), options (dict)
    """

    input_path = Path(input_dir)

    # === 1. Selector.in ===
    selector_file = input_path / "Selector.in"
    config = {}

    with open(selector_file, "r") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    # No Fortran, as leituras eram feitas em sequência, ex:
    # read(30,*) -> linha ignorada
    # read(30,'(a)') Hed -> cabeçalho
    config["header"] = lines[2]  # exemplo: Hed
    config["LUnit"] = lines[4]
    config["TUnit"] = lines[5]
    config["MUnit"] = lines[6]

    # Flags principais (True/False)
    # Aqui simplificado: convertemos "1/0" ou ".true./.false." para bool
    def to_bool(val):
        return str(val).strip().lower() in ["1", "true", ".true."]

    flags = lines[8].split()
    config["lWat"]   = to_bool(flags[0])
    config["lChem"]  = to_bool(flags[1])
    config["lTemp"]  = to_bool(flags[2])
    config["SinkF"]  = to_bool(flags[3])
    config["lRoot"]  = to_bool(flags[4])
    config["ShortO"] = to_bool(flags[5])
    # ... e assim por diante

    # Materiais e camadas
    nm_line = lines[12].split()
    config["NMat"]  = int(nm_line[0])
    config["NLay"]  = int(nm_line[1])
    config["CosAlf"] = float(nm_line[2])

    # Tolerâncias
    tol_line = lines[16].split()
    config["MaxIt"] = int(tol_line[0])
    config["TolTh"] = float(tol_line[1])
    config["TolH"]  = float(tol_line[2])

    # Condições de contorno superiores
    top_line = lines[18].split()
    config["TopInF"] = to_bool(top_line[0])
    config["WLayer"] = to_bool(top_line[1])
    config["KodTop"] = int(top_line[2])
    config["lInitW"] = to_bool(top_line[3])

    # (o resto continua... este é só o esqueleto)

    # === 2. Profile.dat ===
    profile_file = input_path / "Profile.dat"
    profile = {}
    with open(profile_file, "r") as f:
        profile["raw"] = f.read()

    # === 3. Options.in ===
    options_file = input_path / "Options.in"
    options = {}
    if options_file.exists():
        with open(options_file, "r") as f:
            options["raw"] = f.read()

    return config, profile, options
