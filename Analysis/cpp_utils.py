import subprocess
import pathlib as pl
import pandas as pd

import os

def rootPath():
    cur_path = pl.Path(__file__)
    root_path = cur_path

    while root_path.name != "variational-monte-carlo-fys4411":
        root_path = root_path.parent

    return root_path

def vmcPath():
    vmc_path = rootPath() / pl.Path("build/vmc")
    return vmc_path

def dataPath(filename):
    filename_path = rootPath() / pl.Path(f"Data/{filename}")
    return filename_path

def vmcRun(D=3, N=10, logMet=6, logEq=5, omega=1.0, alpha=0.5, stepLength=0.1, importance=False, analytical=True, timing=False, GD=False, filename="test.txt"):
    vmc_path = vmcPath()
    filename_path = dataPath(filename)

    assert vmc_path.exists(), f"I cannot find {vmc_path} :((, are you sure you have compiled?"
    args = [
        vmc_path,
        D,
        N,
        logMet,
        logEq,
        omega,
        alpha,
        stepLength,
        int(importance),
        int(analytical),
        int(timing),
        int(GD),
        filename_path,
    ]

    if not filename:
        args.pop()

    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)

def vmcRunTiming(D=3, N=10, logMet=6, logEq=5, omega=1.0, alpha=0.5, stepLength=0.1, analytical=True, filename="timing.txt"):
    vmc_path = vmcPath()
    filename_path = dataPath(filename)

    assert vmc_path.exists(), f"I cannot find {vmc_path} :((, are you sure you have compiled?"
    args = [
        vmc_path,
        D,
        N,
        logMet,
        logEq,
        omega,
        alpha,
        stepLength,
        int(analytical),
        filename_path,
    ]

    if not filename:
        args.pop()

    args_run = [str(arg) for arg in args]

    subprocess.run(args_run)

def vmcLoad(filename):
    filename_path = dataPath(filename)

    df = pd.read_csv(filename_path, delim_whitespace=True)

    int_cols = ["Dimensions", "Particles" ,"Metro-steps"]
    numeric_cols = [col for col in df.columns if col not in int_cols]
    
    for col in numeric_cols:
        df[col] = df[col].astype(float)

    return df