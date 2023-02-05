import subprocess
import pathlib as pl

def runVMC(*args, **kwargs):
    opt = {
        "filename": ""
    }

    opt.update(kwargs)

    cur_path = pl.Path(__file__)
    root_path = cur_path

    while root_path.name != "variational-monte-carlo-fys4411":
        root_path = root_path.parent

    vmc_path = root_path / pl.Path("build/vmc")
    filename_path = root_path / pl.Path(f"Data/{opt['filename']}")

    assert vmc_path.exists(), f"I cannot find {vmc_path} :((, are you sure you have compiled?"

    subprocess.run([
        vmc_path,
        *args,
        filename_path,
    ])