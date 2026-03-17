import yaml
import os

from src.pipeline import run_pt_bin

with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

for pt_bin in config["pt_bins"]:

    print(f"Running pT bin {pt_bin}")

    run_pt_bin(pt_bin, config)    