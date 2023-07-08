import os
import yaml

my_path = os.path.dirname(__file__)

potcar_md5sum_data = {}
potcar_nelect = {}

potcar_sets = [
    "PBE",
    "PBE_52",
    "PBE_54",
    "PBE_54r",
    "LDA_54r",
    "GGA",
    "USPP_GGA",
    "LDA",
    "LDA_52",
    "LDA_54",
    "USPP_LDA",
]

for potcar_set in potcar_sets:
    with open(f"{os.path.join(my_path, potcar_set)}_md5.yaml", "r") as stream:
        potcar_md5sum_data[potcar_set] = yaml.load(stream, Loader=yaml.SafeLoader)

for potcar_set in potcar_sets:
    with open(f"{os.path.join(my_path, potcar_set)}_nelect.yaml", "r") as stream:
        potcar_nelect[potcar_set] = yaml.load(stream, Loader=yaml.SafeLoader)
