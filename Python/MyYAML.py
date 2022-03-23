import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# python .\Driver.py --yaml=meta.yml '.\ncFiles\067.20201225_1200.20201225_1600.nc'
# pip install pyyaml
def load_meta(fn:str) -> dict:
    try:
        with open(fn, "r") as fp:
            lines = fp.read()
            data = yaml.load(lines, Loader=yaml.SafeLoader)

            meta = xr.Dataset ({
                "frequency": data["meta"]["frequency"],
                "longitude": data["meta"]["deployLongitude"],
                "latitude": data["meta"]["deployLatitude"],
                "depth": data["meta"]["depth"],
                "declination": data["meta"]["declination"]
            })
        print(meta)
        return meta
    except Exception as e:
        raise e

def merge(data):
    """takes a list and merges the elements into a dataset"""
    dataset = xr.DataArray(
        np.array(data, dtype=object),
        dims = ("timeblock", "value"),
        coords = {
            "value": (np.arange(0, len(data[0]))).tolist(),
            "timeblock": (np.arange(0, len(data))).tolist(),
        }
    )
    return dataset 

