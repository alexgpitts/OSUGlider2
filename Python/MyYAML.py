#! /usr/bin/env python3

import yaml
import logging
import xarray as xr
import netCDF4 as nc4
import os
from argparse import ArgumentParser


def load_meta(fn:str) -> dict:
    # try:
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
    # except Exception as e:
    #     raise e
    
# def main():
#     parser = ArgumentParser()
#     parser.add_argument("yaml", nargs="+", metavar="fn.yml", help="YAML file(s) to load")
#     args = parser.parse_args()

#     for fn in args.yaml:
#         load_meta(fn)

# if __name__ == "__main__":
#     main()