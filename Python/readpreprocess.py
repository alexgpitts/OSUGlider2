
"""
Description: Testing file to read a NETCDF file

to run: 
    python .\readNETcdf.py "your file.nc"

Author: Alex Pitts
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr


def process(filename: str, args: ArgumentParser) -> None:

    # normal
    meta_xr = xr.open_dataset(filename, group="Meta")
    print("Meta\n", meta_xr, "\n\n")

    XYZ_xr = xr.open_dataset(filename, group="XYZ")
    print("XYZ\n", XYZ_xr, "\n\n")

    wave_xr = xr.open_dataset(filename, group="Wave")
    print("Wave\n", wave_xr, "\n\n")

    
    
    


def main(raw_args=None):
    # command line stuff
    parser = ArgumentParser()
    grp = parser.add_mutually_exclusive_group()
    
    # required
    parser.add_argument("nc", nargs="+", type=str, help="netCDF file to process") 
    
    args = parser.parse_args(raw_args)
    # for each nc filename passed
    for fn in args.nc:
        process(fn, args)


if __name__ == "__main__":
    main()