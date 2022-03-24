
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
    PSD_xr = xr.open_dataset(filename, group="PSD")
    print("PSD\n", PSD_xr, "\n\n")
    
    # welch
    wPSD_xr = xr.open_dataset(filename, group="WelchPSD")
    print("Welch PSD\n", wPSD_xr, "\n\n")

    # banded
    bPSD_xr = xr.open_dataset(filename, group="BandedPSD")
    print("Banded PSD\n", bPSD_xr, "\n\n")
    


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