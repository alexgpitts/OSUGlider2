
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
    CSD_xr = xr.open_dataset(filename, group="CSD")

    print("PSD\n", PSD_xr, "\n\n")
    print("CSD\n", CSD_xr, "\n\n")

    # welch
    wPSD_xr = xr.open_dataset(filename, group="WelchPSD")
    wCSD_xr = xr.open_dataset(filename, group="WelchCSD")
    
    print("Welch PSD\n", wPSD_xr, "\n\n")
    print("Welch CSD\n", wCSD_xr, "\n\n")

    # banded
    bPSD_xr = xr.open_dataset(filename, group="BandedPSD")
    bCSD_xr = xr.open_dataset(filename, group="BandedCSD")
    
    print("Banded PSD\n", bPSD_xr, "\n\n")
    print("Banded CSD\n", bCSD_xr, "\n\n")


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