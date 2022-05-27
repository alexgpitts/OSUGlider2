
"""
Description: Testing "sanbox" file to read a NETCDF file for testing during development

to run: 
    python .\readNETcdf.py "your file.nc"

Author: Alex Pitts
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr
import getdata as gd


def process(filename: str, args: ArgumentParser) -> None:

    # # normal
    # PSD_xr = xr.open_dataset(filename, group="PSD")
    # print("PSD\n", PSD_xr, "\n\n")
    
    # # welch
    # wPSD_xr = xr.open_dataset(filename, group="WelchPSD")
    # print("Welch PSD\n", wPSD_xr, "\n\n")

    # # banded
    # bPSD_xr = xr.open_dataset(filename, group="BandedPSD")
    # print("Banded PSD\n", bPSD_xr, "\n\n")
    
    XYZ = xr.open_dataset(filename, group="XYZ")
    SR = XYZ.SampleRate.to_numpy()
    
    XYZ_ACC = {
        "t": XYZ.t.to_numpy(),
        "x": XYZ.x.to_numpy(),
        "y": XYZ.y.to_numpy(),
        "z": XYZ.z.to_numpy(),
    }
    print("t=", XYZ_ACC["t"], "\n\n")
    print("x=", gd.calcAcceleration(XYZ_ACC["x"], SR), "\n\n")
    print("y=", gd.calcAcceleration(XYZ_ACC["y"], SR), "\n\n")
    print("z=", gd.calcAcceleration(XYZ_ACC["z"], SR), "\n\n")

    # print("t=", XYZ_ACC["t"], "\n\n")
    # print("x=", XYZ_ACC["x"], "\n\n")
    # print("y=", XYZ_ACC["y"], "\n\n")
    # print("z=", XYZ_ACC["z"], "\n\n")
    


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