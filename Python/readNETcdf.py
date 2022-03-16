import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr
import netCDF4 as nc4
import os

def process(filename: str, args: ArgumentParser) -> None:
    PSD_xr = xr.open_dataset(filename, group="PSD")
    CSD_xr = xr.open_dataset(filename, group="CSD")

    print(PSD_xr)
    print(CSD_xr)



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