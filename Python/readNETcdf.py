import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr


def process(filename: str, args: ArgumentParser) -> None:

    # normal
    PSD_xr = xr.open_dataset(filename, group="PSD")
    CSD_xr = xr.open_dataset(filename, group="CSD")

    print("PSD\n", PSD_xr)
    print("CSD\n", CSD_xr)

    # welch
    wPSD_xr = xr.open_dataset(filename, group="wPSD")
    wCSD_xr = xr.open_dataset(filename, group="wCSD")
    
    print("wPSD\n", wPSD_xr)
    print("wCSD\n", wCSD_xr)

    # banded
    bPSD_xr = xr.open_dataset(filename, group="bPSD")
    bCSD_xr = xr.open_dataset(filename, group="bCSD")
    
    print("bPSD\n", bPSD_xr)
    print("bCSD\n", bCSD_xr)


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