import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr
import netCDF4 as nc4
import os


import processdata as processdata



def process(filename: str, args: ArgumentParser) -> None:
    # construct an output .nc file
    folder_name = os.path.split(args.nc[0])[0]
    output_name = os.path.splitext(os.path.basename(args.nc[0]))[0] + "_output.nc"
    output_dir = os.path.join(folder_name, output_name)
    nc4.Dataset(output_dir, 'w', format='NETCDF4')

    # calculate CSDs and PSDs for each time block 
    CSDs, PSDs = processdata.getPSDs(filename)
    
    print(PSDs)
    print(CSDs)
    

    
    # where I left off

    # PSDs_dataset = xr.Dataset(PSDs)
    # CSDs_dataset = xr.Dataset(CSDs)


    # write PSDs and CSDs to output file in PSD and CSD group
    # CSDs_dataset.to_netcdf(output_dir, mode="w", group="PSD")
    # PSDs_dataset.to_netcdf(output_dir, mode="a", group="CSD")
    
        



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