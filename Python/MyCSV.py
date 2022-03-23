import numpy as np
import pandas as pd
import xarray as xr
import os
import time
from argparse import ArgumentParser
from getdata import Data, calcAcceleration, output

"""
Finds a CSV file in current directory, parses it and stores it in nc file
"""

def parse_csv(fn: str, output: str, group: str) -> None:
    if (not os.path.isfile(fn)): # if nc file isn't in current directory 
        raise Exception(f"Input CSV file, {fn}, not found!") # throw error if file isn't found

    # os.listdir
    # return a list of files in directory thats given
    # glob - look up the module for csv files
    # return all csv files

    # reads csv file, stores it in a dataframe then into a nc file
    data = pd.read_csv(fn)
    ds = xr.Dataset.from_dataframe(data) # format and stores csv file into dataset
    ds.to_netcdf(output, mode="a", group=group)


"""
Finds CSV file to parse
"""
def store_data(args: ArgumentParser) -> None:
    # # if group is entered and csv file is found in the current directory, then parse that file
    if ("Meta" in args.group):
        parse_csv(args.meta, output(args), "Meta") # be a place to add attributes by using dict
    if ("Wave" in args.group):
        parse_csv(args.wave, output(args), "Wave")
    if ("XYZ" in args.group):
        parse_csv(args.xyz, output(args), "XYZ")
    
    
"""
Stores acceleration data in nc file
"""
def acceleration_XYZ(fn: str, args: ArgumentParser) -> None:
    if ("XYZ" in args.group): # prevents a clash if acc XYZ is already written in
        return

    xyz_xr = xr.open_dataset(args.nc[0], group="XYZ")
    frequency = float(xyz_xr.SampleRate)

    data = Data(args.nc[0])
    # # calculates for acceleration and stores in XYZ
    ds = xr.Dataset(
        data_vars=dict(
            x=(("time",), calcAcceleration(data["TXYZ"]["x"], frequency), {"units": "m/s^2"}),
            y=(("time",), calcAcceleration(data["TXYZ"]["y"], frequency), {"units": "m/s^2"}),
            z=(("time",), calcAcceleration(data["TXYZ"]["z"], frequency), {"units": "m/s^2"}),
            ),
        attrs=dict(
            comment="Acceleration Values",
            ),
        )
    ds.to_netcdf(output(args), mode="a", group="XYZ") # writes acceleration data to nc file


