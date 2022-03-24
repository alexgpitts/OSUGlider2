import numpy as np
import pandas as pd
import xarray as xr
import os
import time
import processdata as pr
from argparse import ArgumentParser
from getdata import Data, calcAcceleration, output

"""
Finds a CSV file in current directory, parses it and stores it in nc file
"""

def parse_csv(fn: str, output: str):
    if (not os.path.isfile(fn)): # if nc file isn't in current directory 
        raise Exception(f"Input CSV file, {fn}, not found!") # throw error if file isn't found

    # global attribute
    ds = xr.Dataset(
        attrs=dict(
            date_created=str(np.datetime64(round(1e3*time.time()), "ms")),
            ),
        )

    # This stops is it from resetting the file over and over because of "w"
    if (not os.path.isfile(output)): # if this file is not found in the directory
        ds.to_netcdf(output, 'w', format='NETCDF4') # constructs an nc file

    # reads csv file, stores it in a dataframe then into a nc file
    data = pd.read_csv(fn)
    return xr.Dataset.from_dataframe(data) # format and stores csv file into dataset
     


"""
Finds CSV file to parse
"""
def store_data(args: ArgumentParser) -> None:
    substring = "_output.nc"
    output_name = output(args, substring)
    # # if group is entered and csv file is found in the current directory, then parse that file
    if ("Meta" in args.group):
        meta = parse_csv(args.meta, output_name) # be a place to add attributes by using dict
        meta.to_netcdf(output_name, mode="a", group="Meta", engine="netcdf4")
    if ("Wave" in args.group):
        wave = parse_csv(args.wave, output_name)

        # construct array of timebounds
        tlow = wave.time_lower.to_numpy()
        tup = wave.time_upper.to_numpy()
        TimeBounds = np.vstack((tlow, tup)).T
        wave["TimeBounds"] = pr.merge(TimeBounds)

        #construct frequency bounds and bandwidth
        freq = parse_csv(args.freq, output_name)
        flow = freq.fLower.to_numpy()
        fup = freq.fUpper.to_numpy()
        FreqBounds = np.vstack((flow, fup)).T

        #add frequency bounds to wave
        wave["FreqBounds"] = pr.merge(FreqBounds)
        wave["Bandwidth"] = freq.Bandwidth
        
        wave.to_netcdf(output_name, mode="a", group="Wave", engine="netcdf4")
    if ("XYZ" in args.group):
        xyz = parse_csv(args.xyz, output_name)
        xyz.to_netcdf(output_name, mode="a", group="XYZ", engine="netcdf4")
    

"""
Stores acceleration data in nc file
"""
def acceleration_XYZ(fn: str, args: ArgumentParser) -> None:
    substring = "_output.nc"
    # if ("XYZ" in args.group): # prevents a clash if acc XYZ is already written in
    #     print("Already exists!")
    #     return 

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
    ds.to_netcdf(output(args, substring), mode="a", group="XYZ") # writes acceleration data to nc file


def main(raw_args=None):
    # command line stuff
    parser = ArgumentParser()
    grp = parser.add_mutually_exclusive_group()
    meta = "./ncFiles/067.20201225_1200.20201225_1600_meta.csv"
    wave = "./ncFiles/067.20201225_1200.20201225_1600_wave.csv"
    xyz = "./ncFiles/067.20201225_1200.20201225_1600_xyz.csv"
    freq = "./ncFiles/067.20201225_1200.20201225_1600_freq.csv"
    # required
    parser.add_argument("nc", nargs="+", type=str, help="netCDF file to process") 
    parser.add_argument("--meta", type=str, metavar="meta.csv", default=meta, help="For Meta Group")
    parser.add_argument("--wave", type=str, metavar="wave.csv", default=wave, help="For Wave Group")
    parser.add_argument("--xyz", type=str, metavar="acceleration.csv", default=xyz, help="For XYZ Group")
    parser.add_argument("--freq", type=str, metavar="freq.csv", default=freq, help="For freq bounds")
    parser.add_argument("--group", type=str, action="append", required=False, choices=("Meta", "Wave", "XYZ"), help="Enter Meta, Wave or XYZ")
    args = parser.parse_args(raw_args)

    # for each nc filename passed
    for fn in args.nc:
        # acceleration_XYZ(fn, args)
        store_data(args)
                

if __name__ == "__main__":
    main()