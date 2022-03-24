
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
import getdata as gd
import csv


def process(filename: str, args: ArgumentParser) -> None:

    
    data = gd.Data(filename)
    
    
    wave_name = gd.output(args, "_wave.csv")
    meta_name = gd.output(args, "_meta.csv")
    xyz_name = gd.output(args, "_xyz.csv")
    # write the time bounds from wave group to wave.csv
    with open(wave_name, 'w', newline='') as f: 
        writer = csv.writer(f)
        header = ['time_lower','time_upper']
        writer.writerow(header)
        for i in range(len(data["Wave"]["Timebounds"])):
            row = [data["Wave"]["time_lower"][i], data["Wave"]["time_upper"][i]]
            
            writer.writerow(row)

    # write the time bounds from wave group to wave.csv
    with open(meta_name, 'w', newline='') as f: 
        writer = csv.writer(f)
        header = ['SampleRate','DeployLatitude','DeployLongitude','WaterDepth','Declination']
        writer.writerow(header)
        
        row = [data["Meta"]["frequency"], data["Meta"]["latitude"], data["Meta"]["longitude"], 
            data["Meta"]["depth"], data["Meta"]["declination"]]
            
        writer.writerow(row)

    # write the xyz data from XYZ group to xyz_acc.csv
    with open(xyz_name, 'w', newline='') as f: 
        writer = csv.writer(f)
        header = ['t','x','y','z','SampleRate']
        writer.writerow(header)
        x = gd.calcAcceleration(data["XYZ"]["x"], data["Meta"]["frequency"])
        y = gd.calcAcceleration(data["XYZ"]["x"], data["Meta"]["frequency"])
        z = gd.calcAcceleration(data["XYZ"]["x"], data["Meta"]["frequency"])
        for i in range(len(data["XYZ"]["t"])):
            row = [data["XYZ"]["t"][i], x[i], y[i], z[i], data["Meta"]["frequency"]]
            writer.writerow(row)

    # write the frequency bound data from FreqBounds to freq.csv
    # coming soon...
    # with open(xyz_name, 'w', newline='') as f: 
    #     writer = csv.writer(f)
    #     header = ['t','x','y','z','SampleRate']
    #     writer.writerow(header)
    #     x = gd.calcAcceleration(data["XYZ"]["x"], data["Meta"]["frequency"])
    #     y = gd.calcAcceleration(data["XYZ"]["x"], data["Meta"]["frequency"])
    #     z = gd.calcAcceleration(data["XYZ"]["x"], data["Meta"]["frequency"])
    #     for i in range(len(data["XYZ"]["t"])):
    #         row = [data["XYZ"]["t"][i], x[i], y[i], z[i], data["Meta"]["frequency"]]
    #         writer.writerow(row)
    

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