"""
This file contains all the code for opening our NetCDF file and storing the data in dictionaries

Authors: Alex Pitts, Benjamin Cha, Clayton Surgeon 
"""
from argparse import ArgumentParser
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os

def output(args: ArgumentParser, substring: str) -> str:
    folder_name = os.path.split(args.nc[0])[0]
    output_name = os.path.splitext(os.path.basename(args.nc[0]))[0] + substring
    output_dir = os.path.join(folder_name, output_name)
    return output_dir

def calcAcceleration(x: np.array, fs: float) -> np.array:
    """converts displacement data to acceleration.
    We will not use this in final implementation because 
    the xyz data will already be in acceleration from the glider"""
    dx2 = np.zeros(x.shape)
    dx2[2:] = np.diff(np.diff(x))
    dx2[0:2] = dx2[2]
    return dx2 * fs * fs

def Data(filename: str, cdip: bool) -> dict:
    """Master data reading function. Reads the .nc file.
    The data is stored in dictionary (data), which contains many dictionaries 
    to hold information. Examples include: acceleration data, frequency bounds, etc."""
   
    meta_xr = xr.open_dataset(filename, group="Meta")  # For water depth
    wave_xr = xr.open_dataset(filename, group="Wave")
    xyz_xr = xr.open_dataset(filename, group="XYZ")
    try:
        frequency = float(xyz_xr.SampleRate)
    except:
        frequency = float(xyz_xr.SampleRate[0])

    # read in the meta, xyz, and wave groups
    data = {

        "Meta": {
            "frequency": frequency,
            "latitude": float(meta_xr.DeployLatitude),
            "longitude": float(meta_xr.DeployLongitude),
            "depth": float(meta_xr.WaterDepth),
            "declination": float(meta_xr.Declination)
        },    
          
        "Wave": {
            "Timebounds": wave_xr.TimeBounds.to_numpy(),
            "time_lower": wave_xr.TimeBounds[:, 0].to_numpy(),
            "time_upper": wave_xr.TimeBounds[:, 1].to_numpy(),
            "FreqBounds": wave_xr.FreqBounds.to_numpy(),
            "Bandwidth": wave_xr.Bandwidth.to_numpy()
        },
        

        "Freq": {
            "lower": wave_xr.FreqBounds[:, 0].to_numpy(),
            "upper": wave_xr.FreqBounds[:, 1].to_numpy(),
            "joint": wave_xr.FreqBounds[:, :].to_numpy()
        }
    }

  
    # if we are dealing with CDIP data, then convert from displacement to acceleration
    if cdip:
        data["XYZ"] = {
            "t": xyz_xr.t.to_numpy(),
            "x": calcAcceleration(xyz_xr.x.to_numpy(), frequency),
            "y": calcAcceleration(xyz_xr.y.to_numpy(), frequency),
            "z": calcAcceleration(xyz_xr.z.to_numpy(), frequency)
        }
        print("in CDIP")

    # otherwise just read the acceleration data from the glider
    else:
        data["XYZ"] = {
            "t": xyz_xr.t.to_numpy(),
            "x": xyz_xr.x.to_numpy(),
            "y": xyz_xr.y.to_numpy(),
            "z": xyz_xr.z.to_numpy()
        }
    
    

    return data