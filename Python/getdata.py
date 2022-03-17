"""
This file contains all the code for opening our NetCDF file and storing the data in dictionaries

Authors: 
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

def calcAcceleration(x: np.array, fs: float) -> np.array:
    """converts displacement data to acceleration.
    We will not use this in final implementation because 
    the xyz data will already be in acceleration from the glider"""
    dx2 = np.zeros(x.shape)
    dx2[2:] = np.diff(np.diff(x))
    dx2[0:2] = dx2[2]
    return dx2 * fs * fs


def Data(filename) -> dict:
    """Master data reading function. Reads the .nc file from CDIP.
    The data is stored in dictionary (data), which contains many dictionaries 
    to hold information. Examples include: acceleration data, frequency bounds, etc."""
   

    meta_xr = xr.open_dataset(filename, group="Meta")  # For water depth
    wave_xr = xr.open_dataset(filename, group="Wave")
    xyz_xr = xr.open_dataset(filename, group="XYZ")

    frequency = float(xyz_xr.SampleRate)

    # read in the meta, xyz, and wave groups
    data = {

        "Meta": {
            "frequency": frequency,
            "latitude": float(meta_xr.DeployLatitude),
            "longitude": float(meta_xr.DeployLongitude),
            "depth": float(meta_xr.WaterDepth),
            "declination": float(meta_xr.Declination)
        },
        
        "XYZ": {
            "t": xyz_xr.t.to_numpy(),
            "x": calcAcceleration(xyz_xr.x.to_numpy(), frequency),
            "y": calcAcceleration(xyz_xr.y.to_numpy(), frequency),
            "z": calcAcceleration(xyz_xr.z.to_numpy(), frequency),
        },
        
        "Wave": {
            "Timebounds": wave_xr.TimeBounds,
            "time_lower": wave_xr.TimeBounds[:, 0],
            "time_upper": wave_xr.TimeBounds[:, 1],
            "FreqBounds": wave_xr.FreqBounds,
            "Bandwidth": wave_xr.Bandwidth
        }
    }

    return data