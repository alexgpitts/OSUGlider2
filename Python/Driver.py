"""
Driver file to process a .nc file and store PSDs, CSDs and other calculations in an _output.nc file
Authors: Alex Pitts, Benjamin Cha, Clayton Surgeon 
"""


import numpy as np
import getdata as gd
from argparse import ArgumentParser
from getdata import output
import xarray as xr
import netCDF4 as nc4
import os
import processdata as pr
import MyYAML as meta

import warnings
warnings.filterwarnings("error")

# python .\Driver.py --group=Meta --group=Wave --group=XYZ '.\ncFiles\067.20201225_1200.20201225_1600.nc'
# (at least 1 group has to be used) python .\Driver.py --yaml=meta.yaml --group=Wave '.\ncFiles\067.20201225_1200.20201225_1600.nc'

def __process(fn: str, args: ArgumentParser) -> None:
    # stores meta data from yaml file
    if args.yaml:
        xr.Dataset(meta.load_meta(fn)).to_netcdf(output(args), mode="a", group="Meta")

def process(filename: str, args: ArgumentParser) -> None:
    # A function in getdata has this implemented so it can be pull in other files and reduce repeated code
    substring = "_processed.nc"
    nc4.Dataset(output(args, substring), 'w', format='NETCDF4') # constructs an nc file
    
    data = gd.Data(filename, args.cdip)
    
    
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]

    # lists of upper, lower, and midpoint frequencies for banding
    freq_bounds = data["Wave"]["FreqBounds"]
    freq_lower = data["Freq"]["lower"]
    freq_upper = data["Freq"]["upper"]
    freq_midpoints = data["Freq"]["joint"].mean(axis=1)

    # lists to store PSDs
    PSDs = []       #normal PSD
    wPSDs = []      #Welch method PSD
    bPSDs = []      #banding methods PSDs
    wCalcs = []     #welch method calculations
    bCalcs = []     #banding method calculations

    for i in range(len(timebounds)):
        print(f"window {i}")
        time_lower = data["Wave"]["time_lower"][i]
        time_upper = data["Wave"]["time_upper"][i]

        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )
        
        averaging_window = 2

        # select xyz data from block
        acc = {
            # x is northwards
            "x": pr.Rolling_mean(data["XYZ"]["x"][select], averaging_window),
            # y is eastwards
            "y": pr.Rolling_mean(data["XYZ"]["y"][select], averaging_window),
            # z is upwards
            "z": pr.Rolling_mean(data["XYZ"]["z"][select], averaging_window)
        }

        ################################# 
        # PSDs and CSDs
        #################################

        # preform FFT on block
        FFT = {
            "x": np.fft.rfft(acc["x"], n=acc["z"].size),  # northwards
            "y": np.fft.rfft(acc["y"], n=acc["z"].size),  # eastwards
            "z": np.fft.rfft(acc["z"], n=acc["z"].size),  # upwards
        }

        window_type = "hann"
        # preform FFT on block using welch mothod
        wFFT = {
            "x": pr.wfft(acc["x"], 2**8, window_type),
            "y": pr.wfft(acc["y"], 2**8, window_type),
            "z": pr.wfft(acc["z"], 2**8, window_type),
        }

        # Calculate PSD of data from normal FFT
        PSD = {
            # imaginary part is zero
            "xx": pr.calcPSD(FFT["x"], FFT["x"], data["Meta"]["frequency"], "boxcar").real,
            "yy": pr.calcPSD(FFT["y"], FFT["y"], data["Meta"]["frequency"], "boxcar").real,
            "zz": pr.calcPSD(FFT["z"], FFT["z"], data["Meta"]["frequency"], "boxcar").real,

            "xy": pr.calcPSD(FFT["x"], FFT["y"], data["Meta"]["frequency"], "boxcar"),
            "zx": pr.calcPSD(FFT["z"], FFT["x"], data["Meta"]["frequency"], "boxcar"),
            "zy": pr.calcPSD(FFT["z"], FFT["y"], data["Meta"]["frequency"], "boxcar"),

            "freq_space": np.fft.rfftfreq(acc["z"].size, 1/data["Meta"]["frequency"])
        }
        

        # calculate PSD on output from welch method FFT
        wPSD = {
            "xx": pr.wcalcPSD(wFFT["x"], wFFT["x"], data["Meta"]["frequency"], window_type).real,
            "yy": pr.wcalcPSD(wFFT["y"], wFFT["y"], data["Meta"]["frequency"], window_type).real,
            "zz": pr.wcalcPSD(wFFT["z"], wFFT["z"], data["Meta"]["frequency"], window_type).real,

            "xy": pr.wcalcPSD(wFFT["x"], wFFT["y"], data["Meta"]["frequency"], window_type),
            "zx": pr.wcalcPSD(wFFT["z"], wFFT["x"], data["Meta"]["frequency"], window_type),
            "zy": pr.wcalcPSD(wFFT["z"], wFFT["y"], data["Meta"]["frequency"], window_type),

            "freq_space": np.fft.rfftfreq(wFFT["z"][0].size*2-1, 1/data["Meta"]["frequency"])
        }        

        # bit mask so as to select only within the bounds of one lower:upper range pair
        freq_select = np.logical_and(
            np.less_equal.outer(freq_lower, PSD["freq_space"]),
            np.greater_equal.outer(freq_upper, PSD["freq_space"])
        )
        count = freq_select.sum(axis=1)
        window_type = "boxcar"
        windowing_method = pr.Bias(len(PSD["freq_space"]), window_type)
        

        # calculate PSD on With banded method
        bPSD = {
            "xx": (freq_select * PSD["xx"] * windowing_method).sum(axis=1) / count,
            "yy": (freq_select * PSD["yy"] * windowing_method).sum(axis=1) / count,
            "zz": (freq_select * PSD["zz"] * windowing_method).sum(axis=1) / count,

            "xy": (freq_select * PSD["xy"] * windowing_method).sum(axis=1) / count,
            "zx": (freq_select * PSD["zx"] * windowing_method).sum(axis=1) / count,
            "zy": (freq_select * PSD["zy"] * windowing_method).sum(axis=1) / count,

            "freq_space": freq_midpoints
        }
    

        # append PSDs
        PSDs.append(PSD)
        wPSDs.append(wPSD)
        bPSDs.append(bPSD)

        ################################# 
        # calculation
        #################################

        if(any(i > 100 for i in abs(acc["x"])) or 
        any(i > 100 for i in abs(acc["y"])) or 
        any(i > 100 for i in abs(acc["z"]))):
            print(f"block {i}:  bad data containing extremely large values")
            wCalcs.append(pr.errorCalc(len(wPSD["freq_space"])))     
            bCalcs.append(pr.errorCalc(len(freq_midpoints)))    
        else:
            try:
                # welch method calculations
                wCalcs.append(pr.welchCalc(wPSD, data))
                # banding method calculations
                bCalcs.append(pr.bandedCalc(bPSD, data))
            except:
                print(f"Error on window {i}")
                wCalcs.append(pr.errorCalc(len(wPSD["freq_space"])))     
                bCalcs.append(pr.errorCalc(len(freq_midpoints)))    
                pass
           

    #next step write  PSDs, wPSDs, bPSDs, wCalcs, and bCalcs to netCDF file using some type of custom dictionary merging functionl

    # PSDs and CSDs with normal calculation method
    PSD_Norm = pr.formatPSD(PSDs)
    xr.Dataset(PSD_Norm).to_netcdf(output(args, substring), mode="a", group="PSD")

    # PSDs and CSDs with Welch calculation method
    wPSD_Welch = pr.formatPSD(wPSDs)
    xr.Dataset(wPSD_Welch).to_netcdf(output(args, substring), mode="a", group="WelchPSD")

    # PSDs and CSDs with banding calculation method
    bPSD_Banded = pr.formatPSD(bPSDs)
    xr.Dataset(bPSD_Banded).to_netcdf(output(args, substring), mode="a", group="BandedPSD")
  
    # PSDs and CSDs with banded welch calculation method
    # comming soon...

    # calculations using welchPSD
    calcs_welch = pr.formatCalc(wCalcs)
    xr.Dataset(calcs_welch).to_netcdf(output(args, substring), mode="a", group="WelchWave")
    
    # calculations using BandedPSD
    calcs_banded = pr.formatCalc(bCalcs)
    xr.Dataset(calcs_banded).to_netcdf(output(args, substring), mode="a", group="Wave")

    # calculations using BandedWelchPSD
    # comming soon


def main(raw_args=None):
    # command line stuff
    parser = ArgumentParser()
    grp = parser.add_mutually_exclusive_group()
    
    # required
    parser.add_argument("--cdip", action="store_true", help="to choose if we are working with cdip data")
    parser.add_argument("nc", nargs="+", type=str, help="netCDF file to process") 
    parser.add_argument("--yaml", nargs="+", type=str, help="YAML file(s) to load")
    
    args = parser.parse_args(raw_args)

    # for each nc filename passed
    for fn in args.nc:
        process(fn, args) # then goes to store calculated data 
        
    # cannot merge a yaml and nc file in one for loop because of different file names
    if args.yaml: 
        for fn in args.yaml:
            __process(fn, args) # optional to store meta from yaml file


if __name__ == "__main__":
    main()