"""
Description: compares directional values with CDIPs values

To compare: 
    python .\tests\test_directional_values.py --graph ".\ncFiles\067.20201225_1200.20201225_1600_output.nc" ".\ncFiles\067.20201225_1200.20201225_1600.nc"


Author: Alex Pitts
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr
import helperfunctions as hf


def process(args: ArgumentParser) -> None:

    # normal
    Wave_xr = xr.open_dataset(args.nc[1], group="Wave")
    Wave_xr_Banded = xr.open_dataset(args.nc[0], group="Wave")
    Wave_xr_Welch = xr.open_dataset(args.nc[0], group="WelchWave")

    
    PSD_xr_Banded = xr.open_dataset(args.nc[0], group="BandedPSD")
    PSD_xr_Welch = xr.open_dataset(args.nc[0], group="WelchPSD")
    freq_space_Banded = PSD_xr_Banded.freq_space.to_numpy()
    freq_space_Welch = PSD_xr_Welch.freq_space.to_numpy()

    Wave = {
        "A1": Wave_xr.A1.to_numpy(),
        "B1": Wave_xr.B1.to_numpy(),
        "A2": Wave_xr.A2.to_numpy(),
        "B2": Wave_xr.B2.to_numpy(),
    }
    bWave = {
        "A1": Wave_xr_Banded.A1.to_numpy(),
        "B1": Wave_xr_Banded.B1.to_numpy(),
        "A2": Wave_xr_Banded.A2.to_numpy(),
        "B2": Wave_xr_Banded.B2.to_numpy(),
    }
    wWave = {
        "A1": Wave_xr_Welch.A1.to_numpy(),
        "B1": Wave_xr_Welch.B1.to_numpy(),
        "A2": Wave_xr_Welch.A2.to_numpy(),
        "B2": Wave_xr_Welch.B2.to_numpy(),
    }
    
    if(len(wWave["A1"])!=len(Wave["A1"])):
        print("error: non matching files")
        exit(1)

    # print comparison
    # blocks = []
    # for i in range(len(Wave["A1"])):
    #     blocks.append("Block" + str(i))

    # for i in Wave:
    #     data1 = pd.DataFrame(Wave[i])
    #     data2 = pd.DataFrame(bWave[i])
    #     data3 = pd.DataFrame(wWave[i])
    #     data = {"CDIP": data1, "Banded": data2, "Welch":data3}

    #     print(i, "\n", data,"\n")
        
    # graph comparison
    if args.graph:          
        for i in range(len(Wave["A1"])):
            figure = [
                ["Directional Spectra", "", "A1", [freq_space_Welch[i], freq_space_Banded[i], freq_space_Banded[i]], [wWave["A1"][i], bWave["A1"][i], Wave["A1"][i]]],
                ["", "", "B1", [freq_space_Welch[i], freq_space_Banded[i], freq_space_Banded[i]], [wWave["B1"][i], bWave["B1"][i], Wave["B1"][i]]],
                ["", "", "A2", [freq_space_Welch[i], freq_space_Banded[i], freq_space_Banded[i]], [wWave["A2"][i], bWave["A2"][i], Wave["A2"][i]]],
                ["", "freq (Hz)", "B2", [freq_space_Welch[i], freq_space_Banded[i], freq_space_Banded[i]], [wWave["B2"][i], bWave["B2"][i], Wave["B2"][i]]]
            ]
            fig, axs = plt.subplots(nrows=4, ncols=1)
            plt.suptitle("Time Block " + str(i))
            legend = ['Banded', 'Welch', 'CDIP']
            
            hf.Plotter(fig, axs, figure, legend)
            
           

        plt.show()


def main(raw_args=None):
    # command line stuff
    parser = ArgumentParser()
    grp = parser.add_mutually_exclusive_group()
    parser.add_argument("--graph", action="store_true",
                        help="to graph data")
    # required
    parser.add_argument("nc", nargs="+", type=str, help="netCDF file to process") 
    
    args = parser.parse_args(raw_args)
    process(args)


if __name__ == "__main__":
    main()