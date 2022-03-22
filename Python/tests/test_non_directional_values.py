"""
Description: compares non-directional values with CDIPs values

To compare: 
    python .\tests\test_non_directional_values.py "Your File.nc" "CDIPs file.nc"

To compare with graphs: 
    python .\tests\test_non_directional_values.py --graph "Your File.nc" "CDIPs file.nc"

Author: Alex Pitts
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr


def process(args: ArgumentParser) -> None:

    # normal
    Wave_xr = xr.open_dataset(args.nc[1], group="Wave")
    Wave_xr_Banded = xr.open_dataset(args.nc[0], group="Wave")
    Wave_xr_Welch = xr.open_dataset(args.nc[0], group="WelchWave")

    Wave = {
        "Hs": Wave_xr.Hs.to_numpy(),
        "Ta": Wave_xr.Ta.to_numpy(),
        "Tp": Wave_xr.Tp.to_numpy(),
        "Tz": Wave_xr.Tz.to_numpy(),
        "PeakPSD": Wave_xr.PeakPSD.to_numpy(),
        "Dp": Wave_xr.Dp.to_numpy()
    }
    bWave = {
        "Hs": Wave_xr_Banded.Hs.to_numpy(),
        "Ta": Wave_xr_Banded.Ta.to_numpy(),
        "Tp": Wave_xr_Banded.Tp.to_numpy(),
        "Tz": Wave_xr_Banded.Tz.to_numpy(),
        "PeakPSD": Wave_xr_Banded.PeakPSD.to_numpy(),
        "Dp": Wave_xr_Banded.Dp_true.to_numpy()
    }
    wWave = {
        "Hs": Wave_xr_Welch.Hs.to_numpy(),
        "Ta": Wave_xr_Welch.Ta.to_numpy(),
        "Tp": Wave_xr_Welch.Tp.to_numpy(),
        "Tz": Wave_xr_Welch.Tz.to_numpy(),
        "PeakPSD": Wave_xr_Welch.PeakPSD.to_numpy(),
        "Dp": Wave_xr_Welch.Dp_true.to_numpy()
    }
    
    if(len(wWave["Hs"])!=len(Wave["Hs"])):
        print("error: non matching files")
        exit(1)

    # print comparison
    blocks = []
    for i in range(len(Wave["Hs"])):
        blocks.append("Block" + str(i))

    for i in Wave:
        data = pd.DataFrame({"CDIP": Wave[i], "Banded": bWave[i], "Welch": wWave[i]}, index=blocks)
        print(i, "\n", data,"\n")
    
    # graph comparison
    if args.graph:          
        for i in Wave:
            # set width of bar
            barWidth = 0.25
            fig = plt.subplots(figsize =(12, 8))
            
            # set height of bar
            data = [Wave[i], bWave[i], wWave[i]]
            
            # Set position of bar on X axis
            br1 = np.arange(len(Wave[i]))
            br2 = [x + barWidth for x in br1]
            br3 = [x + barWidth for x in br2]
            
            # Make the plot
            plt.bar(br1, data[0], color ='r', width = barWidth,
                    edgecolor ='grey', label ='CDIP')
            plt.bar(br2, data[1], color ='g', width = barWidth,
                    edgecolor ='grey', label ='Banded')
            plt.bar(br3, data[2], color ='b', width = barWidth,
                    edgecolor ='grey', label ='Welch')
            
            # Adding Xticks
            plt.xlabel(i, fontweight ='bold', fontsize = 15)
            plt.ylabel('', fontweight ='bold', fontsize = 15)
            plt.xticks([r + barWidth for r in range(len(Wave[i]))],
                    np.arange(len(Wave[i])))
            
            plt.legend()

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