"""
Description: Compares all values with CDIP and displays comparison results in histograms for each value

To compare: 
    python .\tests\test_vals_v2.py "Your File.nc" "CDIPs file.nc"


Author: Alex Pitts
In collaboration with Pat Welch
"""

from distutils.log import error
from ftplib import error_perm
from tkinter import E
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
    # ###################################################
    indices = Wave_xr_Banded.index.to_numpy()
    index = np.nan
    for i in indices:
        if np.isnan(i):
            pass
        else:
            index = int(i)
    if np.isnan(index): 
        print("error with index. All time windows failed")
        exit(1) 

    
    bWave = {
        "Hs": Wave_xr_Banded.Hs.to_numpy(),
        "Ta": Wave_xr_Banded.Ta.to_numpy(),
        "Tp": Wave_xr_Banded.Tp.to_numpy(),
        "Tz": Wave_xr_Banded.Tz.to_numpy(),
        "PeakPSD": Wave_xr_Banded.PeakPSD.to_numpy(),
        "Dp": Wave_xr_Banded.Dp_true.to_numpy(), 
        "A1_sigg": Wave_xr_Banded.A1_sigg.to_numpy(),
        "B1_sigg": Wave_xr_Banded.B1_sigg.to_numpy(),
        "A2_sigg": Wave_xr_Banded.A2_sigg.to_numpy(), 
        "B2_sigg": Wave_xr_Banded.B2_sigg.to_numpy()
    }
    

    wWave = {
        "Hs": Wave_xr_Welch.Hs.to_numpy(),
        "Ta": Wave_xr_Welch.Ta.to_numpy(),
        "Tp": Wave_xr_Welch.Tp.to_numpy(),
        "Tz": Wave_xr_Welch.Tz.to_numpy(),
        "PeakPSD": Wave_xr_Welch.PeakPSD.to_numpy(),
        "Dp": Wave_xr_Welch.Dp_true.to_numpy(),
        "A1_sigg": Wave_xr_Welch.A1_sigg.to_numpy(),
        "B1_sigg": Wave_xr_Welch.B1_sigg.to_numpy(),
        "A2_sigg": Wave_xr_Welch.A2_sigg.to_numpy(), 
        "B2_sigg": Wave_xr_Welch.B2_sigg.to_numpy()
    }
   

    A1 = Wave_xr.A1.to_numpy()
    B1 = Wave_xr.B1.to_numpy()
    A2 = Wave_xr.A2.to_numpy()
    B2 = Wave_xr.B2.to_numpy()

    A1_sigg = []
    B1_sigg = []
    A2_sigg = []
    B2_sigg = []

    for i in range(len(A1)):
        A1_sigg.append(A1[i][index])
        B1_sigg.append(B1[i][index])
        A2_sigg.append(A2[i][index])
        B2_sigg.append(B2[i][index])
   

    Wave = {
        "Hs": Wave_xr.Hs.to_numpy(),
        "Ta": Wave_xr.Ta.to_numpy(),
        "Tp": Wave_xr.Tp.to_numpy(),
        "Tz": Wave_xr.Tz.to_numpy(),
        "PeakPSD": Wave_xr.PeakPSD.to_numpy(),
        "Dp": Wave_xr.Dp.to_numpy(),
        "A1_sigg": np.asarray(A1_sigg),
        "B1_sigg": np.asarray(B1_sigg),
        "A2_sigg": np.asarray(A2_sigg), 
        "B2_sigg": np.asarray(B2_sigg)
    }
    

    # ###################################################

    Bdiff_numerical = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [], "A1_sigg": [], "B1_sigg": [], "A2_sigg": [], "B2_sigg": []}
    Bdiff_proportional = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [], "A1_sigg": [], "B1_sigg": [], "A2_sigg": [], "B2_sigg": [] }

    Wdiff_numerical = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [], "A1_sigg": [], "B1_sigg": [], "A2_sigg": [], "B2_sigg": [] }
    Wdiff_proportional = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [], "A1_sigg": [], "B1_sigg": [], "A2_sigg": [], "B2_sigg": [] }

    for j in wWave: 
        # building X-axis
        error_count = 0
        for i in range(len(wWave[j])):
            if np.isnan(wWave[j][i]):
                error_count += 1
                
            else:

                Bdiff_proportional[j].append((bWave[j][i] - Wave[j][i]) / bWave[j][i])
                Bdiff_numerical[j].append(bWave[j][i] - Wave[j][i])

                Wdiff_proportional[j].append((wWave[j][i] - Wave[j][i]) / wWave[j][i])
                Wdiff_numerical[j].append(wWave[j][i] - Wave[j][i])
    

    # graphing results
    num_bins = args.bins
    for i in Wave:
        (fig, axes) = plt.subplots(nrows=3, ncols=4, squeeze=True, sharex=True)
        plt.suptitle(f"{i}: n_windows used = {len(Bdiff_proportional[i])} \nwindows ommetted due to errors = {error_count}")

        # setting up lables
        axes[0][0].set_ylabel("Cnts")
        axes[1][0].set_ylabel("Cnts / max(Cnts)")
        axes[2][0].set_ylabel("Cnts / sum(Cnts)")

        axes[0][0].set_title("banded proportional")
        axes[0][1].set_title("banded numerical")
        axes[0][2].set_title("Welch proportional")
        axes[0][3].set_title("Welch numerical")
 
        # banded proportional
        (cnts, bins) = np.histogram(Bdiff_proportional[i], bins=num_bins, 
            range=None if args.lower is None or args.upper is None else (args.lower, args.upper),
            )
        axes[0][0].stairs(cnts, edges=bins, fill=True)
        axes[1][0].stairs(cnts / cnts.max(), edges=bins, fill=True)
        axes[2][0].stairs(cnts / cnts.sum(), edges=bins, fill=True)

        # banded numerical
        (cnts, bins) = np.histogram(Bdiff_numerical[i], bins=num_bins, 
            range=None if args.lower is None or args.upper is None else (args.lower, args.upper),
            )
        axes[0][1].stairs(cnts, edges=bins, fill=True)
        axes[1][1].stairs(cnts / cnts.max(), edges=bins, fill=True)
        axes[2][1].stairs(cnts / cnts.sum(), edges=bins, fill=True)
        
        # Welch proportional
        (cnts, bins) = np.histogram(Wdiff_proportional[i], bins=num_bins, 
            range=None if args.lower is None or args.upper is None else (args.lower, args.upper),
            )
        axes[0][2].stairs(cnts, edges=bins, fill=True)
        axes[1][2].stairs(cnts / cnts.max(), edges=bins, fill=True)
        axes[2][2].stairs(cnts / cnts.sum(), edges=bins, fill=True)
        
        # Welch numerical
        (cnts, bins) = np.histogram(Wdiff_numerical[i], bins=num_bins, 
            range=None if args.lower is None or args.upper is None else (args.lower, args.upper),
            )
        axes[0][3].stairs(cnts, edges=bins, fill=True)
        axes[1][3].stairs(cnts / cnts.max(), edges=bins, fill=True)
        axes[2][3].stairs(cnts / cnts.sum(), edges=bins, fill=True)
        
    print(error_count, " errors")
        
    plt.tight_layout()
    plt.show()

    
        
            
   
    
   

    
    
    


def main(raw_args=None):
    # command line stuff
    parser = ArgumentParser()
    grp = parser.add_mutually_exclusive_group()
    parser.add_argument("--graph", action="store_true",
                        help="to graph data")
    # required
    parser.add_argument("nc", nargs="+", type=str, help="netCDF file to process") 
    parser.add_argument("--bins", type=int, default=10, help="Number of histogram bins")
    parser.add_argument("--lower", type=float, help="Left hand bin lower edge")
    parser.add_argument("--upper", type=float, help="Right hand bin upper edge")
    
    args = parser.parse_args(raw_args)
    process(args)


if __name__ == "__main__":
    main()