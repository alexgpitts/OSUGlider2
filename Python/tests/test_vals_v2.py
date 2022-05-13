#! /usr/bin/env python3
#
# A sample histogram for Alex Pitts
#
# May-2022, Pat Welch, pat@mousebrains.com

from distutils.log import error
from ftplib import error_perm
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

    Bdiff_numerical = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [] }
    Bdiff_proportional = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [] }

    Wdiff_numerical = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [] }
    Wdiff_proportional = { "Hs": [], "Ta": [], "Tp": [], "Tz": [], "PeakPSD": [], "Dp": [] }

    for j in Wave: 
        # building X-axis
        error_count = 0
        for i in range(len(wWave[j])):
            if np.isnan(wWave[j][i]):
                error_count += 1
                
            else:
                # hs
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