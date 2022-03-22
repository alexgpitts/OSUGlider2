"""
Description:



Author: Alex Pitts
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr


def process(args: ArgumentParser) -> None:
    pass
    


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