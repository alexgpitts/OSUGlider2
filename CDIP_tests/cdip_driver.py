# import Data from CDIP
from helper_functions import (
    Data, Rolling_mean, calcPSD, wcalcPSD, wfft, Plotter, Bias)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import xarray as xr
import netCDF4 as nc4
import os


def process(fn: str, args: ArgumentParser) -> None:

    # output_name = args.output[0]
    folder_name = os.path.split(args.nc[0])[0]
    output_name = os.path.splitext(os.path.basename(args.nc[0]))[
        0] + "_output.nc"

    # windowing perameters
    if(args.hann):
        window_type = "hann"
    else:
        window_type = "boxcar"

    # making sure command line was passed for a calculation method
    if(args.welch or args.banding):
        pass
    else:
        print("please enter a calculation method (--welch or --banding)")
        return

    # call master data function to extract all our data from the .nc file
    data = Data(args.nc[0])
    outputs = []
    datasets = []

    # time bounds holds an array starting and ending times for each analisis block
    time_bounds = data["time-bounds"]
    time = data["time"]

    # lists of upper, lower, and midpoint frequencies for banding
    freq_bounds = data["freq"]["bounds"]
    freq_lower = freq_bounds["lower"]
    freq_upper = freq_bounds["upper"]
    freq_midpoints = freq_bounds["joint"].mean(axis=1)

    # loop runs through every analisis block, displaying output calculations
    # and graphs at the end of each loop run
    for i in range(len(time_bounds["lower"])):
        print("Processing Block {0}".format(i))

        time_lower = time_bounds["lower"][i]
        time_upper = time_bounds["upper"][i]
        # print(time_lower, " - ", time_upper)

        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )

        # use select to filter to select the acc data corresponding to the current block
        # time = data["time"][select]

        averaging_window = 2
        acc = {
            # x is northwards
            "x": Rolling_mean(data["acc"]["x"][select], averaging_window),
            # y is eastwards
            "y": Rolling_mean(data["acc"]["y"][select], averaging_window),
            # z is upwards
            "z": Rolling_mean(data["acc"]["z"][select], averaging_window)

        }
        if(any(i > 500 for i in acc["x"]) or any(i > 500 for i in acc["y"]) or any(i > 500 for i in acc["z"])):
            print("bad data containing extremely large values\n\n")
            outputs.append({})
            outputs[i]["error"] = {
                "Hs": "error",
                "Ta": "error",  # average period
                "Tp": "error",  # peak wave period
                "wave_energy_ratio": "error",
                "Tz": "error",
                "Dp": "error",
                "PeakPSD": "error",
                "te": "error",  # mean energy period
                "dp_true": "error",
                "dp_mag": "error",
                "a1": "error",
                "b1": "error",
                "a2": "error",
                "b2": "error",
            }
            continue

        # preform FFT on block
        FFT = {
            "x": np.fft.rfft(acc["x"], n=acc["z"].size),  # northwards
            "y": np.fft.rfft(acc["y"], n=acc["z"].size),  # eastwards
            "z": np.fft.rfft(acc["z"], n=acc["z"].size),  # upwards
        }

        # preform FFT on block using welch mothod
        wFFT = {
            "x": wfft(acc["x"], 2**8, window_type),
            "y": wfft(acc["y"], 2**8, window_type),
            "z": wfft(acc["z"], 2**8, window_type),
        }

        # Calculate PSD of data from normal FFT
        PSD = {
            # imaginary part is zero
            "xx": calcPSD(FFT["x"], FFT["x"], data["frequency"], "boxcar").real,
            "yy": calcPSD(FFT["y"], FFT["y"], data["frequency"], "boxcar").real,
            "zz": calcPSD(FFT["z"], FFT["z"], data["frequency"], "boxcar").real,

            "xy": calcPSD(FFT["x"], FFT["y"], data["frequency"], "boxcar"),
            "zx": calcPSD(FFT["z"], FFT["x"], data["frequency"], "boxcar"),
            "zy": calcPSD(FFT["z"], FFT["y"], data["frequency"], "boxcar"),
        }

        # calculate PSD on output from welch method FFT
        wPSD = {
            "xx": wcalcPSD(wFFT["x"], wFFT["x"], data["frequency"], window_type).real,
            "yy": wcalcPSD(wFFT["y"], wFFT["y"], data["frequency"], window_type).real,
            "zz": wcalcPSD(wFFT["z"], wFFT["z"], data["frequency"], window_type).real,

            "xy": wcalcPSD(wFFT["x"], wFFT["y"], data["frequency"], window_type),
            "zx": wcalcPSD(wFFT["z"], wFFT["x"], data["frequency"], window_type),
            "zy": wcalcPSD(wFFT["z"], wFFT["y"], data["frequency"], window_type),

            "freq_space": np.fft.rfftfreq(wFFT["z"][0].size*2-1, 1/data["frequency"])
        }

        # frequency space for plotting FFT
        freq_space = np.fft.rfftfreq(acc["z"].size, 1/data["frequency"])

        # bit mask so as to select only within the bounds of one lower:upper range pair
        freq_select = np.logical_and(
            np.less_equal.outer(freq_lower, freq_space),
            np.greater_equal.outer(freq_upper, freq_space)
        )

        count = freq_select.sum(axis=1)

        windowing_method = Bias(len(freq_space), window_type)

        Band = {
            "xx": (freq_select * PSD["xx"] * windowing_method).sum(axis=1) / count,
            "yy": (freq_select * PSD["yy"] * windowing_method).sum(axis=1) / count,
            "zz": (freq_select * PSD["zz"] * windowing_method).sum(axis=1) / count,

            "xy": (freq_select * PSD["xy"] * windowing_method).sum(axis=1) / count,
            "zx": (freq_select * PSD["zx"] * windowing_method).sum(axis=1) / count,

            "zy": (freq_select * PSD["zy"] * windowing_method).sum(axis=1) / count,
        }

        ##########################################
        # Calculations
        ####################### ##################

        outputs.append({})
        ##########################################
        # Calculations using the welch method
        ##########################################

        def welch(run: bool):
            if run == False:
                return

            a0 = wPSD["zz"][1:] / \
                np.square(np.square(2 * np.pi * wPSD["freq_space"][1:]))

            m0 = (a0 * wPSD["freq_space"][1]).sum()
            m1 = (a0*wPSD["freq_space"][1:]*wPSD["freq_space"][1]).sum()

            mm1 = (a0/wPSD["freq_space"][1:]*wPSD["freq_space"][1]).sum()
            te = mm1/m0  # mean energy period

            m2 = (a0*np.square(wPSD["freq_space"][1:])
                  * wPSD["freq_space"][1]).sum()

            tp = 1/wPSD["freq_space"][1:][a0.argmax()]

            denom = np.sqrt(wPSD["zz"] * (wPSD["xx"] + wPSD["yy"]))
            a1 = wPSD["zx"].imag / denom
            b1 = -wPSD["zy"].imag / denom
            denom = wPSD["xx"] + wPSD["yy"]

            dp = np.arctan2(b1[a0.argmax()], a1[a0.argmax()])  # radians

            outputs[i]["welch"] = {
                "Hs": 4 * np.sqrt(m0),
                "Ta": m0/m1,  # average period
                "Tp": tp,  # peak wave period
                "wave_energy_ratio": te/tp,
                "Tz": np.sqrt(m0/m2),
                "Dp": np.arctan2(b1[a0.argmax()], a1[a0.argmax()]),
                "PeakPSD": a0.max(),
                "te": te,  # mean energy period
                "dp_true": np.degrees(dp) % 360,
                "dp_mag": np.degrees(dp+data["declination"]) % 360,
                "a1": a1,
                "b1": b1,
                "a2": (wPSD["xx"] - wPSD["yy"]) / denom,
                "b2": -2 * wPSD["xy"].real / denom,
            }
            datasets.append(xr.Dataset(outputs[i]["welch"]))
            print("Calculated Data using Welch method \"{0}\" window: ".format(
                window_type))
            for j in outputs[i]["welch"]:
                if np.isscalar(outputs[i]["welch"][j]):
                    print(j, "=", outputs[i]["welch"][j])

        welch(args.welch)

        ##########################################
        # Calculations using the banded method
        ##########################################

        def banded(run: bool):
            # Preform Baniding on the PSD. Averages the data withen each bin.\

            if run == False:
                return
            a0 = Band["zz"] / np.square(np.square(2 * np.pi * freq_midpoints))

            m0 = (a0 * data["freq"]["bandwidth"]).sum()
            m1 = (a0*freq_midpoints*data["freq"]["bandwidth"]).sum()
            mm1 = (a0/freq_midpoints*data["freq"]["bandwidth"]).sum()
            te = mm1/m0  # mean energy period

            m2 = (a0*np.square(freq_midpoints)
                  * data["freq"]["bandwidth"]).sum()

            tp = 1/freq_midpoints[a0.argmax()]

            denom = np.sqrt(Band["zz"] * (Band["xx"] + Band["yy"]))
            a1 = Band["zx"].imag / denom
            b1 = -Band["zy"].imag / denom
            denom = Band["xx"] + Band["yy"]

            dp = np.arctan2(
                b1[a0.argmax()], a1[a0.argmax()])  # radians

            outputs[i]["banded"] = {
                "Hs": 4 * np.sqrt(m0),
                "Ta": m0/m1,
                "Tp": tp,  # peak wave period
                "wave_energy_ratio": te/tp,
                "Tz": np.sqrt(m0/m2),
                "Dp": np.arctan2(b1[a0.argmax()], a1[a0.argmax()]),
                "PeakPSD": a0.max(),
                "te": te,
                "dp_true": np.degrees(dp) % 360,
                "dp_mag": np.degrees(dp+data["declination"]) % 360,
                "a1": a1,
                "b1": b1,
                "a2": (Band["xx"] - Band["yy"]) / denom,
                "b2": -2 * Band["xy"].real / denom,
            }
            datasets.append(xr.Dataset(outputs[i]["welch"]))
            print("Calculated Data using Banding and \"{0}\" window: ".format(
                window_type))
            for j in outputs[i]["banded"]:
                if np.isscalar(outputs[i]["banded"][j]):
                    print(j, "=", outputs[i]["banded"][j])

        banded(args.banding)

        print("\n\nCDIP Data: ")
        for j in data["wave"]:
            if np.isscalar(data["wave"][j][i]):
                print(j, "=", data["wave"][j][i])

        ##########################################
        # plotting
        ##########################################

        if(args.raw):
            figure = [
                ["X Acc", "", "m/s^2", data["time"][select], acc["x"]],
                ["Y Acc", "", "m/s^2", data["time"][select], acc["y"]],
                ["Z Acc", "Time (s)", "m/s^2", data["time"][select], acc["z"]]
            ]
            fig, axs = plt.subplots(nrows=3, ncols=1)
            Plotter(fig, axs, figure)

        if(args.norm and args.graph):
            figure = [
                ["X PSD", "", "", freq_space, PSD["xx"]],
                ["Y PSD", "", "", freq_space, PSD["yy"]],
                ["Z PSD", "freq (Hz)", "", freq_space, PSD["zz"]]
            ]
            fig, axs = plt.subplots(nrows=3, ncols=1)
            Plotter(fig, axs, figure)

        if(args.banding and args.graph):
            figure = [
                ["X Banded PSD", "", "", freq_midpoints, Band["xx"]],
                ["Y Banded PSD", "", "", freq_midpoints, Band["yy"]],
                ["Z Banded PSD", "freq (Hz)", "", freq_midpoints, Band["zz"]]
            ]
            fig, axs = plt.subplots(nrows=3, ncols=1)
            Plotter(fig, axs, figure)

        if(args.welch and args.graph):
            figure = [
                ["X Windowed PSD", "", "", wPSD["freq_space"], wPSD["xx"]],
                ["Y Windowed PSD", "", "", wPSD["freq_space"], wPSD["yy"]],
                ["Z Windowed PSD", "freq (Hz)", "",
                 wPSD["freq_space"], wPSD["zz"]]
            ]
            fig, axs = plt.subplots(nrows=3, ncols=1)
            Plotter(fig, axs, figure)

        if(args.ds):
            if args.banding:
                figure = [
                    ["Directional Spectra with banding", "", "A1",
                        freq_midpoints, [outputs[i]["banded"]["a1"], data["wave"]["a1"][i]]],
                    ["", "", "B1", freq_midpoints, [
                        outputs[i]["banded"]["b1"], data["wave"]["b1"][i]]],
                    ["", "", "A2", freq_midpoints, [
                        outputs[i]["banded"]["a2"], data["wave"]["a2"][i]]],
                    ["", "freq (Hz)", "B2", freq_midpoints, [
                        outputs[i]["banded"]["b2"], data["wave"]["b2"][i]]]
                ]
                fig, axs = plt.subplots(nrows=4, ncols=1)
                Plotter(fig, axs, figure)
            elif(args.welch):
                figure = [
                    ["Directional Spectra with welch", "", "A1",
                        [wPSD["freq_space"], freq_midpoints], [outputs[i]["welch"]["a1"], data["wave"]["a1"][i]]],
                    ["", "", "B1", [wPSD["freq_space"], freq_midpoints], [
                        outputs[i]["welch"]["b1"], data["wave"]["b1"][i]]],
                    ["", "", "A2", [wPSD["freq_space"], freq_midpoints], [
                        outputs[i]["welch"]["a2"], data["wave"]["a2"][i]]],
                    ["", "freq (Hz)", "B2", [wPSD["freq_space"], freq_midpoints], [
                        outputs[i]["welch"]["b2"], data["wave"]["b2"][i]]]
                ]
                fig, axs = plt.subplots(nrows=4, ncols=1)
                Plotter(fig, axs, figure)

        if(args.welch or args.banding or args.raw or args.ds or args.norm):
            plt.show()

        print("\n--------------------------\n")

        # exit(0)  # comment out if you want to proccess all the blocks of data
    # print(outputs)
    if(args.output):
        output_dir = os.path.join(args.output, output_name)
    else:
        output_dir = os.path.join(folder_name, output_name)

    nc4.Dataset(output_dir, 'w', format='NETCDF4')
    for i in datasets:
        i.to_netcdf(output_dir, mode="a", group="wave")


# parser.add_argument(
#     "-o", "--output", help="Directs the output to a name of your choice")
# args = parser.parse_args()
# with open(args.output, 'w') as output_file:
#     output_file.write("%s\n" % item)


def main(raw_args=None):

    #######################################
    # command line stuff
    #######################################
    parser = ArgumentParser()
    grp = parser.add_mutually_exclusive_group()

    # calculation options
    grp.add_argument("--welch", action="store_true", help="Welch Method")
    grp.add_argument("--banding", action="store_true", help="Banding Method")

    # optional args
    parser.add_argument("--hann", action="store_true",
                        help="to choose hann windowing method")
    parser.add_argument("--boxcar", action="store_true",
                        help="to choose boxcar windowing method")
    parser.add_argument("--norm", action="store_true", help="Normal FFT PSD")
    parser.add_argument("--ds", action="store_true",
                        help="Directional Spectrum coefficients")
    parser.add_argument("--graph", action="store_true", help="Turns graphs on")
    parser.add_argument("--raw", action="store_true",
                        help="Raw acceleration data")

    parser.add_argument("-o", "--output", type=str,
                        help="netCDF file write out too")  # typed after commands

    # required
    parser.add_argument("nc", nargs="+", type=str,
                        help="netCDF file to process")  # typed after commands

    args = parser.parse_args(raw_args)

    for fn in args.nc:
        process(fn, args)


if __name__ == "__main__":
    main()
