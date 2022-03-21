"""
This file containes all the calculations for calculating PSDs, CSDs, directional values, and non directional values

Authors: 
"""
from asyncio.windows_events import NULL
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
import getdata as gd
from argparse import ArgumentParser

def Rolling_mean(x: np.array, w: np.array) -> np.array:
    """Smoothes the raw acceleration data with a rolling mean. 
    Accepts data to be smoothed and a window width for the moving average. 
    """ 
    df = pd.DataFrame({'a': x.tolist()})
    return df['a'].rolling(window=w, center=True).mean().fillna(0).to_numpy()

def Bias(width: int, window: str = "hann") -> np.array:
    """returns a either a boxcar, or hann window"""
    return np.ones(width) if window == "boxcar" else (
        0.5*(1 - np.cos(2*np.pi*np.arange(width) / (width - 0)))
    )


def wfft(data: np.array, width: int, window: str = "hann") -> list[np.array]:
    """Splits the acceleration data into widows, 
    preforms FFTs on them returning a list of all the windows
    """
    bias = Bias(width, window)

    ffts = []
    for i in range(0, data.size-width+1, width//2):
        w = data[i:i+width]
        ffts.append(np.fft.rfft(w*bias))

    return ffts


def wcalcPSD(
        A_FFT_windows: list[np.array],
        B_FFT_windows: list[np.array],
        fs: float,
        window: str) -> np.array:
    """calculates the PSD of the FFT output preformed with the windowing method.
    After calculateing the PSD of each window, the resulting lists are averaged together"""

    width = A_FFT_windows[0].size
    spectrums = np.complex128(np.zeros(width))
    for i in range(len(A_FFT_windows)):
        A = A_FFT_windows[i]
        B = B_FFT_windows[i]

        spectrum = calcPSD(A, B, fs, window=window)
        spectrums += spectrum
    return spectrums / len(A_FFT_windows)

  
def calcPSD(xFFT: np.array, yFFT: np.array, fs: float, window: str) -> np.array:
    "calculates the PSD on an output of a FFT"
    nfft = xFFT.size
    qOdd = nfft % 2
    n = (nfft - qOdd) * 2  # Number of data points input to FFT
    w = Bias(n, window)  # Get the window used
    wSum = (w * w).sum()
    psd = (xFFT.conjugate() * yFFT) / (fs * wSum)
    if not qOdd:       # Even number of FFT bins
        psd[1:] *= 2   # Real FFT -> double for non-zero freq
    else:              # last point unpaired in Nyquist freq
        psd[1:-1] *= 2  # Real FFT -> double for non-zero freq
    return psd

def merge(data):
    dataset = xr.DataArray(
        np.array(data, dtype=object),
        dims = ("timeblock", "value"),
        coords = {
            "value": (np.arange(0, len(data[0]))).tolist(),
            "timeblock": (np.arange(0, len(data))).tolist(),
        }
    )
    return dataset

def getPSDs(filename: str) -> dict:
    """Calculates the PSD without banding or welch given a filename containing xyz data"""
    # get all data from .nc file
    data = gd.Data(filename)

    frequency = data["Meta"]["frequency"]
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]

    # arrays for storing PSDs and CSDs real and imaginary numbers
    PSDxs = []
    PSDys = []
    PSDzs = []
    PSDxs_imag = []
    PSDys_imag = []
    PSDzs_imag = []
    CSDxys = []
    CSDzxs = []
    CSDzys = []
    CSDxys_imag = []
    CSDzxs_imag = []
    CSDzys_imag = []

    # process all blocks of time
    for i in range(len(timebounds)):
        # grab the upper and lower time window for this timeblock
        time_lower = data["Wave"]["time_lower"].to_numpy()[i]
        time_upper = data["Wave"]["time_upper"].to_numpy()[i]

        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )

        # use "select" to grab the x, y, and z values in question 
        acc = {
            # x is northwards
            "x": data["XYZ"]["x"][select],
            # y is eastwards
            "y": data["XYZ"]["y"][select],
            # z is upwards
            "z": data["XYZ"]["z"][select]
        }

        # calculate the FFTs
        FFT = {
            "x": np.fft.rfft(acc["x"], n=acc["z"].size),  # northwards
            "y": np.fft.rfft(acc["y"], n=acc["z"].size),  # eastwards
            "z": np.fft.rfft(acc["z"], n=acc["z"].size),  # upwards
        }

        # Calculate PSD of data from normal FFT
        PSDx = calcPSD(FFT["x"], FFT["x"], frequency, "boxcar")
        PSDy = calcPSD(FFT["y"], FFT["y"], frequency, "boxcar")
        PSDz = calcPSD(FFT["z"], FFT["z"], frequency, "boxcar")

        # append real and imaginary numbers
        PSDxs.append(PSDx.real)
        PSDys.append(PSDy.real)   
        PSDzs.append(PSDz.real)
        PSDxs_imag.append(PSDx.imag)
        PSDys_imag.append(PSDy.imag)   
        PSDzs_imag.append(PSDz.imag)
        
        # Calculate CSD of data from normal FFT
        CSDxy = calcPSD(FFT["x"], FFT["y"], frequency, "boxcar")
        CSDzx = calcPSD(FFT["z"], FFT["x"], frequency, "boxcar")
        CSDzy = calcPSD(FFT["z"], FFT["y"], frequency, "boxcar")

        # append real and imaginary numbers
        CSDxys.append(CSDxy.real)
        CSDzxs.append(CSDzx.real)
        CSDzys.append(CSDzy.real)
        CSDxys_imag.append(CSDxy.imag)
        CSDzxs_imag.append(CSDzx.imag)
        CSDzys_imag.append(CSDzy.imag)
        
    # merge all PSD lists into a dataset
    PSDss = xr.Dataset({
        "x": merge(PSDxs),
        "y": merge(PSDys),
        "z": merge(PSDzs),
        "x_imag": merge(PSDxs_imag),
        "y_imag": merge(PSDys_imag),
        "z_imag": merge(PSDzs_imag)
    })
    
    # merge all CSD lists into a dataset    
    CSDss = xr.Dataset({
        "xy": merge(CSDxys),
        "zx": merge(CSDzxs),
        "zy": merge(CSDzys),
        "xy_imag": merge(CSDxys_imag),
        "zx_imag": merge(CSDzxs_imag),
        "zy_imag": merge(CSDzys_imag)
    })

    return PSDss, CSDss


def testThis(fn: str) -> dict: 
    data = gd.Data(fn)

    frequency = data["Meta"]["frequency"]
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]

    for i in range(len(timebounds)):
            
        time_lower = data["Wave"]["time_lower"].to_numpy()[i]
        time_upper = data["Wave"]["time_upper"].to_numpy()[i]

        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )

        # use "select" to grab the x, y, and z values in question 
        acc = {
            # x is northwards
            "x": data["XYZ"]["x"][select],
            # y is eastwards
            "y": data["XYZ"]["y"][select],
            # z is upwards
            "z": data["XYZ"]["z"][select]
        }

        wFFT = {
            "x": wfft(acc["x"], 2**8, "hann"),
            "y": wfft(acc["y"], 2**8, "hann"),
            "z": wfft(acc["z"], 2**8, "hann"),
        }

        datas = {
            "wPSDxx": wcalcPSD(wFFT["x"], wFFT["x"], data["Meta"]["frequency"], "hann").real,
            "wPSDyy": wcalcPSD(wFFT["y"], wFFT["y"], data["Meta"]["frequency"], "hann").real,
            "wPSDzz": wcalcPSD(wFFT["z"], wFFT["z"], data["Meta"]["frequency"], "hann").real,
            
            "wCSDxy": wcalcPSD(wFFT["x"], wFFT["y"], data["Meta"]["frequency"], "hann").real,
            "wCSDzx": wcalcPSD(wFFT["z"], wFFT["x"], data["Meta"]["frequency"], "hann").real,
            "wCSDzy": wcalcPSD(wFFT["z"], wFFT["y"], data["Meta"]["frequency"], "hann").real,

            "freq_space": np.fft.rfftfreq(wFFT["z"][0].size*2-1, 1/data["Meta"]["frequency"])
        }
   
    return datas




def getWelchPSDs(filename: str) -> dict:
    """(incomplete) Calculates the PSD with the Welch method given a filename containing xyz data"""
    # get all data from .nc file
    data = gd.Data(filename)

    frequency = data["Meta"]["frequency"]
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]

    wPSDxxs = []
    wPSDyys = []
    wPSDzzs = []
    wPSDxxs_imag = []
    wPSDyys_imag = []
    wPSDzzs_imag = []
    wCSDxys = []
    wCSDzxs = []
    wCSDzys = []
    wCSDxys_imag = []
    wCSDzxs_imag = []
    wCSDzys_imag = []


    # process all blocks of time
    for i in range(len(timebounds)):      
        time_lower = data["Wave"]["time_lower"].to_numpy()[i]
        time_upper = data["Wave"]["time_upper"].to_numpy()[i]

        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )

        # use "select" to grab the x, y, and z values in question 
        acc = {
            # x is northwards
            "x": data["XYZ"]["x"][select],
            # y is eastwards
            "y": data["XYZ"]["y"][select],
            # z is upwards
            "z": data["XYZ"]["z"][select]
        }

        wFFT = {
            "x": wfft(acc["x"], 2**8, "hann"),
            "y": wfft(acc["y"], 2**8, "hann"),
            "z": wfft(acc["z"], 2**8, "hann"),
        }

        wPSDxx = wcalcPSD(wFFT["x"], wFFT["x"], data["Meta"]["frequency"], "hann").real
        wPSDyy = wcalcPSD(wFFT["y"], wFFT["y"], data["Meta"]["frequency"], "hann").real
        wPSDzz = wcalcPSD(wFFT["z"], wFFT["z"], data["Meta"]["frequency"], "hann").real

        # Storing every instance of a real and imag x,y,z value
        wPSDxxs.append(wPSDxx.real)
        wPSDyys.append(wPSDyy.real)
        wPSDzzs.append(wPSDzz.real)
        wPSDxxs_imag.append(wPSDxx.imag)
        wPSDyys_imag.append(wPSDyy.imag)
        wPSDzzs_imag.append(wPSDzz.imag)

        wCSDxy = wcalcPSD(wFFT["x"], wFFT["y"], data["Meta"]["frequency"], "hann").real
        wCSDzx = wcalcPSD(wFFT["z"], wFFT["x"], data["Meta"]["frequency"], "hann").real
        wCSDzy = wcalcPSD(wFFT["z"], wFFT["y"], data["Meta"]["frequency"], "hann").real

        wCSDxys.append(wCSDxy.real)
        wCSDzxs.append(wCSDzx.real)
        wCSDzys.append(wCSDzy.real)
        wCSDxys_imag.append(wCSDxy.imag)
        wCSDzxs_imag.append(wCSDzx.imag)
        wCSDzys_imag.append(wCSDzy.imag)

    wPSDss = xr.Dataset({
        "x": merge(wPSDxxs),
        "y": merge(wPSDyys),
        "z": merge(wPSDzzs),
        "x_imag": merge(wPSDxxs_imag),
        "y_imag": merge(wPSDyys_imag),
        "z_imag": merge(wPSDzzs_imag)
    })
    
    # merge all CSD lists into a dataset    
    wCSDss = xr.Dataset({
        "xy": merge(wCSDxys),
        "zx": merge(wCSDzxs),
        "zy": merge(wCSDzys),
        "xy_imag": merge(wCSDxys_imag),
        "zx_imag": merge(wCSDzxs_imag),
        "zy_imag": merge(wCSDzys_imag)
    })

    return wPSDss, wCSDss

def getBandPSDs(filename: str) -> dict:
    """(incomplete) Calculates the PSD with a banding method given a filename containing xyz data"""
    # get all data from .nc file
    data = gd.Data(filename)

    frequency = data["Meta"]["frequency"]
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]

    freq_upper = data["Freq"]["upper"]
    freq_lower = data["Freq"]["lower"]

    Bandedxxs = []
    Bandedyys = []
    Bandedzzs = []
    Bandedxxs_imag = []
    Bandedyys_imag = []
    Bandedzzs_imag = []
    Bandedxys = []
    Bandedzxs = []
    Bandedzys = []
    Bandedxys_imag = []
    Bandedzxs_imag = []
    Bandedzys_imag = []

    # process all blocks of time
    for i in range(len(timebounds)):
        time_lower = data["Wave"]["time_lower"].to_numpy()[i]
        time_upper = data["Wave"]["time_upper"].to_numpy()[i]
        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )

        # use "select" to grab the x, y, and z values in question 
        acc = {
            # x is northwards
            "x": data["XYZ"]["x"][select],
            # y is eastwards
            "y": data["XYZ"]["y"][select],
            # z is upwards
            "z": data["XYZ"]["z"][select]
        }

        freq_space = np.fft.rfftfreq(acc["z"].size, 1/data["Meta"]["frequency"])
        windowing_method = Bias(len(freq_space), "boxcar")
        freq_select = np.logical_and(
            np.less_equal.outer(freq_lower, freq_space),
            np.greater_equal.outer(freq_upper, freq_space)
        )
        
        count = select.sum(axis=0)

        FFT = {
            "x": np.fft.rfft(acc["x"], n=acc["z"].size),  # northwards
            "y": np.fft.rfft(acc["y"], n=acc["z"].size),  # eastwards
            "z": np.fft.rfft(acc["z"], n=acc["z"].size),  # upwards
        }

        PSD = {
            # imaginary part is zero
            "xx": calcPSD(FFT["x"], FFT["x"], data["Meta"]["frequency"], "boxcar").real,
            "yy": calcPSD(FFT["y"], FFT["y"], data["Meta"]["frequency"], "boxcar").real,
            "zz": calcPSD(FFT["z"], FFT["z"], data["Meta"]["frequency"], "boxcar").real,

            "xy": calcPSD(FFT["x"], FFT["y"], data["Meta"]["frequency"], "boxcar"),
            "zx": calcPSD(FFT["z"], FFT["x"], data["Meta"]["frequency"], "boxcar"),
            "zy": calcPSD(FFT["z"], FFT["y"], data["Meta"]["frequency"], "boxcar"),
        }

        Bandedxx = (freq_select * PSD["xx"] * windowing_method).sum(axis=0) / count
        Bandedyy = (freq_select * PSD["xx"] * windowing_method).sum(axis=0) / count
        Bandedzz = (freq_select * PSD["xx"] * windowing_method).sum(axis=0) / count

        Bandedxxs.append(Bandedxx.real)
        Bandedyys.append(Bandedyy.real)
        Bandedzzs.append(Bandedzz.real)
        Bandedxxs_imag.append(Bandedxx.imag)
        Bandedyys_imag.append(Bandedxx.imag)
        Bandedzzs_imag.append(Bandedxx.imag)

        Bandedxy = (freq_select * PSD["xy"] * windowing_method).sum(axis=0) / count
        Bandedzx = (freq_select * PSD["zx"] * windowing_method).sum(axis=0) / count
        Bandedzy = (freq_select * PSD["zy"] * windowing_method).sum(axis=0) / count

        Bandedxys.append(Bandedxy.real)
        Bandedzxs.append(Bandedzx.real)
        Bandedzys.append(Bandedzy.real)
        Bandedxys_imag.append(Bandedxy.imag)
        Bandedzxs_imag.append(Bandedzx.imag)
        Bandedzys_imag.append(Bandedzy.imag)

    BandedPSDss = xr.Dataset({
        "x": merge(Bandedxxs),
        "y": merge(Bandedyys),
        "z": merge(Bandedzzs),
        "x_imag": merge(Bandedxxs_imag),
        "y_imag": merge(Bandedyys_imag),
        "z_imag": merge(Bandedzzs_imag)
    })
    
    # merge all CSD lists into a dataset    
    BandedCSDss = xr.Dataset({
        "xy": merge(Bandedxys),
        "zx": merge(Bandedzxs),
        "zy": merge(Bandedzys),
        "xy_imag": merge(Bandedxys_imag),
        "zx_imag": merge(Bandedzxs_imag),
        "zy_imag": merge(Bandedzys_imag)
    })

    return BandedPSDss, BandedCSDss
   
def calcWelch(filename: str) -> None:
    data = gd.Data(filename)

    frequency = data["Meta"]["frequency"]
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]

    test = testThis(filename)

    outputs = []
    for i in range(len(timebounds)):
        datasets = []
        outputs.append({})
        time_lower = data["Wave"]["time_lower"].to_numpy()[i]
        time_upper = data["Wave"]["time_upper"].to_numpy()[i]
        # print(time_lower, " - ", time_upper)

        # bit mask so as to select only within the bounds of one lower:upper range pair
        select = np.logical_and(
            time >= time_lower,
            time <= time_upper
        )

        #test["wPSDzz"]
        a0 = a0 = test["wPSDzz"][1:] / \
            np.square(np.square(2 * np.pi * test["freq_space"][1:]))

        m0 = (a0 * test["freq_space"][1]).sum()
        m1 = (a0*test["freq_space"][1:]*test["freq_space"][1]).sum()

        mm1 = (a0/test["freq_space"][1:]*test["freq_space"][1]).sum()
        te = mm1/m0  # mean energy period

        m2 = (a0*np.square(test["freq_space"][1:])
                * test["freq_space"][1]).sum()

        tp = 1/test["freq_space"][1:][a0.argmax()]

        denom = np.sqrt(test["wPSDzz"] * (test["wPSDxx"] + test["wPSDyy"]))
        a1 = test["wCSDzx"].imag / denom
        b1 = test["wCSDzy"].imag / denom
        denom =test["wPSDxx"] + test["wPSDyy"]

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
            "dp_mag": np.degrees(dp+data["Meta"]["declination"]) % 360,
            "a1": a1,
            "b1": b1,
            "a2": (test["wPSDxx"] - test["wPSDyy"]) / denom,
            "b2": -2 * test["wCSDxy"].real / denom,
        }
        datasets.append(xr.Dataset(outputs[i]["welch"]))
        print("Calculated Data using Welch method \"{0}\" window: ".format(
            "hann"))
        for j in outputs[i]["welch"]:
            if np.isscalar(outputs[i]["welch"][j]):
                print(j, "=", outputs[i]["welch"][j])

            
def bandedCalc(fn: str, args: ArgumentParser) -> None:
    return