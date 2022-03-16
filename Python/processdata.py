from asyncio.windows_events import NULL
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
import getdata as gd

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
        dims = ("y", "x"),
        coords = {
            "x": (np.arange(1, len(data[0])+1)).tolist(),
            "y": (np.arange(1, len(data)+1)).tolist(),
        }
    )
    return dataset




def getPSDs(filename: str) -> dict:
    # get all data from .nc file
    data = gd.Data(filename)

    frequency = data["Meta"]["frequency"]
    time = data["XYZ"]["t"]
    timebounds = data["Wave"]["Timebounds"]
    PSDxs = []
    PSDys = []
    PSDzs = []
    CSDxys = []
    CSDzxs = []
    CSDzys = []
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

        PSDxs.append(PSDx)
        PSDys.append(PSDy)    
        PSDzs.append(PSDz)

        # PSDxs.append(xr.Dataset({'real': PSDx.real, 'imag': PSDx.imag}))
        # PSDys.append(xr.Dataset({'real': PSDy.real, 'imag': PSDy.imag}))    
        # PSDzs.append(xr.Dataset({'real': PSDz.real, 'imag': PSDz.imag}))  
        
        CSDxy = calcPSD(FFT["x"], FFT["y"], frequency, "boxcar")
        CSDzx = calcPSD(FFT["z"], FFT["x"], frequency, "boxcar")
        CSDzy = calcPSD(FFT["z"], FFT["y"], frequency, "boxcar")

        CSDxys.append(CSDxy)
        CSDzxs.append(CSDzx)
        CSDzys.append(CSDzy)

        # CSDxys.append(xr.Dataset({'real': CSDxy.real, 'imag': CSDxy.imag}))
        # CSDzxs.append(xr.Dataset({'real': CSDzx.real, 'imag': CSDzx.imag}))
        # CSDzys.append(xr.Dataset({'real': CSDzy.real, 'imag': CSDzy.imag}))
    
    
    PSDs = {
        "x": merge(PSDxs),
        "y": merge(PSDys),
        "z": merge(PSDzs)
    }

    CSDs = {
        "xy": merge(CSDxys),
        "zx": merge(CSDzxs),
        "zy": merge(CSDzys)
    }

    
    return CSDs, PSDs