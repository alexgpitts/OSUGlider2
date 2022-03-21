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

def merge(data):
    """takes a list and merges the elements into a dataset"""
    dataset = xr.DataArray(
        np.array(data, dtype=object),
        dims = ("timeblock", "value"),
        coords = {
            "value": (np.arange(0, len(data[0]))).tolist(),
            "timeblock": (np.arange(0, len(data))).tolist(),
        }
    )
    return dataset 

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

def errorCalc() -> dict:
    output = {
        "Hs": np.NAN,
        "Ta": np.NAN,  # average period
        "Tp": np.NAN,  # peak wave period
        "wave_energy_ratio": np.NAN,
        "Tz": np.NAN,
        "PeakPSD": np.NAN,
        "Te": np.NAN,
        "Dp": np.NAN,
        "Dp_mag": np.NAN,
        "Dp_true": np.NAN,
        "A1": [np.NAN],
        "B1": [np.NAN],
        "A2": [np.NAN],
        "B2": [np.NAN]
    }
    return output

def welchCalc(PSD: dict, Data: dict) -> dict:
    
    # non directional
    a0 = PSD["zz"][1:] / np.square(np.square(2 * np.pi * PSD["freq_space"][1:]))
    m0 = (a0 * PSD["freq_space"][1]).sum()
    m1 = (a0*PSD["freq_space"][1:]*PSD["freq_space"][1]).sum()
    mm1 = (a0/PSD["freq_space"][1:]*PSD["freq_space"][1]).sum()
    te = mm1/m0  # mean energy period
    m2 = (a0*np.square(PSD["freq_space"][1:]) * PSD["freq_space"][1]).sum()
    tp = 1/PSD["freq_space"][1:][a0.argmax()]

    # directional
    denom = np.sqrt(PSD["zz"] * (PSD["xx"] + PSD["yy"]))
    a1 = PSD["zx"].imag / denom
    b1 = -PSD["zy"].imag / denom
    denom = PSD["xx"] + PSD["yy"]
    dp = np.arctan2(b1[a0.argmax()], a1[a0.argmax()])  # radians

    output = {
        "Hs": 4 * np.sqrt(m0),
        "Ta": m0/m1,  # average period
        "Tp": tp,  # peak wave period
        "wave_energy_ratio": te/tp,
        "Tz": np.sqrt(m0/m2),
        "PeakPSD": a0.max(),
        "Te": te,  # mean energy period
        "Dp": np.arctan2(b1[a0.argmax()], a1[a0.argmax()]),
        "Dp_mag": np.degrees(dp+Data["Meta"]["declination"]) % 360,
        "Dp_true": np.degrees(dp) % 360,
        "A1": a1,
        "B1": b1,
        "A2": (PSD["xx"] - PSD["yy"]) / denom,
        "B2": -2 * PSD["xy"].real / denom,
    }

    return output



def bandedCalc(PSD: dict, Data: dict):

     # non directional
    a0 = PSD["zz"] / np.square(np.square(2 * np.pi * PSD["freq_space"]))
    m0 = (a0 * Data["Wave"]["Bandwidth"]).sum()
    m1 = (a0*PSD["freq_space"]*Data["Wave"]["Bandwidth"]).sum()
    mm1 = (a0/PSD["freq_space"]*Data["Wave"]["Bandwidth"]).sum()
    te = mm1/m0  # mean energy period
    m2 = (a0 * np.square(PSD["freq_space"]) * Data["Wave"]["Bandwidth"]).sum()
    tp = 1/PSD["freq_space"][a0.argmax()]
    
    # directional
    denom = np.sqrt(PSD["zz"] * (PSD["xx"] + PSD["yy"]))
    a1 = PSD["zx"].imag / denom
    b1 = -PSD["zy"].imag / denom
    denom = PSD["xx"] + PSD["yy"]

    dp = np.arctan2(b1[a0.argmax()], a1[a0.argmax()])  # radians

    output = {
        "Hs": 4 * np.sqrt(m0),
        "Ta": m0/m1,
        "Tp": tp,  # peak wave period
        "wave_energy_ratio": te/tp,
        "Tz": np.sqrt(m0/m2),
        "PeakPSD": a0.max(),
        "Te": te,
        "Dp": np.arctan2(b1[a0.argmax()], a1[a0.argmax()]),
        "Dp_mag": np.degrees(dp+Data["Meta"]["declination"]) % 360,
        "Dp_true": np.degrees(dp) % 360,
        "A1": a1,
        "B1": b1,
        "A2": (PSD["xx"] - PSD["yy"]) / denom,
        "B2": -2 * PSD["xy"].real / denom,
    }
    
    return output