import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from matplotlib.widgets import Slider



def readoutnoise_calculator(biasfiles, **kwargs):

    if "gain" in kwargs.keys():
        gain = kwargs["gain"] # Readout noise is in Electron fluctuations
    else:
        gain = 1 # Readout noise is in ADU fluctuations

    if len(biasfiles[1])==0:
        hdu = fits.open(biasfiles[0])
        bias = hdu[0].data
        readoutnoise = np.sqrt(np.mean((bias - np.mean(bias)) ** 2))

    else:
        biasfile1,biasfile2 = biasfiles
        hdu1 = fits.open(biasfile1)
        bias1 = hdu1[0].data
        hdu2 = fits.open(biasfile2)
        bias2 = hdu2[0].data
        bias = bias1 - bias2
        readoutnoise = np.std(bias) / np.sqrt(2)

    return readoutnoise*gain

