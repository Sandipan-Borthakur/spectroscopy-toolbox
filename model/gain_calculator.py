from astropy.io import fits
import numpy as np

def gain_calculator(flatfiles,biasfiles,**kwargs):
    flatfile1,flatfile2 = flatfiles
    biasfile1,biasfile2 = biasfiles
    hdu1 = fits.open(flatfile1)
    flat1 = hdu1[0].data
    hdu2 = fits.open(flatfile2)
    flat2 = hdu2[0].data

    hdu1 = fits.open(biasfile1)
    bias1 = hdu1[0].data
    hdu2 = fits.open(biasfile2)
    bias2 = hdu2[0].data

    if "cropinds" in kwargs.keys():
        minx, maxx, miny, maxy = kwargs["cropinds"]
        flat1,flat2 = flat1[minx:maxx,miny:maxy],flat2[minx:maxx,miny:maxy]
        bias1,bias2 = bias1[minx:maxx,miny:maxy],bias2[minx:maxx,miny:maxy]

    flat_diff = flat1 - flat2
    bias_diff = bias1 - bias2
    gain = (flat1 + flat2 - (bias1+bias2))/(np.std(flat_diff)**2 + np.std(bias_diff)**2)

    return gain