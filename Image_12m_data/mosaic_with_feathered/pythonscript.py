import numpy
import bottleneck
from astropy.io import fits
import numpy as np
from glob import glob

fitsfiles = glob('*.fits')
datas = np.array([fits.open(w)[0].data for w in fitsfiles])
newarr = bottleneck.nanmean(datas, axis=0)

hdul_out = fits.open(fitsfiles[0])
hdul_out[0].data = newarr
hdul_out.writeto('mosaic.fits',overwrite=True)

