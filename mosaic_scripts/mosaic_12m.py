import os
import numpy as np
from astropy.stats import mad_std
from astropy import wcs
from astropy.io import fits
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from FITS_tools import regrid_cube
from FITS_tools.downsample import *

# Script to regrid the field cubes to a common velocity axis
# Generates both pbcor cubes and noise cubes
# Also generates a mosaic template for the mosaic_cube script

imtype    = '12meter'
sline     = ['12CO']
prefix    = './images/30DOR_F'
basename  = '_12m7mTP_feather+nomodel.'+sline[0]
baseint   = '_12m+7m+TP_CLEAN_smallvelo+nomodel'
outfile   = '30Dor_feather_mosaic_'
dofields  = ['1', '2', '3', '4', '5']
binfactor = 2
pbcut = 0.2
velregrid = False

# Only if input cubes have inconsistent velocity axes
if velregrid:
    vstart = 220.1
    vend   = 285
    delv   = 0.1
    naxis3 = int(round((vend-vstart)/delv + 1))

# Initialize arrays
hdu_slic = [None] * len(dofields)
hdu_cube = [None] * len(dofields)
var_slic = [None] * len(dofields)
var_cube = [None] * len(dofields)

# Main loop over spectral lines to be processed
for line in sline:

    # Loop over fields to produce individual pbcor and variance cubes
    for i, field in enumerate(dofields):

        # --- Regrid the noise-flattened cube and calculate rms
        flatfile = prefix + field + basename + '.image.fits'
        flatimg, hdr = fits.getdata(flatfile, header=True)
        newhdr = hdr.copy()
        if velregrid:
            newhdr['crpix3'] = 1.
            newhdr['cdelt3'] = delv*1000
            newhdr['crval3'] = vstart*1000
            newhdr['naxis3'] = naxis3
            flatimg = regrid_cube(flatimg, hdr, newhdr)
        noisech = np.r_[0:5, flatimg.shape[0]-5:flatimg.shape[0]]
        rms = mad_std(flatimg[noisech,:,:], axis=None, ignore_nan=True)

        # --- Regrid the pb cube and apply pb correction
        if imtype != 'TP':
            pbfile = prefix + field + baseint + '.flux.pbcoverage.fits'
            pbimage, pbhd = fits.getdata(pbfile, header=True)
            pbimage[pbimage<pbcut] = np.nan
            if velregrid:
                pbimage = regrid_cube(pbimage, pbhd, newhdr)
            # Currently we regenerate the pbcor image - option to use existing?
            pbcoimg = flatimg / pbimage
            varimg  = (rms / pbimage)**2
        else:
            pbcoimg = flatimg
            varimg = pbcoimg*0 + rms**2

        # --- Write out the pb corrected and variance cubes
        newhdr['datamin'] = np.nanmin(pbcoimg)
        newhdr['datamax'] = np.nanmax(pbcoimg)
        pbcofile = prefix + field + basename + '.rg_pbcor.fits'
        fits.writeto(pbcofile, 
                     pbcoimg.astype(np.float32), newhdr, overwrite=True)
        newhdr['datamin'] = 0.9*np.nanmin(varimg)
        newhdr['datamax'] = 1.1*np.nanmax(varimg)
        varfile = prefix + field + basename + '.rg_var.fits'
        fits.writeto(varfile, 
                     varimg.astype(np.float32), newhdr, overwrite=True)

        # --- Save the WCS of each field
        if line == sline[0]:
            hdu_slic[i] = fits.PrimaryHDU(data=pbcoimg[0], header=newhdr)
            hdu_slic[i].header['WCSAXES'] = 2
            for key in ['CRVAL3', 'CTYPE3', 'CRPIX3', 'CDELT3', 'CUNIT3',
                        'PC3_1', 'PC3_2', 'PC1_3', 'PC2_3', 'PC3_3']:
                if key in hdu_slic[i].header.keys():
                    del hdu_slic[i].header[key]
            print('\nWCS for field',field)
            print(repr(wcs.WCS(hdu_slic[i])))
        print('\nFinished regridding field',field,'for',line)

    # --- Calculate the mosaic WCS
    if line == sline[0]:
        wcs_out, shape_out = find_optimal_celestial_wcs(hdu_slic)
        hd2d = wcs_out.to_header()
        print('\nOutput header:')
        print(repr(hd2d))
        print('Output shape is',shape_out)
        arr, foot = reproject_and_coadd(hdu_slic, wcs_out, 
                    shape_out=shape_out, reproject_function=reproject_interp)
        # Downsample in RA and DEC if requested
        if binfactor > 1:
            print('Downsampling spatially by a factor of', binfactor)
            arr1 = downsample_axis(arr, binfactor, axis=1)
            hdr1 = downsample_header(hd2d, binfactor, axis=1)
            arr  = downsample_axis(arr1, binfactor, axis=0)
            hd2d = downsample_header(hdr1, binfactor, axis=2)
            print('\nOutput header:')
            print(repr(hd2d))
        fits.writeto('template_'+line+'_'+imtype+'.fits', arr, hd2d,
                     overwrite=True)

    # --- Get geometry of output array
    naxis3 = newhdr['naxis3']
    hd2d = fits.getheader('template_'+sline[0]+'_'+imtype+'.fits')
    naxis2 = hd2d['naxis2']
    naxis1 = hd2d['naxis1']
    print('\nDimensions of output array:',naxis3,naxis2,naxis1)
    #print(repr(hd2d))

    # --- Allocate the output cubes and load the individual fields
    mcube = np.zeros(shape=(naxis3,naxis2,naxis1))
    foot  = np.zeros(shape=(naxis3,naxis2,naxis1))
    wtnse = np.zeros(shape=(naxis3,naxis2,naxis1))
    for i, field in enumerate(dofields):
        hdu_cube[i] = fits.open(prefix+field+basename+'.rg_pbcor.fits')[0]
        var_cube[i] = fits.open(prefix+field+basename+'.rg_var.fits')[0]

    # --- Loop over velocity channels to generate cube slices
    for ich in range(naxis3):
        variance_wts = []
        for i, field in enumerate(dofields):
            hdu_slic[i] = fits.PrimaryHDU(data=hdu_cube[i].data[ich],
                                          header=hdu_cube[i].header)
            hdu_slic[i].header['WCSAXES'] = 2
            for key in ['CRVAL3', 'CTYPE3', 'CRPIX3', 'CDELT3', 'CUNIT3',
                        'PC3_1', 'PC3_2', 'PC1_3', 'PC2_3', 'PC3_3']:
                if key in hdu_slic[i].header.keys():
                    del hdu_slic[i].header[key]
            variance_wts.append(1/var_cube[i].data[ich])
        if (ich % 10 == 1):
            print('Working on channel',ich)
        mcube[ich], foot[ich] = reproject_and_coadd(hdu_slic, hd2d, 
                input_weights=variance_wts, reproject_function=reproject_interp)

    # --- Assemble the mosaic cubes
    hd3d = hdu_cube[0].header
    for key in ['CRPIX1', 'CDELT1', 'CTYPE1', 'CRVAL1', 'CRPIX2', 'CDELT2', 
                'CTYPE2', 'CRVAL2', 'LONPOLE', 'LATPOLE']:
        hd3d[key] = hd2d[key]
    hd3d['datamin'] = np.nanmin(mcube)
    hd3d['datamax'] = np.nanmax(mcube)
    fits.writeto(outfile+line+'_'+imtype+'.pbcor.fits', 
                 mcube.astype(np.float32), hd3d, overwrite=True)
    hd3d['datamin'] = np.nanmin(foot)
    hd3d['datamax'] = np.nanmax(foot)
    fits.writeto(outfile+line+'_'+imtype+'.foot.fits', 
               foot.astype(np.float32), hd3d, overwrite=True)
    foot[foot==0] = np.nan
    wtnse = np.nanmean(1/np.sqrt(foot), axis=0, keepdims=True)
    hd3d['datamin'] = np.min(wtnse[np.isfinite(wtnse)])
    hd3d['datamax'] = np.max(wtnse[np.isfinite(wtnse)])
    fits.writeto(outfile+line+'_'+imtype+'.rms.fits', 
                 wtnse.astype(np.float32), hd3d, overwrite=True)

    # --- Clean up
    for field in dofields:
        os.system("rm -rf "+prefix+field+basename+'.rg_pbcor.fits')
        os.system("rm -rf "+prefix+field+basename+'.rg_var.fits')
