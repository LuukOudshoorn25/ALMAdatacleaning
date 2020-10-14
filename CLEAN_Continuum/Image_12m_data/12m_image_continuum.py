# step2: Imaging

######################################################
# Use original MS sets


# Field 2: X2180
vis12m_F1 = "../../2019.1.00843.S/science_goal.uid___A001_X1465_X2180/group.uid___A001_X1465_X2181/member.uid___A001_X1465_X2182/calibrated/working/uid___A002_Xe3da01_X18fa.ms"
# Field 2: X2188
vis12m_F2 = "../../2019.1.00843.S/science_goal.uid___A001_X1465_X2188/group.uid___A001_X1465_X2189/member.uid___A001_X1465_X218a/calibrated/working/uid___A002_Xe48598_Xafeb.ms"
# Field 3: X21a0
vis12m_F3 = ["../../2019.1.00843.S/science_goal.uid___A001_X1465_X21a0/group.uid___A001_X1465_X21a1/member.uid___A001_X1465_X21a2/calibrated/working/uid___A002_Xe3f5bb_X21c9.ms",
            "../../2019.1.00843.S/science_goal.uid___A001_X1465_X21a0/group.uid___A001_X1465_X21a1/member.uid___A001_X1465_X21a2/calibrated/working/uid___A002_Xe407cf_X1887b.ms"]
# Field 4: X2198
vis12m_F4 = "../../2019.1.00843.S/science_goal.uid___A001_X1465_X2198/group.uid___A001_X1465_X2199/member.uid___A001_X1465_X219a/calibrated/working/uid___A002_Xe407cf_X15ea.ms"
# Field 5: X2190
vis12m_F5 = "../../2019.1.00843.S/science_goal.uid___A001_X1465_X2190/group.uid___A001_X1465_X2191/member.uid___A001_X1465_X2192/calibrated/working/uid___A002_Xe45e29_X9839.ms"


# Identify line free channels. (done "live")
plotms(vis12m_F1,iteraxis='spw',xaxis='chan',yaxis='amp',avgtime='1e8')
######################################################
# Parameters:

Nchans = 1
NumIter = 15000
cleanthres = '0.5mJy'
cellsize = '0.25arcsec'
mapsize = [1000,1000,1000,1000,1000]
# We image each field seperately and have these phasecenters
phasecenters = ['J2000 05h38m32.0 -69d02m18.0','J2000 05h38m34.0 -69d04m38.0','J2000 05h39m0.0 -69d02m45.0','J2000 05h39m0.0 -69d04m42.0','J2000 05h38m36.0 -69d07m00.0']

###########################################################


spw = "31:20~200,31:350~430,29:20~550,29:1600~1900,27:0~400,27:1400~1900,25:10~550,25:1500~1900,33:10~50,33:400~450,35:10~50,35:400~450,37:10~50,37:400~450"
visdata = [vis12m_F1,vis12m_F2,vis12m_F3,vis12m_F4,vis12m_F5]
imagenames = ['30DOR_F'+str(w)+"_cont_all_spw_taper0arcsec" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]



for i in [0,1,2,3,4]:
    if i==2:
        spws = 2*[spw]
    else:
        spws = spw
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    tclean(vis=visdata[i],
        imagename = out_file,
        field = '30_Doradus',
        spw = spws,
        specmode = 'mfs',
        outframe = 'LSRK',
        niter = NumIter,
        threshold = cleanthres,
        deconvolver = 'hogbom',
        gridder = 'mosaic',
        imsize = mapsize[i],
        cell = cellsize,
        weighting = 'briggs',
        robust = 0.5,
        restoringbeam=["1.75arcsec","1.75arcsec","0deg"],
        phasecenter = phasecenters[i],
        #uvtaper=['10arcsec','10arcsec','0deg'],
        parallel=True,
        pbcor = False)



spw = "31:20~200,31:350~430,33:10~200,33:400~450,35:10~200,35:400~450,37:10~200,37:400~450"
visdata = [vis12m_F1,vis12m_F2,vis12m_F3,vis12m_F4,vis12m_F5]
imagenames = ['30DOR_F'+str(w)+"_cont_spw_31_33_35_37_spw_taper0arcsec" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]



for i in [1,2,3,4,0]:
    if i==2:
        spws = 2*[spw]
    else:
        spws = spw
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    tclean(vis=visdata[i],
        imagename = out_file,
        field = '30_Doradus',
        spw = spws,
        specmode = 'mfs',
        outframe = 'LSRK',
        niter = NumIter,
        threshold = cleanthres,
        deconvolver = 'hogbom',
        gridder = 'mosaic',
        imsize = mapsize[i],
        cell = cellsize,
        weighting = 'briggs',
        robust = 0.5,
        restoringbeam=["1.75arcsec","1.75arcsec","0deg"],
        phasecenter = phasecenters[i],
        #uvtaper=['10arcsec','10arcsec','0deg'],
        parallel=True,
        pbcor = False)



# Only H30a window:27
spw = "27:10~700,27:1300~1910"
visdata = [vis12m_F1,vis12m_F2,vis12m_F3,vis12m_F4,vis12m_F5]
imagenames = ['30DOR_F'+str(w)+"_cont_spw_27_spw_taper0arcsec" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]
for i in [1,2,3,4,0]:
    if i==2:
        spws = 2*[spw]
    else:
        spws = spw
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    tclean(vis=visdata[i],
        imagename = out_file,
        field = '30_Doradus',
        spw = spws,
        specmode = 'mfs',
        outframe = 'LSRK',
        niter = NumIter,
        threshold = cleanthres,
        deconvolver = 'hogbom',
        gridder = 'mosaic',
        imsize = mapsize[i],
        cell = cellsize,
        weighting = 'briggs',
        robust = 0.5,
        restoringbeam=["1.75arcsec","1.75arcsec","0deg"],
        phasecenter = phasecenters[i],
        #uvtaper=['10arcsec','10arcsec','0deg'],
        parallel=True,
        pbcor = False)












# MOSAICS
# Weighted mosaics using python
import numpy as np
from astropy.stats import mad_std
from astropy import wcs
from astropy.io import fits
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from FITS_tools import regrid_cube
from FITS_tools.downsample import *
from glob import glob

imlist = glob('../field?/*cont_all_spw_taper10arcsec.image')
for im in imlist:
    outfile = './fitsfiles/'+im.split('/')[-1].replace('.image','.fits')
    if not os.path.exists(outfile):        
        exportfits(im, outfile, velocity=True,dropdeg=True,overwrite=True)

imlist = glob('../field?/*cont_all_spw_taper10arcsec.pb')
for im in imlist:
    outfile = './fitsfiles/'+im.split('/')[-1].replace('.pb','.pb.fits')
    if not os.path.exists(outfile):        
        exportfits(im, outfile, velocity=True,dropdeg=True,overwrite=True)




# Python: get var maps
def get_var_map(im):
    flatimg, hdr = fits.getdata(im,header=True)
    rms = mad_std(flatimg, axis=None, ignore_nan=True)
    pbfile = im.replace('.fits','.pb.fits')
    pbimage, pbhd = fits.getdata(pbfile,header=True)
    pbcorimg = flatimg / pbimage
    varimg  = (rms / pbimage)**2
    fits.writeto(im.replace('.fits','.var.fits'),
                     varimg.astype(np.float32), hdr, overwrite=True)

imlist = glob('./fitsfiles/*arcsec.fits')
for im in imlist:
    get_var_map(im,)


hd2d = fits.getheader('./12m_continuum.model.fits')
naxis2 = hd2d['naxis2']
naxis1 = hd2d['naxis1']


naxis3 = 1
imlist = glob('./fitsfiles/*arcsec.fits')
outputnames='30Dor_continuum'

mcube = np.zeros(shape=(naxis2,naxis1))
foot  = np.zeros(shape=(naxis2,naxis1))
wtnse = np.zeros(shape=(naxis2,naxis1))

varlist = [w.replace('.fits','.var.fits') for w in imlist]
N_ims = len(imlist)
hdu_cube = [None] * N_ims
hdu_slic = [None] * N_ims
var_cube = [None] * N_ims
var_slic = [None] * N_ims
wcs_obj = [None] * N_ims
for i, im in enumerate(imlist):
    hdu_cube[i] = fits.open(imlist[i])[0]
    var_cube[i] = fits.open(varlist[i])[0]



variance_wts = []
for i, im in enumerate(imlist):
    hdu_slic[i] = fits.PrimaryHDU(data=hdu_cube[i].data, header=hdu_cube[i].header)
    var_slic[i] = fits.PrimaryHDU(data=var_cube[i].data, header=var_cube[i].header)
    variance_wts.append(1/var_slic[i].data)
    wcs_obj[i] = wcs.WCS(hdu_slic[i])
data_tuple = [None]*N_ims
for i, im in enumerate(imlist):
    data_tuple[i] = (hdu_slic[i].data,wcs_obj[i])
mcube, foot = reproject_and_coadd(data_tuple, hd2d,input_weights=variance_wts,
				        reproject_function=reproject_interp)

fits.writeto(outputnames+'.mosaic.fits',mcube.astype(np.float32), hd2d, overwrite=True)
fits.writeto(outputnames+'.foot.fits',foot.astype(np.float32), hd2d, overwrite=True)
foot_sqrt = np.sqrt(foot)



foot_sqrt[foot_sqrt==0] = np.nan
wtnse = np.nanmean(1/foot_sqrt, axis=0, keepdims=True)
fits.writeto(outputnames+'.rms.fits',wtnse.astype(np.float32), hd2d, overwrite=True)



