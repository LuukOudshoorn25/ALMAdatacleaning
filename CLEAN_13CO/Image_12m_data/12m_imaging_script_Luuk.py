# step2: Imaging

######################################################
# Visibilities
vis12m_F1 = ["../12m_visdata/uid___A002_Xe3da01_X18fa.ms.split29.contsub"]
# Field 2: X2188
vis12m_F2 = ["../12m_visdata/uid___A002_Xe48598_Xafeb.ms.split29.contsub"]
# Field 3: X21a0
vis12m_F3 = ["../12m_visdata/uid___A002_Xe3f5bb_X21c9.ms.split29.contsub",
             "../12m_visdata/uid___A002_Xe407cf_X1887b.ms.split29.contsub"]
# Field 4: X2198
vis12m_F4 = ["../12m_visdata/uid___A002_Xe407cf_X15ea.ms.split29.contsub"]
# Field 5: X2190
vis12m_F5 = ["../12m_visdata/uid___A002_Xe45e29_X9839.ms.split29.contsub"]

# Field 2 +5 are imaged simultaneous
vis12m_F2_5 = ["../12m_visdata/uid___A002_Xe48598_Xafeb.ms.split29.contsub",
               "../12m_visdata/uid___A002_Xe45e29_X9839.ms.split29.contsub"]

# Zero-spacing data (either TP or TP+7m)


######################################################
# Parameters:
# 13CO(2-1)
linefreq = '220.39870060GHz'
Vstart = '282.715km/s'
chanwidth = '-250m/s'
Nchans = 300
NumIter = 8000
cleanthres = '0.18mJy'
cellsize = '0.25arcsec'
# Clean field 1,3,4
mapsize = [1000,1000,1000,1000,1000]
phasecenters = ['J2000 05h38m32.0 -69d02m18.0','J2000 05h38m34.0 -69d04m38.0','J2000 05h39m0.0 -69d02m45.0','J2000 05h39m0.0 -69d04m42.0','J2000 05h38m36.0 -69d07m00.0']

###########################################################

def process(out_file,zerospacing):
    TP_bmaj = str(imhead(zerospacing,hdkey='bmaj',mode='get')['value'])+'arcsec'
    TP_bmin = str(imhead(zerospacing,hdkey='bmin',mode='get')['value'])+'arcsec'
    TP_bpa   = str(imhead(zerospacing,hdkey='bpa',mode='get')['value'])+'deg'
    print('=============================================================')
    print(' PB correction....')
    # PB correction
    os.system("rm -rf "+out_file+".pbcor.image")
    impbcor(imagename=out_file+".image",pbimage=out_file+".flux.pbcoverage", outfile=out_file+".pbcor.image")

    ###########################################################
    print('=============================================================')
    print(' Running Feathering....')
    # Feather
    os.system("rm -rf "+out_file+".pbcor.feather.image")
    feather(imagename=out_file+".pbcor.feather.image",
	    highres=out_file+".pbcor.image",
	    lowres=zerospacing)
    os.system("rm -rf "+out_file+".feather.image")
    feather(imagename=out_file+".feather.image",
	    highres=out_file+".image",
	    lowres=zerospacing)

    ###########################################################
    # Convolve to same resolution
    print('=============================================================')
    print(' Running convolutions....')
    os.system("rm -rf "+out_file+".pbcor.image.TPres")
    imsmooth(imagename=out_file+".pbcor.image",outfile=out_file+".pbcor.image.TPres",
	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".image.TPres")
    imsmooth(imagename=out_file+".image",outfile=out_file+".image.TPres",
	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".feather.image.TPres")
    imsmooth(imagename=out_file+".feather.image",outfile=out_file+".feather.image.TPres",
	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+out_file+".pbcor.feather.image.TPres")
    imsmooth(imagename=out_file+".pbcor.feather.image",outfile=out_file+".pbcor.feather.image.TPres",
	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)

    imregrid(out_file+".image",'12m.modelimage',out_file+".FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".pbcor.image",'12m.modelimage',out_file+".pbcor.FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".pbcor.feather.image",'12m.modelimage',out_file+".pbcor.feather.FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".feather.image",'12m.modelimage',out_file+".feather.FullMapRegrid.image",overwrite=True)

    exportfits(out_file+".FullMapRegrid.image",out_file+".FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".pbcor.FullMapRegrid.image",out_file+".pbcor.FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".pbcor.feather.FullMapRegrid.image",out_file+".pbcor.feather.FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".feather.FullMapRegrid.image",out_file+".feather.FullMapRegrid.fits",overwrite=True)


###########################################################
visdata = [vis12m_F1,vis12m_F2,vis12m_F3,vis12m_F4,vis12m_F5]
imagenames = ['30DOR_F'+str(w)+"_12m+7m+TP_CLEAN+nomodel" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]

for i in [4]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    #if not os.path.exists(out_file+'.image'):
    zerospacing = '../Image_7m_data/field'+str(i+1)+'/30DOR_F'+str(i+1)+'_7m+TP_CLEAN+nomodel.feather.image'
    clean(vis=visdata[i],
        imagename = out_file,
        field = '30_Doradus',
        spw = '*',
        mode = 'velocity',
        nchan = Nchans,
        interpolation='linear',
        start = Vstart,
        width = chanwidth,
        restfreq = linefreq,
        outframe = 'LSRK',
        niter = NumIter,
        threshold = cleanthres,
        psfmode = 'hogbom',
        imagermode = 'mosaic',
        interactive = False,
        imsize = mapsize[i],
        cell = cellsize,
        weighting = 'briggs',
        robust = 0.5,
        #modelimage=zerospacing,
        restoringbeam=["1.75arcsec","1.75arcsec","0deg"],
        phasecenter = phasecenters[i],
        pbcor = False)
    process(out_file,zerospacing)



imagenames = ['30DOR_F'+str(w)+"_12m+7m+TP_CLEAN+model" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]

for i in [0,2,3]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    #if not os.path.exists(out_file+'.image'):
    zerospacing = '../Image_7m_data/field'+str(i+1)+'/30DOR_F'+str(i+1)+'_7m+TP_CLEAN+nomodel.feather.image'
    clean(vis=visdata[i],
        imagename = out_file,
        field = '30_Doradus',
        spw = '*',
        mode = 'velocity',
        nchan = Nchans,
        interpolation='linear',
        start = Vstart,
        width = chanwidth,
        restfreq = linefreq,
        outframe = 'LSRK',
        niter = NumIter,
        threshold = cleanthres,
        psfmode = 'hogbom',
        imagermode = 'mosaic',
        interactive = False,
        imsize = mapsize[i],
        cell = cellsize,
        weighting = 'briggs',
        robust = 0.5,
        modelimage=zerospacing,
        restoringbeam=["1.75arcsec","1.75arcsec","0deg"],
        phasecenter = phasecenters[i],
        pbcor = False)
    process(out_file,zerospacing)



####################################
# Image field 2 and 5 since there are strong emission features on the edges between the two
os.system('mkdir ./field2_5_new/')
importfits('../Image_7m_data/mosaics/7m+TP+nomodel.feather.fits','./field2_5/7m+TP+nomodel.feather.image')
zerospacing = './field2_5/7m+TP+nomodel.feather.image'
linefreq = '220.39870060GHz'
Vstart = '282.715km/s'
chanwidth = '-250m/s'
Nchans = 300
NumIter = 14000
cleanthres = '0.18mJy'
cellsize = '0.25arcsec'
mapsize = [1000,1600]
phasecenter = 'J2000 05h38m34.61 -69d05m42.8'

imagename = './field2_5/30DOR_F2+5_12m+7m+TP_CLEAN+nomodel'
out_file = imagename
clean(vis=vis12m_F2_5,
	imagename = out_file,
	field = '30_Doradus',
	spw = '*',
	mode = 'velocity',
	nchan = Nchans,
	interpolation='linear',
	start = Vstart,
	width = chanwidth,
	restfreq = linefreq,
	outframe = 'LSRK',
	niter = NumIter,
	threshold = cleanthres,
	psfmode = 'hogbom',
	imagermode = 'mosaic',
	interactive = False,
	imsize = mapsize,
	cell = cellsize,
	weighting = 'briggs',
	robust = 0.5,
	#modelimage=zerospacing,
	restoringbeam=["1.75arcsec","1.75arcsec","0deg"],
	phasecenter = phasecenter,
	pbcor = False)
process(out_file,zerospacing)



# MOSAICS
# Mosaic with weigthing
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
# make and go to folder weighted_mosaic
imlist = ['../field1/30DOR_F1_12m+7m+TP_CLEAN+nomodel.feather.image',
          '../field2_5/30DOR_F2+5_12m+7m+TP_CLEAN+nomodel.feather.image',
          '../field3/30DOR_F3_12m+7m+TP_CLEAN+nomodel.feather.image',
          '../field4/30DOR_F4_12m+7m+TP_CLEAN+nomodel.feather.image',

          '../field1/30DOR_F1_12m+7m+TP_CLEAN+nomodel.flux.pbcoverage',
          '../field2_5/30DOR_F2+5_12m+7m+TP_CLEAN+nomodel.flux.pbcoverage',
          '../field3/30DOR_F3_12m+7m+TP_CLEAN+nomodel.flux.pbcoverage',
          '../field4/30DOR_F4_12m+7m+TP_CLEAN+nomodel.flux.pbcoverage',

          '../field1/30DOR_F1_12m+7m+TP_CLEAN+nomodel.pbcor.feather.image',
          '../field2_5/30DOR_F2+5_12m+7m+TP_CLEAN+nomodel.pbcor.feather.image',
          '../field3/30DOR_F3_12m+7m+TP_CLEAN+nomodel.pbcor.feather.image',
          '../field4/30DOR_F4_12m+7m+TP_CLEAN+nomodel.pbcor.feather.image',

          '../field1/30DOR_F1_12m+7m+TP_CLEAN+model.feather.image',
          '../field2/30DOR_F2_12m+7m+TP_CLEAN+model.feather.image',
          '../field3/30DOR_F3_12m+7m+TP_CLEAN+model.feather.image',
          '../field4/30DOR_F4_12m+7m+TP_CLEAN+model.feather.image',
          '../field5/30DOR_F5_12m+7m+TP_CLEAN+model.feather.image',

          '../field1/30DOR_F1_12m+7m+TP_CLEAN+model.flux.pbcoverage',
          '../field2/30DOR_F2_12m+7m+TP_CLEAN+model.flux.pbcoverage',
          '../field3/30DOR_F3_12m+7m+TP_CLEAN+model.flux.pbcoverage',
          '../field4/30DOR_F4_12m+7m+TP_CLEAN+model.flux.pbcoverage',
          '../field5/30DOR_F5_12m+7m+TP_CLEAN+model.flux.pbcoverage',

          '../field1/30DOR_F1_12m+7m+TP_CLEAN+model.pbcor.feather.image',
          '../field2/30DOR_F2_12m+7m+TP_CLEAN+model.pbcor.feather.image',
          '../field3/30DOR_F3_12m+7m+TP_CLEAN+model.pbcor.feather.image',
          '../field4/30DOR_F4_12m+7m+TP_CLEAN+model.pbcor.feather.image',
          '../field5/30DOR_F5_12m+7m+TP_CLEAN+model.pbcor.feather.image']
os.system('mkdir fitsfiles)
# CASA: Export all fitsfiles
for im in imlist:
    if 'image' in im:
        outfile = './fitsfiles/'+im.split('/')[-1].replace('.image','.fits')
        if not os.path.exists(outfile):        
            exportfits(im, outfile, velocity=True,dropdeg=True,overwrite=True)
    elif 'pbco' in im:
        outfile = './fitsfiles/'+im.split('/')[-1].replace('.flux.pbcoverage','.pb.fits')
        exportfits(im, outfile, velocity=True,dropdeg=True,overwrite=True)




# Python: get var maps
def get_var_map(im,pb_replace):
    flatimg, hdr = fits.getdata(im,header=True)
    rms = mad_std(flatimg, axis=None, ignore_nan=True)
    pbfile = im.replace(pb_replace,'.pb')
    pbimage, pbhd = fits.getdata(pbfile,header=True)
    pbcorimg = flatimg / pbimage
    varimg  = (rms / pbimage)**2
    fits.writeto(im.replace('.feather','.var'),
                     varimg.astype(np.float32), hdr, overwrite=True)

imlist = glob('./fitsfiles/*pbcor.feather.fits')
for im in imlist:
    get_var_map(im,'.pbcor.feather')
imlist = glob('./fitsfiles/*model.feather.fits')
for im in imlist:
    get_var_map(im,'.feather')


hd2d = fits.getheader('template_30Dor.fits')

head0 = fits.getheader(imlist[0])
naxis3 = head0['naxis3']
naxis2 = hd2d['naxis2']
naxis1 = hd2d['naxis1']

# To Do: nomodel: normal + pbcor
#        model:   normal + pbcor

def process_coaddition(imlist, outputnames):
    mcube = np.zeros(shape=(naxis3,naxis2,naxis1))
    foot  = np.zeros(shape=(naxis3,naxis2,naxis1))
    wtnse = np.zeros(shape=(naxis3,naxis2,naxis1))
    #
    varlist = [w.replace('.feather','.var') for w in imlist]
    N_ims = len(imlist)
    hdu_cube = [None] * N_ims
    hdu_slic = [None] * N_ims
    var_cube = [None] * N_ims
    var_slic = [None] * N_ims
    wcs_obj = [None] * N_ims
    for i, im in enumerate(imlist):
        hdu_cube[i] = fits.open(imlist[i])[0]
        var_cube[i] = fits.open(varlist[i])[0]
    # --- Loop over velocity channels
    for ich in range(naxis3):
        variance_wts = []
        for i, im in enumerate(imlist):
            hdu_slic[i] = fits.PrimaryHDU(data=hdu_cube[i].data[ich], header=hdu_cube[i].header)
            hdu_slic[i].header['WCSAXES'] = 2
            var_slic[i] = fits.PrimaryHDU(data=var_cube[i].data[ich], header=var_cube[i].header)
            var_slic[i].header['WCSAXES'] = 2
            for key in ['CRVAL3', 'CTYPE3', 'CRPIX3', 'CDELT3', 'CUNIT3']:
                del hdu_slic[i].header[key]
                del var_slic[i].header[key]
            variance_wts.append(1/var_slic[i].data)
            wcs_obj[i] = wcs.WCS(hdu_slic[i]).dropaxis(2)
        print('Working on channel',ich,'of cube ',outputnames,end='\r')
        data_tuple = [None]*N_ims
        for i, im in enumerate(imlist):
            data_tuple[i] = (hdu_slic[i].data,wcs_obj[i])
        mcube[ich], foot[ich] = reproject_and_coadd(data_tuple, hd2d,input_weights=variance_wts,
						    reproject_function=reproject_interp)
    hd3d = hdu_cube[0].header
    for key in ['CRPIX1', 'CDELT1', 'CTYPE1', 'CRVAL1', 'CRPIX2', 
		'CDELT2', 'CTYPE2', 'CRVAL2', 'LONPOLE', 'LATPOLE']:
        hd3d[key] = hd2d[key]
    fits.writeto(outputnames+'.cube.fits',mcube.astype(np.float32), hd3d, overwrite=True)
    wtnse = np.nanmean(1/np.sqrt(foot), axis=0, keepdims=True)
    fits.writeto(outputnames+'.rms.fits',wtnse.astype(np.float32), hd3d, overwrite=True)


imlist = glob('./fitsfiles/*+nomodel.feather.fits')
process_coaddition(imlist,'30Dor_13CO')

imlist = glob('./fitsfiles/*+nomodel.pbcor.feather.fits')
process_coaddition(imlist,'30Dor_13CO_pbcor')

imlist = glob('./fitsfiles/*+newmodel.feather.fits')
process_coaddition(imlist,'30Dor_13CO_hybrid')

imlist = glob('./fitsfiles/*+newmodel.pbcor.feather.fits')
process_coaddition(imlist,'30Dor_13CO_hybrid_pbcor')


























