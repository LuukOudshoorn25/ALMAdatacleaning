# step2: Imaging

######################################################
# Visibilities
vis12m_F1 = ["../12m_visdata/uid___A002_Xe3da01_X18fa.ms.split25.contsub"]
# Field 2: X2188
vis12m_F2 = ["../12m_visdata/uid___A002_Xe48598_Xafeb.ms.split25.contsub"]
# Field 3: X21a0
vis12m_F3 = ["../12m_visdata/uid___A002_Xe3f5bb_X21c9.ms.split25.contsub",
             "../12m_visdata/uid___A002_Xe407cf_X1887b.ms.split25.contsub"]
# Field 4: X2198
vis12m_F4 = ["../12m_visdata/uid___A002_Xe407cf_X15ea.ms.split25.contsub"]
# Field 5: X2190
vis12m_F5 = ["../12m_visdata/uid___A002_Xe45e29_X9839.ms.split25.contsub"]


# Field 2 +5 are imaged simultaneous
vis12m_F2_5 = ["../12m_visdata/uid___A002_Xe48598_Xafeb.ms.split25.contsub",
               "../12m_visdata/uid___A002_Xe45e29_X9839.ms.split25.contsub"]
# Zero-spacing data (either TP or TP+7m)


######################################################
# Parameters:
# 12CO(2-1)
linefreq = '230.5380GHz'
Vstart = '282.715km/s'
chanwidth = '-250m/s'
Nchans = 300
NumIter = 8000
cleanthres = '0.18mJy'
cellsize = '0.25arcsec'
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
    #os.system("rm -rf "+out_file+".pbcor.image.TPres")
    #imsmooth(imagename=out_file+".pbcor.image",outfile=out_file+".pbcor.image.TPres",
#	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
#    os.system("rm -rf "+out_file+".image.TPres")
#    imsmooth(imagename=out_file+".image",outfile=out_file+".image.TPres",
#	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".feather.image.TPres")
    imsmooth(imagename=out_file+".feather.image",outfile=out_file+".feather.image.TPres",
	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".pbcor.feather.image.TPres")
    imsmooth(imagename=out_file+".pbcor.feather.image",outfile=out_file+".pbcor.feather.image.TPres",
	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)

    #imregrid(out_file+".image",'12m.modelimage',out_file+".FullMapRegrid.image",overwrite=True)
    #imregrid(out_file+".pbcor.image",'12m.modelimage',out_file+".pbcor.FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".pbcor.feather.image",'12m.modelimage',out_file+".pbcor.feather.FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".feather.image",'12m.modelimage',out_file+".feather.FullMapRegrid.image",overwrite=True)

    #exportfits(out_file+".FullMapRegrid.image",out_file+".FullMapRegrid.fits",overwrite=True)
    #exportfits(out_file+".pbcor.FullMapRegrid.image",out_file+".pbcor.FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".pbcor.feather.FullMapRegrid.image",out_file+".pbcor.feather.FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".feather.FullMapRegrid.image",out_file+".feather.FullMapRegrid.fits",overwrite=True)


###########################################################
visdata = [vis12m_F1,vis12m_F2,vis12m_F3,vis12m_F4,vis12m_F5]
imagenames = ['30DOR_F'+str(w)+"_12m+7m+TP_CLEAN+nomodel" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]

for i in [0,2,3]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    zerospacing = '../Image_7m_data/field'+str(i+1)+'/30DOR_F'+str(i+1)+'_7m+TP_CLEAN+nomodel.feather.image'
    """clean(vis=visdata[i],
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
        pbcor = False)"""
    process(out_file,zerospacing)



imagenames = ['30DOR_F'+str(w)+"_12m+7m+TP_CLEAN+model" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]

for i in [0,1,2,3,4]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
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
os.system('mkdir ./field2_5/')
importfits('../Image_7m_data/mosaics/7m+TP+nomodel.feather.fits','./field2_5/7m+TP+nomodel.feather.image')
zerospacing = './field2_5/7m+TP+nomodel.feather.image'
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
os.system('mkdir mosaics')
os.system('cp -r ./field*/*F1*+nomodel*Full*Regrid*fits ./mosaics')
os.system('cp -r ./field*/*F3*+nomodel*Full*Regrid*fits ./mosaics')
os.system('cp -r ./field*/*F4*+nomodel*Full*Regrid*fits ./mosaics')
os.system('cp -r ./field2_5/30DOR_F2+5_12*+nomodel*Full*Regrid*fits ./mosaics')
# Start ipython!
import numpy
import os
import bottleneck
from astropy.io import fits
import numpy as np
from glob import glob
import time

def run(regex,outname):
    fitsfiles = glob(regex)
    datas = np.array([fits.open(w)[0].data for w in fitsfiles])
    newarr = bottleneck.nanmedian(datas, axis=0)

    hdul_out = fits.open(fitsfiles[0])
    hdul_out[0].data = newarr
    hdul_out.writeto(outname,overwrite=True)
maps = ['30DOR_F*_*7m+TP_CLEAN+nomodel.feather.FullMapRegrid.fits',
        '30DOR_F*_*7m+TP_CLEAN+nomodel.FullMapRegrid.fits',
        '30DOR_F*_*7m+TP_CLEAN+nomodel.pbcor.feather.FullMapRegrid.fits',
        '30DOR_F*_*7m+TP_CLEAN+nomodel.pbcor.FullMapRegrid.fits',
        #'30DOR_F*_*7m+TP_CLEAN+model.feather.FullMapRegrid.fits',
        #'30DOR_F*_*7m+TP_CLEAN+model.FullMapRegrid.fits',
        #'30DOR_F*_*7m+TP_CLEAN+model.pbcor.feather.FullMapRegrid.fits',
        #'30DOR_F*_*7m+TP_CLEAN+model.pbcor.FullMapRegrid.fits',
]
names = ['12+7m+TP+nomodel.feather.fits',
         '12+7m+TP+nomodel.fits',   
         '12+7m+TP+nomodel.pbcor.feather.fits',
         '12+7m+TP+nomodel.pbcor.fits',
         #'12+7m+TP+model.feather.fits',
         #'12+7m+TP+model.fits',   
         #'12+7m+TP+model.pbcor.feather.fits',
         #'12+7m+TP+model.pbcor.fits',
]
for i in range(0,len(maps)):
    run(maps[i],names[i])


# run casa in mosaics folder

# run casa in mosaics folder
# Cut out the region with data
box='41,77,1261,1727'
# Drop first few channels since there is no data. 
# Drop last few since there is no emission

chans='3~283'
imsubimage(imagename='12+7m+TP+nomodel.feather.fits',box=box,chans=chans,outfile='30Dor_12CO.image/')
imsubimage(imagename='12+7m+TP+nomodel.pbcor.feather.fits',box=box,chans=chans,outfile='30Dor_12CO_pbcor.image/')

# Export v_cubes and freq_cubes
exportfits('30Dor_12CO.image/','30Dor_12CO.fits')
exportfits('30Dor_12CO.image/','30Dor_12CO_vcube.fits', velocity=True)

exportfits('30Dor_12CO_pbcor.image/','30Dor_12CO_pbcor.fits')
exportfits('30Dor_12CO_pbcor.image/','30Dor_12CO_pbcor_vcube.fits', velocity=True)




















