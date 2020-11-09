# step2: Imaging
######################################################
# Visibilities
# 7m Data
# Field 1: X2180
vis7m_F1 = ["../7m_visdata/uid___A002_Xe1baa0_X8097.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe220f2_X4ea1.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe230a1_X28f6.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe230a1_X32eb.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe220f2_X5842.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe247d0_Xeb.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe220f2_X62c1.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe27761_X172c.ms.split16.contsub"]

# Field 2: X2188
vis7m_F2 = ["../7m_visdata/uid___A002_Xe31981_Xf7ea.ms.split.cal.split16.contsub",
	"../7m_visdata/uid___A002_Xe32bed_Xdce2.ms.split.cal.split16.contsub",
	"../7m_visdata/uid___A002_Xe32bed_Xe889.ms.split.cal.split16.contsub",
	"../7m_visdata/uid___A002_Xe37224_X3836.ms.split.cal.split16.contsub",
	"../7m_visdata/uid___A002_Xe37224_X461a.ms.split.cal.split16.contsub",
	"../7m_visdata/uid___A002_Xe37224_Xdfe9.ms.split.cal.split16.contsub",
	"../7m_visdata/uid___A002_Xe37224_Xe70a.ms.split.cal.split16.contsub"
	]
# Field 3: X21a0
vis7m_F3 = ["../7m_visdata/uid___A002_Xe1a561_X2fad.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe1baa0_X477e.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe1baa0_X7737.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe1baa0_X7d43.ms.split16.contsub"]
# Field 4: X2198
vis7m_F4 = ["../7m_visdata/uid___A002_Xe2ada9_X18144.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe2ada9_X18ea1.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe31981_X3cee.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe31981_Xebf5.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe2ada9_X187a8.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe31981_X30e4.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe31981_X47c7.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe32bed_X3f50.ms.split16.contsub"]
# Field 5: X2190
vis7m_F5 = ["../7m_visdata/uid___A002_Xe3da01_X305b.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe407cf_X9d8c.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe3da01_X16c5.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe407cf_X17cba.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe407cf_Xbb9d.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe3da01_X2352.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe407cf_X20de.ms.split16.contsub",
	"../7m_visdata/uid___A002_Xe407cf_Xf987.ms.split16.contsub"]
# Zero-spacing data (either TP or TP+7m)


# Regrid to the model image frame
#imregrid(imagename='../TotalPower/TP_12CO.line.subim',template='7m.modelimage/',output='TP_12CO_regrid.image')
zerospacing = '../TotalPower/12CO_TP_29sept.fits'

######################################################
# Parameters:
# 12CO(2-1)
# Rest freq of CO2-1 line
linefreq = '230.5380GHz'

Vstart = '290km/s'
chanwidth = '-83.33333333m/s'

# Channels could be a bit smaller, but this seems to give more than enough resolution
#Nchans = 300
Nchans = 1100
NumIter = 5000
cleanthres = '0.35mJy'
cellsize = '1.5arcsec'
mapsize = [180,180,180,180,180]
#mapsize_large = 600
#phasecenter_large="J2000 05h38m40.7 -69d04m29.47" # for complete field
phasecenters = ['J2000 05h38m32.0 -69d02m18.0','J2000 05h38m34.0 -69d04m38.0','J2000 05h39m0.0 -69d02m45.0','J2000 05h39m0.0 -69d04m42.0','J2000 05h38m36.0 -69d07m00.0']
SDFAC = 1.

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
    # We feather with the lower resolution (7m+TP data). Important step. 
    # For field 1,3,4 these are the seperate 7m+TP fields. For 2+5 this is 
    # the mosaic of the 7m+TP data from the individual fields

    print('=============================================================')
    print(' Running Feathering....')
    # Feather
    os.system("rm -rf "+out_file+".pbcor.feather.image")
    feather(imagename=out_file+".pbcor.feather.image",
	    highres=out_file+".pbcor.image",
	    lowres=zerospacing,sdfactor=SDFAC)
    os.system("rm -rf "+out_file+".feather.image")
    feather(imagename=out_file+".feather.image",
	    highres=out_file+".image",
	    lowres=zerospacing,sdfactor=SDFAC)

    ###########################################################
    # Convolve to same resolution
    # Smoothen everything to resolution of TP array to compare 
    # Since the combined image should result in the same flux as 
    # The TP alone data
    #print('=============================================================')
    #print(' Running convolutions....')
    os.system("rm -rf "+out_file+".pbcor.image.TPres")
    #imsmooth(imagename=out_file+".pbcor.image",outfile=out_file+".pbcor.image.TPres",
#	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".image.TPres")
    #imsmooth(imagename=out_file+".image",outfile=out_file+".image.TPres",
#	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".feather.image.TPres")
    #imsmooth(imagename=out_file+".feather.image",outfile=out_file+".feather.image.TPres",
#	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    os.system("rm -rf "+out_file+".pbcor.feather.image.TPres")
    #imsmooth(imagename=out_file+".pbcor.feather.image",outfile=out_file+".pbcor.feather.image.TPres",
#	    major=TP_bmaj,minor=TP_bmin,pa=TP_bpa,targetres=True)
    
    os.system('rm -rf '+out_file+'*FullMap*')

    # Regrid everything to the model image such that we can easily make a mosaic
    #imregrid(out_file+".image",'7m.modelimage',out_file+".FullMapRegrid.image",overwrite=True)
    #imregrid(out_file+".pbcor.image",'7m.modelimage',out_file+".pbcor.FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".pbcor.feather.image",'7m.modelimage',out_file+".pbcor.feather.FullMapRegrid.image",overwrite=True)
    imregrid(out_file+".feather.image",'7m.modelimage',out_file+".feather.FullMapRegrid.image",overwrite=True)

    #exportfits(out_file+".FullMapRegrid.image",out_file+".FullMapRegrid.fits",overwrite=True)
    #exportfits(out_file+".pbcor.FullMapRegrid.image",out_file+".pbcor.FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".pbcor.feather.FullMapRegrid.image",out_file+".pbcor.feather.FullMapRegrid.fits",overwrite=True)
    exportfits(out_file+".feather.FullMapRegrid.image",out_file+".feather.FullMapRegrid.fits",overwrite=True)


###########################################################
# Image everything without model first. This is the most natural
# way. Only afterwards we feather with the TP data. 
#os.system('mkdir field1 field2 field3 field4 field5')
visdata = [vis7m_F1,vis7m_F2,vis7m_F3,vis7m_F4,vis7m_F5]
imagenames = ['30DOR_F'+str(w)+"_7m+TP_CLEAN_smallvelo_largepix2+nomodel" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]
for i in [0,1,2,3,4]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
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
        restoringbeam=["7.0arcsec","7.0arcsec","0.0deg"],
        phasecenter = phasecenters[i],
        pbcor = False)
    #process(out_file,zerospacing)

# Now image with modelimage (the TP model). This is called hybrid imaging, since
# we also feather the two images afterwards
visdata = [vis7m_F1,vis7m_F2,vis7m_F3,vis7m_F4,vis7m_F5]
imagenames = ['30DOR_F'+str(w)+"_7m+TP_CLEAN_smallvelo_largepix+model" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]
for i in [0,1,2,3,4]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
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
        restoringbeam=["7.0arcsec","7.0arcsec","0.0deg"],
        phasecenter = phasecenters[i],
        pbcor = False)
    #process(out_file,zerospacing)


# MOSAICS
import os
os.system('mkdir mosaics')
os.system('cp -r ./field?/*F*model*Full*Regrid*fits ./mosaics')

# Start ipython in the mosaics folder!

import numpy
import bottleneck
import time
from astropy.io import fits
import numpy as np
from glob import glob

def run(regex,outname):
    fitsfiles = glob(regex)
    datas = np.array([fits.open(w)[0].data for w in fitsfiles])
    newarr = bottleneck.nanmedian(datas, axis=0)

    hdul_out = fits.open(fitsfiles[0])
    hdul_out[0].data = newarr
    hdul_out.writeto(outname,overwrite=True)
maps = ['30DOR_F?_7m+TP_CLEAN+model.feather.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+model.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+model.pbcor.feather.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+model.pbcor.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+nomodel.feather.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+nomodel.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+nomodel.pbcor.feather.FullMapRegrid.fits',
        '30DOR_F?_7m+TP_CLEAN+nomodel.pbcor.FullMapRegrid.fits']
names = ['7m+TP+model.feather.fits',
         '7m+TP+model.fits',   
         '7m+TP+model.pbcor.feather.fits',
         '7m+TP+model.pbcor.fits',
         '7m+TP+nomodel.feather.fits',
         '7m+TP+nomodel.fits',   
         '7m+TP+nomodel.pbcor.feather.fits',
         '7m+TP+nomodel.pbcor.fits']
for i in range(0,len(maps)):
    run(maps[i],names[i])



