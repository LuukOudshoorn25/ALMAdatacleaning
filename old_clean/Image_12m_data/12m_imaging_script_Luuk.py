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

# Zero-spacing data (either TP or TP+7m)


######################################################
# Parameters:
# 12CO(2-1)
linefreq = '230.5380GHz'
Vstart = '260.715km/s'
chanwidth = '-200m/s'
Nchans = 40
NumIter = 5000
cleanthres = '0.25mJy'
cellsize = '0.333333arcsec'
mapsize =[900,900,900,900,900]
#mapsize_large = 600
#phasecenter_large="J2000 05h38m40.7 -69d04m29.47" # for complete field
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

    imregrid(out_file+".image",zerospacing,out_file+".TPregrid.image")
    imregrid(out_file+".pbcor.image",zerospacing,out_file+".pbcor.TPregrid.image")
    imregrid(out_file+".pbcor.feather.image",zerospacing,out_file+".pbcor.feather.TPregrid.image")

    exportfits(out_file+".TPregrid.image",out_file+".TPregrid.fits")
    exportfits(out_file+".pbcor.TPregrid.image",out_file+".pbcor.TPregrid.fits")
    exportfits(out_file+".pbcor.feather.TPregrid.image",out_file+".pbcor.feather.TPregrid.fits")


###########################################################


visdata = [vis12m_F1,vis12m_F2,vis12m_F3,vis12m_F4,vis12m_F5]
imagenames = ['30DOR_F'+str(w)+"_12m+7m+TP_CLEAN+model" for w in range(1,6)]
#imagenames_nomodel = ['30DOR_F'+str(w)+"_7m+TP_CLEAN_nomodel" for w in range(1,6)]
#imagenames_oldmodel = ['30DOR_F'+str(w)+"_7m+TP_CLEAN_oldmodel" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]
for i in [2]:
    target_dir = targetdirs[i]
    out_file = target_dir + imagenames[i]
    #if not os.path.exists(out_file+'.image'):
    zerospacing = '../Image_7m_data/field'+str(i+1)+'/30DOR_F'+str(i+1)+'_7m+TP_CLEAN+model.feather.image'
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

