linefreq = '230.5380GHz'
Vstart = '282.715km/s'
chanwidth = '-80m/s'
Nchans = 701

NumIter = 150
cleanthres = '0.2mJy'
cellsize = '0.3arcsec'
mapsize = [1500, 1890]
phasecenter = "J2000 05h38m40.7 -69d04m29.47"
zerospacing = '../Image_7m_data/mosaic_with_feathered/mosaic_using_feathers.image.smooth/' 
SDFAC = 1.05


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
	    lowres=zerospacing,sdfactor=SDFAC)
    os.system("rm -rf "+out_file+".feather.image")
    feather(imagename=out_file+".feather.image",
	    highres=out_file+".image",
	    lowres=zerospacing,sdfactor=SDFAC)

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
target_dir = './AllFields_deeper/'
imagename = '12m_alldata_CO+model_allchans'
out_file = target_dir + imagename
clean(vis='../12m_visdata/12m.concat25.contsub/',
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
    modelimage=zerospacing,
    phasecenter = phasecenter,
    pbcor = False)
process(out_file,zerospacing)
