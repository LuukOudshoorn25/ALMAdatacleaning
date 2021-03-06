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
zerospacing = '../TotalPower/TP_concat.image.line.oldmodel/' 

######################################################
# Parameters:
# 12CO(2-1)
linefreq = '230.5380GHz'
Vstart = '282.715km/s'
chanwidth = '-80m/s'
Nchans = 701
NumIter = 5000
cleanthres = '0.25mJy'
cellsize = '1.0arcsec'
mapsize = [300,300,300,300,300]
#mapsize_large = 600
#phasecenter_large="J2000 05h38m40.7 -69d04m29.47" # for complete field
phasecenters = ['J2000 05h38m32.0 -69d02m18.0','J2000 05h38m34.0 -69d04m38.0','J2000 05h39m0.0 -69d02m45.0','J2000 05h39m0.0 -69d04m42.0','J2000 05h38m36.0 -69d07m00.0']
SDFAC = 1.51
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


visdata = [vis7m_F1,vis7m_F2,vis7m_F3,vis7m_F4,vis7m_F5]
imagenames = ['30DOR_F'+str(w)+"_7m+TP_CLEAN+model" for w in range(1,6)]
targetdirs = ['./field'+str(w)+'/' for w in range(1,6)]
for i in [1,2,3,4]:
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
    process(out_file,zerospacing)
