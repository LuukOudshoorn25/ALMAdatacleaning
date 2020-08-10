from glob import glob

imlist = glob('../field?/30*odel.feather.image')
for i,f in enumerate(imlist):
    if not os.path.exists('Field{}.regrid'.format(i+1)):
        imregrid(f,'12m_mosaic.model','Field{}.regrid'.format(i+1))

imlist = glob('./*.regrid')

for f in imlist:
    if not os.path.exists(f+'.fits'):
        exportfits(f, f+'.fits')


