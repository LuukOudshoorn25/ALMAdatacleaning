#
# Combine 7m and TP data
#
# This script regrids the 250 x 250 x 1100 TP data
# (2.5" x 2.5" x 64 kHz pixels) to match the 600 x 600 x 300
# 7m single-field images (0.5" x 0.5" x 0.25 km/s pixels)
# and feathers it with the 7m images.

#import os

dofields  = ['1', '2', '3', '4', '5']

prefix12 = 'PCC_12m_'
#prefixtp = '../TP/30DOR_TP_'
prefixtp = '../TP/TP_29sept_'
prefixint = '30DOR_F'
baseint = '_7m+TP_CLEAN+nomodel'
baseout   = '_7mTP_feather+nomodel.'
linename = ['12CO']

# Use TP_29sept_12CO.fits as the TP mosaic

for i in range(0,len(linename)):
    if (os.path.exists(prefixtp+linename[i]+'.image') == False):
        importfits(fitsimage=prefixtp+linename[i]+'.fits', 
            imagename=prefixtp+linename[i]+'.image',overwrite=True,
            defaultaxes=True, defaultaxesvalues=['', '', '', 'I'])
    for field in dofields:
        # Regrid TP data to match 7m data:
        imregrid(imagename=prefixtp+linename[i]+'.image',
             template=prefixint+field+baseint+'.image',
             output=prefixint+field+'.'+linename[i]+'.tpregrid', overwrite=True)
        # Apply 7m sensitivity to TP data:
        rmtables(prefixint+field+'.'+linename[i]+'.tpregrid.pb')
        immath(imagename=[prefixint+field+'.'+linename[i]+'.tpregrid',
                     prefixint+field+baseint+'.flux.pbcoverage'],
            expr='IM0*IM1',
            outfile=prefixint+field+'.'+linename[i]+'.tpregrid.pb')
        exportfits(imagename=prefixint+field+'.'+linename[i]+'.tpregrid.pb', 
            fitsimage=prefixint+field+'.'+linename[i]+'.tpregrid.pb.fits',
            dropdeg=True, velocity=True, overwrite=True)
        rmtables(prefixint+field+'.'+linename[i]+'.tpregrid')
        # Also write out 7m sensitivity and image for mosaic script:
        exportfits(imagename=prefixint+field+baseint+'.flux.pbcoverage', 
            fitsimage=prefixint+field+baseint+'.flux.pbcoverage.fits',
            dropdeg=True, velocity=True, overwrite=True)
        exportfits(imagename=prefixint+field+baseint+'.image', 
            fitsimage=prefixint+field+baseint+'.image.fits',
            dropdeg=True, velocity=True, overwrite=True)
        # Do the feathering:
        rmtables(prefixint+field+baseout+linename[i]+'.image')
        feather(imagename=prefixint+field+baseout+linename[i]+'.image',
                highres=prefixint+field+baseint+'.image',
                lowres=prefixint+field+'.'+linename[i]+'.tpregrid.pb')
        # Export to FITS:
        exportfits(imagename=prefixint+field+baseout+linename[i]+'.image', 
            fitsimage=prefixint+field+baseout+linename[i]+'.image.fits',
            dropdeg=True, velocity=True, overwrite=True)
        # Primary beam correction:
        impbcor(imagename=prefixint+field+baseout+linename[i]+'.image',
                pbimage=prefixint+field+baseint+'.flux.pbcoverage',
                outfile=prefixint+field+baseout+linename[i]+'.pbcor',
                overwrite=True)
        # Export to FITS:
        exportfits(imagename=prefixint+field+baseout+linename[i]+'.pbcor', 
            fitsimage=prefixint+field+baseout+linename[i]+'.pbcor.fits',
            dropdeg=True, velocity=True, overwrite=True)

