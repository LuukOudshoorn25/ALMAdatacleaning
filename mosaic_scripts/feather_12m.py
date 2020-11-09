#
# Combine 12m and 7m+TP data
#

#import os

dofields  = ['1','2','3','4','5']

os.system('mkdir 7meter')
os.system('cp -r ../../Image_7m_data/mosaic_largepix/images/*pbcor ./7meter')
os.system('mkdir images')
os.system('cp -r ../field*/30DOR_F?_12m+7m+TP_CLEAN_smallvelo+nomodel.image ./images')
os.system('cp -r ../field*/30DOR_F?_12m+7m+TP_CLEAN_smallvelo+nomodel.flux.pbcoverage ./images')

prefixint = './images/30DOR_F'
dir7tp = './7meter/'
prefix7 = './7meter/30DOR_F'
base7tp = '_7m+TP_feather_smallvelo_largepix+nomodel.'
base12m = '_12m+7m+TP_CLEAN_smallvelo+nomodel'
#base12m = '_12m+7m+TP_CLEAN+nomodel.2as'
baseout   = '_12m7mTP_feather+nomodel.'
linename = ['12CO']

for i in range(0,len(linename)):
    for field in dofields:

#         if (os.path.exists(prefixint+field+baseint+'.flux3d') == False):
#             imsubimage(imagename=prefixint+field+baseint+'.flux.pbcoverage',
#                 outfile=prefixint+field+baseint+'.flux3d', dropdeg=True,
#                 overwrite=True)
#         if (os.path.exists(prefixint+field+baseint+'.image3d') == False):
#             imsubimage(imagename=prefixint+field+baseint+'.image',
#                 outfile=prefixint+field+baseint+'.image3d', dropdeg=True,
#                 overwrite=True)
#             imsmooth(imagename=prefix12+linename[i]+'.image',
#                 outfile=prefix12+linename[i]+'.'+bmtyp+'.smooth',overwrite=True,
#                 targetres=True, major='2arcsec', minor='2arcsec', pa='0deg')
 
        # Regrid 7m+TP data to match 12m data:
        imregrid(imagename=prefix7+field+base7tp+linename[i]+'.pbcor',
             template=prefixint+field+base12m+'.image',
             output=prefixint+field+'.'+linename[i]+'.7tprgd', overwrite=True)
        # Apply 12m sensitivity to 7m+TP data:
        rmtables(prefixint+field+'.'+linename[i]+'.7tprgd.pb')
        immath(imagename=[prefixint+field+'.'+linename[i]+'.7tprgd',
                     prefixint+field+base12m+'.flux.pbcoverage'],
            expr='IM0*IM1',
            outfile=prefixint+field+'.'+linename[i]+'.7tprgd.pb')
        exportfits(imagename=prefixint+field+'.'+linename[i]+'.7tprgd.pb', 
            fitsimage=prefixint+field+'.'+linename[i]+'.7tprgd.pb.fits',
            dropdeg=True, velocity=True, overwrite=True)
        rmtables(prefixint+field+'.'+linename[i]+'.7tprgd')
        # Also write out 12m sensitivity and image for mosaic script:
        exportfits(imagename=prefixint+field+base12m+'.flux.pbcoverage', 
            fitsimage=prefixint+field+base12m+'.flux.pbcoverage.fits',
            dropdeg=True, velocity=True, overwrite=True)
        exportfits(imagename=prefixint+field+base12m+'.image', 
            fitsimage=prefixint+field+base12m+'.image.fits',
            dropdeg=True, velocity=True, overwrite=True)
        # Do the feathering:
        rmtables(prefixint+field+baseout+linename[i]+'.image')
        feather(imagename=prefixint+field+baseout+linename[i]+'.image',
                highres=prefixint+field+base12m+'.image',
                lowres=prefixint+field+'.'+linename[i]+'.7tprgd.pb')
        # Export to FITS:
        exportfits(imagename=prefixint+field+baseout+linename[i]+'.image', 
            fitsimage=prefixint+field+baseout+linename[i]+'.image.fits',
            dropdeg=True, velocity=True, overwrite=True)
        # Primary beam correction:
        impbcor(imagename=prefixint+field+baseout+linename[i]+'.image',
                pbimage=prefixint+field+base12m+'.flux.pbcoverage',
                outfile=prefixint+field+baseout+linename[i]+'.pbcor',
                overwrite=True)
        # Export to FITS:
        exportfits(imagename=prefixint+field+baseout+linename[i]+'.pbcor', 
            fitsimage=prefixint+field+baseout+linename[i]+'.pbcor.fits',
            dropdeg=True, velocity=True, overwrite=True)

