#########################################################
## Image Quality Assessment (IQA) scripts for CASA
##
## Created: Feb. 2020
## Last modified: Aug. 2020
## Created: A. Hacar, IQA group
##
#########################################################
## Functions

## keywords
##
## FITSfile = path to FITS file
## convo_file = path to CASA image to be convolved
## beam_final = final beamsize (in arcsec)
## target_image = target image (e.g., Interferometric image)
## ref_image = reference image (e.g., TP image)
##

##-------------------------------------------------------
## Additional libraries

from scipy import fftpack
from astropy.io import fits
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import os
import copy
from numpy import inf

##-------------------------------------------------------
## Image manipulation

## Convert FITS files into CASA images
def fits2CASA(FITSfile):
	print FITSfile
	os.system("rm -rf "+  FITSfile+".image")
	importfits(fitsimage=FITSfile,imagename=FITSfile+'.image')

## same as fits2CASA but for a list of FITS
def fitslist2CASA(FITSfile):
	for i in FITSfile:
		print i
		os.system("rm -rf "+ i+".image")
		importfits(fitsimage=i,imagename=i+'.image')

## Convolve CASA image with a final resolution (beam_final)
def get_convo(convo_file,beam_final):
	imsmooth(imagename= convo_file,
		outfile= convo_file + '_conv' + str(beam_final),
		kernel='gauss',
		major=str(beam_final)+'arcsec',
		minor=str(beam_final)+'arcsec',
		pa='0deg',
		targetres=True)

## same as get_convo but for FITS files
def getFITS_convo(FITSfile,beam_final):
	# FITS into CASA
	convo_file = FITSfile+'.image'
	os.system("rm -rf "+convo_file)
	importfits(fitsimage=FITSfile,imagename=convo_file)
	# convolution
	imsmooth(imagename= convo_file,
		outfile= convo_file + '_conv' + str(beam_final),
		kernel='gauss',
		major=str(beam_final)+'arcsec',
		minor=str(beam_final)+'arcsec',
		pa='0deg',
		targetres=True)

## Convolve CASA image with a final resolution (beam_final)
def get_convo2target(convo_file,ref_image):
	# Get beam info from refence image
	hdr = imhead(ref_image,mode='summary')
	beam_major = hdr['restoringbeam']['major']
	beam_minor = hdr['restoringbeam']['minor']
	beam_PA = hdr['restoringbeam']['positionangle']
	# Convolution into the same beam as the reference
	os.system("rm -rf " + convo_file + '_conv' + str(round(beam_major['value'])))
	imsmooth(imagename= convo_file,
		outfile= convo_file + '_conv' + str(round(beam_major['value'])),
		kernel='gauss',
		major=beam_major,
		minor=beam_minor,
		pa=beam_PA,
		targetres=True)

## same as get_convo2target but for FITS
def getFITS_convo2target(convo_file,ref_image):
	# FITS into CASA
	convo_file = FITSfile+'.image'
	os.system("rm -rf "+convo_file)
	importfits(fitsimage=FITSfile,imagename=convo_file)
	# Get beam info from refence image
	hdr = imhead(ref_image,mode='summary')
	beam_major = hdr['restoringbeam']['major']
	beam_minor = hdr['restoringbeam']['minor']
	beam_PA = hdr['restoringbeam']['positionangle']
	# Convolution into the same beam as the reference
	os.system("rm -rf " + convo_file + '_conv' + str(round(beam_major['value'])))
	imsmooth(imagename= convo_file,
		outfile= convo_file + '_conv' + str(round(beam_major['value'])),
		kernel='gauss',
		major=beam_major,
		minor=beam_minor,
		pa=beam_PA,
		targetres=True)

## Convolve CASA image with a final resolution (beam_final)
def get_convo2target(convo_file,ref_image):
	# Get beam info from refence image
	hdr = imhead(ref_image,mode='summary')
	beam_major = hdr['restoringbeam']['major']
	beam_minor = hdr['restoringbeam']['minor']
	beam_PA = hdr['restoringbeam']['positionangle']
	# Convolution into the same beam as the reference
	os.system("rm -rf convo2ref")
	imsmooth(imagename= convo_file,
		outfile= "convo2ref",
		kernel='gauss',
		major=beam_major,
		minor=beam_minor,
		pa=beam_PA,
		targetres=True)

## Resample 
def resample_velaxis(image,template):
	os.system('rm -rf '+image+'_resvel')
	imregrid(imagename= image,
		 template= template,
		 axes=[2],
		 output= image+'_resvel')

# Mask data (typically reference image)
def mask_image(myimage,threshold=0.0,relative=False):
	# If Relative = False, threshold is taken as absolute value (see image flux units)
	# If relative = True, threshold value is measured as rms fraction (aka sigmas)
	print "========================================="
	print " mask_image(): masking image "
	print "========================================="
	print " Original image : " + str(myimage)
	# Create a copy of your image
	os.system('rm -rf masked.tmp')
	os.system('cp -r ' + str(myimage) + ' masked.tmp')
	# Create your mask
	ia.open('masked.tmp')
	if (relative == False):
		ia.calcmask(mask= 'masked.tmp >= '+str(threshold), name='mymask')
	if (relative == True):
		ima_sigma = imstat(myimage)["rms"][0]
		ia.calcmask(mask= 'masked.tmp >= '+str(threshold*ima_sigma), name='mymask')
	ia.done()
	os.system('mv masked.tmp ' + str(myimage) + '_masked')
	#makemask(mode='copy',inpimage='masked.tmp',inpmask=['masked.tmp:mymask'],output=str(myimage) + '_masked',overwrite=True)
	print " New masked image : " + str(myimage) + '_masked'
	print "-----------------------------------------"
	print " mask_image(): DONE "
	print "========================================="

# Check if images have the same axis order as reference
def check_axis(ref_image,target_im=[]):
	print "========================================="
	print " check_axis(): checking axis consistency "
	print "========================================="
	# reference: check axis 
	axis_ref = imhead(ref_image).get("axisnames")
	print " Reference image: " + str(ref_image)
	print " Axis : (",
	for j in np.arange(np.shape(axis_ref)[0]):
		print axis_ref[j] + " ",
	print ")"
	print "-----------------------------------------"
	# Targets: check axis
	n = 0
	for i in target_im:
		n = n+1
		axis_target = imhead(i).get("axisnames")
		print " Target image #" + str(n) + ": " + str(i)
		print " Axis : (",
		transpose = -1
		for j in np.arange(np.shape(axis_target)[0]):
			print axis_target[j] + " ",
			if (axis_ref[j] != axis_target[j]):
				transpose = j
		print ")"
		if (transpose != -1):
			print " WARNING: Axis do not match reference -> transpose OR drop axis"
			print "          (see also drop_axis function)"
		print "-----------------------------------------"
	print " check_axis(): DONE "
	print "========================================="



def drop_axis(ref_image):
	print "================================================="
	print " drop_axis(): drop additional axis (e.g. Stokes) "
	print "================================================="
	# reference: check axis 
	os.system("rm -rf " + ref_image + "_subimage")
	imsubimage(imagename=ref_image,outfile=ref_image + "_subimage",dropdeg=True)
	print " Reference image: " + str(ref_image)
	print " New image: " + str(ref_image) +"_subimage"
	print "-----------------------------------------"
	print " drop_axis(): DONE "
	print "========================================="

##-------------------------------------------------------
## Quality estimators
## see a detailed discussion in: https://library.nrao.edu/public/memos/ngvla/NGVLA_67.pdf

## Calculate Image Accuracy parameter (Apar)
## Apar = (image-reference)/abs(reference)
## image & reference have to have the smae resolution
## image will be resampled into reference frame
def image_Apar(image,ref_image):
	# Resampling
	os.system('rm -rf tmp_resampled')
	imregrid(imagename= image,
		 template= ref_image,
		 #axes=[0, 1],
		 output= 'tmp_resampled')
	os.system('rm -rf ' + image + '_Apar')
	# Q parameter
	immath(imagename=['tmp_resampled',ref_image], 
		outfile= image + '_Apar',
		expr='(IM0-IM1)/abs(IM1)')

# Same as image_Apar but for a list of images
# imagelist = [im1,im2,...,imN]  
# e.g., list of Int. images
def imagelist_Apar(imagelist,ref_image):
	# Loop over the list
	for image in imagelist:
		# Resampling
		os.system('rm -rf tmp_resampled')
		imregrid(imagename= image,
			 template= ref_image,
			 #axes=[0, 1],
			 output= 'tmp_resampled')
		os.system('rm -rf ' + image + '_Apar')
		# Q parameter
		immath(imagename=['tmp_resampled',ref_image], 
			outfile= image + '_Apar',
			expr='(IM0-IM1)/abs(IM1)')


## Calculate image Fidelity
## Fidelity = abs(reference)/abs(reference-image)
## image & reference have to have the smae resolution
## image will be resampled into reference frame
def image_Fidelity(image,ref_image):
	# Resampling
	os.system('rm -rf tmp_resampled')
	imregrid(imagename= image,
		 template= ref_image,
		 #axes=[0, 1],
		 output= 'tmp_resampled')
	# Fidelity parameter
	os.system('rm -rf ' + image + '_Fidelity')
	immath(imagename=['tmp_resampled',ref_image], 
		outfile= image + '_Fidelity',
		expr='abs(IM1)/abs(IM1-IM0)')

# Same as image_Fidelity but for a list of images
# imagelist = [im1,im2,...,imN]  
# e.g., list of Int. images
def imagelist_Fidelity(imagelist,ref_image):
	# Loop over the list
	for image in imagelist:
		# Resampling
		os.system('rm -rf tmp_resampled')
		imregrid(imagename= image,
			 template= ref_image,
			 #axes=[0, 1],
			 output= 'tmp_resampled')
		# Fidelity parameter
		os.system('rm -rf ' + image + '_Fidelity')
		immath(imagename=['tmp_resampled',ref_image], 
			outfile= image + '_Fidelity',
			expr='abs(IM1)/abs(IM1-IM0)')



##-------------------------------------------------------
## Wrappers

# IQA methods: Accuracy, Fidelity, etc...
def get_IQA(ref_image = '',target_image=['']):
	# Reference image
	print "============================================="
	print " get_IQA(): Obtain IQA estimators"
	print " Reference : "+ str(ref_image)
	print " Depending on the image/cube size, this process may take a while..."
	print "---------------------------------------------"
	# Target images
	for j in np.arange(0,np.shape(target_image)[0],1):
		# print file
		print " Target image " + str(j+1) + " : " + str(target_image[j])
		# Convolve data into reference resolution
		get_convo2target(target_image[j],ref_image)
		os.system("rm -rf " + target_image[j] + "_convo2ref")
		os.system("mv convo2ref " + target_image[j] + "_convo2ref")
		# Get Apar, Fidelity, etc... images
		image_Apar(target_image[j] + "_convo2ref",ref_image)
		image_Fidelity(target_image[j] + "_convo2ref",ref_image)
		# export into FITS
		os.system('rm -rf '+target_image[j] + "_convo2ref*.fits")
		exportfits(imagename=target_image[j] + "_convo2ref",fitsimage=target_image[j] + "_convo2ref.fits",dropdeg=True)
		exportfits(imagename=target_image[j] + "_convo2ref_Apar",fitsimage=target_image[j] + "_convo2ref_Apar.fits",dropdeg=True)
		exportfits(imagename=target_image[j] + "_convo2ref_Fidelity",fitsimage=target_image[j] + "_convo2ref_Fidelity.fits",dropdeg=True)		
		print "---------------------------------------------"
	print " IQA estimators... DONE"	
	print "============================================="

##-------------------------------------------------------
## Plot functions:

# Show Accuracy parameter comparisons
def Compare_Apar_cubes(ref_image = '',target_image=['']):
	# Reference image
	print "============================================="
	print " Display Accuracy parameter:"
	print " Reference : "+str(ref_image)
	print "---------------------------------------------"
	# Number of plots
	Nplots = np.shape(target_image)[0]
	# Generate figure
	plt.show()
	ysizeplots = 4.
	fig = plt.figure(figsize=(15,Nplots*ysizeplots))
	grid = plt.GridSpec(ncols=3,nrows=Nplots, wspace=0.1, hspace=0.2)
	# Create Q-plots
	for j in np.arange(0,Nplots,1):
		xvalues, yvalues = plot_Apar(image2plot=target_image[j]+"_convo2ref_Apar.fits",Nplots=Nplots,Ny=j,title="Target image"+str(j+1))
		# Store results for comparison
		if (j == 0):
			results = yvalues
		else:
			results = np.row_stack([results,yvalues])
		#
		print " Target image " + str(j+1) + " : " + str(target_image[j])
		print "---------------------------------------------"
	plt.savefig("Apar_tmp.png")
	# Global comparisons (all channels are considered together)
	plt.figure(figsize=(7,7))
	for m in np.arange(Nplots):
		print " Target image " + str(m+1) + " : " + str(target_image[m])
		nchans, b, mids, h = get_ALLvalues(FITSfile=target_image[m]+"_convo2ref_Apar.fits",xmin=-1.5,xmax=1.5,xstep=0.1)
		plt.plot(mids,h,label="Target"+str(m+1),linewidth=3)
		meanvalue = np.round(np.average(mids,weights=h),2)
		sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=h)),2)
		print " A-parameter = " + str(meanvalue) + " +/- " + str(sigmavalue)
	plt.vlines(0.,1.,np.max(h),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
	plt.xlim(-1.5,1.5)
	plt.yscale('log')	# Make y axis in log scale
	plt.ylim(1,)
	plt.legend()
	plt.title(" Accuracy parameter (global)",fontsize=25)
	plt.xlabel("Accuracy",fontsize=25)
	plt.ylabel(r'# pixels$',fontsize=20)
	plt.savefig("AparALL_tmp.png")
	print "---------------------------------------------"
	print " See results: Apar_tmp.png, AparALL_tmp.png"
	print " Display Accuracy parameter... DONE"
	print "============================================="

# Accuracy parameter comparisons
def Compare_Apar_continuum(ref_image = '',target_image=['']):
	# Reference image
	print "============================================="
	print " Accuracy parameter comparisons:"
	print " Reference : "+str(ref_image)
	print "---------------------------------------------"
	# Number of plots
	Nplots = np.shape(target_image)[0]
	# Global comparisons 
	plt.figure(figsize=(7,7))
	for m in np.arange(Nplots):
		print " Target image " + str(m+1) + " : " + str(target_image[m])
		nchans, b, mids, h = get_ALLvalues(FITSfile=target_image[m]+"_convo2ref_Apar.fits",xmin=-1.5,xmax=1.5,xstep=0.1)
		plt.plot(mids,h,label="Target"+str(m+1),linewidth=3)
		meanvalue = np.round(np.average(mids,weights=h),2)
		sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=h)),2)
		print " A-parameter = " + str(meanvalue) + " +/- " + str(sigmavalue) 
	plt.vlines(0.,1.,np.max(h),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
	plt.xlim(-1.5,1.5)
	plt.yscale('log')	# Make y axis in log scale
	plt.ylim(1,)
	plt.legend()
	plt.xlabel("A-par",fontsize=25)
	plt.ylabel(r'# pixels',fontsize=20)
	plt.title("Accuracy Parameter Comparisons")
	plt.savefig("AparALL_tmp.png")
	print "---------------------------------------------"
	print " See results: AparALL_tmp.png"
	print " Accuracy parameter comparisons... DONE"
	print "============================================="

# Show Fidelity comparisons (cubes)
def Compare_Fidelity_cubes(ref_image = '',target_image=['']):
	# Reference image
	print "============================================="
	print " Display Fidelity:"
	print " Reference : "+str(ref_image)
	print "---------------------------------------------"
	# Number of plots
	Nplots = np.shape(target_image)[0]
	# Generate figure
	plt.show()
	ysizeplots = 4.
	fig = plt.figure(figsize=(15,Nplots*ysizeplots))
	grid = plt.GridSpec(ncols=3,nrows=Nplots, wspace=0.1, hspace=0.2)
	# Create Q-plots
	for j in np.arange(0,Nplots,1):
		xvalues, yvalues = plot_Fidelity(image2plot=target_image[j]+"_convo2ref_Fidelity.fits",Nplots=Nplots,Ny=j,title="Target image"+str(j+1))
		# Store results for comparison
		if (j == 0):
			results = yvalues
		else:
			results = np.row_stack([results,yvalues])
		#
		print " Target image " + str(j+1) + " : " + str(target_image[j])
		print "---------------------------------------------"
	plt.savefig("Fidelity_tmp.png")
	# Global comparisons (all channels are considered together)
	plt.figure(figsize=(7,7))
	for m in np.arange(Nplots):
		nchans, b, mids, h = get_ALLvalues(FITSfile=target_image[m]+"_convo2ref_Fidelity.fits",xmin=0,xmax=100,xstep=0.5)
		plt.plot(mids,h,label="Target"+str(m+1),linewidth=3)
	plt.vlines(0.,0.,np.max(results),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
	plt.xlim(0,)
	plt.yscale('log')	# Make y axis in log scale
	plt.ylim(1,)
	plt.legend()
	plt.xlabel("Fidelity",fontsize=25)
	plt.ylabel(r'log$_{10}(\# pixels)$',fontsize=20)
	plt.savefig("FidelityALL_tmp.png")
	print "---------------------------------------------"
	print " See results: Fidelity_tmp.png, FidelityALL_tmp.png"
	print " Display Fidelity... DONE"
	print "============================================="

# Fidelity comparisons
def Compare_Fidelity_continuum(ref_image = '',target_image=['']):
	# Reference image
	print "============================================="
	print " Fidelity comparisons:"
	print " Reference : "+str(ref_image)
	print "---------------------------------------------"
	# Number of plots
	Nplots = np.shape(target_image)[0]
	# Global comparisons 
	plt.figure(figsize=(7,7))
	for m in np.arange(Nplots):
		print " Target image " + str(m+1) + " : " + str(target_image[m])
		nchans, b, mids, h = get_ALLvalues(FITSfile=target_image[m]+"_convo2ref_Fidelity.fits",xmin=0.,xmax=100.,xstep=0.5)
		plt.plot(mids,h,label="Target"+str(m+1),linewidth=3)
		meanvalue = np.round(np.average(mids,weights=h),2)
		print " Fidelity (mean only) = " + str(meanvalue) 
	plt.xlim(0,100.)
	plt.yscale('log')	# Make y axis in log scale
	plt.ylim(1,)
	plt.legend()
	plt.xlabel("Fidelity",fontsize=25)
	plt.ylabel(r'# pixels',fontsize=20)
	plt.title("Fidelity Comparisons")
	plt.savefig("FidelityALL_tmp.png")
	print "---------------------------------------------"
	print " See results: FidelityALL_tmp.png"
	print " Fidelity comparisons... DONE"
	print "============================================="

##-------------------------------------------------------
## Plot functions

# get values from FITS (ALL)
def get_ALLvalues(FITSfile,xmin,xmax,xstep):
	# FITS file
	image = fits.open(FITSfile)
	# Histogram
	bins_histo = np.arange(xmin,xmax,xstep)
	bins_mids = bins_histo[1:]-xstep/2.
	subimage = image[0].data.flatten()
	idxs = np.isfinite(subimage)
	hist, bin_edges = np.histogram(subimage[idxs],bins=bins_histo)
	## values in log scale
	##values_log = np.log10(hist.T)
	# return
	return 0. , bins_histo, bins_mids, hist.T

# get values from FITS for a given channel
def get_CHANvalues(FITSfile,xmin,xmax,xstep,channel):
	# FITS file
	image = fits.open(FITSfile)
	# Histogram
	bins_histo = np.arange(xmin,xmax,xstep)
	bins_mids = bins_histo[1:]-xstep/2.
	subimage = image[0].data[channel].flatten()
	idxs = np.isfinite(subimage)
	hist, bin_edges = np.histogram(subimage[idxs],bins=bins_histo)
	## values in log scale
	##values_log = np.log10(hist.T)
	# return
	return 0. , bins_histo, bins_mids, hist.T

# get values from FITS (per channel)
def get_values(FITSfile,xmin,xmax,xstep):
	# FITS file
	image = fits.open(FITSfile)
	# Histogram per channel
	nchan = image[0].shape[0]
	bins_histo = np.arange(xmin,xmax,xstep)
	bins_mids = bins_histo[1:]-xstep/2.	
    	for i in np.arange(0,nchan):
        	subimage = image[0].data[i,:,:].flatten()
		idxs = np.isfinite(subimage)
		hist, bin_edges = np.histogram(subimage[idxs],bins=bins_histo)
		# values in log scale
		values_log = np.log10(hist.T)
        	if (i == 0):
			Histogram = copy.deepcopy(values_log)	
        	else:
        		Histogram = np.column_stack((Histogram,values_log))
    	return nchan, bins_histo, bins_mids, Histogram




def plot_Apar(image2plot,Nplots,Ny,title):
	grid = plt.GridSpec(ncols=3,nrows=Nplots, wspace=0.1, hspace=0.2)
	# A-par plot
	xminplot = -1.5; xmaxplot = 1.5
	nchans, b, mids, h = get_values(FITSfile=image2plot,xmin=xminplot,xmax=xmaxplot,xstep=0.1)
	h[h == -inf] = np.nan
	# 2D plot per channel
	ax1 = plt.subplot(grid[Ny, :2])
	plt.title(title+" - Accuracy",fontsize=20,x=0,ha="left",position=(0.1,0.8))
	#plt.title("Q parameter",fontsize=10,x=1,ha="right",color='grey', style='italic')
	plt.imshow(h, extent =(-0.5,nchans-0.5,xminplot,xmaxplot), aspect='auto',vmin=0,cmap="jet", interpolation='none',origin='lower')
	plt.hlines(0,-0.5,nchans-0.5,linestyle="--",color="black",linewidth=3,label="Goal")	
	plt.xlabel("Channel number",fontsize=15)
	plt.ylabel("Accuracy",fontsize=18)
	h[np.isnan(h)] = 0.0
	cbar = plt.colorbar()
	cbar.set_label(r'log$_{10}(\# pixels)$',fontsize=15)
	# Histogram
	ax2 = plt.subplot(grid[Ny, 2])
	plt.hlines(0,0.1,max(h.flat),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
	# Get mean values (only orientative)
	for i in np.arange(0,nchans):
		plt.step(h[:,i],mids,color="grey",linewidth=1, alpha=0.1)
	##h[h == 0.0] = np.nan
	# We use linear weights rather than log10
	plt.step(np.log10(np.average(10.**h,axis=1)),mids,color="red",linewidth=3,label="mean")
	meanvalue = np.round(np.average(mids,weights=np.log10(np.average(10.**h,axis=1))),2)
	sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=np.log10(np.average(10.**h,axis=1)))),2)
	# The alternative would be
	#plt.step(np.nanmean(h,axis=1)[:],mids,color="red",linewidth=3,label="mean")
	#meanvalue = np.round(np.average(mids,weights=np.nanmean(h,axis=1)[:]),2)
	#sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=np.nanmean(h,axis=1)[:])),2)
	plt.title('{:.2f}'.format(meanvalue)+"+/-"+'{:.2f}'.format(sigmavalue),fontsize=15,x=0,ha="left",position=(0.6,0.1),color="blue")
	plt.hlines(meanvalue,0.1,max(h.flat),color="blue",linewidth=3,label="average",alpha=0.5,zorder=-2)
	# Plot parameters
	plt.xlim(0,max(h.flat))
	plt.ylabel("Accuracy",fontsize=18)
	plt.xlabel(r'log$_{10}(\# pixels)$',fontsize=15)
	plt.legend(loc="upper left",prop={"size":10})
	ax2.tick_params(labelbottom=True, labelleft=False, labelright=True,bottom=True, top=True, left=True, right=True)
	ax2.yaxis.set_label_position("right")
	# return for comparisons
	return mids, np.nanmean(h,axis=1)[:]

def plot_Fidelity(image2plot,Nplots,Ny,title):
	grid = plt.GridSpec(ncols=3,nrows=Nplots, wspace=0.1, hspace=0.2)
	# Fidelity plot
	xminplot = -1; xmaxplot = 100.
	nchans, b, mids, h = get_values(FITSfile=image2plot,xmin=xminplot,xmax=xmaxplot,xstep=1)
	h[h == -inf] = np.nan
	# 2D plot per channel
	ax1 = plt.subplot(grid[Ny, :2])
	plt.title(title+" - Fidelity",fontsize=20,x=0,ha="left",position=(0.1,0.8))
	#plt.title("Q parameter",fontsize=10,x=1,ha="right",color='grey', style='italic')
	plt.imshow(h, extent =(-0.5,nchans-0.5,xminplot,xmaxplot), aspect='auto',vmin=0,cmap="jet", interpolation='none',origin='lower')
	#plt.hlines(0,-0.5,nchans-0.5,linestyle="--",color="black",linewidth=3,label="Goal")	
	plt.xlabel("Channel number",fontsize=15)
	plt.ylabel("Fidelity",fontsize=18)
	h[np.isnan(h)] = 0.0
	cbar = plt.colorbar()
	cbar.set_label(r'log$_{10}(\# pixels)$',fontsize=15)
	# Histogram
	ax2 = plt.subplot(grid[Ny, 2])
	plt.hlines(0,0.1,max(h.flat),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
	# Get mean values (only orientative)
	for i in np.arange(0,nchans):
		plt.step(h[:,i],mids,color="grey",linewidth=1, alpha=0.1)
	##h[h == 0.0] = np.nan
	# We use linear weights rather than log10
	plt.step(np.log10(np.average(10.**h,axis=1)),mids,color="red",linewidth=3,label="mean")
	meanvalue = np.round(np.average(mids,weights=np.log10(np.average(10.**h,axis=1))),2)
	sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=np.log10(np.average(10.**h,axis=1)))),2)
	# The alternative would be
	#plt.step(np.nanmean(h,axis=1)[:],mids,color="red",linewidth=3,label="mean")
	#meanvalue = np.round(np.average(mids,weights=np.nanmean(h,axis=1)[:]),2)
	#sigmavalue = np.round(np.sqrt(np.cov(mids, aweights=np.nanmean(h,axis=1)[:])),2)
	plt.title('{:.2f}'.format(meanvalue)+"+/-"+'{:.2f}'.format(sigmavalue),fontsize=15,x=0,ha="left",position=(0.6,0.1),color="blue")
	plt.hlines(meanvalue,0.1,max(h.flat),color="blue",linewidth=3,label="average",alpha=0.5,zorder=-2)
	# Plot parameters
	plt.xlim(0,max(h.flat))
	plt.ylabel("Fidelity",fontsize=18)
	plt.xlabel(r'log$_{10}(\# pixels)$',fontsize=15)
	plt.legend(loc="upper left",prop={"size":10})
	ax2.tick_params(labelbottom=True, labelleft=False, labelright=True,bottom=True, top=True, left=True, right=True)
	ax2.yaxis.set_label_position("right")
	# return for comparisons
	return mids, np.nanmean(h,axis=1)[:]

def show_Apar_map(target_image,channel=0):
	plt.figure(figsize=(15,7))
	grid = plt.GridSpec(ncols=2,nrows=1, wspace=0.2, hspace=1.0)
	# Map
	ax1 = plt.subplot(grid[0, 0])
	image = fits.open(target_image+"_convo2ref_Apar.fits")
	# Number of axis = Dimentions (Cont vs cubes)
	Ndims = np.shape(np.shape(image[0].data))
	if (Ndims[0] == 2):
		# Continuum
		im = ax1.imshow(image[0].data,vmin=-1.5,vmax=1.5,cmap='bwr_r')
	else:
		# Cubes
		im = ax1.imshow(image[0].data[channel],vmin=-1.5,vmax=1.5,cmap='bwr_r')
	cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
	cbar.ax.set_ylabel('Accuracy parameter', fontsize=15)
	plt.gca().invert_yaxis()
	plt.show()
	plt.xlabel("X (pixel units)",fontsize=15)
	plt.ylabel("Y (pixel units)",fontsize=15)
	plt.title(" Accuracy map (Channel # " + str(channel) + ")")
	# Histogram
	ax2 = plt.subplot(grid[0, 1])
	nchans, b, mids, h = get_ALLvalues(FITSfile=target_image+"_convo2ref_Apar.fits",xmin=-1.5,xmax=1.5,xstep=0.1)
	plt.plot(mids,h,label="ALL channels",linewidth=3,c="red")
	if (Ndims[0] > 2): # cubes only
		nchans_chan, b_chan, mids_chan, h_chan = get_CHANvalues(FITSfile=target_image+"_convo2ref_Apar.fits",xmin=-1.5,xmax=1.5,xstep=0.1,channel=channel)
		plt.plot(mids_chan,h_chan,label="Channel # " +str(channel),c="blue",linewidth=3,linestyle="dotted")
	plt.vlines(0.,1.,np.max(h),linestyle="--",color="black",linewidth=3,label="Goal",alpha=1.,zorder=-2)
	plt.xlim(-1.5,1.5)
	plt.yscale('log')	# Make y axis in log scale
	plt.ylim(1,)
	plt.legend()
	plt.xlabel("Accuracy parameter",fontsize=20)
	plt.ylabel(r'Number of pixels pixels',fontsize=20)

def show_Fidelity_map(target_image,channel=0):
	plt.figure(figsize=(15,7))
	grid = plt.GridSpec(ncols=2,nrows=1, wspace=0.2, hspace=1.0)
	# Map
	ax1 = plt.subplot(grid[0, 0])
	image = fits.open(target_image+"_convo2ref_Fidelity.fits")
	# Number of axis = Dimentions (Cont vs cubes)
	Ndims = np.shape(np.shape(image[0].data))
	if (Ndims[0] == 2):
		# Continuum
		im = ax1.imshow(image[0].data,vmin=0.,vmax=100.,cmap='jet')
	else:
		# Cubes
		im = ax1.imshow(image[0].data[channel],vmin=0.,vmax=100.,cmap='jet')
	cbar = plt.colorbar(im, ax=ax1,orientation='vertical')
	cbar.ax.set_ylabel('Fidelity', fontsize=15)
	plt.gca().invert_yaxis()
	plt.show()
	plt.xlabel("X (pixel units)",fontsize=15)
	plt.ylabel("Y (pixel units)",fontsize=15)
	plt.title(" Fidelity map (Channel # " + str(channel) + ")")
	# Histogram
	ax2 = plt.subplot(grid[0, 1])
	nchans, b, mids, h = get_ALLvalues(FITSfile=target_image+"_convo2ref_Fidelity.fits",xmin=0.,xmax=100.,xstep=0.5)
	plt.plot(mids,h,label="ALL channels",linewidth=3,c="red")
	if (Ndims[0] > 2): # cubes only
		nchans_chan, b_chan, mids_chan, h_chan = get_CHANvalues(FITSfile=target_image+"_convo2ref_Fidelity.fits",xmin=0.,xmax=100.,xstep=0.5,channel=channel)
		plt.plot(mids_chan,h_chan,label="Channel # " +str(channel),c="blue",linewidth=3,linestyle="dotted")
	plt.xlim(0.,100.)
	plt.yscale('log')	# Make y axis in log scale
	plt.ylim(1,)
	plt.legend()
	plt.xlabel("Fidelity",fontsize=20)
	plt.ylabel(r'Number of pixels pixels',fontsize=20)



##-------------------------------------------------------
## Analysis
def get_data():
	os.system('rm -rf *.image')

	os.system("cp -r /allegro1/allegro/data/projects/aGqnEW3z/analysis/oudshoorn/CLEAN_13CO/Image_7m_data/TP_13CO_regrid.image ./30Dor_TPalone.image")
	TP_image = "30Dor_TPalone.image"
	#importfits(fitsimage='',imagename='')

	importfits(fitsimage='../30Dor_13CO.cube.fits',imagename="30Dor_12m_feather.image")
	Feather_image = "30Dor_12m_feather.image"

	importfits(fitsimage='../30Dor_13CO_pbcor.cube.fits',imagename="30Dor_12m_feather+pbcor.image")
	Feather_pbco_image = "30Dor_12m_feather+pbcor.image"

	importfits(fitsimage='../30Dor_13CO_hybrid.cube.fits',imagename="30Dor_12m_hybrid.image")
	Hybrid_image = "30Dor_12m_hybrid.image"

	importfits(fitsimage='../30Dor_13CO_hybrid_pbcor.cube.fits',imagename="30Dor_12m_hybrid+pbcor.image")
	Hybrid_pbco_image = "30Dor_12m_hybrid+pbcor.image"

	return TP_image, Feather_image, Feather_pbco_image, Hybrid_image, Hybrid_pbco_image

TP_image, Feather_image, Feather_pbco_image, Hybrid_image, Hybrid_pbco_image = get_data()

## Data: paths to images
## TP_image = "myTP.image"
## Feather_image = "ALMAalone.image"
## Feather_pbco_image = "ALMAfeather.image"
## Hybrid_image = "Hybrid.image"

## STEP 1a: Check axis
## Example:
#check_axis(ref_image=TP_image,target_im=[Hybrid_pbco_image,Hybrid_image])

## STEP 1b: Drop image degeneracies (optional)
## Example:
drop_axis(ref_image=TP_image)
drop_axis(ref_image=Feather_image)
drop_axis(ref_image=Feather_pbco_image)
drop_axis(ref_image=Hybrid_image)
drop_axis(ref_image=Hybrid_pbco_image)
## check again
#check_axis(ref_image=TP_image+"_subimage",target_im=[Feather_image+"_subimage",Feather_pbco_image+"_subimage"])

# OR
## STEP 1c: Transpose cubes in case of axis mismatch (optional)
## Example:
#os.system("rm -rf " + str(TP_image)+ "_T")
#imtrans(imagename=TP_image,outfile=TP_image+"_T",order="0132")

## STEP 2: Mask reference image
## Example:
mask_image(str(TP_image)+"_subimage",threshold=5.0,relative=False)

## STEP 2: get ALL Image Quality Assesments
## a) Make all images in the same resolution + grid + beam as the reference
## b) calculates all the IQAs: Accuracy, Fidelity, etc...
## IMPORTANT: some images may give errors if they were opened in viewer in the same CASA session (restart CASA if necessary)
## This may take a while...
## Execute only once! (it can be skipped if images were produced before)
## Example:
get_IQA(ref_image = TP_image+"_subimage_masked",target_image=[Feather_image+"_subimage",Feather_pbco_image+"_subimage",Hybrid_pbco_image+"_subimage",Hybrid_image+"_subimage"])

## STEP 3: Display IQA parameters
## (1) Accuracy parameter, (2) Fidelity, ...
## Example:
Compare_Apar_cubes(ref_image = TP_image+"_subimage_masked",target_image=[Feather_image+"_subimage",Feather_pbco_image+"_subimage",Hybrid_pbco_image+"_subimage",Hybrid_image+"_subimage"])



