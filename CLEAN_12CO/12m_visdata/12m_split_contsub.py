from glob import glob

vis = glob('uid*.ms')
# Split visibilities + continuum subtraction
# For 12m data, where spw 25 is used
# Only for 12 CO 2-1
for i in vis:
	os.system("rm -rf "+i+".split25*")
	split(vis=i,spw="25",field="30_Dor*",outputvis=i+".split25")
	uvcontsub(vis=i+".split25",field="30_Dor*",fitspw="0:230.266~230.321GHz;230.346~230386GHz")


concat(glob('*.contsub'),'12m.concat25.contsub')
