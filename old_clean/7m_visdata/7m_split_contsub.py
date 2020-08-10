from glob import glob

vis = glob('uid*.ms')
# Split visibilities + continuum subtraction
# For 7  where spw 16 is used
# Only for 12 CO 2-1
for i in vis:
	os.system("rm -rf "+i+".split16*")
	split(vis=i,spw="16",field="30_Dor*",outputvis=i+".split16")
	uvcontsub(vis=i+".split16",field="30_Dor*",fitspw="0:230.266~230.321GHz;230.346~230386GHz")


concat(glob('*.contsub'),'7m.contsub.concat.ms')
