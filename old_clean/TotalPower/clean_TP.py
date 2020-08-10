from glob import glob

infiles            =  glob('*ms')
# Split visibilities into spw17, which contains the 12CO 
# J=2-1 for the TP data
for vis in infiles:
    split(vis,vis+'.split17',field='30*',spw='17')
concat(glob('*split17'),concatvis='TP.split17.ms')



os.system('rm -rf *.image')
sdimaging(
taskname           = "sdimaging"
infiles            =  'TP.split17.ms'
outfile            =  "TP.image"
overwrite          =  True
field              =  "30_Doradus"
spw                =  "*"
antenna            =  ""
scan               =  ""
intent             =  "OBSERVE_TARGET#ON_SOURCE"
mode               =  "velocity"
nchan              =  701
start              =  "282.715km/s"
width              =  "-80m/s"
veltype            =  "radio"
outframe           =  "lsrk"
gridfunction       =  "SF"
convsupport        =  6
truncate           =  -1
gwidth             =  -1
jwidth             =  -1
imsize             =  [200, 200]
cell               =  "3arcsec"
phasecenter        =  'J2000 05h38m32.0 -69d02m18.0'
projection         =  "SIN"
ephemsrcname       =  ""
pointingcolumn     =  "direction"
restfreq           =  "230.5380GHz"
stokes             =  "I"
minweight          =  0.1
brightnessunit     =  "Jy/beam"
clipminmax         =  False
)

imcontsub("TP.image",chans="0~30;670~699",linefile='TP.line',contfile='TP.final.image')
