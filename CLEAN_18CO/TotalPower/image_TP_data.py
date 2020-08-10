os.system('rm -rf *.image* *log *last')
os.system('mv *image* ../to_remove')
sdimaging(
infiles            =  'TP_18CO.ms',
outfile            =  "TP_18CO.image",
overwrite          =  True,
field              =  "30_Doradus",
spw                =  "*",
antenna            =  "",
scan               =  "",
intent             =  "OBSERVE_TARGET#ON_SOURCE",
mode               =  "velocity",
nchan              =  150,
start              =  "285km/s",
width              =  "-500m/s",
veltype            =  "radio",
outframe           =  "lsrk",
gridfunction       =  "SF",
convsupport        =  6,
truncate           =  -1,
gwidth             =  -1,
jwidth             =  -1,
imsize             =  [180, 220],
cell               =  "3arcsec",
phasecenter        =  'J2000 05h38m43.0 -69d04m40.0',
projection         =  "SIN",
ephemsrcname       =  "",
pointingcolumn     =  "direction",
restfreq           =  "219.56035410GHz",
stokes             =  "I",
minweight          =  0.1,
brightnessunit     =  "Jy/beam",
clipminmax         =  False
)

imcontsub("TP_18CO.image",chans="0~15;135~149",linefile='TP_18CO.line.image',contfile='TP_18CO.res.image')


