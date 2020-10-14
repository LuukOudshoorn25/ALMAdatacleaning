### 13 CO

os.system('rm -rf *.image *last *log')
sdimaging(
infiles            =  ['ms_data_F1', 'ms_data_F2', 'ms_data_F3', 'ms_data_F4', 'ms_data_F5'],
outfile            =  "TP_13CO_notconcat.image",
overwrite          =  True,
field              =  "30_Doradus",
spw                =  "*",
antenna            =  "",
scan               =  "",
intent             =  "OBSERVE_TARGET#ON_SOURCE",
mode               =  "velocity",
nchan              =  1100,
start              =  "290km/s",
width              =  "-83.33333333m/s",
veltype            =  "radio",
outframe           =  "lsrk",
gridfunction       =  "SF",
convsupport        =  6,
truncate           =  -1,
gwidth             =  -1,
jwidth             =  -1,
imsize             =  [250, 250],
cell               =  "2.5arcsec",
phasecenter        =  'J2000 05h38m46.839 -69d04m08.354',
projection         =  "SIN",
ephemsrcname       =  "",
pointingcolumn     =  "direction",
restfreq           =  "220.39870060GHz",
stokes             =  "I",
minweight          =  0.1,
brightnessunit     =  "Jy/beam",
clipminmax         =  False
)

imcontsub("TP_13CO_notconcat.image",chans="30~150;950~1080",linefile='TP_13CO_notconcat.line.image',contfile='TP_13CO_notconcat.res.image',fitorder=1)
