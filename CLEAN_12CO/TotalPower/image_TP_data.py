os.system('rm -rf *.image')
sdimaging(
infiles            =  ['ms_data_F1','ms_data_F2','ms_data_F3','ms_data_F4','ms_data_F5'],
outfile            =  "TP_12CO_new.image",
overwrite          =  True,
field              =  "30_Doradus",
spw                =  "*",
antenna            =  "",
scan               =  "",
intent             =  "OBSERVE_TARGET#ON_SOURCE",
mode               =  "velocity",
nchan              =  300,
start              =  "282.715km/s",
width              =  "-250m/s",
veltype            =  "radio",
outframe           =  "lsrk",
gridfunction       =  "SF",
convsupport        =  6,
truncate           =  -1,
gwidth             =  -1,
jwidth             =  -1,
imsize             =  [280, 320],
cell               =  "3arcsec",
phasecenter        =  'J2000 05h38m32.0 -69d02m18.0',
projection         =  "SIN",
ephemsrcname       =  "",
pointingcolumn     =  "direction",
restfreq           =  "230.5380GHz",
stokes             =  "I",
minweight          =  0.1,
brightnessunit     =  "Jy/beam",
clipminmax         =  False
)

imcontsub("TP_12CO_new.image",chans="0~15;260~299",linefile='TP_12CO_new.line.image',contfile='TP_12CO_new.res.image')
imsubimage(imagename='TP_12CO_new.line.image',box="45,22,180,202",outfile="TP_12CO_new.line.subim")

