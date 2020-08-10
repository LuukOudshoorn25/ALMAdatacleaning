from glob import glob
import shutil
# To be executed in 7m_vis folder

mslist = ['uid___A002_Xe3da01_X18fa.ms',
          "uid___A002_Xe48598_Xafeb.ms",
	  "uid___A002_Xe3f5bb_X21c9.ms",
 	  "uid___A002_Xe407cf_X1887b.ms",
	  "uid___A002_Xe407cf_X15ea.ms",
	  "uid___A002_Xe45e29_X9839.ms"]

for ms in mslist:
    print(ms)
    if not os.path.exists(ms+".split20"):
        fpath = glob('../../2019*/science*/*/*/calibrated/'+ms+'*')[0]
        print('Splitting',ms)
        split(vis=fpath,spw="29",field="30_Dor*",outputvis=ms+".split29")
        print('UVcontsub',ms)
        uvcontsub(vis=ms+".split29",field="30_Dor*",fitspw="0:220.126~220.182GHz;220.207~220.247GHz")

concat(glob('*.contsub'),'13CO_12m.contsub.concat.ms')


