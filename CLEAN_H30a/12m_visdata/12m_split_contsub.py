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
        #fpath = glob('../../2019*/science*/*/*/calibrated/'+ms+'*')[0]
        #print('Splitting',ms)
        #split(vis=fpath,spw="27",field="30_Dor*",outputvis=ms+".split27")
        #print('UVcontsub',ms)
    uvcontsub(vis=ms+".split27",field="30_Dor*",fitspw="*:230.688~231.5GHz;231.95~232.687GHz",fitorder=1)








#concat(glob('*.contsub'),'H30a_12m.contsub.concat.ms')

for ms in glob('*.split27'):
    print('UVcontsub',ms)
    split(vis=ms,field="30_Dor*",spw='0:231.47~232.55GHz',outputvis=ms+'.splitchan',datacolumn='data')


for ms in glob('*.split27.splitchan'):
    print('UVcontsub',ms)
    uvcontsub(vis=ms,field="30_Dor*",fitspw="*:231.47~231.6GHz;231.95~232.5GHz",fitorder=0)
