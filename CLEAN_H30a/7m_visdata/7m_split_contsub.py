from glob import glob
import shutil
# To be executed in 7m_vis folder

mslist = ['uid___A002_Xe1baa0_X8097.ms',
'uid___A002_Xe220f2_X4ea1.ms',
'uid___A002_Xe220f2_X5842.ms',
'uid___A002_Xe220f2_X62c1.ms',
'uid___A002_Xe220f2_X6b1e.ms',
'uid___A002_Xe230a1_X28f6.ms',
'uid___A002_Xe230a1_X32eb.ms',
'uid___A002_Xe247d0_Xeb.ms',
'uid___A002_Xe27761_X172c.ms',
'uid___A002_Xe1a561_X2fad.ms',
'uid___A002_Xe1baa0_X477e.ms',
'uid___A002_Xe1baa0_X7737.ms',
'uid___A002_Xe1baa0_X7d43.ms',
'uid___A002_Xe2ada9_X18144.ms',
'uid___A002_Xe2ada9_X187a8.ms',
'uid___A002_Xe2ada9_X18ea1.ms',
'uid___A002_Xe31981_X30e4.ms',
'uid___A002_Xe31981_X3cee.ms',
'uid___A002_Xe31981_X47c7.ms',
'uid___A002_Xe31981_Xebf5.ms',
'uid___A002_Xe32bed_X3f50.ms',
'uid___A002_Xe3da01_X16c5.ms',
'uid___A002_Xe3da01_X2352.ms',
'uid___A002_Xe3da01_X305b.ms',
'uid___A002_Xe31981_Xf7ea.ms.split.cal',
'uid___A002_Xe32bed_Xdce2.ms.split.cal',
'uid___A002_Xe32bed_Xe889.ms.split.cal',
'uid___A002_Xe37224_X3836.ms.split.cal',
'uid___A002_Xe37224_X461a.ms.split.cal',
'uid___A002_Xe37224_Xdfe9.ms.split.cal',
'uid___A002_Xe37224_Xe70a.ms.split.cal',
'uid___A002_Xe407cf_X17cba.ms',
'uid___A002_Xe407cf_X20de.ms',
'uid___A002_Xe407cf_X9d8c.ms',
'uid___A002_Xe407cf_Xbb9d.ms',
'uid___A002_Xe407cf_Xf987.ms']

for ms in mslist:
    print(ms)
    if not os.path.exists(ms+".split18"):
        fpath = glob('../../2019*/science*/*/*/calibrated/'+ms+'*')[0]
        print('Splitting',ms)
        split(vis=fpath,spw="18",field="30_Dor*",outputvis=ms+".split18")
for ms in mslist:
    if not os.path.exists(ms+".split18"):
        print(ms)    
        fpath = glob('../../2019*/science*/*/*/calibrated/'+ms+'*')[0]
        print('Splitting',ms)
        split(vis=fpath,spw="18",field="30_Dor*",outputvis=ms+".split18",datacolumn='data')
     

for ms in glob('*.split18'):
    print('UVcontsub',ms)
    split(vis=ms,field="30_Dor*",spw='0:231.47~232.55GHz',outputvis=ms+'.splitchan',datacolumn='data')


for ms in glob('*.split18.splitchan'):
    print('UVcontsub',ms)
    uvcontsub(vis=ms,field="30_Dor*",fitspw="*:231.47~231.6GHz;231.95~232.5GHz",fitorder=0)

# emission range: 231.677-231.742





