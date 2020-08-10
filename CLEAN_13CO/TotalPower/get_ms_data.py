spw = 21 # 13 CO line
#### Copy and split all MS data for the 13CO line
from glob import glob
science_goals = ['X2180','X2188','X2190','X2198','X21a0']
member_ids    = ['X2186','X218e','X2196','X219e','X21a6']
import os
for iter_ in range(len(science_goals)):
    sc_goal = science_goals[iter_]
    member_id = member_ids[iter_]
    vislist = glob('../../2019*/science*'+sc_goal+'/*/*'+member_id+'/calibrated/*ms*')
    outputdir = './ms_data_F'+str(iter_+1)+'/'
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    for vis in vislist:
        outputname=vis.split('/')[-1]+'.split'+str(spw)
        #print(vis,outputdir+outputname)
        split(vis,outputdir+outputname,spw=spw, field='30*',datacolumn='data')

# Merge all data into one MS per field
dirs = ['ms_data_F'+str(w) for w in range(1,6)]
for iter_, dir_ in enumerate(dirs):
    os.chdir(dir_)
    vislist = glob('*ms*')
    concat(vislist,'concatvis_F'+str(iter_+1))
    os.chdir('../')

# Merge all
vislist = glob('./ms_data*/*ms*split*')
concat(vislist,'TP_13CO.ms')
