from definitions import *

from lineage import timeNamed,timesNamed

#RENAME FUSED DATA
for t in range(begin,end+1):
	os.system('mv '+fuse_Path+'fused_t'+str(t).zfill(3)+'.inr '+timeNamed(fused_files,t))

#RENAME SEGMENTED DATA
for t in range(begin,end+1):
	os.system('mv '+segmented_Path+'Segmented_t'+str(t).zfill(3)+'.inr '+timeNamed(segmentation_files,t))



