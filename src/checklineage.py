from definitions import *

from lineage import read_lineage_tree
lin_tree_information=read_lineage_tree(lineage_tree_filename) # Read the lineage tree (in case it was previously created)

for x in lin_tree_information:
	print "Lineage Tree contains "+x
	cellat={}
	for y in lin_tree_information[x]:
		t=y/10**4
		if t not in cellat:
			cellat[t]=1
		else:
			cellat[t]+=1
	for t in cellat:
		print ' --> '+str(cellat[t])+" cells at "+str(t)
