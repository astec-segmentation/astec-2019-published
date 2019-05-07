from definitions import *

print 'Check Lineage File '+named_lineage_tree_filename
from lineage import read_lineage_tree


lin_tree_information_named=read_lineage_tree(named_lineage_tree_filename) # Read the lineage tree (in case it was previously created)
lin_tree_named=lin_tree_information_named['lin_tree']

lin_tree_information_post=read_lineage_tree(post_lineage_tree_filename) 
lin_tree_post=lin_tree_information_post['lin_tree']


for x in lin_tree_named:
	if x not in lin_tree_post:
		print 'Miss '+str(x)+' in post '
	elif lin_tree_named[x]!=lin_tree_post[x]:
		print 'For '+str(x)+'  '+str(lin_tree_named[x])+'!='+str(lin_tree_post[x])

for x in lin_tree_post:
	if x not in lin_tree_named:
		print 'Miss '+str(x)+' in named '
