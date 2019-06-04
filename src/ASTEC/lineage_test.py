import sys
import os
sys.path.append(os.path.dirname(__file__)+'/CommunFunctions/')

from ImageHandling import imread
from morpheme_lineage import Morpheme_lineage
from lineage import read_lineage_tree
from pkl_converter import pklDictToTimeOrderedDict


def extractTimeFromName(path, extension=['.inr','.inr.gz']):
	'''
	Returns the integer corresponding to the time associated to the image according to file name.
	Expected file name format is :
	<path>/<a>_<list>_<of>_<things>_t([0-9]+)[_<etc>].<extension>
	'''
	import re
	sub_path=path.split('/')[-1]
	for ext in extension:
		sub_path=sub_path.rstrip('.'+ext.lstrip('.'))
	fields=sub_path.split('_')

	pattern = re.compile("^t([0-9]+)+$")

	for f in fields:
		if pattern.match(f):
			return int(f[1:])
	return None

def labelsInImage(image, background=[0,1]):
	'''
	Returns a list which corresponds to the listing of labels that are contained in an image.
	This listing excludes the background label(s) (int or list or None).
	background value is set to [0, 1] by default.
	/!\ works for images in 8 or 16 bits ! 
	'''
	from numpy import histogram
	if type(background) != list:
		background = [background]

	bin_counts, bin_edges = histogram(image, bins=range(0,2**16))
	idx = [i for i, x in enumerate(bin_counts > 0) if x]

	listing=bin_edges[idx]

	return [x for x in listing.tolist() if not background.count(x)] 




def imageDict(path='./*'):
	'''
	Returns a dictionnary  whose keys are the time-points of image names that matches the specified path, 
	and corresponding values are the image names themeselves.
	'''	
	import glob
	dico={}
	for fileName in glob.glob(path):
		time = extractTimeFromName(fileName)
		dico[time] = fileName
	return dico


def pkl_lineage_test(pkl_lineage, image_dict, file_out=None, background=[0,1], light=0, verbose=1, title=None):
	'''
	'''
	return lineage_test(Morpheme_lineage("tmp", pklDictToTimeOrderedDict(pkl_lineage)), image_dict, file_out, background, light, verbose, title) 

def lineage_test(lineage, image_dict, file_out=None, background=[0,1], light=0, verbose=1, title=None):
	'''

	Morpheme_lineage( 'Ralph_seg_post_named', pklDictToTimeOrderedDict( read_lineage_tree(path_lineage_ralph_seg_post_named) ) )
	'''

	def subfinder(mylist, pattern, reverse=1):
	    pattern = set(pattern)
	    if reverse:
	    	return [x for x in mylist if x not in pattern]	
	    return [x for x in mylist if x in pattern]


	NL='\n'
	f=None
	if file_out:
		f = open(file_out, 'w')
	else:
		f = open('/dev/null','w')

	if title:
		print title
		f.write(title+NL)


	f.write(lineage_labels_unicity_test(lineage, file_out=None, verbose=verbose, title=None))

	report={}

	is_ok=True

	time_images=image_dict.keys()
	time_lineage=lineage.timepoints()

	report['missing lineage times']=subfinder(time_images,time_lineage) # time-points d'images inexistants dans le lignage
	report['missing image times']=subfinder(time_lineage,time_images) # time-points de lignage inexistants dans les images

	if report['missing lineage times']:
		if verbose:
			print "Missing lineage times"
		f.write("\n  Missing lineage times:\n")
		is_ok=False
		if verbose:
			print report['missing lineage times']
		for time in report['missing lineage times']:
		    f.write('    t%03d\n' % time)

	if report['missing image times']:
		if verbose:
			print "Missing image times"
		f.write("\n  Missing image times:\n")
		is_ok=False
		if verbose:
			print report['missing image times']
		for time in report['missing image times']:
		    f.write('    t%03d\n' % time)

	if not light:
		from pprint import pprint
		report['inconsistent number of cells (lineage,images)']={}
		report['missing lineage labels']={}
		report['missing image labels']={}


		for t,fname in image_dict.iteritems():
			if lineage.has_time(t):
				labels_image=labelsInImage(imread(fname), background=background)
				ncells_lineage = lineage.ncells(t)
				labels_lineage=lineage.labels_at_time(t)
				missing_image_labels=subfinder(labels_lineage, labels_image)
				missing_lineage_labels=subfinder(labels_image, labels_lineage)
				if missing_image_labels:
					report['missing image labels'][t]=missing_image_labels
				if missing_lineage_labels:
					report['missing lineage labels'][t]=missing_lineage_labels

				if ncells_lineage != len(labels_image):
					report['inconsistent number of cells (lineage,images)'][t]=(ncells_lineage,len(labels_image))

		#for t in report['missing image times']:
		#	report['inconsistent number of cells (lineage,images)'][t]=(lineage.ncells(t), None)


		if report['inconsistent number of cells (lineage,images)']:
			if verbose:
				print "Inconsistent number of cells (lineage,images)"
			f.write("\n  Inconsistent number of cells (lineage,images):\n")
			is_ok=False
			if verbose:
				pprint( report['inconsistent number of cells (lineage,images)'])
			for time, v in report['inconsistent number of cells (lineage,images)'].iteritems():
			    f.write('    t%03d: %s\n' % (time, str(v)))


		if report['missing lineage labels']:
			if verbose:
				print "Missing lineage labels (time: label_list)"
			f.write("\n  Missing lineage labels (time: label_list):\n")
			is_ok=False
			if verbose:
				pprint( report['missing lineage labels'])
			for time, v in report['missing lineage labels'].iteritems():
			    f.write('    t%03d: %s\n' % (time, str(v)))		
		if report['missing image labels']:
			if verbose:
				print "Missing image labels (time: label_list)"
			f.write("\n  Missing image labels (time: label_list):\n")
			is_ok=False
			if verbose:
				pprint( report['missing image labels'])
			for time, v in report['missing image labels'].iteritems():
			    f.write('    t%03d: %s\n' % (time, str(v)))		
		#if report['missing image labels']:

	if is_ok:
		if verbose:
			print "Given lineage and image list are consistent."
		f.write("\n  Given lineage and image list are consistent.\n")

	if f:
		f.close()

	return report



def pkl_lineage_labels_unicity_test(pkl_lineage, file_out=None, verbose=1, title=None):
	'''
	'''
	return lineage_labels_unicity_test(Morpheme_lineage("tmp", pklDictToTimeOrderedDict(pkl_lineage)), file_out, verbose, title) 

def lineage_labels_unicity_test(lineage, file_out=None, verbose=1, title=None):
	'''
	Controle de l'unicite de la mere de chaque cellule dans l'arbre de lignage (hypothese biologique qu'il n'y a pas de fusions de cellules).
	'''

	NL='\n'
	f=None
	if file_out:
		f = open(file_out, 'w')
	else:
		f = open('/dev/null','w')

	if title:
		print title
		f.write(title+NL)

	report=""

	is_ok=True

	# STUFF

	time_lineage=lineage.timepoints()

	inv_lineage={}
	fusing_labels=[]
	for time in time_lineage:
		lin_at_time=lineage.lineage[time]
		#keys=[]
		#values=[]
		for mother,daughters in lin_at_time.iteritems():
			if len(daughters):
				for daughter in daughters:
					t_daughter, l_daughter = daughter[0], daughter[1]
					if not inv_lineage.has_key(daughter):
						inv_lineage[daughter]=[]
					else:
						fusing_labels.append(daughter)
						is_ok = False
					inv_lineage[daughter].append((time,mother))

	if is_ok:
		if verbose:
			print "Given lineage is consistent for labels mother unicity."
		report=NL+"Given lineage is consistent for labels mother unicity."+NL
		f.write("\n  Given lineage is consistent for labels mother unicity.\n")
	else:
		report=NL+"Given lineage is unconsistent for labels mother unicity:"+NL
		for daughter in fusing_labels:
			if verbose:
				print "Label %s has non-unique mother: %s"%(str(daughter),str(inv_lineage[daughter]))
			report=report+("Label %s has non-unique mother: %s"%(str(daughter),str(inv_lineage[daughter])))+NL
			f.write("\n  Label %s has non-unique mother: %s \n"%(str(daughter),str(inv_lineage[daughter])))

	if f:
		f.close()

	return report
