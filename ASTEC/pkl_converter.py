
from lineage import read_lineage_tree
from morpheme_lineage import versionWriter, cmp_cell_tuplet, getEmbryoKey, getTimeKey, getIDKey, writeTextLineage, control_persistency

#version='1.0'

#def version():
#	'''
#	Version du convertisseur de format
#	'''
#	return '1.0'

def getLongIDFromTimeAndLabel(time, label):
	'''
	Usage : 
		getLongIDFromTimeAndLabel(14, 112) returns 140112
	The LongID format an integer begining with the time, and followed by the label encoded with four numbers.
	'''
	return int(label + time*10000)

def getIDFromLongID( longID):
	'''
	Extrait d'un longID la composante label selon la convention que 
	longID=<composante temps><composante label> avec <composante label> ecrit sur 4 chiffres.
	'''
	t = round(longID / 10000);
	return int(longID - t * 10000);
def getTimeFromLongID ( longID):
	'''
	Extrait d'un longID la composante temps selon la convention que 
	longID=<composante temps><composante label> avec <composante label> ecrit sur 4 chiffres.
	'''
	return  int(round(longID / 10000));


def getTimeKeyFromLongID(longID):
	"""
	Usage : 
		getTimeKeyFromLongID(1560209)	returns 'T 156'
	"""
	return getTimeKey(getTimeFromLongID(longID))
	#return 'T ' + str(getTimeFromLongID(longID))

def getIDKeyFromLongID(longID):
	"""
	Usage : 
		getIDKeyFromLongID(1560209)	returns 'L 209'
	"""
	return getIDKey(getIDFromLongID(longID))
	#return 'L' + str(getIDFromLongID(longID))



def getKeyFromLongID(longID,embryoName=''):
	"""
	Usage : 
		getKeyFromLongID(1560209)	returns 'T 156, L 209'
		getKeyFromLongID(1560209,'161213_lucie')	returns 'E 161213_lucie, T 156, L 209'
	"""
	key=""
	sub_key=getEmbryoKey(embryoName)
	if sub_key:
		if key:
			key+=','+sub_key
		else:
			key=sub_key

	sub_key=getTimeKeyFromLongID(longID)
	if sub_key:
		if key:
			key+=','+sub_key
		else:
			key=sub_key

	sub_key=getIDKeyFromLongID(longID)
	if sub_key:
		if key:
			key+=','+sub_key
		else:
			key=sub_key

	return '['+key+']'





def pklDictToTimeOrderedDict(lineageDict):
	'''
	Convertisseur de dictionnaire de lineages en dictionnaire de timepoint:[mother1:[(time11,daughter11), (time12,daughter12)], mother2:[(time2,daughter2)]), ...] 
	avec en convention :
		Soit (m,D) tel que lineageDict.has_key(m), alors 
			out.has_key(getTimeFromLongID(m))
			out[getTimeFromLongID(m)].has_key(getIDFromLongID(m))
			type(out[getTimeFromLongID(m)][getIDFromLongID(m)]) == 'list'	(liste de tuplets de la forme (time_daughter, ID_daughter))		
			Pour tout tuplet tup dans D, out[getTimeFromLongID(m)][getIDFromLongID(m)].count(tup) == 1
		Soit t, id, (t1,d1), (t2,d2) tel que (t1,d1) et (t2,d2) existent dans out[t][id]
			Alors (t1,d1) < (t2,d2) (au sens de la fonction "cmp_cell_tuplet") <=> out[t][id].index((t1,d1)) < out[t][id].index((t2,d2))
			ie les enfants d'un element sont tries selon l'instant qui leur est associe puis selon le label (ou identifiant) qui leur est associe
	#Proprietes verifiees en sortie :
	#	Si (mother, [...]) est dans dict[timepoint], alors getTimeFromLongID(mother)=timepoint 
	#	Il existe au plus une occurence de mother dans l'ensemble des elements (mother, ~) contenus dans le dictionnaire (ie <out>[timepoint].index(mother)<=1 pour tout mother et tout timepoint)
	#	Soient mother1, mother2 tels que getTimeFromLongID(mother1)==getTimeFromLongID(mother2). Alors 
	#la liste associee est ordonnee selon l'ID de la mere 
	'''
	D={}

	if type(lineageDict) == list: # 140317-Patrick-St8_fuse_seg_post_lin_tree_named_byleo.pkl
		lin_tree=lineageDict[0]
		lineageDict = lineageDict[1]
	elif lineageDict.has_key('lin_tree'): # 140317-Patrick-St8_fuse_seg_post_lin_tree_named.pkl
		lin_tree=lineageDict['lin_tree']
	else:
		lin_tree=lineageDict

	dico={}


	for mother,daughters in lin_tree.iteritems():
		time=getTimeFromLongID(mother)
		if not dico.has_key(time):
			dico[time]={}
		ID=getIDFromLongID(mother)
		if not dico[time].has_key(ID):
			dico[time][ID]=[]
		for daughter in daughters:
			dico[time][ID].append((getTimeFromLongID(daughter),getIDFromLongID(daughter)))
			d_time=getTimeFromLongID(daughter)
			d_ID=getIDFromLongID(daughter)
			if not dico.has_key(d_time):
				dico[d_time]={}
			if not dico[d_time].has_key(d_ID):
				dico[d_time][d_ID]=[]
		dico[time][ID].sort(cmp_cell_tuplet)

	D['lin_tree']=dico


	if type(lineageDict)==dict:
		import numbers
		for key in lineageDict.keys():
			pkl_dico=lineageDict[key]
			if type(pkl_dico)==list and len(pkl_dico) == 2: # cas particulier embryon Patrick...
				pkl_dico=pkl_dico[0]
			dico={}
			if not key == 'lin_tree':
				#if key == 'cell_cell_surface_information': #/!\ Embryon Ralph : problemes d'identifiants ?
				#	dico = pkl_dico
				#else:

				for k,v in pkl_dico.iteritems():
					if isinstance(k, numbers.Number):
						time=getTimeFromLongID(k)
						if not dico.has_key(time):
							dico[time]={}
						ID=getIDFromLongID(k)
						if not dico[time].has_key(ID):
							dico[time][ID]={}
						if type(v)==dict:
								# e.g. case for key == 'surface'
							for k2, v2 in v.iteritems():
								if isinstance(k2, numbers.Number):
									try:
										assert getTimeFromLongID(k2) == time
									except:
										print key + ' ' + str(time) + ' ' + str( ID) + ' : unable to get time from long id from ' + str(k2)
										return
									ID2=getIDFromLongID(k2)
									try:
										dico[time][ID][ID2]=v2
									except:
										print key + ' ' + str(time) + ' ' + str( ID) + ' ' + str(k2) + ' ' + str(ID2)
										return None
								else:
									dico[time][ID][k2]=v2
						else:
							# other cases
							dico[time][ID] = v
					else:
						dico[k]=v
				D[key]=dico
	'''
	volumes={}

	if lineageDict.has_key('volumes_information'):
		lin_volumes=lineageDict['volumes_information']

		for cell_key, v in lin_volumes.iteritems():
			time=getTimeFromLongID(cell_key)
			if not volumes.has_key(time):
				volumes[time]={}
			ID=getIDFromLongID(cell_key)
			if not volumes[time].has_key(ID):
				volumes[time][ID]=v

		D['volumes_information']=volumes
	'''
	#if not control_persistency(D):
	#	print "Warning: given lineage from pkl format may contain inconsistent information."

	return D




def pklToTextFile(pklFileNameIn, textFileNameOut, embryoName='', version=versionWriter):
	''' 
	Convertisseur de fichiers de lineages cellulaires du format "pkl" a un format texte respectant les conventions morpheme 
	'''
	
	lin_tree=read_lineage_tree(pklFileNameIn)

	#if version == '0.0':
	#	writeTextLineage(lin_tree, textFileNameOut, embryoName,version)
	if version == '1.0':
		writeTextLineage(pklDictToTimeOrderedDict(lin_tree),textFileNameOut, embryoName, version)
	else:
		print 'Version '+ version +' unavailable or not implemented: returning\n'
		return 


	'''
	lineageText="# Version " + version + "\n"
	
	clique_opener='{'
	clique_closer='}'
	element_list_separator=','
	end_correspondance='\n'

	if version=='0.0':
		for mother, daughters in lin_tree.iteritems():
			lineageText+= clique_opener + getKeyFromLongID(mother,embryoName=embryoName) + clique_closer + ' - ' + clique_opener
			is_first_daughter=True
			for daughter in daughters:
				if is_first_daughter:
					is_first_daughter=False
				else:
					lineageText+=element_list_separator
				lineageText+=getKeyFromLongID(daughter,embryoName=embryoName)
			lineageText+=clique_closer + end_correspondance
	elif version=='1.0':
		lineage=pklDictToTimeOrderedDict(lin_tree)
		for time,IDs in lineage.iteritems():
			for ID,daughters in IDs.iteritems():
				lineageText+= clique_opener + getKey(time=time, ID=ID, embryoName=embryoName) + clique_closer + ' - ' + clique_opener
				is_first_daughter=True
				for daughter in daughters[:]:
					if is_first_daughter:
						is_first_daughter=False
					else:
						lineageText+=element_list_separator
					lineageText+=getKey(time=daughter[0],ID=daughter[1],embryoName=embryoName)
				lineageText+=clique_closer + end_correspondance
	else:
		lineageText+='Version unavailable or not implemented\n'			
	
	file=open(textFileNameOut, "w")
	file.write(lineageText)
	file.close()

	'''