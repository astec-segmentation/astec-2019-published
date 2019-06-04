
import os



###### DEFINE LINEAGE FILE SYNTAXES ######
versionWriter='2.0'
clique_opener='{'
clique_closer='}'
cliques_separator=' - '
element_list_separator=','
end_correspondance='\n'
key_separator=','
key_opener='['
key_closer=']'
key_time='T'
key_ID='L'
key_embryoName='E'
key_cellName='N'
key_space=' '
comment="#"



class Morpheme_lineage:
	"""Morpheme lineage class"""
	def __init__(self, name="", lineage={}):
		'''
		name : embryo name
		lineage : ordered lineage dictionnary given after pkl conversion or path to pkl file
		
		E.g:

		ml_ralph=Morpheme_lineage( 'Ralph', pklDictToTimeOrderedDict( read_lineage_tree('160707-Ralph-St8_fuse_seg_post_lin_tree_named.pkl') ) )
		
		'''
		if type(lineage) == str:
			from lineage import read_lineage_tree
			from pkl_converter import pklDictToTimeOrderedDict 
			lineage=pklDictToTimeOrderedDict( read_lineage_tree(lineage) )

		if type(lineage) == dict:
			if not lineage.has_key('lin_tree'):
				print "Lineage expected to have a 'lin_tree' key."
				return

			self.info = lineage
			self.lineage = self.info['lin_tree']
			#if self.info.has_key('volumes_information'):
			#	self.volumes=lineage['volumes_information']
		else:
			print "Lineage type not yield yet."
			return

		self.lineage_relabelling, self.luts = lineageRelabelling(self.lineage)
		self.reverse_luts = reverse_luts(self.luts)
		self.labels = to_label_classification(self.lineage_relabelling)

		self.name = name

	#@name.getter
	def get_name(self):
		return self.name


	#@lineage.getter
	def get_lineage(self):
		return self.lineage

	#@lineage_relabelling.getter
	def get_lineage_relabelling(self):
		return self.lineage_relabelling

	#@luts.getter
	def get_luts(self):
		return self.luts

	#@labels.getter
	def get_labels(self):
		return self.labels

	#@info.getter
	def get_info(self):
		return self.info

	#@name.setter
	def set_name(self, name):
		self.name=name


	#@lineage.setter
	def set_lineage(self, lineage):
		"""
		/!\ Modifies a priori the following attributes : 'lineage', 'lineage_relabelling', 'luts', 'labels' /!\
		"""
		self.lineage=copy.deepcopy(lineage)
		self.lineage_relabelling, self.luts = lineageRelabelling(lineage)
		self.labels = to_label_classification(self.lineage_relabelling)
		self.reverse_luts = reverse_luts(self.luts)


	def copy(self):
		'''
		Returns a deep copy of the object.
		'''
		import copy
		copie_lineage=Morpheme_lineage()
		copie_lineage.name=self.name
		copie_lineage.lineage=copy.deepcopy(self.lineage)
		copie_lineage.lineage_relabelling=copy.deepcopy(self.lineage_relabelling)
		copie_lineage.labels=copy.deepcopy(self.labels)
		copie_lineage.luts=copy.deepcopy(self.luts)
		copie_lineage.reverse_luts=copy.deepcopy(self.reverse_luts)
		copie_lineage.info=copy.deepcopy(self.info)
		return  copie_lineage


	def ncells(self, time=None):
		'''
		Returns lineage tree number of labels at specified time(s).
		time can be a integer or a vector of integers
		'''
		Time=[]
		res=[]
		if time==None:
			Time=self.lineage.keys()
		elif type(time) == int:
			Time=[time]
		elif type(time)==list:
			Time=time
		else:
			print "Type d'entree incorrect : seulement du 'int' ou du 'list' ou None"

		for t in Time:
			if not self.lineage.has_key(t):
				res.append(None)
			else :
				res.append(len(self.lineage[t].keys()))

		if len(res)>1:
			return res
		elif len(res)==1:
			return res[0]
		return None

	def timepoints(self):
		'''
		Returns the list of time-points of the lineage
		'''
		return self.lineage.keys()

	def tp(self, index):
		'''
		Returns the lineage time-point at the corresponding index
		'''
		return self.lineage.keys()[index]

	def labels_at_time(self, time):
		'''
		Returns the list of existing labels at given time.
		'''
		return self.lineage[time].keys()

	def has_time(self, time):
		'''
		Returns True if the time exists in the lineage, False otherwise.
		'''
		return self.lineage.has_key(time)

	def has_label_at_time(self, label, time):
		'''
		Returns True if the label exists at time in the attribute lineage, False otherwise.
		'''
		return self.has_time(time) and self.lineage[time].has_key(label)

	def has_relabelled_id(self, label):
		'''
		Returns True if the label exists in the attribute lineage_relabelled (equivalently in labels), False otherwise.
		'''
		return self.labels.has_key(label)

	def label_to_relabelled_id(self, label, time):
		"""
		Returns the relabelled id corresponding to label 'label' at time 'time' (ie returns self.luts[time][label]) if existing, None otherwise.
		"""
		if self.has_time(time) and self.lineage[time].has_key(label):
			return self.luts[time][label] # must exist by construction
		return None

	def relabelled_to_labels(self, label):
		"""
		Returns a dictionnary whose keys are times for which the relabelled_id 'label' is defined and whose associated values are the labels before relabelling.

		"""
		if self.has_relabelled_id(label):
			labels_at_time={}
			for time in range(self.relabelled_birth(label), self.relabelled_death(label)+1):
				labels_at_time[time]=self.reverse_luts[time][label]
				#for key, value in self.reverse_luts[time].iteritems():
				#	labels_at_time[time][value]=key
			return labels_at_time
		return None

	def relabelled_to_label_at_time(self, label, time):
		"""
		Returns a dictionnary whose keys are times for which the relabelled_id 'label' is defined and whose associated values are the labels before relabelling.

		"""
		if self.has_relabelled_id(label) and self.relabelled_birth(label)<=time and self.relabelled_death(label)>=time:
			return self.reverse_luts[time][label]
		return None

	def birth(self, label, time=None):
		'''
		Returns the birth time of the cell of the given initial label at the given time if time is given.
		Returns the birth time of the relabelled cell if time==None.
		'''
		if time==None:
			return self.relabelled_birth(label)
		assert self.exists(label, time), "Attempted to access to label %d at time %d but these do not exist in this lineage."%(label, time)
		return self.relabelled_birth(self.label_to_relabelled_id(label, time))
		#return self.labels[self.label_to_relabelled_id(label, time)]['birth']

	def death(self, label, time=None):
		'''
		Returns the death time of the cell of the given initial label at the given time if time is given.
		Returns the death time of the relabelled cell if time==None.
		'''
		if time==None:
			return self.relabelled_death(label)
		assert self.exists(label, time), "Attempted to access to label %d at time %d but these do not exist in this lineage."%(label, time)
		return self.labels[self.label_to_relabelled_id(label, time)]['death']


	def age(self,label,time):
		''' Returns the age of the cell corresponding to the relabelled key <<label>> at the specified time (may be negative value if time < birth(label)) '''
		return time - self.birth(label, time)

	def normalized_age(self,label,time):
		''' Returns the age of the label (before relabelling) at the specified time, normalized by the life duration of the label.
		Definition:
			if birth(label) == death(label), then 
				normalized_age(label,birth(label)) =  0.5
				normalized_age(label,time) =  (time - birth(time)) / 0.5
			else:
				normalized_age(label,birth(label)) = 0
				normalized_age(label,death(label)) = 1
				normalized_age(label,time) =  (time - birth(time)) / (death(time)-birth(time))
		Properties:
			normalized_age(label,time)<0 if time < birth(label)
			normalized_age(label,time)>1 if time > death(label)
		 '''
		return self.relabelled_normalized_age(self.label_to_relabelled_id(label, time) ,time)


	def relabelled_lineage_lifespan(self, label):
		"""
		Returns the lifespan (tuple of two elements) of the complete lineage the given relabelled cell belongs to.
		"""
		def recursive_lineage_birth(self, l):
			#Birth:
			m=self.relabelled_mother(l)
			if not m:
				return self.relabelled_birth(l)
			return recursive_lineage_birth(self,m)
		def recursive_lineage_death(self,l):
			#Death:
			D=self.relabelled_daughters(l)
			if not D:
				return self.relabelled_death(l)
			v=[]
			for d in D:
				v.append(recursive_lineage_death(self,d))
			return max(v)
		return (recursive_lineage_birth(self,label), recursive_lineage_death(self,label))



	def relabelled_cells_at_time(self,time):
		"""
		Returns the sorted list of relabelled cells that are defined at given time 
		"""
		L=[self.label_to_relabelled_id(l,time) for l in self.labels_at_time(time)]
		L.sort()
		return L

	def relabelled_daughters(self,label):
		"""
		Returns the list of relabelled daughters of given label ([] if no daughters)
		"""
		return self.labels[label]['daughters']

	def relabelled_mother(self,label):
		"""
		Returns the relabelled mother of given label (None if no mother)		
		"""
		return self.labels[label]['mother']

	def relabelled_sisters(self, label):
		'''
		Returns the list of relabelled sisters of given label (including itself)
		'''
		mother=self.relabelled_mother(label)
		if not mother:
			return [label]
		return self.relabelled_daughters(mother)

	def relabelled_sister(self, label):
		"""
		Returns the relabelled sister of given label (None if no sister)
		"""
		from copy import copy
		sisters=copy(self.relabelled_sisters(label))
		n_sisters = len(sisters) 
		if n_sisters < 2:
			return None
		assert n_sisters == 2
		l=sisters.pop()
		if l == label :
			l=sisters.pop()
		return l

	def relabelled_ancestors(self, label):
		"""
		Returns the list of relabelled ancestors of given label (empty list if no ancestor)
		Index i corresponds to the (i+1)th ancestor where 1st ancestor is the mother.
		"""
		mother=self.relabelled_mother(label)
		if not mother:
			return []
		return [mother]+self.relabelled_ancestors(mother)

	def relabelled_descendents(self, labels):
		"""
		Returns the list of relabelled ancestors of given label (empty list if no ancestor)
		Index i corresponds to the (i+1)th ancestor where 1st ancestor is the mother.
		"""

		if not labels:
			return []
		
		descendents=[]

		if type(labels)!=list:
			labels=[labels]

		for label in labels:
			daughters=self.relabelled_daughters(label)
			descendents += daughters + self.relabelled_descendents(daughters)

		return descendents

	def relabelled_birth(self,label):
		'''
		Returns the birth time of the cell corresponding to the relabelled key <<label>>
		'''
		assert self.exists(label)
		return self.labels[label]['birth']

	def relabelled_death(self,label):
		'''
		Returns the death time of the cell corresponding to the relabelled key <<label>>
		'''
		assert self.exists(label)
		return self.labels[label]['death']

	def relabelled_age(self,label,time):
		''' Returns the age of the cell corresponding to the relabelled key <<label>> at the specified time (may be negative value if time < birth(label)) '''
		return time - self.relabelled_birth(label)

		
		
	def relabelled_time_points(self, label):
		'''
		Returns a list of times at which given re-labelled label exists.
			self.relabelled_to_labels(label).keys()
		'''
		if self.exists(label):
			return self.relabelled_to_labels(label).keys()
		return []

	def relabelled_lineage_at_time(self, label, time):
		'''
		Returns a list of relabellings that belong to given label lineage and that exist at specified time. 
		If the lineage does not exist at specified time, returns [].
		'''
		assert(self.exists(label))
		assert(time in self.timepoints())
		if self.relabelled_birth(label) <= time :
			if self.relabelled_death(label) >= time:
				return [label]
			daughters = self.relabelled_daughters(label)
			out=[]
			for l in daughters:
				out += self.relabelled_lineage_at_time(l, time)
			return out
		mother = self.relabelled_mother(label)
		if mother == None:
			return []
		return self.relabelled_lineage_at_time(mother, time)


	def relabelled_lineage(self, label):
		'''
		Returns the (sorted) list of re-labellings that belong to specified label's lineage.
		'''
		assert(self.exists(label))
		return sorted(self.relabelled_ancestors(label)+[label]+self.relabelled_descendents(label))

	def relabelled_time_points(self, label):
		'''
		Returns a list of times at which given re-labelled label exists.
			self.relabelled_to_labels(label).keys()
		'''
		if self.exists(label):

			return self.relabelled_to_labels(label).keys()
		return []

	def relabelled_volume_at_time(self, label, at_time, field_volume='real_volume'):
		'''
		Returns the volume of the cell defined by its (relabel) information at the at_time specified. 
		'''
		#To do :
		#	If at that at_time, the cell has divided or was not born yet, then it returns the volume of its descendants/ancestor. 
		#	If at that at_time, there is no ancestor neither descendant for this cell, then it returns an None value.
		
		assert self.exists(label)
		assert at_time>=self.relabelled_birth(label) and  at_time<=self.relabelled_death(label) 
		if not self.info.has_key(field_volume):
			print "Warning: the asked volume field '"+field_volume+"' was not found. Trying 'volume_information' instead, but processus may fail..."
			field_volume='volume_information'
		assert(self.info[field_volume].has_key(at_time))
		assert(self.reverse_luts.has_key(at_time))
		assert(self.reverse_luts[at_time].has_key(label))
		return self.extract_cell_volume(self.reverse_luts[at_time][label], at_time, field_volume=field_volume)

	def relabelled_barycenter_at_time(self, label, at_time, field_barycenter='real_barycenter', field_volume='real_volume'):
		'''
		Returns a tuple of the barycenter (3 coordinates) and volume of the cell defined by its (relabel) information at the at_time specified. 
		If at that at_time, the cell has divided or was not born yet, then it returns the barycenter of its descendants/ancestor. 
		If at that at_time, there is no ancestor neither descendant for this cell, then it returns an None value.
		'''
		assert(at_time in self.timepoints())
		if not self.info.has_key(field_barycenter):
			print "Warning: the asked barycenter field '"+field_barycenter+"' was not found. Trying 'barycenter' instead, but processus may fail..."
			field_barycenter='barycenter'
		assert(self.info.has_key(field_barycenter))
		assert(self.info[field_barycenter].has_key(at_time))
		time=self.relabelled_birth(label)
		assert(self.reverse_luts.has_key(time))
		assert(self.reverse_luts[time].has_key(label))
		return self.extract_cell_barycenter_at_time(self.reverse_luts[time][label], time, at_time, field_barycenter=field_barycenter, field_volume=field_volume)

	def relabelled_barycenters_at_time(self, at_time, field_barycenter='real_barycenter', field_volume='real_volume'):
		"""
		Returns a dictionnary of labels existing at specified time as keys and (x, y, z) corresponding barycenter coordinates as values
		"""

		assert(at_time in self.timepoints())
		if not self.info.has_key(field_barycenter):
			print "Warning: the asked barycenter field '"+field_barycenter+"' was not found. Trying 'barycenter' instead, but processus may fail..."
			field_barycenter='barycenter'
		assert(self.info.has_key(field_barycenter))
		assert(self.info[field_barycenter].has_key(at_time))
		dico={}
		list_of_relabels=self.reverse_luts[at_time].keys()
		for relabel in list_of_relabels:
			value=self.relabelled_barycenter_at_time(relabel, at_time, field_barycenter=field_barycenter, field_volume=field_volume)
			if value:
				dico[relabel]=value
		return dico

	def relabelled_normalized_age(self,label,time):
		''' Returns the age of the label (after relabelling) at the specified time, normalized by the life duration of the label.
		Definition:
			if birth(label) == death(label), then 
				relabelled_normalized_age(label,birth(label)) =  0.5
				relabelled_normalized_age(label,time) =  (time - birth(time)) / 0.5
			else:
				relabelled_normalized_age(label,birth(label)) = 0
				relabelled_normalized_age(label,death(label)) = 1
				relabelled_normalized_age(label,time) =  (time - birth(time)) / (death(time)-birth(time))
		Properties:
			relabelled_normalized_age(label,time)<0 if time < birth(label)
			relabelled_normalized_age(label,time)>1 if time > death(label)
		'''
		if self.relabelled_birth(label) == self.relabelled_death(label):
			if time == self.relabelled_birth(label):
				return 0.5
			return (time - self.relabelled_birth(label)) * 2 # on divise par 0.5
		else:
			return float(time - self.relabelled_birth(label)) / (self.relabelled_death( label) - self.relabelled_birth(label))

	def exists(self, label, time=None):
		'''
		Returns True if the (label,time) or (relabel) cell exists in the lineage or lineage_relabelling, False otherwise. 
		'''
		if time==None:
			return self.labels.has_key(label)
		return (self.has_time(time) and self.lineage[time].has_key(label))


	def extract_cell_segment(self, label, time):
		'''
		Returns a dict of time keys at which given label at given time exists, and of values corresponding to labels at times.
			self.relabelled_to_labels(self.label_to_relabelled_id(label, time))
		'''
		return self.relabelled_to_labels(self.label_to_relabelled_id(label, time))

	def extract_cell_mother(self, label, time):
		'''
		Returns a dict of time keys at which the mother of given label at given time exists, and of values corresponding to labels at times.
		'''
		return self.relabelled_to_labels(self.relabelled_mother(self.label_to_relabelled_id(label, time)))

	def extract_cell_daughters(self, label, time):
		'''
		Returns a list of dicts of time keys at which the daughters of given label at given time exist, and of values corresponding to labels at times.
		'''
		return [self.relabelled_to_labels(i) for i in self.relabelled_daughters(self.label_to_relabelled_id(label, time))]

	def extract_cell_sister(self, label, time):
		'''
		Returns the label (or list of labels) of the sister (or sister's descendents) corresponding to the given cell at given time.
		If not existing, returns None.
		'''
		'''
		relabelled_cell = self.label_to_relabelled_id(label, time)
		relabelled_mother = self.relabelled_mother(relabelled_cell)
		if not relabelled_mother:
			return None
		'''
		relabelled_sister = self.relabelled_sister(self.label_to_relabelled_id(label, time))
		if relabelled_sister:
			if self.relabelled_death(relabelled_sister)>=time:
				return self.relabelled_to_label_at_time(relabelled_sister, time)
			time_sister, label_sister = self.relabelled_to_labels(relabelled_sister).items()[0]
			descendents=self.extract_cell_descendents(label_sister, time_sister)
			if descendents.has_key(time):
				return descendents[time]
		return None

	def extract_cell_ancestors(self, label, time):
		'''
		Returns dictionnary of ancestors of given label at given time
		'''

		relabelled_ancestors=self.relabelled_ancestors(self.label_to_relabelled_id(label, time))
		cell_ancestors={}
		for i in relabelled_ancestors:
			cell_ancestors=merge_two_dicts(cell_ancestors, self.relabelled_to_labels(i))
		return cell_ancestors
	
	def extract_cell_descendents(self, label, time):
		"""

		"""
		relabelled_descendents=self.relabelled_descendents(self.label_to_relabelled_id(label, time))
		cell_descendents={}
		for i in relabelled_descendents:
			cell_descendents=cumul_two_dicts(cell_descendents, self.relabelled_to_labels(i))

		return cell_descendents

	#def extract_cell_lineage():
	#	'''
	#	
	#	'''
	#
	#	return None


	def extract_cell_info(self, label, time):
		'''
		Returns all the available information for given cell (taken from self.info)
		'''
		info={}
		for key in self.info.keys():
			if self.info[key].has_key(time) and self.info[key][time].has_key(label):
				info[key] = self.info[key][time][label]

		return info


	def extract_cell_lineage_at_time(self, label, time, at_time):
		"""
		
		"""

		if at_time >= self.birth(label, time) :
			if at_time <= self.death(label, time):
				return self.extract_cell_segment(label, time)[at_time]

			cell_descendents = self.extract_cell_descendents(label, time)
			if cell_descendents.has_key(at_time):
				return cell_descendents[at_time]
			return None

		cell_ancestors=self.extract_cell_ancestors(label, time)
		if cell_ancestors.has_key(at_time):
			return cell_ancestors[at_time]

		return None

	def extract_cell_barycenter_at_time(self, label, time, at_time, field_barycenter='real_barycenter', field_volume='real_volume'):
		"""
		Returns a tuple of the barycenter (3 coordinates) and volume of the cell defined by its (label,time) information at the at_time specified. 
		If at that at_time, the cell has divided or was not born yet, then it returns the barycenter of its descendants/ancestor. 
		If at that at_time, there is no ancestor neither descendant for this cell, then it returns an None value.
		"""
		assert(time in self.timepoints())
		assert(at_time in self.timepoints())
		if not self.info.has_key(field_barycenter):
			print "Warning: the asked barycenter field '"+field_barycenter+"' was not found. Trying 'barycenter' instead, but processus may fail..."
			field_barycenter='barycenter'
		assert(self.info.has_key(field_barycenter))
		if not self.info.has_key(field_volume):
			print "Warning: the asked volume field '"+field_volume+"' was not found. Trying 'volumes' instead, but processus may fail..."
			field_volume='volumes_information'
		assert(self.info.has_key(field_volume))
		labels_at_time=self.extract_cell_lineage_at_time(label, time, at_time)
		if not labels_at_time:
			return None
		x=0.0
		y=0.0
		z=0.0
		v=0.0
		label_barycenters_at_time = self.extract_barycenters_at_time(at_time, field_barycenter=field_barycenter)
		if type(labels_at_time)==int:
			labels_at_time=[labels_at_time]
		for l in labels_at_time:
			if not label_barycenters_at_time.has_key(l):
				print "Warning: unexpected empty field '"+field_barycenter+"' at time "+str(at_time)+" for label "+str(l)+"..."
			else:
				(l_x,l_y,l_z)=label_barycenters_at_time[l]
				l_v=self.extract_cell_volume(l, at_time, field_volume=field_volume)
				if not l_v:
					l_v=1.0
					print "Warning: volume not found for label "+str(l)+" at time "+str(at_time)+"..."
				x = x + l_x*l_v 
				y = y + l_y*l_v 
				z = z + l_z*l_v 
				v = v + l_v
		x = x/v
		y = y/v
		z = z/v
		return x,y,z,v

	def next_death(self, time):
		"""
		Returns (labels_at_time, label_death_time).
		labels_at_time can be a list with one to N elements. 
		"""
		assert(time in self.timepoints())
		labels_at_time=[]
		next_death_time=None
		for label in self.labels_at_time(time):
			if next_death_time == None :
				labels_at_time = [label]
				next_death_time = self.death(label, time)
			else:
				if next_death_time >= self.death(label, time):
					if next_death_time > self.death(label, time):
						labels_at_time = [label]
						next_death_time = self.death(label, time)
					else:
						labels_at_time.append(label)
		return labels_at_time, next_death_time

	def relabelled_next_death(self, time):
		"""
		Returns (relabels_at_time, relabel_death_time).
		relabels_at_time can be a list with one to N elements. 
		"""
		assert(time in self.timepoints())
		labels_at_time=[]
		next_death_time=None
		for label in self.reverse_luts[time].keys():
			if next_death_time == None :
				labels_at_time = [label]
				next_death_time = self.relabelled_death(label)
			else:
				if next_death_time >= self.relabelled_death(label):
					if next_death_time > self.relabelled_death(label):
						labels_at_time = [label]
						next_death_time = self.relabelled_death(label)
					else:
						labels_at_time.append(label)
		return labels_at_time, next_death_time

	def extract_barycenters_at_time(self, time, field_barycenter='real_barycenter'):
		"""
		Returns a dictionnary of labels existing at specified time as keys and (x, y, z) corresponding barycenter coordinates as values
		"""
		assert(time in self.timepoints())
		if not self.info.has_key(field_barycenter):
			print "Warning: the asked barycenter field '"+field_barycenter+"' was not found. Trying 'barycenter' instead, but processus may fail..."
			field_barycenter='barycenter'
		assert(self.info.has_key(field_barycenter))
		assert(self.info[field_barycenter].has_key(time))
		return self.info[field_barycenter][time]

	def extract_field_at_time(self, field, time):
		"""
		Returns a dictionnary of labels existing at specified time as keys and corresponding information contained in the asked field
		"""
		assert(time in self.timepoints())
		assert(self.info.has_key(field))
		assert(self.info[field].has_key(time))
		return self.info[field][time]

	def extract_cell_volume(self, label, time, field_volume='real_volume'):
		"""
		Returns the volume of given cell label at given time with respect to the specified volume field.
		"""
		assert(time in self.timepoints())
		if not self.info.has_key(field_volume):
			print "Warning: the asked volume field '"+field_volume+"' was not found. Trying 'volumes_information' instead, but processus may fail..."
			field_volumes='volumes_information'
		#assert(self.info.has_key(field_volume))
		if self.info.has_key(field_volume):
			if self.info[field_volume].has_key(time):
				if self.info[field_volume][time].has_key(label):
					return self.info[field_volume][time][label]
		return None

	def relabelled_volumes(self,label, t1=None, t2=None, field_volume="real_volume"):
		'''
		Returns the volume variation of the specified relabelled cell from its birth to its death.
		If specified, from t1 to t2 instead.
		Parameter field_volume corresponds to the method extract_cell_volume volume field (default="real_volume").
		'''
		if t1==None:
			t1=self.relabelled_birth(label)
		if t2==None:
			t2=self.relabelled_death(label)
		D={}
		for t in [tmp for tmp in self.timepoints() if (tmp >=t1 and tmp<=t2)]:
			D[t]=self.relabelled_volume_at_time(label,t, field_volume=field_volume)
		return D

	def relabelled_growth_ratio(self, label, t1=None, t2=None, field_volume='real_volume'):
		"""
		Returns the ratio of volumes volume_death/volume_birth of the specified relabelled cell between its birth and its death.
		If specified, returns the ratio of volumes volume_t2/volume_t1 of the specified relabelled cell between t1 and t2.
		Parameter field_volume corresponds to the method extract_cell_volume volume field (default="real_volume").
		"""
		if t1==None:
			t1=self.relabelled_birth(label)
		if t2==None:
			t2=self.relabelled_death(label)
		volumebirth=self.relabelled_volume_at_time(label,t1, field_volume=field_volume)
		volumedeath=self.relabelled_volume_at_time(label,t2, field_volume=field_volume)
		if not volumebirth or volumedeath == None:
			return None
		return volumedeath/volumebirth

	def relabelled_growth_ratios(self, label, t1=None, t2=None, field_volume='real_volume'):
		"""
		Returns the dictionnary of ratios of volumes volume_t/volume_ of the specified relabelled cell between its birth and its death.
		Keys are the time-points and values are the corresponding volume ratios.
		If specified, returns the ratio of volumes volume_t2/volume_t1 of the specified relabelled cell between t1 and t2.
		Parameter field_volume corresponds to the method extract_cell_volume volume field (default="real_volume").
		"""
		if t1==None:
			t1=self.relabelled_birth(label)
		if t2==None:
			t2=self.relabelled_death(label)
		volumet1=self.relabelled_volume_at_time(label,t1, field_volume=field_volume)
		D={}
		if volumet1 :
			D[t1]=1.0
		else:
			return {tmp:None for tmp in self.timepoints() if (tmp >=t1 and tmp<=t2)}
		for t in [tmp for tmp in self.timepoints() if (tmp >t1 and tmp<=t2)]:
			volumet=self.relabelled_volume_at_time(label,t, field_volume=field_volume)
			if volumet == None:
				D[t]=None
			else:
				D[t]=volumet/volumet1
		return D

	def get_fields(self):
		"""
		Returns the list of keys contained in the info dictionnary.
		(those which can be used in the method "Morpheme_lineage.extract_field_at_time" as fields)
		"""
		return self.info.keys()



def print_label_barycenters_at_time(morpheme_lineage, time, filename):
	"""
	Prints in the specified filename the list of label barycenters at specified time in the given morpheme_lineage
	"""
	assert(not os.path.dirname(filename) or os.path.exists(os.path.dirname(filename)))
	dico=morpheme_lineage.extract_label_barycenters_at_time(time)
	f = open(filename, 'w')
	for k,v in dico.iteritems():
		print >>f, k, v[0], v[1], v[2]
	f.close()		

def control_persistency(lineageDict,verbose=False):
	'''
	Controls that every key of the lineageDict is a dict with the same keys (time-points) and sub-keys (image labels).
	If true, returns True. Otherwise, returns False. Verbose mode prints detailed information.
	'''

	key_list=lineageDict.keys()
	if len(key_list)<=1:
		return True


	return True


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def cumul_two_dicts(x, y):
    """Given two dicts, cumul them into a new dict as a shallow copy."""
    z = {}

    for k,v in x.iteritems():
    	z[k]=v

    for k,v in y.iteritems():
	    if z.has_key(k):
	    	z[k].append(v)
	    else:
	    	z[k]=[v]

    return z

def reverse_luts(luts):
	'''
	Returns a "reversed" luts dictionnary :
	if luts[time][id] == value   then   rev_luts[time][value] == id
	'''
	rev_luts={}
	for time in luts.keys():
		rev_luts[time]={}

		for label, relabel in luts[time].iteritems():
			rev_luts[time][relabel]=label

	return rev_luts

def to_label_classification(lineage):
	'''
	Returns a dictionnary whose keys are the labels contained in "lineage_relabelling" and whose elements are dictionnaries relatively to the label.
	/!\ CAUTION /!\ "lineage" must be a lineage so that one cell has the same id during all its life and this id must be unique.
	If l is a label in "lineage_relabelling" that lives from time t1 to t2, has the mother-label ml and the daughter-labels [dl1, dl2], then :
	out[l] exists and out[l] = ( 'birth':t1, 'death':t2, 'mother':[ml], 'daughters':[dl1, dl2] )
	'''
	labels={}
	if type(lineage).__name__ == 'instance' and lineage.__class__.__name__ == 'Morpheme_lineage'  : # cas ou l'argument est un objet de la classe Morpheme_lineage_relabelled
		lineage=lineage.lineage_relabelling
	times=lineage.keys()
	

	for t in times:
		t_labels=lineage[t].keys()
		for l in t_labels:
			if not labels.has_key(l):
				#labels[l]=(t, t, [], [])
				labels[l]={'birth': t, 'death':t, 'mother': None, 'daughters': []}
			labels[l]['death']=t
			daughters=lineage[t][l]
			for time_daughter, label_daughter in daughters:
				if label_daughter==l:
					labels[l]['death']=time_daughter
				else:
					if not labels.has_key(label_daughter):
						labels[label_daughter]={'birth':time_daughter, 'death':time_daughter, 'mother':l, 'daughters': []}
					else:
						labels[label_daughter]['mother']=l
					labels[l]['daughters'].append(label_daughter)

	return labels


def cmp_cell_tuplet(a,b):
    '''
    if a[0] > b[0] :
        return 1
    elif a[0] < b[0]:
        return -1
    else:
    	if a[1] > b[1]:
    		return 1
    	elif a[1] < b[1]:
    		return -1
    return 0
    '''
    if a[0] > b[0] :
        return 1
    elif a[0] < b[0]:
        return -1
    else:
    	if a[1] > b[1]:
    		return 1
    	elif a[1] < b[1]:
    		return -1
    return 0



def getTimeKey(time):
	"""
	Usage : 
		getTimeKey(156)	returns 'T 156'
	"""
	return key_time + key_space + str(time)

def getIDKey(ID):
	"""
	Usage : 
		getIDKey(209)	returns 'L 209'
	"""
	return key_ID + key_space + str(ID)


def getEmbryoKey(embryoName=''):
	"""
	Usage : 
		getEmbryoKey('161213_lucie')	returns 'E 161213_lucie'
		getEmbryoKey('')	returns ''
	"""
	if embryoName:
		return key_embryoName + key_space + embryoName
	else:
		return ''



def getKey(time='', ID='', embryoName=''):
	'''
	Usage : 
		getKey(156,209)	returns '[T 156, L 209]'
		getKey(time=156,ID=209,'161213_lucie')	returns '[E 161213_lucie, T 156, L 209]'
	'''
	key=""
	sub_key=getEmbryoKey(embryoName)
	if sub_key:
		if key:
			key += key_separator+ ' ' +sub_key
		else:
			key = sub_key

	sub_key=str(time)
	if sub_key:
		sub_key=key_time+ key_space +sub_key
		if key:
			key += key_separator+ ' ' +sub_key
		else:
			key = sub_key

	sub_key=str(ID)
	if sub_key:
		sub_key=key_ID+ key_space +sub_key
		if key:
			key += key_separator+ ' ' +sub_key
		else:
			key = sub_key

	return key_opener+key+key_closer

def readKeyTime(key):
	k_list=key.lstrip(key_opener).rstrip(key_closer).split(key_separator)
	for k in k_list:
		k_split=k.strip().split(key_space,1)
		if k_split[0] == key_time:
			if len(k_split)==2:
				return int(k_split[1])
			else :
				return None
	return None
	 
def readKeyID(key):
	k_list=key.lstrip(key_opener).rstrip(key_closer).split(key_separator)
	for k in k_list:
		k_split=k.strip().split(key_space,1)
		if k_split[0] == key_ID:
			if len(k_split)==2:
				return int(k_split[1])
			else :
				return None
	return None

def readKeyEmbryoName(key):
	k_list=key.lstrip(key_opener).rstrip(key_closer).split(key_separator)
	for k in k_list:
		k_split=k.strip().split(key_space,1)
		if k_split[0] == key_embryoName:
			if len(k_split)==2:
				return k_split[1]
			else :
				return None
	return None
	 
def readKey(key):
	'''
	Retourne un dictionnaire D avec les cles suivantes : 
	'time', 'ID', 'embryoName'
	Exemple : 
	D = readKey('[T 156,L 209]') 
	print D['time'] -> 156
	print D['ID'] -> 209
	print D['embryoName'] -> None
	'''
	dic={'time':None, 'ID':None, 'embryoName':None}
	if key.startswith(key_opener) & key.endswith(key_closer) & (key.count(key_opener) + key.count(key_closer) == 2):
		dic['time']=readKeyTime(key)
		dic['ID']=readKeyID(key)
		dic['embryoName']=readKeyEmbryoName(key)
	return dic


def timeToString(time, timeFormat='%03d'):
	'''
	Retourne la chaine time en string selon le format specifie (return timeFormat % time).
	'''
	if not timeFormat:
		return str(time)
	return timeFormat % time


def lineageRelabelling(orderedLineage):
	'''
	Re-ecriture d'un dictionnaire de lignees cellulaires temporellement ordonnees en un nouveau dictionnaire verifiant qu'un objet conserve le meme label tout au long de sa "vie" 
	et construction d'un dictionnaire de correspondances entre les anciens et les nouveaux labels de chaque time-point.
	Autrement dit, si orderedLineage[time][label]=[(other_time, other_label)] (ie ce label a l'instant time a une unique correspondance au temps other_time)
	alors renumberedLineage[time][new_label]=[(other_time, new_label)] et lut[time][label]=lut[other_time][other_label]=new_label.
	'''

	timepoints=orderedLineage.keys()

	last_new_label=orderedLineage[timepoints[0]].keys()[0]

	renumberedLineage={}
	lut={}

	for time in timepoints[:]:
		if not lut.has_key(time):
			lut[time] = {}
		if not renumberedLineage.has_key(time):
			renumberedLineage[time] = {}
		for label_at_time in orderedLineage[time].keys()[:]:
			if not lut[time].has_key(label_at_time):
				lut[time][label_at_time] = last_new_label
				last_new_label += 1

		for label_at_time,daughters_at_time in orderedLineage[time].iteritems():
			if not lut[time].has_key(label_at_time):
				lut[time][label_at_time] = last_new_label
				last_new_label += 1

			new_label_at_time = lut[time][label_at_time]
			# On considere le label (time, label_at_time) du lineage d'origine qui a ete renumorote en new_label_at_time
			if not renumberedLineage[time].has_key(new_label_at_time):
				renumberedLineage[time][new_label_at_time]=[]

			# Les enfants de label_at_time sont dans la liste daughters_at_time provenant du lineage d'origine
			# Si 0 enfant: on ne fait rien	
			# Si 1 enfant: on "propage" le nouveau label au time-point suivant (dans lut) et on y ajoute au nouveau lineage
			# Si 2 enfants ou plus: on attribue des nouveaux labels aux enfants (dans lut) et on y ajoute au nouveau lineage
			if len(daughters_at_time) == 1:
				daughter_at_time=daughters_at_time[0]
				time_daughter=daughter_at_time[0]
				label_daughter=daughter_at_time[1]
				if not lut.has_key(time_daughter):
					lut[time_daughter]={}
				lut[time_daughter][label_daughter]=new_label_at_time
				renumberedLineage[time][new_label_at_time].append((time_daughter,lut[time_daughter][label_daughter]))
			elif len(daughters_at_time) > 1:
				for daughter_at_time in daughters_at_time[:]:
					time_daughter=daughter_at_time[0]
					label_daughter=daughter_at_time[1]
					if not lut.has_key(time_daughter):
						lut[time_daughter]={}
					lut[time_daughter][label_daughter] = last_new_label
					last_new_label += 1
					renumberedLineage[time][new_label_at_time].append((time_daughter,lut[time_daughter][label_daughter]))
	return renumberedLineage, lut



def writeTextLineage(lin_tree, textFileNameOut, embryoName='', version=versionWriter):
	''' 
	Ecriture de fichiers de lineages cellulaires au format texte respectant les conventions morpheme 
	'''
	lineageText="# Version " + version + "\n"
	
	if type(lin_tree).__name__ == 'instance' and lin_tree.__class__.__name__ == 'Morpheme_lineage' : # cas ou l'argument est un objet de la classe Morpheme_lineage
		lin_tree=lin_tree.lineage

	'''if version=='0.0':
		for mother, daughters in lin_tree.iteritems():
			lineageText+= clique_opener + pkl_converter.getKeyFromLongID(mother,embryoName=embryoName) + clique_closer + cliques_separator + clique_opener
			is_first_daughter=True
			for daughter in daughters:
				if is_first_daughter:
					is_first_daughter=False
				else:
					lineageText+=element_list_separator
				lineageText+=getKeyFromLongID(daughter,embryoName=embryoName)
			lineageText+=clique_closer + end_correspondance
	'''
	if version=='1.0':
		for time,IDs in lin_tree.iteritems():
			for ID,daughters in IDs.iteritems():
				lineageText+= clique_opener + getKey(time=time, ID=ID, embryoName=embryoName) + clique_closer + cliques_separator + clique_opener
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



def writeTextLineageRelabelled(lin_tree, textFileNameOut, embryoName='', version=versionWriter):
	''' 
	Ecriture de fichiers de lineages cellulaires relabellises au format texte respectant les conventions morpheme 
	'''
	
	if type(lin_tree).__name__ == 'instance' and lin_tree.__class__.__name__ == 'Morpheme_lineage' : # cas ou l'argument est un objet de la classe Morpheme_lineage
		writeTextLineage(lin_tree.lineage_relabelling, textFileNameOut, embryoName, version)



def writeTextLUT(lut, fileNamePrefix='lut_t', fileNameSuffix='', fileNameExtension='.txt', timeFormat='%03d'):
	'''
	Ecriture de fichiers de look up table pour les renumotations de labels de sequences d'images a partir de lineages
	Format : une ligne par association de label ancien / nouveau au format : str(old_label) + ' ' + str(new_label) + '\n'
	Fichiers generes de la forme fileNameOut=fileNamePrefix + timeToString(time, timeFormat='%03d') + fileNameSuffix + fileNameExtension
	'''
	if type(lut).__name__ == 'instance' and lut.__class__.__name__ == 'Morpheme_lineage'  : # cas ou l'argument est un objet de la classe Morpheme_lineage_relabelled
		lut=lut.luts

	for time, lut_at_time in lut.iteritems():
		lutText=''
		for old_label, new_label in lut_at_time.iteritems():
			lutText+=str(old_label) + ' ' + str(new_label) +'\n'
		fileNameOut=fileNamePrefix + timeToString(time, timeFormat) + fileNameSuffix + fileNameExtension
		file=open(fileNameOut, "w")
		file.write(lutText)
		file.close()



def readTextLineage(filename,t=None):
	'''
	Return a lineage tree dictionnary (morpheme lineage) from a text file if exist until time t other return an empty lineage tree.
	'''
	lin_tree={}
	if os.path.exists(filename):
		f=open(filename)
		for line in f:
			li=line.strip()
			if not li.startswith(comment):
				info=li.split(cliques_separator.strip())
				if len(info) == 2: # info[0] : mother ; info[1] : daughters 
					mothers=info[0].strip()
					daughters=info[1].strip()
					if (mothers.startswith(clique_opener) & mothers.endswith(clique_closer) & daughters.startswith(clique_opener) & daughters.endswith(clique_closer)): # respecte la syntaxe attendue des ouvertures - fermetures de cliques
						mos=mothers.lstrip(clique_opener).rstrip(clique_closer) # String des mothers sans les opener-closer de cliques
						das=daughters.lstrip(clique_opener).rstrip(clique_closer) # String des daughters sans les opener-closer de cliques
						# Parcours de mothers 
						parsed_mothers_list_of_keys=[]
						mos_split = mos.split(element_list_separator)
						for i in range(len(mos_split)):
							m_split=mos_split[i]
							if m_split.startswith(key_opener):
								begin=i
							if m_split.endswith(key_closer):
								mother_key=key_separator.join(mos_split[begin:i+1]).strip()
								#mother_dict=readKey(mother_key)
								parsed_mothers_list_of_keys.append(mother_key)
								
						# Parcours de daughters
						parsed_daughters_list_of_keys=[]
						das_split = das.split(element_list_separator)
						for i in range(len(das_split)):
							d_split=das_split[i]
							if d_split.startswith(key_opener):
								begin=i
							if d_split.endswith(key_closer):
								parsed_daughters_list_of_keys.append(key_separator.join(das_split[begin:i+1]).strip())

						# Remplissage du dictionnaire
						#for mother,daughters in lineageDict.iteritems():
						for mother_key in parsed_mothers_list_of_keys :
							mother_dict=readKey(mother_key)
							time=mother_dict['time']
							if not lin_tree.has_key(time):
								lin_tree[time]={}
							ID=mother_dict['ID']
							if not lin_tree[time].has_key(ID):
								lin_tree[time][ID]=[]
							for daughter_key in parsed_daughters_list_of_keys :
								daughter_dict=readKey(daughter_key)
								lin_tree[time][ID].append((daughter_dict['time'],daughter_dict['ID']))
							lin_tree[time][ID].sort(cmp_cell_tuplet)
		f.close()
	return lin_tree;


def readPointCloud(file):
	'''
	Return a dictionnary of point-cloud of format {label:(x, y, z), ...} or {label:(x, y, z, volume), ...}.
	'''
	cloud={}
	
	coef_x=1.0
	coef_y=1.0
	coef_z=1.0

	if os.path.exists(file):
		f=open(file)
		for line in f:
			li=line.strip()
			if not li.startswith('#'):
				info=li.split()
				if info[0]=="Voxelsize":
					assert(len(info)==4)
					coef_x=float(info[1])
					coef_y=float(info[2])
					coef_z=float(info[3])
				else:
					if len(info) == 4:
						if not cloud.has_key(int(info[0])):
							cloud[int(info[0])]=None
						cloud[int(info[0])]=(float(info[1]), float(info[2]), float(info[3]))
					if len(info) == 5:
						if not cloud.has_key(int(info[0])):
							cloud[int(info[0])]=None
						cloud[int(info[0])]=(float(info[1]), float(info[2]), float(info[3]), float(info[4]))
	
	if coef_x != 1.0 or coef_y != 1.0 or coef_z != 1.0:
		for k in cloud.keys():
			if len(cloud[k])==3:
				cloud[k] = (cloud[k][0]*coef_x, cloud[k][1]*coef_y, cloud[k][2]*coef_z)
			if len(cloud[k])==4:
				cloud[k] = (cloud[k][0]*coef_x, cloud[k][1]*coef_y, cloud[k][2]*coef_z, cloud[k][3]*coef_x*coef_y*coef_z)
	return cloud

def label_couple_cost(lineage_in, label_in, time_in, lineage_ext, label_ext, time_ext, weight = 1.0):
	'''
	Returns the cost associated to a pair of labels (original labels) depending on the absolute difference between their normalized ages in their corresponding lineages
	'''
	from math import fabs
	cost = fabs( lineage_in.normalized_age(label_in,time_in) - lineage_ext.normalized_age(label_ext,time_ext) )
	return weight * cost

def relabelled_couple_cost(lineage_in, label_in, time_in, lineage_ext, label_ext, time_ext, weight = 1.0):
	'''
	Returns the cost associated to a pair of labels (relabelled) depending on the absolute difference between their normalized ages in their corresponding lineages
	'''
	from math import fabs
	cost = fabs( lineage_in.relabelled_normalized_age(label_in,time_in) - lineage_ext.relabelled_normalized_age(label_ext,time_ext) )
	return weight * cost


def write_barycenters(barycenters_dict, filename, ID_key="ID", voxelsize=None):
	'''
	Writes in filename a text file with information contained in barycenters_dict at adapted format for the morpheme vt program pointCloudRegistration.	
	'''
	text="# "+ID_key+" X Y Z <VOLUME>"
	text = text+"\n"
	if voxelsize:
		assert(type(voxelsize)==tuple or type(voxelsize)==list)
		assert(len(voxelsize)==3)
		text = "Voxelsize "+str(voxelsize).lstrip('(').lstrip('[').rstrip(')').rstrip(']').replace(',','')+"\n"+text
	keys=barycenters_dict.keys()
	keys.sort()
	for k in keys:
		line=str(k)+" "+str(barycenters_dict[k]).lstrip('(').rstrip(')').replace(',','')
		text = text + line + "\n"
	text_file = open(filename, "w")
	text_file.write(text)
	text_file.close()

def write_correspondences(correspondences_dict, filename, ID_key="LABEL"):
	'''
	Writes in filename a text file with information contained in correspondences_dict at adapted format for the morpheme vt point cloud registration programs.
	'''
	text="# "+ID_key+"-ref "+ID_key+"-flo\n"
	keys=correspondences_dict.keys()
	keys.sort()
	for k in keys:
		text=text+str(k)+' '+str(correspondences_dict[k])+'\n'
	text_file = open(filename, "w")
	text_file.write(text)
	text_file.close()

def lineage_relabelled_point_cloud_at_time(lineage, time, filename=None, field_barycenter='real_barycenter', field_volume='real_volume', voxelsize=None):
	"""
	Returns at dict format (or write in filename if specified at adapted format for the morpheme vt program pointCloudRegistration) 
	the list of relabelled points with barycenter coordinates and eventually volume information and eventually voxelsize at 
	specified time-point.
	"""
	text="# RELABELLED_ID X Y Z"
	if field_volume:
		text = text+" VOLUME"
	text = text+"\n"
	if voxelsize:
		assert(type(voxelsize)==tuple or type(voxelsize)==list)
		assert(len(voxelsize)==3)
		text = "Voxelsize "+str(voxelsize).lstrip('(').lstrip('[').rstrip(')').rstrip(']').replace(',','')+"\n"+text
	assert(lineage.info.has_key(field_barycenter))
	if field_volume:
		assert(lineage.info.has_key(field_volume))
	assert(lineage.has_time(time))
	assert(lineage.luts.has_key(time))
	dico={}
	for label in lineage.luts[time].keys():
		dico[lineage.luts[time][label]]=lineage.extract_cell_info(label,time)[field_barycenter]
		#lineage.luts[time][label]
		#line=str(lineage.luts[time][label])+" "+str(lineage.extract_cell_info(label,time)[field_barycenter]).lstrip('(').rstrip(')').replace(',','')
		if field_volume:
			dico[lineage.luts[time][label]] += (lineage.extract_cell_info(label,time)[field_volume],)
			#line = line+" "+str(lineage.extract_cell_info(label,time)[field_volume])
		#text = text + line + "\n"
	if not filename:
		return dico
	for k,v in dico.items():
		line=str(k)+" "+str(v).lstrip('(').rstrip(')').replace(',','')
		text = text + line + "\n"
	text_file = open(filename, "w")
	text_file.write(text)
	text_file.close()
    
def lineage_relabelled_point_cloud_at_time_from_list(lineage, time, relabelled_list, filename=None, field_barycenter='real_barycenter', field_volume='real_volume', voxelsize=None):
	"""
	Returns at dict format (or write in filename if specified at adapted format for the morpheme vt program pointCloudRegistration) 
	the list of relabelled points with barycenter coordinates and eventually volume information and eventually voxelsize at 
	specified time-point with respect to the given relabelled_list.
	If needed, a barycenter propagation is processed when specified relabelling does not exist anymore (or yet) at given time.
	"""
	text="# RELABELLED_ID X Y Z"
	if field_volume:
		text = text+" VOLUME"
	text = text+"\n"
	if voxelsize:
		assert(type(voxelsize)==tuple or type(voxelsize)==list)
		assert(len(voxelsize)==3)
		text = "Voxelsize "+str(voxelsize).lstrip('(').lstrip('[').rstrip(')').rstrip(']').replace(',','')+"\n"+text
	assert(lineage.info.has_key(field_barycenter))
	if field_volume:
		assert(lineage.info.has_key(field_volume))
	assert(lineage.has_time(time))
	dico={}
	for label in relabelled_list:
		assert(lineage.exists(label))
		bary=lineage.relabelled_barycenter_at_time(label, time, field_barycenter=field_barycenter, field_volume=field_volume)
		if bary != None:
			dico[label]=bary
	if not filename:
		return dico
	for k,v in dico.items():
		line=str(k)+" "+str(v).lstrip('(').rstrip(')').replace(',','')
		text = text + line + "\n"
	text_file = open(filename, "w")
	text_file.write(text)
	text_file.close()

def get_relabelled_lists_at_time(label_correspondences, lineage_ref, lineage_flo, time_ref, time_flo, keep_divided_cells=True):
	'''
	Returns a pair of lists that correspond to relabel equivalences between lineage_ref and lineage_flo at time_ref and time_flo. 
	Useful function for cell correspondence propagation algorithm.
	with respect to label_correspondences dictionnary of known correspondences between both lineages.
		label_correspondences : dictionnary with lineage_ref relabels as keys and lineage_flo corresponding relabels as values 
		lineage_ref : Morpheme_lineage reference object
		lineage_flo : Morpheme_lineage floating object
		time_ref : lineage_ref time-point of interest
		time_flo : lineage_flo time-point of interest
		keep_divided_cells : True if the user is expecting to keep relabels associations even if one (and only one) of the relabels undergoes some divisions (default) ; False otherwise
	'''
	assert(lineage_ref.has_time(time_ref))
	assert(lineage_flo.has_time(time_flo))

	list_at_time_ref=[]
	list_at_time_flo=[]

	labels_at_time_ref = lineage_ref.reverse_luts[time_ref].keys()
	labels_at_time_flo = lineage_flo.reverse_luts[time_flo].keys()

	label_correspondances_ref = label_correspondences.keys()
	label_correspondances_flo = label_correspondences.values()

	for l_ref in labels_at_time_ref:
		if l_ref in label_correspondances_ref:
			l_flo_associated = label_correspondences[l_ref]
			if l_flo_associated in labels_at_time_flo: # Easy case: both corresponding labels exist in respective lineages at specified times
				list_at_time_ref.append(l_ref)
				list_at_time_flo.append(l_flo_associated)
			else:
				# On ajoute la paire (l_ref, l_flo) ssi le flag keep_divided_cells est a True et la cellule correspondante a l_ref dans le lineage_flo est deja divisee a time_flo
				if keep_divided_cells:
					assert(lineage_flo.exists(l_flo_associated))
					l_flo_associated_descendants = (set(lineage_flo.relabelled_descendents(l_flo_associated)) & set(labels_at_time_flo))
					if l_flo_associated_descendants: 
						# La cellule correspondante a l_ref dans le lineage_flo est deja divisee a time_flo
						list_at_time_ref.append(l_ref)
						list_at_time_flo.append(l_flo_associated)
					#else : nothing.
	if keep_divided_cells: # The easy case is already treaten with the first loop
		inv_label_correspondences = { v: k for k,v in label_correspondences.iteritems() }
		for l_flo in labels_at_time_flo:
			if l_flo in label_correspondances_flo:
				l_ref_associated = inv_label_correspondences[l_flo]
				if not l_ref_associated in labels_at_time_ref:
					assert(lineage_ref.exists(l_ref_associated))
					l_ref_associated_descendants = (set(lineage_ref.relabelled_descendents(l_ref_associated)) & set(labels_at_time_ref))
					if l_ref_associated_descendants:
						# La cellule correspondante a l_flo dans le lineage_ref est deja divisee au propagation_time_ref
						list_at_time_ref.append(l_ref_associated)
						list_at_time_flo.append(l_flo)

	return list_at_time_ref, list_at_time_flo

def labels_at_time_to_relabelled_correspondences(label_correspondences, lineage_ref, lineage_flo, time_ref, time_flo):
	"""
	Returns the dictionnary of relabelled correspondences given the original label correspondences at given times.
	"""
	relabelled_correspondences={}
	assert(lineage_ref.has_time(time_ref))
	assert(lineage_flo.has_time(time_flo))
	if type(label_correspondences)==dict:
		relabelled_correspondences = { lineage_ref.luts[time_ref][k] : lineage_flo.luts[time_flo][v] for k, v in label_correspondences.iteritems() }
	else:
		if type(label_correspondences)==list:
			relabelled_correspondences = { lineage_ref.luts[time_ref][label_correspondences[i]] : lineage_flo.luts[time_flo][label_correspondences[i+1]] for i in range(0,len(label_correspondences),2) }
		else:
			print "Unexpected label_correspondences type. Returning empty dictionnary..."
	return relabelled_correspondences

def relabelled_to_labels_at_time_correspondences(relabelled_correspondences, lineage_ref, lineage_flo, time_ref, time_flo):
	"""
	Returns the dictionnary of original label correspondences at given time considering the dictionnary of relabelled correspondences.
	"""
	labels_at_time={}
	assert(lineage_ref.has_time(time_ref))
	assert(lineage_flo.has_time(time_flo))
	if type(relabelled_correspondences)==dict:
		labels_at_time = { lineage_ref.reverse_luts[time_ref][k] : lineage_flo.reverse_luts[time_flo][v] for k, v in relabelled_correspondences.iteritems() }
	else:
		print "Unexpected relabelled_correspondences type. Returning empty dictionnary..."
	return labels_at_time


def lineages_relabelled_correspondences_propagation(lineage_ref, lineage_flo, relabelled_correspondences, label_ref, field_barycenter='real_barycenter', field_volume='real_volume', voxelsize=None, keep_divided_cells=True):
	"""
	Function for extraction of point-clouds of cell-to-cell correspondences after propagation of the given label_ref (in "relabelling mode").
	Parameters:
		lineage_ref/flo : reference and floating lineages
		relabelled_correspondences : one-to-one cell established cell mapping between lineages ref and flo (relabelled only) (type dict)
		label_ref : the one from reference that has to be propagated (as a relabelled one). Must belong to relabelled_correspondences.
	Optional parameters:
		field_barycenter : parameter for lineage_relabelled_point_cloud_at_time_from_list function 
		field_volume(='real_barycenter') : parameter for lineage_relabelled_point_cloud_at_time_from_list function
		voxelsize(='real_volume') : parameter for lineage_relabelled_point_cloud_at_time_from_list function
		keep_divided_cells(=True) : True if the user is expecting to keep relabels associations even if one (and only one) of the relabels undergoes some divisions (default) ; False otherwise
	Outputs:
		point_cloud_ref, point_cloud_flo, daughters_ref, daughters_flo
	The outputs daughters_ref and daughters_flo are corresponding to lists of two elements which are the 
	daughter cells of respectively label_ref and label_ref's corresponding label in lineage_flo with respect to relabelled_correspondences.
	"""
	assert type(label_ref)==int, "label_ref parameter is expected to be of type int, but is instead of type %s." % type(label_ref)
	assert lineage_ref.exists(label_ref), "Unexpected problem : label_ref = %d not found in lineage_ref." % int(label_ref)
	assert label_ref in relabelled_correspondences.keys(), "Unexpected problem : label_ref = %d not found in relabelled_correspondences keys." % int(label_ref)
	
	time_death_ref = lineage_ref.relabelled_death(label_ref)
	label_flo = relabelled_correspondences[label_ref]

	assert lineage_flo.exists(label_flo), "Unexpected problem : corresponding label_flo = %d computed from relabelled_correspondences not found in lineage_flo." % int(label_flo)

	
	assert label_flo in relabelled_correspondences.values(), "Unexpected problem : label_flo = %d not found in relabelled_correspondences values (wiered case)." % int(label_flo)

	time_death_flo = lineage_flo.relabelled_death(label_flo)	

	propagation_time_ref=time_death_ref+1
	propagation_time_flo=time_death_flo+1

	if not (lineage_ref.has_time(propagation_time_ref) and lineage_flo.has_time(propagation_time_flo)):
		#return barys_ref, barys_flo, established_correspondences, [], time_death_ref, time_death_flo
		return None, None, None, None


	relabelled_list_ref, relabelled_list_flo = get_relabelled_lists_at_time(relabelled_correspondences, lineage_ref, lineage_flo, propagation_time_ref, propagation_time_flo, keep_divided_cells=keep_divided_cells)
	assert not label_ref in relabelled_list_ref
	assert not label_flo in relabelled_list_flo

	daughters_ref = lineage_ref.relabelled_daughters(label_ref)
	daughters_flo = lineage_flo.relabelled_daughters(label_flo)

	relabelled_list_ref += daughters_ref
	relabelled_list_flo += daughters_flo

	point_cloud_ref=lineage_relabelled_point_cloud_at_time_from_list(lineage_ref, propagation_time_ref, relabelled_list_ref, filename=None, field_barycenter=field_barycenter, field_volume=field_volume, voxelsize=voxelsize)
	point_cloud_flo=lineage_relabelled_point_cloud_at_time_from_list(lineage_flo, propagation_time_flo, relabelled_list_flo, filename=None, field_barycenter=field_barycenter, field_volume=field_volume, voxelsize=voxelsize)

	return point_cloud_ref, point_cloud_flo, daughters_ref, daughters_flo


