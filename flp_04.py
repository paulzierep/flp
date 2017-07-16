########################
# written by Paul Zierep
########################

import rdkit
from rdkit import Chem
import numpy as np
from rdkit.Chem import Draw
from PIL import Image
import ast
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from pks_path import *

'''classes'''

class pks_mol():
	'''pks_mol class - the minimum parameter required is a smiles, then the object can be used to create long paths trough the 
	molecule, identify staring points for the parts, possbile outputs are shown in the examples below.'''

	def __init__(self, input):
		'''initiation of the object'''
		try:
			#creat the smiles and mol object of the input
			self.mol = Chem.MolFromSmiles(input)
			self.smiles = Chem.MolToSmiles(Chem.MolFromSmiles(input))
		except:
			print 'input must be smiles(str)'

	def name(self, name):
		'''optional, can be useful for sorting'''
		self.name = name

	def draw(self, size = (300,300)):
		'''get an image of the molecule'''
		return(Draw.MolToImage(self.mol, size))

	def draw_and_show(self, size = (300,300)):
		'''get an image of the molecule and direct output'''
		(Draw.MolToImage(self.mol, size)).show()

	def get_starting_matches(self, smarts_list = []):
		#smarts_list are a list of set that has all the needed information for a starting unit (name, smarts_pattern, idx of start atom), example: [(lactone, '[#6][#6](=O)[#8][#6]', 1), (amide, '[#6][#6]([#7])=O', 1), ...]

		'''this function gets all the needed starter units for the pks path algorithm, for self written smiles patterns use the custom function'''

		if smarts_list == []:

			lactone = Chem.MolFromSmarts('[#6][#6](=O)[#8][#6]') 
			amide = Chem.MolFromSmarts('[#6][#6]([#7])=O')	
			carboxyl = Chem.MolFromSmarts('[#6][#6]([OH])=O')	
			carboxylate = Chem.MolFromSmarts('[#6]C([O-])=O') 	
			#aldehyde = Chem.MolFromSmarts('C[!N&D2]=O')	#so far the best def. for an aldehyde which does not match other subs. (expl: an carbon, followed by any atom which 
														#has two explicit bonds but not an nitrogen, tested for ~500 compounds no false match so far)

			aldehyde = Chem.MolFromSmarts('[CX3H1:1]=[OX1;!$(O=C~[!#1!#6]):2]') #better definition from https://www.chemaxon.com/forum/ftopic3049.html

			match_names = ['lactone','amide','carboxyl','carboxylate','aldehyde'] #name of the starting points
			smart_list = [lactone,amide,carboxyl,carboxylate,aldehyde] #smarts
			start_idx = [1,1,1,1,0]
			
			#(Draw.MolToImage(aldehyde, size = (300,300))).show()

		'''this part allows to use custom designed smarts patterns to find the longest paths'''

		if smarts_list != []:

			match_names = [] #name of the starting points
			smart_list = [] #smarts
			start_idx = []

			for sets in smarts_list:
				match_names.append(sets[0])
				smart_list.append(Chem.MolFromSmarts(sets[1]))
				start_idx.append(sets[2])


		match = {}
		idx = 0
		for smarts in smart_list:

			substruct = self.mol.GetSubstructMatches(smarts)
			match[match_names[idx]] = substruct, start_idx[idx]
			idx += 1 

		return(match)


	def find_starting_atom(self):

		#function which returns the atom idx in the mol of the starting atom for a given starting unit.

		match = self.get_starting_matches()

		'''Takes the output of find_starting_matches and returns a list with the indices of the carbon atoms from the substructure starting points, 
		ex.: 'amide': ((3, 4, 6, 5), (10, 8, 6, 9), (19, 17, 16, 18)) ---> [4, 8, 17]'''
		start_atom_list = []
		for name in match:
			for sub in match[name][0]:
				start_atom_list.append(sub[match[name][1]]) #is based on the idx supplied in get_starting_matches

		return(start_atom_list)

	def path_recursion(self, no_peptide = True, **kwargs):

		'''this function has a more specific design as the rdkit FindAllPathsOfLengthN() function, but
		it only gives the unique longest paths instead of all paths of the length N, it is also reasonable fast due to the 
		recursive implemention. Still large ring systems can take w long time to process.''' 

		def list_change(liste):
			'''changes a list with direction marker (+) to a normal list, which is needed for plotting or check up'''
			news = []
			for k in liste:
				if type(k) == type(()):
					k = k[0]
				news.append(k)
			return(news)

		def recu(dics):
			'''recursive function, this function is supposed to run until the final end condition is found, how this is achieved is explained where the
			function is called'''

			unique_paths = []	#for every run the unique paths which can not be further explored are stored in this list.
			new_dics = {}		# this dic will contain the last chain atom as key and the already explored path as value.

			for k in dics:	# go through the atoms

				#ring system, back conversion
				k1 = k
				#print 'k:', k
				if type(k) == type(()):
					#print 'hi'
					k1 = k[0]

				neig = self.mol.GetAtomWithIdx(k1).GetNeighbors() # get the neighbors
				#print dics[k]

				#counter = False
				counter = 0
				for n in neig:	#go through the neighbors

					'''the implementation of the tuple, makes it possible to go along a path from both sides of a circle,
					but it would not work for a third dimension'''
					'''list change, changes the tuple annotation back to normal, see function, this is used for the
					comparison and the output as unique path (*)''' 
					'''now special conditions regarding the following atom can be defined, such as, what type of 
					atom it can be, or even its own neighbors'''
					
					'''possible exclusion elements can be inserted directly in the recu function, speed enhancement (*),
					if the element is entered as string, the path will not be continued ones is appears a neighbor'''

					scip = False

					if kwargs:
						if 'exclude' in kwargs.keys():

							exclude = kwargs['exclude']
							'''the given list of elements to exclude is checked, if the symbol matches it is skipped'''
							for elem in exclude:
								if n.GetSymbol() == elem:
									scip = True


					'''the following should cover peptide bonds. But also some strange cases which do not necessarily count as peptide bond.
					Everything with a -N-C=O structure '''


					if no_peptide == True:
						if n.GetSymbol() == 'N':
							for n2 in n.GetNeighbors():
								if n2.GetSymbol() == 'C' and \
								str((self.mol.GetBondBetweenAtoms(n.GetIdx(), n2.GetIdx())).GetBondType()) == 'SINGLE':
									for n3 in n2.GetNeighbors():
										if n3.GetSymbol() == 'O' and \
										str((self.mol.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())).GetBondType()) == 'DOUBLE':
											scip = True


					if scip == False:
						'''*'''
						'''the introduction of the tuples with the +, ++ and +++ are needed to give distinctive keys for paths which meet at the same time at one loop,
						as the valance for C atoms is max. 4, this three special cases are enough to take care of that matter, therefore even very complex cyclic systems should 
						be covered. The system was tested with the following smiles: C(=O)CCCC12(CC(C)(C2)C1), which yielded in all possible unique paths''' 
						if n.GetIdx() not in list_change(dics[k]): # (*) check if the neighbor was already in the old_atoms list
							if n.GetIdx() in new_dics.keys():	#here the possibility is checked, that both paths connect
								if (n.GetIdx(),'+') in new_dics.keys():
									if (n.GetIdx(),'++') in new_dics.keys():
										new_list = dics[k] + [k] #concatenate start atom to the old atom list
										new_dics[n.GetIdx(),'+++'] = new_list #dummy index is introduced for ring system
									else :
										new_list = dics[k] + [k] #concatenate start atom to the old atom list
										new_dics[n.GetIdx(),'++'] = new_list #dummy index is introduced for ring system
								else:
							#is the case for ring systems
									new_list = dics[k] + [k] #concatenate start atom to the old atom list
									new_dics[n.GetIdx(),'+'] = new_list #dummy index is introduced for ring system

							else:

								new_list = dics[k] + [k] #concatenate start atom to the old atom list
								new_dics[n.GetIdx()] = new_list #create new dict, with new atoms as keys and the new old atom list as values
							#print 'n', new_dics

						else:
							counter += 1
					else:
						counter += 1

					if counter == len(neig): #if all the neighbors have been already in the old list, this path comes to an end and can therefore be added to the unique path list
						unique_paths.append(list_change(dics[k] + [k])) #(*)'
						#print 'unique:', unique_paths
						#print 'up:', dics[k] + [k] 

			return(new_dics, unique_paths)


		#Pos 1
		start = self.find_starting_atom()	#get a list of start indices, to initiate the recursive function

		longest_paths = {}	# this will be the dic where for each starting atom, a list of lists will be stored which contains 
							# all the unique atom chains starting from this atom

		for sta in start:	# iter trough all the starting atoms.
			dics = {} 		# the first dict contains only the start atom as key and an empty list as value
			dics[sta] = [] 

			unique_p = []	#

			while dics != {}: #when all the paths are run through the recu function will return an empty dict.

				dics, unique_paths = recu(dics) #run the function
				unique_p = unique_p + unique_paths #as the unique path list, gets refreshed every time the recu function runs, the 
				#lists must be stored outside the function.
				#print dics

			longest_paths[sta] = unique_p #the recu function is used for every starting atom, each will be stored in this list


		return(longest_paths)

	def path2matrix_3(self, path):

		'''this path takes into account the side paths of the chain, up to the second order, best way for comparison.'''

		'''takes a path (basically a list with the indices of the atoms in the paths) and a Mol object which needs to contain the indices from the path and 
		 creates the so called matrix, a dict with up-counting sequence of ints as keys and atomic symbols and bonds as values. 
		 Important: The keys do not correspond to
		 the indices of the atoms in the mol object, but only number the sequence of the path (1 - first atom in path, 2 - second atom in path ...),
		 this is by purpose so that independent molecules can be compared later on. The values are the atomic symbol, a list with symbols of the atomic neighbors of the key atom, 
		 the bonds in the second list are the bonds between the key atom and its neighbors, ex.: 0: ('C', ['C'], ['SINGLE'])  '''

		mol = self.mol
		
		matrix = {}
		idx = 0
		#print path
		#exit()
		length = np.linspace(1, len(path), len(path))
		#print length
		#print length
		side_atom_idx = []

		for index in length:
			index = int(index)-1

			if path[index] != path[-1]: #if the atom is not the last in the path, there is a bond 
				current_atom = mol.GetAtomWithIdx(path[index])

				next_atom = mol.GetAtomWithIdx(path[index+1])

				bond = str((mol.GetBondBetweenAtoms(current_atom.GetIdx(), next_atom.GetIdx())).GetBondType())

			else:	#if the atom is the last in path: 
				current_atom = mol.GetAtomWithIdx(path[index])

				bond = 'END'

			neighbors = list(current_atom.GetNeighbors()) #get the neighbors
			# for atom in neighbors: 
			# 	print atom.GetSymbol(), atom.GetIdx()
			# 	continue
			# exit()
			neighbor_list = []

			#print path
			for atom in neighbors: 
				# print atom.GetSymbol(), atom.GetIdx()
				# break
				if atom.GetIdx() in path: #if the atom is already in the path it is not needed a neighbors information (ring systems)
					# print atom.GetSymbol()
					continue
				else:
					side_atom_idx.append(atom.GetIdx())
					#print "in", atom.GetSymbol(), atom.GetIdx()
					# break
					k = Chem.FindAllPathsOfLengthN(mol, 2, useBonds=False, rootedAtAtom=atom.GetIdx())

					tuples = []

					for tuple1 in k:

						s_main = []
						for subtuble1 in tuple1:

							s_main.append(subtuble1)

						if len(s_main) > len(set(s_main)):	#there seems to be a bug in FindAllPathsOfLengthN or at least I am not content with the result, ring atoms are taken twice, which is not good, therefore this removes non unique lists.
					 		pass
					 	else:
					 		tuples.append(s_main)

					#print tuples 
					neighbour_bond = (str(mol.GetBondBetweenAtoms(atom.GetIdx(), current_atom.GetIdx()).GetBondType()))

					second_n = []
					for k in tuples:
						for atom1 in k:
							if atom1 != atom.GetIdx():
								second_n.append(atom1)


					second_n_inf =[]
					for atom2 in second_n:
						if atom2 in path:
							continue
						else:
							side_atom_idx.append(atom2)
							inf = []
							inf.append(str(mol.GetBondBetweenAtoms(atom2, atom.GetIdx()).GetBondType()))
							inf.append(mol.GetAtomWithIdx(atom2).GetSymbol())
							second_n_inf.append(inf)
							
					
					neighbor_list_it = [neighbour_bond, atom.GetSymbol(), second_n_inf]
					neighbor_list.append(neighbor_list_it)

			matrix[index] = [current_atom.GetSymbol(), bond, neighbor_list]


		side_atom_idx = list(set(side_atom_idx))
		return(matrix, side_atom_idx)

	def get_pks_paths(self, path_list):
		'''produces for all the paths the path objects'''
		pks_path_list = []
		#print 'hi'
		for key in path_list:
			for path in path_list[key]:
				#print 'hi'
				matrix, sidepath = self.path2matrix_3(path)
				#print matrix
				path_cls = pks_path(self.name, self.mol, matrix, key, path, sidepath)
				#print path_cls
				pks_path_list.append(path_cls)
				#matrix_list.append(matrix_tuple)

		#matrix_dic = self.matrix_self_filter_02(matrix_dic) #for me it seems reasonable to apply this filter all the time !!
		return(pks_path_list)

'''independent functions'''

def draw_mol_and_path(mol, path):

	Draw.DrawingOptions.elemDict={} #elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
	color_dict={}
	#atom_list = []
	'''the following colors the path step-by-step from red to green'''
	color_change = (np.linspace(0, 1, len(path)))
	red_change = list(color_change)
	green_change = list(color_change[::-1])

	index = 0
	for atom in path:
		color_dict[atom] = [red_change[index],green_change[index],0]
		# print color_dict
		index += 1
	# for atom in side_atom_idx:
	# 	color_dict[atom] = [0,0,1]

	img = Draw.MolToImage(mol, size=(600,600), highlightMap = color_dict, type='svg')
	return(img)

def draw_mol_and_path_save(mol, path, file_path):
	'''drawing function, creates a figure with the molecule and the path inside,
	 should be modified, so that more options are available, ex. begin of path'''
	# Draw.DrawingOptions.elemDict={}#elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
	# color_dict={}
	#atom_list = []
	'''the following colors the path step-by-step from red to green'''
	#color_change = (np.linspace(0, 1, len(path)))
	#red_change = list(color_change)
	#green_change = list(color_change[::-1])
	color_dict = {}
	
	index = 0
	for atom in path:
		if index == 0:
			color_dict[atom] = [0,0,1]
		else:
			color_dict[atom] = [0,1,0]
	# print color_dict
		index += 1

	Draw.MolToFile(mol, file_path, size=(600,600), highlightMap = color_dict, type = 'svg')



##debug and explanatio### 

# #uncomment step by step to learn more about the class.
# #all the major functionalists are demonstrated.

# avermectin_smi = 'CC[C@H](C)[C@@H]1[C@H](C=C[C@@]2(O1)C[C@@H]3C[C@H](O2)C/C=C(/[C@H]([C@H](/C=C/C=C/4\CO[C@H]5[C@@]4([C@@H](C=C([C@H]5O)C)C(=O)O3)O)C)O[C@H]6C[C@@H]([C@H]([C@@H](O6)C)O[C@H]7C[C@@H]([C@H]([C@@H](O7)C)O)OC)OC)\C)C'
# avermectin_smi = 'C1C(CC2(CC3))C(CC2(CC3))CCC1(CCCC=O)'

# #initiate the object
# object1 = pks_mol(avermectin_smi)

# #show it
# object1.draw_and_show(size = (500,500))

# # starting units 
# print object1.get_starting_matches()

# # starting atoms
# print object1.find_starting_atom()

# # all paths for this object (O and S are not in the path)
# k = object1.path_recursion(exclude = ['O','S'])
# print k

# # # draw an example
# for key in k:
# 	longst_p = sorted(k[key], key = lambda x: len(x), reverse=True)
# 	longst_p = longst_p[:3]

# 	for lis in longst_p:
# 		draw_mol_and_path(object1.mol, lis).show()

	

# # #matrix, returns the matrix description of a path which can be used to compare paths to each other (like a fingerprint), the 
# # #documentation for the comparison synthax and the matrix will follow soon,
# # #also the side indices will be returned in order to allow to color them different for example
# matrix, side_idx = object1.path2matrix_3(k[34][-1])
# print matrix
# img_with_path = draw_mol_and_path(object1.mol, side_idx)

# img_with_path.show() #unfortunately rdkit forgets to color single c atoms at the moment.


# # # the paths created can be transformed into path objects, which allows more operations (see the pks_path class)
# path_objects = object1.get_pks_paths(k)
# print path_objects #a list of all the pks_path objects for this pks_mol object

# path_object = path_objects[-1]

# print	path_object.name #name of the object, not assigned in this case
# print	path_object.master_mol #mol object of the molecule where the path is assigned to
# print	path_object.matrix #matrix notation of the path
# print	path_object.path_start #idx of the starting atom
# print	path_object.main_path #idx of the main path trough the molecule
# print	path_object.side_path #idx of the side paths
# print	path_object.mol #straight oriented mol object of the path
# print	path_object.mol_main_path #mol object of the core path, which is used for the orientation


# img1 = Draw.MolToImage(path_object.mol, size=(600,600), type='svg')
# img2 = Draw.MolToImage(path_object.mol_main_path, size=(600,600), type='svg')

# img1.show()
# img2.show()

