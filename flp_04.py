import rdkit
from rdkit import Chem
import numpy as np
from rdkit.Chem import Draw
import Image
import numpy as np
import ast
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from pks_path import *

'''classes'''

'''whenever, there is something already existing in the class, such as the mol file, the smiles ....then the new function will be callable inside the 
function, but in case the function does not need to be part of the class, it will be independent'''

class pks_mol():
	'''pks_mol class - the minimum parameter required is a smiles, then the object can be used to create long paths trough the 
	molecule, identify staring points for the parts ...'''

	def __init__(self, input):
		'''initiation of the object'''
		try:
			self.mol = Chem.MolFromSmiles(input)
			self.smiles = Chem.MolToSmiles(Chem.MolFromSmiles(input))
		except:
			print 'input must be smiles(str)'

	def name(self, name):
		'''optional, can be useful for sorting'''
		self.name = name

	def draw(self):
		'''get an image of the molecule'''
		return(Draw.MolToImage(self.mol))	

	def draw_and_show(self):
		'''get an image of the molecule and direct output'''
		(Draw.MolToImage(self.mol)).show()

	def get_starting_matches(self, smarts_dicts = {}):

		'''this function gets all the needed starter units for the pks path algorithm, for self written smiles patterns use the custom function'''

		if smarts_dicts == {}:

			lactone = Chem.MolFromSmarts('[#6][#6](=O)[#8][#6]') 
			amide = Chem.MolFromSmarts('[#6][#6]([#7])=O')	
			carboxyl = Chem.MolFromSmarts('[#6][#6]([OH])=O')	
			carboxylate = Chem.MolFromSmarts('[#6]C([O-])=O') 	
			aldehyde = Chem.MolFromSmarts('C[!N&D2]=O')	#so far the best def. for an aldehyde which does not match other subs. (expl: an carbon, followed by any atom which 
			#has two explicit bonds but not an nitrogen, tested for ~500 compounds no false match so far)

			match_names = ['lactone','amide','carboxyl','carboxylate','aldehyde'] #name of the starting points
			smarts_list = [lactone,amide,carboxyl,carboxylate,aldehyde] #smarts
			#start_carbon_list = [1, 1, 1, 1, 1]	#number of carbon atom which shell start the path (numbering is oriented on the smarts patter ex.: aldehyde: 1-C-2-[D2]-3-=O), all 1 therefore became useless :)


		'''this part allows to use custom designed smarts patterns to find the longest paths'''

		if smarts_dicts != {}:

			match_names = [] #name of the starting points
			smarts_list = [] #smarts

			for key in smarts_dicts:
				match_names.append(key)
				smarts_list.append(smarts_dicts[key])



		match = {}
		idx = 0
		for smarts in smarts_list:
			substruct = self.mol.GetSubstructMatches(smarts)
			match[match_names[idx]] = substruct
			idx += 1 

		return(match)




		# self.matches = match
	def find_starting_C_atom(self):

		'''maybe implement an automatic algorithm, to fin the correct C atom, must be done for the custom function'''

		match = self.get_starting_matches()


		'''Takes the output of find_starting_matches and returns a list with the indices of the carbon atoms from the substructure starting points, 
		ex.: 'amide': ((3, 4, 6, 5), (10, 8, 6, 9), (19, 17, 16, 18)) ---> [4, 8, 17]'''
		C_list = []
		for key in match:
			for sub in match[key]:
				C_list.append(sub[1])	#can be changed if there is a substructure where the starting C-Atom is at a different position

		return(C_list)

	def path_recursion(self, no_peptide = True, **kwargs):

		'''this function basically does the same as the rdkit function FindAllPathsOfLengthN(), but
		it only gives the unique long paths instead of all paths of the length N, this will hopefully improve 
		the time management, as only a very limited number of paths is directly calculated by a very straight forward approach, 
		instead of first calculation all the paths and then setting conditions...''' 

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
							if n.GetIdx() in new_dics.keys():	#here the possibility is checked, that both paths connect, which
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


		start = self.find_starting_C_atom()	#get a list of start indices, to initiate the recursive function

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
		# print matrix
		# 	#exit()
		# 		#print matrix
		side_atom_idx = list(set(side_atom_idx))
		return(matrix, side_atom_idx)

	def get_pks_paths(self, path_list):
		'''produces for all the paths the matrices'''
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


	def matrix_self_filter_01(self, matrix_list):
		'''makes sure the same matrix is not multiple times created, could be the case for ring system'''
		for matrix in matrix_list:
			new_matrix_list = matrix_list
			counter = matrix_list.count(matrix)
			while counter > 1:
				new_matrix_list.remove(matrix)
				counter -= 1 

		return(matrix_list)

	def matrix_self_filter_02(self, matrix_dic):
		'''makes sure there are no matrices, that code for the same smiles'''
		smiles_dict = {} 

		for matrix in matrix_dic:
			mol, main_path = streight_path(matrix)
			smiles = Chem.MolToSmiles(mol, canonical=True)
			#print smiles
			smiles_dict[smiles] = matrix 

		# for key in smiles_dict:
		# 	print key
		matrix_list = smiles_dict.values() #as the smiles is the key, there can only be one value each. 

		return(matrix_list)


	def matrix_match(self, matrix):
		mol = self.mol
	#index = 0
	#for matrix in matrix_list:
		#if index < 10:
		matrix_mol, main_path = streight_path(matrix)
		mol_list = [mol, matrix_mol]
		comp = rdFMCS.FindMCS(mol_list)
		smart_match = comp.smartsString
		atom_match = comp.numAtoms
		bond_match = comp.numBonds
		Match = mol.GetSubstructMatches(Chem.MolFromSmarts(smart_match), useChirality=False)

		# img = draw_mol_and_path(smi1.mol, list(Match[0]))
		# img.show()	
		# 	#Draw.MolToImage(main_path).show()
		# index += 1
		return(list(Match[0]), atom_match, bond_match)


	# def matrix_self_filter_02(self, matrix_list, min_length):
	# 	'''makes sure the same matrix is not multiple times created, could be the case for ring systems'''
	# 	for matrix in matrix_list:
	# 		new_matrix_list = matrix_list
	# 		counter = matrix_list.count(matrix)
	# 		while counter > 1:
	# 			new_matrix_list.remove(matrix)
	# 			counter -= 1 

	# 	return(matrix_list)





'''independent functions'''

def pks_paths_only_matrix(pks_path_list):
	'''extracts the matrices from the pks_path_list, which are needed for the scoring_algo'''
	matrix_list = []
	for k in pks_path_list:
		matrix_list.append(k.matrix)

	return(matrix_list)


def no_sub_smiles(smiles_list_comp):
	'''takes a list of smiles and compares them with each other, when one smile is a subsmile to a bigger smile in the list, the small smiles will be deleted. 
	Subsequently only main smiles remain which are not subsmiles themselves, Example:
	List_init = [CCC(C), CCC, CC(O)C] ==> List_after = [CCC(C), CC(O)C]'''

	unique_smiles = []
	for smiles in smiles_list_comp:
		copy = smiles_list_comp
		#copy.remove(smiles)
		mol1 = Chem.MolFromSmiles(smiles)
		match_count = 0
		#print copy
		for smi in copy:
			mol2 = Chem.MolFromSmiles(smi)
			comp_liste = [mol1, mol2] 
			comp = rdFMCS.FindMCS(comp_liste)
			smart_match = comp.smartsString
			comp_match = Chem.MolToSmarts(Chem.MolFromSmiles(smiles))
			if comp_match == smart_match:
				match_count += 1 #the first match is the compound itself, if it finds another match besides itself, it is not the biggest unique path !!

		if match_count == 1:
			unique_smiles.append(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))

	unique_smiles = list(set(unique_smiles))
	return(unique_smiles)


def streight_path(dic_path):
	'''Takes a path in the matrix format and creates a mol object, like the inverse of the smile2matrix program, but side paths which would be involved in a ring for example 
	get shown in the way the comparison algo looks at it , ex.: 
		C-O-C		smile2matrix														streight_path		  O-C   C-O
		|	|																									|	|
	C-C-C-C-C-C-C		---->		{0: ['C', 'SINGLE', []], 1: ['C', 'SINGLE' []]...... 	---->			C-C-C-C-C-C-C

	Additionally the mol object is linearised, meaning, that the main_path is taken as a straight template and the entire path is oriented the same way:
	See.: http://www.rdkit.org/docs/GettingStartedInPython.html; Working with 2D molecules: Generating Depictions
	'''

	#dic_path = ast.literal_eval(str(matrix))

	for key in dic_path:
		# print key

		if int(key) == 0:
			last_atom_idx = 0
			last_atom_idx2 = 0
			molecule = Chem.MolFromSmiles(str(dic_path[key][0]))
			molecule = Chem.RWMol(molecule)
			main_path = Chem.RWMol(molecule)
			#print molecule.GetAtomWithIdx(0).GetIdx()
			# first atom index = 0

			#check side paths
			if dic_path[key][2] != []:
				side_paths = dic_path[key][2]
				#print side_paths
				for path in side_paths:
					k = molecule.AddAtom(Chem.Atom(str(path[1])))
					molecule.AddBond(0, k, rdkit_bond_from_str(str(path[0])))
					if path[2] != []:
						#print 'k'
						side_sidepath = path[2]
						for path2 in side_sidepath:
							k2 =  molecule.AddAtom(Chem.Atom(str(path2[1])))
							molecule.AddBond(k, k2, rdkit_bond_from_str(str(path2[0])))



		else:
			index_atom = molecule.AddAtom(Chem.Atom(str(dic_path[key][0])))
			molecule.AddBond(last_atom_idx, index_atom, rdkit_bond_from_str(str(dic_path[int(key)-1][1])))

			index_atom2 = main_path.AddAtom(Chem.Atom(str(dic_path[key][0])))
			# print last_atom_idx
			# print index_atom2
			# exit()
			main_path.AddBond(last_atom_idx2, index_atom2, rdkit_bond_from_str(str(dic_path[int(key)-1][1]))) 


			if dic_path[key][2] != []:
				side_paths = dic_path[key][2]
				#print side_paths
				for path in side_paths:
					k = molecule.AddAtom(Chem.Atom(str(path[1])))
					molecule.AddBond(index_atom, k, rdkit_bond_from_str(str(path[0])))
					if path[2] != []:
						#print 'k'
						side_sidepath = path[2]
						for path2 in side_sidepath:
							k2 =  molecule.AddAtom(Chem.Atom(str(path2[1])))
							molecule.AddBond(k, k2, rdkit_bond_from_str(str(path2[0])))

			last_atom_idx = index_atom
			last_atom_idx2 = index_atom2

	
	mol =  Chem.RemoveHs(molecule.GetMol())
	main_path = Chem.RemoveHs(main_path.GetMol())
	AllChem.Compute2DCoords(main_path)
	AllChem.GenerateDepictionMatching2DStructure(mol ,main_path)

	return(mol, main_path)


def streight_path_all(liste, subsmiles = None):
	mol_list = []
	for elem in liste:
		mol, s_path = streight_path(elem)
		#smiles = Chem.MolToSmiles(mol)
		mol_list.append(mol)

	#if subsmiles == None:
	#smiles_list = no_sub_smiles(smiles_list)
		
	return(mol_list)

def streight_path_all_mol(liste):
	mol_list = []
	for elem in liste:
		mol, s_path = streight_path(elem)
		#smiles = Chem.MolToSmiles(mol)
		mol_list.append(mol)
	return(mol_list)


def rdkit_bond_from_str(string):
	'''as the Chem.BondTypes can not be initialized in python, see rdkit website, http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Bond-class.html:
				__init__(...)
			(Constructor)
				 

			Raises an exception
			This class cannot be instantiated from Python

			Overrides: object.__init__
		this little program helps to create them anyway !! More bond types can be added here, but for my purpose these where enough.'''

	if string == 'SINGLE':
		bondtype = Chem.BondType.SINGLE
	if string == 'DOUBLE':
		bondtype = Chem.BondType.DOUBLE
	if string == 'AROMATIC':
		bondtype = Chem.BondType.UNSPECIFIED #no idea what to do with the aromatic bonds at the moment 
	if string == 'TRIPLE':
		bondtype = Chem.BondType.TRIPLE
	return(bondtype)


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

	# 	# print color_dict
	# 	index += 1
	# for atom in side_atom_idx:
	# 	color_dict[atom] = [0,0,1]

		#atom_list.append(atom)
	#print color_dict
	# print mol
	# print color_dict

	Draw.MolToFile(mol, file_path, size=(600,600), highlightMap = color_dict, type = 'svg')


def list_change(liste):
	'''changes a list with direction marker (+) to a normal list, which is needed for plotting or check up'''
	news = []
	for k in liste:
		if type(k) == type(()):
			k = k[0]
		news.append(k)
	return(news)

def atom_from_smiles(smiles):
	'''creates an rdkit atom object with idx. 0 from a one latter smiles'''
	mol1 = Chem.MolFromSmiles(smiles)
	atom1 = mol1.GetAtomWithIdx(0)
	return(atom1)

def list2mol(mol, liste):
	'''takes a mol object and a list with atom indices contained in the mol object and creates a mol for only the list atoms'''
	atom_list = []
	bond_list = []
	
	index = 0
	for num in liste:
		current_atom = mol.GetAtomWithIdx(liste[index])
		if index > 0:
			last_atom = mol.GetAtomWithIdx(liste[index-1])
			bond_list.append(mol.GetBondBetweenAtoms(current_atom.GetIdx(), last_atom.GetIdx()).GetBondType())
		index += 1
		atom_list.append(current_atom.GetSymbol())

	'''init the RWMol object'''
	mol = Chem.Mol()
	emol = Chem.RWMol(mol)

	'''fill RWMol object with atoms and bonds and save the idx in a dic for later'''
	atom_dic = {}

	#get the atoms
	for ato in atom_list:
		ato = atom_from_smiles(ato)
		idx1 = emol.AddAtom(ato)
		atom_dic[idx1] = ato
	#get the bonds
	index = 0
	for bon in bond_list:

		emol.AddBond(index, (index + 1), bon)
		index += 1

	mol = emol.GetMol()
	return(mol, atom_dic)


###debug### 

#uncomment step by step to learn more about the class.
#all the major functionalists are demonstrated.

# smiles = 'CO[C@H]1C[C@H](O[C@H]2[C@@H](C)/C=C/C=C/3\\CO[C@H]4[C@]3(O)[C@@H](C=C([C@H]4O)C)C(=O)O[C@H]3C[C@@H](C/C=C/2\\C)O[C@]2(C3)C=C[C@@H]([C@H](O2)[C@H](CC)C)C)O[C@H]([C@@H]1O[C@H]1C[C@H](OC)[C@H]([C@@H](O1)C)O)C'

# #initiate the object
# object1 = pks_mol(smiles)

# #show it
# object1.draw_and_show()

# #all the paths with starting idx as key (O and S are not in the path):
# k = object1.path_recursion(exclude = ['O','S'])

# #create a figure with the path inside, can also be donw with iteration
# img_with_path = draw_mol_and_path(object1.mol, k[24][-1])

# #show the figure
# img_with_path.show()

# #matrix, returns the matrix description of a path which can be used to compare paths to each other (like a fingerprint), the 
# #documentation for the comparison synthax and the matrix will follow soon,
# #also the side indices will be returned in order to allow to color them different for example
# matrix, side_idx = object1.path2matrix_3(k[24][-1])
# print matrix
# print side_idx
# img_with_path = draw_mol_and_path(object1.mol, side_idx)
# img_with_path.show() #unfortunately rdkit forgets to color single c atoms at the moment.


# # the paths created can be transformed into path objects, which allows more operations (see the pks_path class)
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

