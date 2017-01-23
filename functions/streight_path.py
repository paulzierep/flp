from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np

def rdkit_bond_from_str(string):
	'''as the Chem.BondTypes can not be initialized in python, see rdkit webside, http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Bond-class.html:
				__init__(...)
			(Constructor)
				 

			Raises an exception
			This class cannot be instantiated from Python

			Overrides: object.__init__
		this little programm helps to create them anyway !! More bond typen can be added here, but for my purpose these where enough.'''

	if string == 'SINGLE':
		bondtype = Chem.BondType.SINGLE
	if string == 'DOUBLE':
		bondtype = Chem.BondType.DOUBLE
	if string == 'AROMATIC':
		bondtype = Chem.BondType.UNSPECIFIED #no idea what to do with the aromaric bonds at the moment 
	if string == 'TRIPLE':
		bondtype = Chem.BondType.TRIPLE
	return(bondtype)


def streight_path(dic_path):
	'''Takes a path in the matrix format and creates a mol object, like the inverse of the smile2matrix program, but side paths which would be involved in a ring for example 
	get shown in the way the comparison algo looks at it , ex.: 
		C-O-C		smile2matrix														streight_path		  O-C   C-O
		|	|																									|	|
	C-C-C-C-C-C-C		---->		{0: ['C', 'SINGLE', []], 1: ['C', 'SINGLE' []]...... 	---->			C-C-C-C-C-C-C

	Additionally the mol object is linearized, meaning, that the main_path is taken as a streight tamplate and the entire path is oriented the same way:
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
