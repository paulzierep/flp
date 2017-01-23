from rdkit import Chem
from rdkit.Chem import RWMol
from rdkit.Chem import Atom
from rdkit.Chem import Draw

# class Atom(Atom):
# 	def __init___(self):
# 		super(Atom, self).__init__()
# 		self.old_idx = 0

def get_mol_from_query(mol, query):
	'''
	as there is often only an output of indices (ex. substructurematch), this function creates an new mol object based on the 
	indices inside the original mol object, if there is a bond to an atom which is not inside the query it is removed.
	The dict was useless as output, as the indices are automatically changed by the rdkit GetMol() synthax,
	therefore I created new property instances for the atoms and bonds stroed in the new molecule, by using GetProp('old_idx') those 
	can be extracted.
	'''

	# smi = Chem.MolToSmiles(mol)
	# print smi

	# mol = Chem.MolFromSmiles(smi)
	#Draw.MolToImage(mol, size = (1000,200)).show()

	mol_for_RW = Chem.MolFromSmiles('') 
	new_mol = RWMol(mol_for_RW)
	atom_index_translater = {}

	'''tranfers all the atoms from the query into the new mol object'''
	for atom in mol.GetAtoms():
		if atom.GetIdx() in query:
			atom.SetProp('old_idx', str(atom.GetIdx()))
			RW_idx = new_mol.AddAtom(atom) #creates the new mol object with identical information about the atom, the new index must be save,
			# as this will be different for the new molecule, the correlation ins saved in atom_index_translater
			atom_index_translater[atom.GetIdx()] = RW_idx #this dict correlates the index of the original molecule with the new made mol object

	#print atom_index_translater
	'''add the correct bonds to the new mol object'''
	#bond_index_translater = {}
	#print atom_index_translater
	for old_idx in atom_index_translater: #go through all the idxs of the old mold
		neighbors = mol.GetAtomWithIdx(old_idx).GetNeighbors() #get the neighbors of the atom
		for neigh in neighbors: #go through all the neighbors
			if neigh.GetIdx() in atom_index_translater.keys(): #if the neighbor is also in the new mol it needs a bond
				bond_type = mol.GetBondBetweenAtoms(old_idx, neigh.GetIdx()).GetBondType() #get the bond information
				#print bond_type 	
				new_idx_atom = atom_index_translater[old_idx] #get the idx where the bond belongs
				#print new_idx_atom
				new_idx_neigh = atom_index_translater[neigh.GetIdx()] #a bond always need two atoms (logic)
				#print new_idx_neigh
				if new_mol.GetBondBetweenAtoms(new_idx_atom, new_idx_neigh) == None: #better then try function !
					old_bond_idx = mol.GetBondBetweenAtoms(old_idx, neigh.GetIdx()).GetIdx()
					#bond.SetProp('old_idx', str(old_bond_idx))
					new_mol.AddBond(new_idx_atom, new_idx_neigh, bond_type) #add the bond to the atoms, unfortunatly it is not possible to jsut add a 
					#bondtype object, therefore in the next code part, the bond is again adressed in order to add the new property.
					idx = 0
					for bond in new_mol.GetBonds():
						if bond.GetBeginAtomIdx() == new_idx_atom and bond.GetEndAtomIdx() == new_idx_neigh:
							new_mol.GetBonds()[idx].SetProp('old_idx', str(old_bond_idx))
						idx += 1

				else:
					pass

	new_mol = new_mol.GetMol() #create the new mol object from the RWMol object, new idxs stay
	return(new_mol) #atom_index_translater, bond_index_translater)

# '''for debugging'''
# from rdkit import Chem
# from rdkit.Chem import Draw
# from rdkit.Chem import rdmolops


# #smiles_temp = 'O=C(O)[C@@H]3[C@@H](O)C[C@@]2(O)C[C@@H](O)C[C@@H](O)[C@H](O)CC[C@@H](O)C[C@@H](O)CC(=O)O[C@@H](C)[C@H](C)[C@H](O)[C@@H](C)C=CC=CC=CC=CC=CC=CC=C[C@H](O[C@@H]1O[C@H](C)[C@@H](O)[C@H](N)[C@@H]1O)C[C@@H]3O2'
# smiles_temp = 'COCC=C(CC)CCC'
# #smiles_temp = 'C=CCOCC=C'

# mol = Chem.MolFromSmiles(smiles_temp)
# rdmolops.RemoveStereochemistry(mol)
# query = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

# new_mol = get_mol_from_query(mol, query)

# for atom in new_mol.GetAtoms():
# 	print atom.old_idx
# exit()


# Draw.MolToImage(new_mol, size = (1000, 1000)).show()
# print atom_index_translater