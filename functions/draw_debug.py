from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
import rdkit

'''only for the comp_with_origin function'''
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdFMCS

'''manipulatable mol class'''
from rdkit.Chem import RWMol

'''self made functions'''
from mol_from_query import get_mol_from_query
from change_bond_type import change_bond_type

save_elem_dict = Draw.DrawingOptions.elemDict

def draw_mol_and_path2(mol, path, file_path, size):

    Draw.DrawingOptions.elemDict= {} #elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
    color_dict={}
    #atom_list = []
    '''the following colors the path step-by-step from red to green'''

    index = 0
    for atom in path:
        if index == 0:
            color_dict[atom] = [0,0,1]
        else:
            color_dict[atom] = [0,1,0]
		# print color_dict
        index += 1
	# for atom in side_atom_idx:
	# 	color_dict[atom] = [0,0,1]

	#print color_dict
	Draw.MolToFile(mol, file_path, size=(size[0],size[1]), highlightMap = color_dict, type='svg')
	# img = Draw.MolToImage(mol, size=(size[0],size[1]), highlightMap = color_dict, type='svg', fitImage=True)
	# return(img)


def draw_mol_and_path(mol, path, file_path, size):

    Draw.DrawingOptions.elemDict={}#elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
    color_dict={}
    #atom_list = []
    '''the following colors the path step-by-step from red to green'''
    color_change = (np.linspace(0, 1, len(path)))
    red_change = list(color_change)
    green_change = list(color_change[::-1])


    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            color_dict[atom.GetIdx()] = [1,0,0]
        if atom.GetSymbol() == 'N':
            color_dict[atom.GetIdx()] = [0,0,1]
        if atom.GetSymbol() == 'S':
            color_dict[atom.GetIdx()] = [1,1,0]

    index = 0
    for atom in path:
        if index == 0:
            color_dict[atom] = [0,0,1]
        else:
            color_dict[atom] = [0,1,0]
        # print color_dict
        index += 1
    # for atom in side_atom_idx:
    #   color_dict[atom] = [0,0,1]

    #print color_dict
    Draw.MolToFile(mol, file_path, size=(size[0],size[1]), highlightMap = color_dict, type='svg')
    # img = Draw.MolToImage(mol, size=(size[0],size[1]), highlightMap = color_dict, type='svg', fitImage=True)
    # return(img)

def draw_mol_and_path3(mol, path, file_path, size):

    Draw.DrawingOptions.elemDict={}#elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
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
    #   color_dict[atom] = [0,0,1]

    #print color_dict
    Draw.MolToFile(mol, file_path, size=(size[0],size[1]), highlightMap = color_dict, type='svg')
    # img = Draw.MolToImage(mol, size=(size[0],size[1]), highlightMap = color_dict, type='svg', fitImage=True)
    # return(img)



def draw_mol(mol, path, size):

	Draw.DrawingOptions.elemDict = save_elem_dict

	Draw.MolToFile(mol, path, size=(size[0],size[1]), type='svg')



# def draw_mol_and_path(mol, path, size):

# 	Draw.DrawingOptions.elemDict={}#elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
# 	color_dict={}
# 	#atom_list = []
# 	'''the following colors the path step-by-step from red to green'''
# 	color_change = (np.linspace(0, 1, len(path)))
# 	red_change = list(color_change)
# 	green_change = list(color_change[::-1])

# 	index = 0
# 	for atom in path:
# 		color_dict[atom] = [red_change[index],green_change[index],0]
# 		# print color_dict
# 		index += 1
# 	# for atom in side_atom_idx:
# 	# 	color_dict[atom] = [0,0,1]

# 	#print color_dict
# 	img = Draw.MolToImage(mol, size=(size[0],size[1]), highlightMap = color_dict, type='svg', fitImage=True)
# 	return(img)

def draw_mol_and_path_save(mol, file_path, size):
	'''drawing function, creates a figure with the molecule and the path inside,
	 should be modified, so that more options are available, ex. begin of path'''
	Draw.DrawingOptions.elemDict={}#elemDict must be set to empty dict, otherwise the drawing is wrong. PZ.
	color_dict={}
	#atom_list = []
	'''the following colors the path step-by-step from red to green'''
	for atom in mol.GetAtoms():
		if atom.GetSymbol() == 'O':
			color_dict[atom.GetIdx()] = [1,0,0]
		if atom.GetSymbol() == 'N':
			color_dict[atom.GetIdx()] = [0,0,1]
		if atom.GetSymbol() == 'S':
			color_dict[atom.GetIdx()] = [1,1,0]
		# print color_dict

	#opts = Draw.DrawingOptions()
	#opts.bgColor = None

	Draw.MolToFile(mol, file_path,  size=size, highlightMap = color_dict, options=opts)


'''
The following function uses the C++ implementation to create the SVG drawing, as described here http://rdkit.blogspot.de/2015/02/new-drawing-code.html
This code is much faster and can create better looking figures then the python implementation
'''

def moltosvg(mol, size=(400,400) , kekulize=True, highlight_atoms = [], highlight_bonds = []): #I just pass all the usual arguments forward to the drawing options
    molSize = size
    mc = Chem.Mol(mol.ToBinary())


    # kekulize should always be inside a try or some other error avoiding sequence, as the rdkit module can be very picky here.
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())

    #not shure about this part, guess there should not be multiple possible confomeres
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)


    #the actual new drawing code
    # mol = Chem.SanitizeMol(mol)
    # mol = Chem.Kekulize(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])

    options = drawer.drawOptions()
    #options.highlightAtoms = (1,1	,0) 
    #options.circleAtoms = False
    #options.continuousHighlight = False
    
    #set your own color values
    a_colors = {}
    b_colors = {}

    color_code = (0.55,1,0.1)

    for a in highlight_atoms:
    	a_colors[a] = color_code
    for b in highlight_bonds:
    	b_colors[b] = color_code


    drawer.DrawMolecule(mc, highlightAtoms = highlight_atoms, highlightBonds = highlight_bonds, highlightAtomColors=a_colors, highlightBondColors=b_colors)# highlightColor = [])#, highlightColor = [])

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return svg.replace('svg:','')

def mol2idx(mol):
	indice_list = []
	for atom in mol.GetAtoms():
		idx = atom.GetIdx()
		indice_list.append(idx)
	return(indice_list)

def mol2idx_bonds(mol):
    indice_list = []
    for bond in mol.GetBonds():
        idx = bond.GetIdx()
        indice_list.append(idx)
    return(indice_list)

def comp_with_origin(template_smiles, comp_smiles, size, save_path):

    # print 'temp', Chem.MolToSmiles(template_smiles)
    # print 'comp', Chem.MolToSmiles(comp_smiles)


    # for some reason the following code does not work with the oriented mol object.
    template_smiles = Chem.MolFromSmiles(Chem.MolToSmiles(template_smiles)) 
    comp_smiles = Chem.MolFromSmiles(Chem.MolToSmiles(comp_smiles))

    '''important: This does not work for ring systems, as the number of bonds can be different'''

    smiles1 = template_smiles #this are the real mol objects, with unspecified bonds
    smiles2 = comp_smiles

    #here the bonds get changed becouse the rdkit bondCompare cannot work with undefined bonds
    #the unspecified bonds, which are usually come from aromatic structures of the original comp are changed to singe, which will then always match.
    nu_smiles1 = Chem.MolFromSmiles(Chem.MolToSmiles(change_bond_type(template_smiles, old_type = 'UNSPECIFIED', new_type = 'SINGLE')))
    nu_smiles2 = Chem.MolFromSmiles(Chem.MolToSmiles(change_bond_type(comp_smiles, old_type = 'UNSPECIFIED', new_type = 'SINGLE')))


    # Draw.MolToImage(smiles1, size = (1000,200)).show()
    # Draw.MolToImage(smiles2, size = (1000,200)).show()

    mols = (nu_smiles1,nu_smiles2) #compare the mol objects, the fact, that the bonds are not correct it not important,
    #as in this step any BondType is matched.

    common_smarts = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareAny).smartsString #molecule comparison
    #print common_smarts
    patt = Chem.MolFromSmarts(common_smarts) #get the commen smarts 

    #Draw.MolToImage(patt, size = (1000,200)).show()
    #exit()
    query1 = nu_smiles1.GetSubstructMatch(patt) #get the match for the smarts in the original mol, this is without considering the bonds 
    query2 = nu_smiles2.GetSubstructMatch(patt) #

    # print query1
    # print query2
    # exit()

    new_mol1 = get_mol_from_query(smiles1,query1) #self made function which gives back a mol object from a query
    # only needed becouse the GetSubstructMatch function from rdkit give back a query instead of a mol object.
    new_mol2 = get_mol_from_query(smiles2,query2)

    Draw.MolToImage(new_mol1, size = (1000,200)).show()
    Draw.MolToImage(new_mol2, size = (1000,200)).show()

    exit()

    diff_bond_idx1 = []
    diff_bond_idx2 = []
    if len(new_mol1.GetBonds()) == len(new_mol2.GetBonds()): #safety first
        '''
        as the molecules must be identical at this stage, except for the bonds the indices should be identical as well (maybe make some tests),
        therefore by just comparing the bond query and extrancting the old indices, the ones to mark should appear


        nooooot true !!
        as the molecules are not identical by idx, the matching pattern must be found by one self. Using the flp_algorithm.
        '''



        #the loop gets the bond indices to highlight, becouse they are different in each mol, this time the real mols must be used.


        # for bond1, bond2 in zip(new_mol1.GetBonds(),new_mol2.GetBonds()):
        #     if bond1.GetBondType() != bond2.GetBondType():
        #         if bond1.GetIdx() == bond2.GetIdx():



        #             #print bond1.GetBondType(), bond2.GetBondType(), bond1.GetBondType(), bond2.GetBondType()

        #             diff_bond_idx1.append(int(bond1.GetProp('old_idx'))) #the indices of the bonds to mark in the first mol object
        #             diff_bond_idx2.append(int(bond2.GetProp('old_idx'))) #the indices of the bonds to mark in the second mol object
        #         else:
        #             print 'different idx'
        #     else:
        #         pass
    else:
        print 'bond querys must have same lenght'

    print diff_bond_idx1
    print diff_bond_idx2

    exit()


    #####
    #this part basically cuts out the areas of the molecule which differ from one another
    #####


    not_highlight = nu_smiles2.GetSubstructMatch(patt) #those are the atoms not to highlight in the molecule

    not_highlight_bond = []
    for bond in new_mol2.GetBonds():
        not_highlight_bond.append(int(bond.GetProp('old_idx'))) 
    #get the bonds not to highlight becouse they are part of the commen structure, the ones which are part of the commen structure but differ will be treated differently

    bnd_idx = mol2idx_bonds(smiles2)
    bond_highlight = [idx for idx in bnd_idx if idx not in not_highlight_bond] #vice verca, if they are not part of the common strucutre, they should be hghlighted

    smi_idx = mol2idx(smiles2) #get the indices for the atoms in the mol object
    atom_highlight = [idx for idx in smi_idx if idx not in not_highlight] #sort out the not to highligh atoms (list comprehensions)

    #####
    #till here
    #####

    bond_highligh = bond_highlight + diff_bond_idx2 #the bonds not in the commen structure and the bonds which differ in type will be highlighted

    smiles2 = change_bond_type(smiles2)

    #bond highlight translation
    #when the bond type is changed the bond indices change as well, this is bad, as the bond_highlight pattern in done for the original mol,
    #therefore the change_bond_type function keeps track of the old idx, as new propertie, similar to the get_mol_from_query function, 
    #this can be used as translation for the highlight map. 

    b_trans_dict = {}

    for bond in smiles2.GetBonds():
        b_trans_dict[int(bond.GetProp('old_idx'))] = bond.GetIdx()  #create the translation dict.

    bond_highligh = [b_trans_dict[idx] for idx in bond_highligh ] #translate the list.

    # print ori
    # print smiles2

    # Draw.MolToImage(ori, size = (1000,200)).show()
    # Draw.MolToImage(smiles2, size = (1000,200)).show()

    #AllChem.GenerateDepictionMatching2DStructure(smiles2 ,ori) #return to the streight path of the mol

    svg1 = (moltosvg(smiles2, highlight_atoms = atom_highlight, highlight_bonds = bond_highligh, size = size)) #pass the lists on to the drawing function


    svg1 = svg1.replace('xmlns:svg=','xmlns=') #the created svg is just fine, but when trying to include into html it does not work,
    # by try and error I found, that this part makes all the difference, amazing ! 
    # print repr(svg1)

    svg_file = open(save_path, 'w')
    svg_file.write(svg1)
    svg_file.close()


#################
#debug
#################


# t_smiles = 'CCC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O'
# c_smiles = 'C~CC=CC(O)C~C(O)C~C(=O)C(C)C=C(C)C(=O)CC=O'

# t_smiles = 'CC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O'
# c_smiles = 'CCCC(O)CC(=O)C(~C(~C)O)~C(~C~C)C(=O)CC(O)C(O)CC(=O)O'

t_smiles = 'CC(O)C(C)C(O)C(C)C=CC=CCCC=CC=CC=CC=CC(O)CC(O)C(C)C(O)CC(=O)CC(O)CCCC(O)CC(O)CC(O)CC=O'
c_smiles = 'COC(C=CC=CC=CC=CCCC=CC=CC(C)C(O)C(C)C(C)O)CC(O)C(C(=O)O)C(O)CC(O)(O)CC(O)C(O)CCC(O)CC(O)CC(O)CC(=O)O'

# temp CC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O
# comp CCC(C)C(C)(OC)C(OC)C(C(O)C(O)C(=O)C(C)(O)CC(C)C(=O)C(C(=O)OC)=C(C)O)C(C)(C)O

# CC(O)C(C)C(O)C(C)C=CC=CCCC=CC=CC=CC=CC(O)CC(O)C(C)C(O)CC(=O)CC(O)CCCC(O)CC(O)CC(O)CC=O
# comp COC(C=CC=CC=CC=CCCC=CC=CC(C)C(O)C(C)C(C)O)CC(O)C(C(=O)O)C(O)CC(O)(O)CC(O)C(O)CCC(O)CC(O)CC(O)CC(=O)O


# temp CC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O
# comp CC(CC(C)(O)C(=O)C(O)C(O)C(C(C)O)C(O)(CO)CC(C)C(C)(C)O)C(=O)C(=CO)C(=O)O

# smiles2 = Chem.MolFromSmiles(c_smiles)

# c_smiles = change_bond_type(smiles2)

#Draw.MolToImage(c_smiles, size = (1000,200)).show()

size = (1000, 200)
save_path1 = './test1.svg'
save_path2 = './test2.svg'



# c_smiles = Chem.MolToSmiles(change_bond_type(Chem.MolFromSmiles(c_smiles)))
#c_smiles = Chem.MolFromSmiles(Chem.MolToSmiles(change_bond_type(Chem.MolFromSmiles(c_smiles))))
#Draw.MolToImage((c_smiles), size = (1000,200)).show()

# svg1 = (moltosvg(c_smiles, size = size)) #pass the lists on to the drawing function


# svg1 = svg1.replace('xmlns:svg=','xmlns=') #the created svg is just fine, but when trying to include into html it does not work,
# # by try and error I found, that this part makes all the difference, amazing ! 
# # print repr(svg1)

# svg_file = open(save_path1, 'w')
# svg_file.write(svg1)
# svg_file.close()

# Draw.MolToImage(Chem.MolFromSmiles(t_smiles), size = (1000,200)).show()
#Draw.MolToImage((c_smiles), size = (1000,200)).show()

#c_smiles = 'CCCC(O)CC(=O)C(~C(~C)O)~C(~C~C)C(=O)CC(O)C(O)CC(=O)O'

# mol = Chem.MolFromSmiles(c_smiles)

# Draw.MolToImage(mol, size = (400,400)).show()

# for bond in mol.GetBonds():
#     print str(bond.GetBondType()), bond.GetIdx()



# for bond in mol.GetBonds():
#     print str(bond.GetBondType()), bond.GetIdx()

#Draw.MolToImage(mol, size = (1000,200)).show()


# size = (1000, 200)
# save_path1 = './test1.svg'
# save_path2 = './test2.svg'

# # mol = Chem.MolFromSmiles(c_smiles)
# # mol1 = change_bond_type(mol)
# #Draw.MolToImage(mol1, size = (1000,200)).show()
# #exit()
comp_with_origin(Chem.MolFromSmiles(t_smiles), Chem.MolFromSmiles(c_smiles), size, save_path1)
comp_with_origin(Chem.MolFromSmiles(c_smiles), Chem.MolFromSmiles(t_smiles), size, save_path2)




'''
temp: CC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O
comp: CC(CC(C)(O)C(=O)C(O)C(O)C(C(C)O)C(O)(CO)CC(C)C(C)(C)O)C(=O)C(=CO)C(=O)O
temp: CC(CC(C)(O)C(=O)C(O)C(O)C(C(C)O)C(O)(CO)CC(C)C(C)(C)O)C(=O)C(=CO)C(=O)O
comp: CC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O
temp: CC=CC(O)CC(O)CC(=O)C(C)C=C(C)C(=O)CC=O
comp: CCCC(O)CC(=O)C(~C(~C)O)~C(~C~C)C(=O)CC(O)C(O)CC(=O)O
[11:47:07] SMARTS Parse Error: syntax error for input: [#6](=[#8])-[#6]-[#6](=,-[#8])-[#6]=,-[#6]-[#6]-[#6]-,[#6]-[#6](-,=[#8])-[#6]-[#6](-[#8])-[#6]=,-[#6]-[#6]
'''