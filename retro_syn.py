########################
# written by Paul Zierep
########################

import os 

#from flp_04 import pks_mol
from graphviz import Digraph
#check if rdkit is installed with all its env paths

for env in ['LD_LIBRARY_PATH','PYTHONPATH','RDBASE']:
	if env not in os.environ:
		print(env + ' must be set with rdkit folder')

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

class reaction():
	'''this class adds some functionality to the rdkit reaction class, like a name'''
	def __init__(self, rxn):
		self.name = rxn[0]
		self.reaction_string = rxn[1]
		self.rdkit_rxn = AllChem.ReactionFromSmarts(rxn[1])
		educt_product = rxn[1].split('>>')
		products = educt_product[1].split('.')
		self.product_num = len(products)

class reaction_node():
	'''this node gets created for each new product, which becomes the new main mol, very handy for 
	iterative functions, stores its educt as mother_mol, can be used to create a tree base on the
	mother_mol indices
	'''

	def __init__(self, mol, mother_mol, node_index, rxn, rest):
		if mother_mol == 'init':
			self.mother_mol = Chem.MolFromSmiles('C')
			self.mother_mol_idx = -1
		elif mother_mol != None:
			self.mother_mol = mother_mol.mol
			self.mother_mol_idx = mother_mol.node_index
		else:
			self.mother_mol = None
			self.mother_mol_idx = None

		self.mol = mol
		self.rxn = rxn
		self.rest = rest
		self.node_index = node_index

#help(rdkit)

# print reaction_class_list

def trim_function(main_mol):
	main_mol = Chem.RemoveHs(main_mol) #new products often have strange proton numbers, even if they are otherwise fine, is hereby avoided
	return(main_mol)

def recursive_reaction_tree(nodes, reaction_class_list, limit = 5):
	'''taks a node object from the nodes list, creates multiple new node objects, based on the simulated reaction and possible products, 
	appnds the newly created nodes to the nodes list, the list is walked trough and each node is processed again, until new nodes 
	can not be created, which is the case if all the reactions fail to produce products, this ends the iteration'''

	test_index = 0

	node_index = 0
	job_index = 0
	for node in nodes:

		#print [len(Chem.MolToSmiles(node.mol)) for node in nodes if node.mol != None]
		test_index += 1
		if test_index >= limit:
			break

		if node.node_index == job_index: #each object in the list gets processed as a job

			for reaction in reaction_class_list:

				# if node.mol == None: #if there is no more mol in the node the iteration stops
				# 	pass 

				# else:
				ps = reaction.rdkit_rxn.RunReactants([node.mol]) #try reaction
				for products in ps:

					# print products[0]

					node_index += 1 #each node gets an individual index, which can be used for backtracking, as a node keeps also track of its mother nodes index

					# print Chem.MolToSmiles(products[0])
					# print Chem.MolToSmiles(node.mother_mol)

					if len(products) == 0: #no more products --> no more nodes
						print 'no product'
						# main_mol = None
						# site_mol = None
						#print node.mother_mol

					elif Chem.MolToSmiles(products[0]) == Chem.MolToSmiles(node.mol): #some reactions can lead to the creation of an identical product, this would
					#lead to an infinit loop, is hereby avoided, this reactions should still be cureated.
						print 'educt == product'
						# main_mol = None
						# site_mol = None

					else:
						main_mol = trim_function(products[0]) #specific function which can be applyed on the product, eg remove water, which would lead to a doublication, as this
						#is not recognised by rdkid as a dublicate -- could also be due to an not perfect reaction string !

						site_mol = products[1:]
						#print reaction.rdkit_rxn.GetReactingAtoms()

						new_node = reaction_node(main_mol, node, node_index, reaction.name, site_mol)

						nodes.append(new_node)

		job_index += 1

	return(nodes)

#test molecule
avermectin_smi = 'CC[C@H](C)[C@@H]1[C@H](C=C[C@@]2(O1)C[C@@H]3C[C@H](O2)C/C=C(/[C@H]([C@H](/C=C/C=C/4\CO[C@H]5[C@@]4([C@@H](C=C([C@H]5O)C)C(=O)O3)O)C)O[C@H]6C[C@@H]([C@H]([C@@H](O6)C)O[C@H]7C[C@@H]([C@H]([C@@H](O7)C)O)OC)OC)\C)C'
a_mol = Chem.MolFromSmiles(avermectin_smi)

# reaction set
reactions = [

['esterfication intramolecular', "([C:1]-O-C(=O)-[C:2])>>([C:1]-O.O-C(=O)-[C:2])"],
['glycosylation', "C1([O:1]-[C:2])OC(C)C(O)C(OC)C1>>[C:2]-O.C1([O:1])OC(C)C(O)C(OC)C1"],
#['AveE_special_modi', "([C:1]1[C:2][C:3]([O:6])[C:4]O1)>>([C:1][C:2][C:3]([O:6])[C:4])"],
#['Ring_Ave_special_modi', "([C:1][C:2](~[C:5]~[C:6]~[C:7]~[C:8]1)[C:3]1[C:4])>>([C:1][C:2](~[C:5]~[C:6]~[C:7]~[C:8]1).[C:3]1[C:4])"], #todo


]

reaction_class_list = [reaction(rxn) for rxn in reactions]


#Draw.MolToImage(a_mol, size = (500,500)).show()

init_reaction_node = [reaction_node(a_mol, 'init', 0, None, ())]

nodes = recursive_reaction_tree(init_reaction_node, reaction_class_list, limit = 100)

#build a tree based on the nodes


#todo, create paths for individual reaction
# take care, that multiple reactions at ones are avoided
# combine with path algorithm

dot = Digraph()

path = './output'

for node in nodes:
	#img_s = Draw.MolToImage(node.mol, size = (1000,1000))
	img_path = os.path.join(path, str(node.node_index)+'.png')
	Draw.MolToFile(node.mol, img_path, size = (500,500), type='png')
	for rest in node.rest:
		rest_path = os.path.join(path, str(node.node_index)+'_rest'+'.png')
		Draw.MolToFile(rest, rest_path, size = (500,500), type='png')
		dot.node(str(node.node_index) + '_rest', image= rest_path, label = '')
		if node.mother_mol_idx != -1:
			dot.edge(str(node.mother_mol_idx), str(node.node_index) + '_rest')
	#print img_path

	dot.node(str(node.node_index), image= img_path, label = '')

	if node.mother_mol_idx != -1:
		dot.edge(str(node.mother_mol_idx), str(node.node_index))

print dot.source

dot.render('test-output/round-table.gv', view=True)
# print nodes

# for node in nodes:
# 	if node.mother_mol_idx == last_node
# 	last_node = node

#print [(node.node_index, node.mother_mol_idx, ) for node in nodes]




#################
# debug
#################

# path = './output'

# for node in nodes:

# 	#img_m = Draw.MolToImage(node.mother_mol, size = (500,500))
# 	img_s = Draw.MolToImage(node.mol, size = (500,500))
# 	img_s.save(os.path.join(path, str(node.node_index)+str(node.mother_mol_idx) + '.jpg'))
# 	for rest in node.rest:
# 		#print rest
# 		img_s = Draw.MolToImage(rest , size = (500,500))
# 		img_s.save(os.path.join(path, str(node.node_index)+str(node.mother_mol_idx)+'_rest' + '.jpg'))




# print nodes




# exit()
# for liste in mol_list:
# 	print liste

# exit()



# rxn3 = ["[C:2]-O-[C:1]>>[C:2]-O.O-[C:1]", 'glycosylation']

# rdkit_rxn = AllChem.ReactionFromSmarts(rxn3[0])

# def reaction_classifier(reaction):
# 	reaction_formula = reaction[0]
# 	rdkit_rxn = AllChem.ReactionFromSmarts(reaction_formula)

# 	reaction_name = reaction[1]
# 	reaction_string = eaction_formula.split('>>')
# 	educt = reaction_string[0]
# 	products = reaction_string[0].split('.')
# 	reaction_dict = {reaction_name: }






# 'C1=CCOCCCCCCC1'
#rxn = AllChem.ReactionFromSmarts("[*:1]=[*:2]>>([C:1]=[C;H2].[C:2]=[C;H2])")
#rxn = AllChem.ReactionFromSmarts('(C(=O)O.OCC)>>C(=O)OCC.O')

# ps = rxn.RunReactants([a_mol])


# #Chem.MolToSmiles(ps[0][0])
# print ps

# left_over = {}

# f_product = [prod[0] for prod in ps]
# side_product = [prod[1] for prod in ps]

# f_product_smi = [(Chem.MolToSmiles(mol, True)) for mol in f_product]

# print f_product_smi


# ps = rdkit_rxn.RunReactants([a_mol])
# for conf in ps:
# 	for mol in conf:
# 		print 'g2'
# 		# print mol
# 		print Chem.MolToSmiles(mol, True)
# 		mol = Chem.RemoveHs(mol)
# 	# for atom in a_mol[0].GetAtoms():
# 	# 	#print atom.GetIdx()
# 		img = Draw.MolToImage(mol, size = (500,500))
# 		img.show()
# 	# 	break