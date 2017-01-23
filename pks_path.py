from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np
from functions.streight_path import streight_path
from functions.draw import draw_mol_and_path, draw_mol_and_path_save, comp_with_origin, draw_mol

class pks_path():

	def __init__(self, name, master_mol, matrix, path_start, path, side_path):
		#collect al the information of the path
		self.name = name
		self.master_mol = master_mol # original molecule the path is extracted from
		self.matrix = matrix # matrix writing
		self.path_start = path_start 
		self.main_path = path
		self.side_path = side_path
		self.mol, self.mol_main_path = streight_path(matrix)


	def draw_m_and_p(self, file_path, size):
		#draws a figure of the main molecule with the path. 
		draw_mol_and_path(self.master_mol, self.main_path, file_path, size)
		# return img

	def draw_path(self, file_path, size = (400,400)):
		#saves the main path to the given file_path
		draw_mol(self.mol, path, size)

	def draw_comp_with_origin(self, template_smiles, size, path, ori):
		#not finished jet.
		comp_with_origin((self.mol), template_smiles, size, path, ori)

		#comp_with_origin(Chem.MolToSmiles(template_smiles), Chem.MolToSmiles(self.mol), size, path)

	def draw_path_debug(self):
		#only for debugging
		Draw.MolToImage(self.mol, size = (400,400)).show()



	# def draw_path_return(self):
	# 	size = (2000, 500)
	# 	#type1 = 'svg'
	# 	img = Draw.MolToImage(self.mol, size=size)
	# 	return img

	# # def draw_path(self, file_path):
	# # 	size = (2000, 500)
	# # 	type1 = 'svg'
	# # 	draw_mol_and_path_save(self.mol, file_path, size, type1)

	# # def draw_path_save(self, file_path):
	# # 	size = (2500,1500)
	# # 	type1 = 'svg'
	# # 	draw_mol_and_path_save(self.mol, file_path, size, type1)

	# def draw_main_path(self):
	# 	Draw.MolToImage(self.mol_main_path).show()

	# def draw_comp_path(self): #....
	# 	pass






