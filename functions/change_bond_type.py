from rdkit.Chem import RWMol
from rdkit import Chem 


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
        bondtype = Chem.BondType.AROMATIC #no idea what to do with the aromaric bonds at the moment 
    if string == 'TRIPLE':
        bondtype = Chem.BondType.TRIPLE

    return(bondtype)

def change_bond_type(mol_object, old_type = 'UNSPECIFIED', new_type = 'AROMATIC'):
    '''
    changes all the bonds from old to new bond type
    '''
    new_type = rdkit_bond_from_str(new_type) #create rdkit type bond
   
    change_mol = RWMol(mol_object) #copy which will be the object to change the bond

    for bond in mol_object.GetBonds():
        #bond.SetProp('old_idx', str(bond.GetIdx()))
        #print str(bond.GetBondType()), bond.GetIdx()

        A = bond.GetBeginAtomIdx()
        B = bond.GetEndAtomIdx()

        if str(bond.GetBondType()) == old_type:

            change_mol.RemoveBond(A,B)
            change_mol.AddBond(A,B,new_type)


        change_mol.GetBondBetweenAtoms(A,B).SetProp('old_idx', str(bond.GetIdx()))


    change_mol = change_mol.GetMol()

    return(change_mol)
