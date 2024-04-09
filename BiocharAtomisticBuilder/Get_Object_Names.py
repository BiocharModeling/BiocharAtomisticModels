import pymol
from pymol import stored
#from chempy import atomic_mass
import sys


#https://pymol.org/pymol-command-ref.html : pymol command reference
#https://pymol.org/dokuwiki/doku.php?id=api : pymol api documments

molecule_names = pymol.cmd.get_object_list()
molecule_names2 = []
for i, name in enumerate(molecule_names):
    if name[:2] != 'p_':
        molecule_names2.append(name)

molecule_names2 = sorted( molecule_names2, key=lambda x: pymol.util.compute_mass(x) )

f = open('object_name_mw_dbe.txt','w')

for i, mol_name in enumerate(molecule_names2):
    f.write('%d %s %f\n' % (i,mol_name,pymol.util.compute_mass(mol_name)))

f.close()


    
