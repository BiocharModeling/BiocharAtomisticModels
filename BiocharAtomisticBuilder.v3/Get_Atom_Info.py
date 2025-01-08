import pymol
from pymol import stored
#from chempy import atomic_mass
import sys

#https://pymol.org/pymol-command-ref.html : pymol command reference
#https://pymol.org/dokuwiki/doku.php?id=api : pymol api documments

stored.atom_index = []
stored.atom_element = []
stored.atom_name = []
pymol.cmd.iterate('mol1', 'stored.atom_index.append(index)')
pymol.cmd.iterate('mol1', 'stored.atom_element.append(elem)')
pymol.cmd.iterate('mol1', 'stored.atom_name.append(name)')

f = open('atom_info.txt','w')
for i in range(len(stored.atom_index)):
    atom_idx = stored.atom_index[i]
    atom_elem = stored.atom_element[i]
    atom_name = stored.atom_name[i]
    f.write('%d %s %s %s\n' % (i,atom_idx,atom_elem,atom_name))

f.close()


    
