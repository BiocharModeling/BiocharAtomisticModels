import pymol
from pymol import stored
#from chempy import atomic_mass
import sys

#https://pymol.org/pymol-command-ref.html : pymol command reference
#https://pymol.org/dokuwiki/doku.php?id=api : pymol api documments

stored.atom_index = []
pymol.cmd.iterate('mol1', 'stored.atom_index.append(index)')

f = open('atom_index.txt','w')
for i, atom_index in enumerate(stored.atom_index):
    f.write('%d %s\n' % (i,atom_index))

f.close()


    
