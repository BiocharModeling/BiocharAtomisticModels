import pymol
from pymol import stored
#from chempy import atomic_mass
import sys


#https://pymol.org/pymol-command-ref.html : pymol command reference
#https://pymol.org/dokuwiki/doku.php?id=api : pymol api documments


stored.atom_names = []
pymol.cmd.iterate('mol1', 'stored.atom_names.append(name)')

f = open('atom_names.txt','w')

for i, atom_name in enumerate(stored.atom_names):
    f.write('%d %s\n' % (i,atom_name))

f.close()


    
