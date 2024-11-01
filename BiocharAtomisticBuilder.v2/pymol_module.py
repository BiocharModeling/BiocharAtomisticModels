import glob, os, sys, subprocess, random
from copy import deepcopy
import numpy as np
import xmlrpc.client as xmlrpclib
from collections import Counter
import time
from rdkit.Chem import AllChem, Draw
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import pymol

def main():

    pdb_files_dir = './pdbfiles'
    files_list = glob.glob(pdb_files_dir+'/*.pdb')
    object_names = []

    for i, pdb_file in enumerate(files_list):
        pymol.cmd.load(pdb_file)
        object_names.append(pdb_file.split('/')[-1].split('.')[0])

    pymol.cmd.orient()
    pymol.cmd.set('sphere_scale', 0.2, '(all)')
    pymol.cmd.set_bond('stick_radius',0.14,'(all)','(all)')
    pymol.cmd.color('red','elem o')
    pymol.cmd.color('grey','elem c')
    pymol.cmd.color('white','elem h')
    pymol.cmd.color('blue','elem n')
    pymol.cmd.color('pink','elem si')
    pymol.cmd.color('yellow','elem s')
    pymol.cmd.show('stick')
    pymol.cmd.show('spheres')


    #Translate molecules as Grid space
    y = 5
    x = int(len(object_names) / y) + 1   
    pos = []

    for x_i in range(x):
        for y_i in range(y):
            x_p = x_i * 30
            y_p = y_i * 30
            pos.append([x_p,y_p,0])

    for i, obj in enumerate(object_names):
        pymol.cmd.translate([pos[i][0],pos[i][1],pos[i][2]],obj)
        pymol.cmd.reset()

    #Find the functional groups
    for i, obj in enumerate(object_names):
        check = examine_h(obj)
        print(check)
        break

def examine_h(object):
    pymol.cmd.zoom('%s' % object)
    pymol.cmd.save('./test.mol','%s' % object,format='mol')
    #refresh_model()
    mol = AllChem.MolFromPNGFile('./test.mol')
    print(mol)
    h_atom_idx    = [atom.GetIdx()+1    for atom in mol.GetAtoms() if atom.GetSymbol() == 'H']
    h_atom_nn_idx = [[n.GetIdx()+1 for n in atom.GetNeighbors()] for atom in mol.GetAtoms() if atom.GetSymbol() == 'H']
    h_atom_nn_sym = [[n.GetSymbol() for n in atom.GetNeighbors()] for atom in mol.GetAtoms() if atom.GetSymbol() == 'H']
    h_atom_nn = [[n for n in atom.GetNeighbors()] for atom in mol.GetAtoms() if atom.GetSymbol() == 'H']

    types = {'Aromatic_C_inRing':[], 'Aliphatic_C_inRing':[],
                'Aromatic_N_inRing':[], 'Aliphatic_N_inRing':[],
                'Aromatic_O_inRing':[], 'Aliphatic_O_inRing':[],
                'Aromatic_S_inRing':[], 'Aliphatic_S_inRing':[],
                'Aromatic_Si_inRing':[], 'Aliphatic_Si_inRing':[],
                'Aromatic_C':[], 'Aliphatic_C':[],
                'Aromatic_N':[], 'Aliphatic_N':[],
                'Aromatic_O':[], 'Aliphatic_O':[],
                'Aromatic_S':[], 'Aliphatic_S':[],
            }
    for i, nns in enumerate(h_atom_nn):
        for n in nns:
            if n.GetIsAromatic() and n.GetSymbol() == 'C' and n.IsInRing():
                types['Aromatic_C_inRing'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'C' and n.IsInRing():
                types['Aliphatic_C_inRing'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'N' and n.IsInRing():
                types['Aromatic_N_inRing'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'N' and n.IsInRing():
                types['Aliphatic_N_inRing'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'Si' and n.IsInRing():
                types['Aromatic_Si_inRing'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'Si' and n.IsInRing():
                types['Aliphatic_Si_inRing'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'O' and n.IsInRing():
                types['Aromatic_O_inRing'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'O' and n.IsInRing():
                types['Aliphatic_O_inRing'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'S' and n.IsInRing():
                types['Aromatic_S_inRing'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'S' and n.IsInRing():
                types['Aliphatic_S_inRing'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'C' and not n.IsInRing():
                types['Aromatic_C'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'C' and not n.IsInRing():
                types['Aliphatic_C'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'N' and not n.IsInRing():
                types['Aromatic_N'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'N' and not n.IsInRing():
                types['Aliphatic_N'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'O' and not n.IsInRing():
                types['Aromatic_O'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'O' and not n.IsInRing():
                types['Aliphatic_O'].append(h_atom_idx[i])
            elif n.GetIsAromatic() and n.GetSymbol() == 'S' and not n.IsInRing():
                types['Aromatic_S'].append(h_atom_idx[i])
            elif not n.GetIsAromatic() and n.GetSymbol() == 'S' and not n.IsInRing():
                types['Aliphatic_S'].append(h_atom_idx[i]) 
    return types

def refresh_model():
    pymol.cmd.orient()
    pymol.cmd.set('sphere_scale', 0.2, '(all)')
    pymol.cmd.set_bond('stick_radius',0.14,'(all)','(all)')
    pymol.cmd.color('red','elem o')
    pymol.cmd.color('grey','elem c')
    pymol.cmd.color('white','elem h')
    pymol.cmd.color('blue','elem n')
    pymol.cmd.color('pink','elem si')
    pymol.cmd.color('yellow','elem s')
    pymol.cmd.show('stick')
    pymol.cmd.show('spheres')


# Examine CNMR function
##### Get the type of carbons ######
def mol_faH(mol):
    num_faH = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    ringinfo = mol.GetRingInfo()
    aring = ringinfo.AtomRings()
    faH_idx = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() =='C' and atom.GetTotalNumHs() == 1 and atom.GetIsAromatic() \
            and atom.IsInRing() and len([r for r in aring if atom.GetIdx() in r and len(r) <= 6]) == 1]
    num_faH = len(faH_idx)
    return num_faH, faH_idx

def mol_faP(mol):
    num_faP = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ringinfo = mol.GetRingInfo()
    aring = [tmp for tmp in ringinfo.AtomRings() if len(tmp) <= 6]

    faP_idx1 = [b.GetBeginAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and (b.GetEndAtom().GetSymbol() == 'O' or b.GetEndAtom().GetSymbol() == 'S') \
        and b.GetBeginAtom().GetIsAromatic() and (len([r for r in aring if b.GetEndAtomIdx() in r ]) == 0) and b.GetEndAtom().GetTotalNumHs() == 1 \
        and b.GetBondType() == Chem.rdchem.BondType.SINGLE] + \
        [b.GetEndAtomIdx() for b in bonds \
        if (b.GetBeginAtom().GetSymbol() == 'O' or b.GetBeginAtom().GetSymbol() == 'S') and b.GetEndAtom().GetSymbol() == 'C' \
        and b.GetEndAtom().GetIsAromatic() and (len([r for r in aring if b.GetBeginAtomIdx() in r ]) == 0) and b.GetBeginAtom().GetTotalNumHs() == 1 \
        and b.GetBondType() == Chem.rdchem.BondType.SINGLE]

    faP_idx2 = [b.GetBeginAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and (b.GetEndAtom().GetSymbol() == 'O' or b.GetEndAtom().GetSymbol() == 'S') \
        and b.GetBeginAtom().GetIsAromatic() and (len([r for r in aring if b.GetEndAtomIdx() in r ]) == 0) and b.GetEndAtom().GetTotalNumHs() == 0 \
        and b.GetBondType() == Chem.rdchem.BondType.SINGLE] + \
        [b.GetEndAtomIdx() for b in bonds \
        if (b.GetBeginAtom().GetSymbol() == 'O' or b.GetBeginAtom().GetSymbol() == 'S') and b.GetEndAtom().GetSymbol() == 'C' \
        and b.GetEndAtom().GetIsAromatic() and (len([r for r in aring if b.GetBeginAtomIdx() in r ]) == 0) and b.GetBeginAtom().GetTotalNumHs() == 0 \
        and b.GetBondType() == Chem.rdchem.BondType.SINGLE]
        
    faP_idx3 = [b.GetBeginAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and (b.GetEndAtom().GetSymbol() == 'O' or b.GetEndAtom().GetSymbol() == 'S') \
        and b.GetBeginAtom().GetIsAromatic() and not b.GetEndAtom().GetIsAromatic() and b.GetEndAtom().GetTotalNumHs() == 0 \
        and b.GetBondType() == Chem.rdchem.BondType.SINGLE] + \
        [b.GetEndAtomIdx() for b in bonds \
        if (b.GetBeginAtom().GetSymbol() == 'O' or b.GetBeginAtom().GetSymbol() == 'S') and b.GetEndAtom().GetSymbol() == 'C' \
        and not b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic() and b.GetBeginAtom().GetTotalNumHs() == 0 \
        and b.GetBondType() == Chem.rdchem.BondType.SINGLE]

    faP_idx4 = [b.GetBeginAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and (b.GetEndAtom().GetSymbol() == 'O' or b.GetEndAtom().GetSymbol() == 'S') \
        and b.GetBeginAtom().GetIsAromatic() and not b.GetEndAtom().GetIsAromatic() and b.GetEndAtom().GetTotalNumHs() == 0 \
        and b.GetBondType() == Chem.rdchem.BondType.DOUBLE] + \
        [b.GetEndAtomIdx() for b in bonds \
        if (b.GetBeginAtom().GetSymbol() == 'O' or b.GetBeginAtom().GetSymbol() == 'S') and b.GetEndAtom().GetSymbol() == 'C' \
        and not b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic() and b.GetBeginAtom().GetTotalNumHs() == 0 \
        and b.GetBondType() == Chem.rdchem.BondType.DOUBLE]

    faP_idx = list(set(faP_idx1+faP_idx2+faP_idx3))
    num_faP = len(faP_idx)
    return num_faP, faP_idx

def mol_faS(mol):
    num_faS = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ringinfo = mol.GetRingInfo()

    aring = [tmp for tmp in ringinfo.AtomRings() if len(tmp) <= 6]

    faS_idx1 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() =='C' and atom.GetTotalNumHs() == 0 and atom.GetIsAromatic() \
            and len([n for n in atom.GetNeighbors() if (n.GetSymbol() == 'S') and (len([r for r in aring if n.GetIdx() in r ]) == 0) and not n.GetIsAromatic() ])==0 \
            and len([n for n in atom.GetNeighbors() if (n.GetSymbol() == 'O') and (len([r for r in aring if n.GetIdx() in r ]) == 0) and not n.GetIsAromatic() ])==0 \
            and len([r for r in aring if atom.GetIdx() in r ]) == 1 ]

    faS_idx2 = [b.GetBeginAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C' \
            and b.GetBeginAtom().GetIsAromatic() and not b.GetEndAtom().GetIsAromatic()] + \
        [b.GetEndAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'C' \
            and not b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic()]

    faS_idx3 = [b.GetBeginAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'C' and b.GetEndAtom().GetSymbol() == 'N' \
        and b.GetBeginAtom().GetIsAromatic() and not b.GetEndAtom().GetIsAromatic() ] + \
        [b.GetEndAtomIdx() for b in bonds \
        if b.GetBeginAtom().GetSymbol() == 'N' and b.GetEndAtom().GetSymbol() == 'C' \
        and not b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic() ]

    faS_idx4 = [b.GetBeginAtomIdx() for b in bonds if b.GetBeginAtom().GetSymbol() == 'C' \
        and (b.GetEndAtom().GetSymbol() != 'O') and (b.GetEndAtom().GetSymbol() != 'S')\
        and b.GetBeginAtom().GetIsAromatic() and not b.GetEndAtom().GetIsAromatic() and b.GetEndAtom().GetTotalNumHs() == 1] + \
        [b.GetEndAtomIdx() for b in bonds if b.GetEndAtom().GetSymbol() == 'C' \
        and (b.GetBeginAtom().GetSymbol() != 'O') and (b.GetBeginAtom().GetSymbol() != 'S') \
        and not b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic() and b.GetBeginAtom().GetTotalNumHs() == 1]
    

    faS_idx = list(set(faS_idx1+faS_idx2+faS_idx3+faS_idx4))
    num_faS = len(faS_idx)
    return num_faS, faS_idx

def mol_faB(mol):
    num_faB = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    ringinfo = mol.GetRingInfo()
    aring = [tmp for tmp in ringinfo.AtomRings() if len(tmp) <= 6]

    faB_idx = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() =='C' and atom.GetIsAromatic() and atom.GetTotalNumHs() == 0 \
        and len([n.GetIdx() for n in atom.GetNeighbors() if n.GetIsAromatic()]) == 3 \
        #and sum([1 if b.GetIsAromatic() else 0 for b in atom.GetBonds()]) >= 2 \
        and len([r for r in aring if atom.GetIdx() in r]) >= 2 ]

    num_faB = len(faB_idx)
    return num_faB, faB_idx

def mol_faC(mol):
    num_faC = 0
    mol = AllChem.RemoveHs(mol)
    bonds = mol.GetBonds()
    atoms = mol.GetAtoms()

    #Aldehyde
    faC_idx1 = [bond.GetBeginAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE \
        and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'O' \
        and not bond.GetBeginAtom().IsInRing() and not bond.GetEndAtom().IsInRing() \
        and bond.GetBeginAtom().GetTotalNumHs() == 1 \
        and len([n.GetSymbol() for n in bond.GetBeginAtom().GetNeighbors() if n.GetSymbol() == 'C' or n.GetSymbol() == 'N' or n.GetSymbol() == 'S']) == 1]
    faC_idx2 = [bond.GetEndAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE \
        and bond.GetEndAtom().GetSymbol() == 'C' and bond.GetBeginAtom().GetSymbol() == 'O' \
        and not bond.GetEndAtom().IsInRing() and not bond.GetBeginAtom().IsInRing() \
        and bond.GetEndAtom().GetTotalNumHs() == 1 \
        and len([n.GetSymbol() for n in bond.GetEndAtom().GetNeighbors() if n.GetSymbol() == 'C' or n.GetSymbol() == 'N' or n.GetSymbol() == 'S']) == 1]

    #Ketone
    faC_idx3 = [bond.GetBeginAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE \
        and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'O' \
        and not bond.GetEndAtom().IsInRing() \
        and bond.GetBeginAtom().GetTotalNumHs() == 0 \
        and len([n.GetSymbol() for n in bond.GetBeginAtom().GetNeighbors() if n.GetSymbol() == 'C' or n.GetSymbol() == 'N' or n.GetSymbol() == 'S']) == 2]
    faC_idx4 = [bond.GetEndAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE \
        and bond.GetEndAtom().GetSymbol() == 'C' and bond.GetBeginAtom().GetSymbol() == 'O' \
        and not bond.GetBeginAtom().IsInRing() \
        and bond.GetEndAtom().GetTotalNumHs() == 0 \
        and len([n.GetSymbol() for n in bond.GetEndAtom().GetNeighbors() if n.GetSymbol() == 'C' or n.GetSymbol() == 'N' or n.GetSymbol() == 'S']) == 2]
    
    #Carboxylic Acid
    faC_idx5 = [bond.GetBeginAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE \
        and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'O' \
        and not bond.GetBeginAtom().GetIsAromatic() and not bond.GetEndAtom().IsInRing() \
        and len([n.GetSymbol() for n in bond.GetBeginAtom().GetNeighbors() if n.GetSymbol() == 'O']) == 2]
    faC_idx6 = [bond.GetEndAtom() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE \
        and bond.GetBeginAtom().GetSymbol() == 'O' and bond.GetEndAtom().GetSymbol() == 'C' \
        and not bond.GetEndAtom().GetIsAromatic() and not bond.GetBeginAtom().IsInRing() \
        and len([n.GetSymbol() for n in bond.GetEndAtom().GetNeighbors() if n.GetSymbol() == 'O']) == 2]

    faC_idx = list(set(faC_idx1+faC_idx2+faC_idx3+faC_idx4+faC_idx5+faC_idx6))
    num_faC = len(faC_idx)
    return num_faC, faC_idx

def mol_fal(mol):
    num_fal = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ringinfo = mol.GetRingInfo()
    aring = ringinfo.AtomRings()

    fal_idx1 = [bond.GetBeginAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE \
        and bond.GetBeginAtom().GetSymbol() == 'C' and not bond.GetBeginAtom().GetIsAromatic() and not bond.GetBeginAtom().IsInRing() \
        and bond.GetBeginAtom().GetTotalNumHs() == 1 \
        and bond.GetEndAtom().GetSymbol() == 'O' and not bond.GetEndAtom().GetIsAromatic() and not bond.GetEndAtom().IsInRing() \
        and bond.GetEndAtom().GetTotalNumHs() > 0 \
        and len([n.GetSymbol() for n in bond.GetBeginAtom().GetNeighbors() if n.GetSymbol() == 'C']) == 1] \
        + [bond.GetEndAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE \
        and bond.GetEndAtom().GetSymbol() == 'C' and not bond.GetEndAtom().GetIsAromatic() and not bond.GetEndAtom().IsInRing() \
        and bond.GetEndAtom().GetTotalNumHs() == 1 \
        and bond.GetBeginAtom().GetSymbol() == 'O' and not bond.GetBeginAtom().GetIsAromatic() and not bond.GetBeginAtom().IsInRing() \
        and bond.GetBeginAtom().GetTotalNumHs() > 0 \
        and len([n.GetSymbol() for n in bond.GetEndAtom().GetNeighbors() if n.GetSymbol() == 'C']) == 1]

    fal_idx2 = [bond.GetBeginAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE \
        and bond.GetBeginAtom().GetSymbol() == 'C' and not bond.GetBeginAtom().GetIsAromatic() and not bond.GetBeginAtom().IsInRing()\
        and bond.GetBeginAtom().GetTotalNumHs() >= 2 ]\
        + [bond.GetEndAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE \
        and bond.GetEndAtom().GetSymbol() == 'C' and not bond.GetEndAtom().GetIsAromatic() and not bond.GetEndAtom().IsInRing()\
        and bond.GetEndAtom().GetTotalNumHs() >= 2 ]

    fal_idx3 = [bond.GetBeginAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE \
        and bond.GetBeginAtom().GetSymbol() == 'C' and not bond.GetBeginAtom().GetIsAromatic() and bond.GetBeginAtom().IsInRing()\
        and bond.GetBeginAtom().GetTotalNumHs() >= 1 ]\
        + [bond.GetEndAtomIdx() for bond in bonds \
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE \
        and bond.GetEndAtom().GetSymbol() == 'C' and not bond.GetEndAtom().GetIsAromatic() and bond.GetEndAtom().IsInRing()\
        and bond.GetEndAtom().GetTotalNumHs() >= 1 ]

    fal_idx4 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.IsInRing() and atom.GetTotalNumHs() == 1 \
        and 'N' in [n.GetSymbol() for n in atom.GetNeighbors()] and len([n.GetSymbol() for n in atom.GetNeighbors()]) == 3] \
        + [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.IsInRing() and atom.GetTotalNumHs() <= 2 \
        and 'O' in [n.GetSymbol() for n in atom.GetNeighbors() if n.IsInRing()] ]

    fal_idx5 =[atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() =='C' and atom.IsInRing() and not atom.GetIsAromatic() and atom.GetTotalNumHs() == 0 \
        and sum([1 if b.GetIsAromatic() else 0 for b in atom.GetBonds()]) < 3 \
        and len([r for r in aring if atom.GetIdx() in r]) > 1 ]

    fal_idx6 = [atom.GetIdx() for atom in atoms \
        if not atom.GetIsAromatic() and atom.GetSymbol() == 'C' and atom.GetTotalNumHs()==0 \
        and len([n for n in atom.GetNeighbors() if (n.GetSymbol() =='C' or n.GetSymbol() =='S') and not n.IsInRing()]) == 1 \
        and len([n for n in atom.GetNeighbors() if (n.GetSymbol() =='O') and (not n.IsInRing()) and (mol.GetBondBetweenAtoms(atom.GetIdx(),n.GetIdx()).GetBondType()==Chem.rdchem.BondType.DOUBLE)]) == 0]

    fal_idx7 = [atom.GetIdx() for atom in atoms \
        if not atom.GetIsAromatic() and atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 \
        and len([n for n in atom.GetNeighbors() if n.GetSymbol() =='O' and not n.IsInRing() and len(n.GetNeighbors())<=2]) == 1 \
        and len([b.GetBeginAtom() for b in atom.GetBonds() if (b.GetBeginAtom().GetSymbol()=='O' or b.GetEndAtom().GetSymbol()=='O') and b.GetBondTypeAsDouble() == 1])==1]
    
    fal_idx8 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 1 \
        and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol()=='C']) == 2]

    fal_idx9 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.IsInRing() and atom.GetTotalNumHs() == 0 \
        and len([n.GetSymbol() for n in atom.GetNeighbors()]) == 4]
    
    fal_idx10 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() <= 1 \
            and len([n for n in atom.GetNeighbors()]) == 3]

    fal_idx11 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 2 \
            and len([n for n in atom.GetNeighbors()]) == 1]

    fal_idx = fal_idx1 + fal_idx2 + fal_idx3 + fal_idx4 + fal_idx5 + fal_idx6 + fal_idx7 + fal_idx8 + fal_idx9 + fal_idx10 + fal_idx11
    fal_idx = list(set(fal_idx))
    num_faC, faC_idx = mol_faC(mol)
    fal_idx = [idx for idx in fal_idx if idx not in faC_idx]

    num_fal = len(fal_idx)
    return num_fal, fal_idx


### Need to be checked
def mol_falas(mol): #CH3 or nonprotonated
    num_fal = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ringinfo = mol.GetRingInfo()
    aring = ringinfo.AtomRings()

    fal_idx1 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 3] 
    fal_idx2 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 0 \
            and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) < 2]
    fal_idx3 =[atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() =='C' and atom.IsInRing() and not atom.GetIsAromatic() and atom.GetTotalNumHs() == 0 \
            and sum([1 if b.GetIsAromatic() else 0 for b in atom.GetBonds()]) < 3 \
            and len([r for r in aring if atom.GetIdx() in r]) > 1 ]
    fal_idx4 = [atom.GetIdx() for atom in atoms \
        if not atom.GetIsAromatic() and atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 \
            and len([n for n in atom.GetNeighbors() if n.GetSymbol() =='O' and not n.IsInRing() and len(n.GetNeighbors())<=2]) == 1 \
            and len([b.GetBeginAtom() for b in atom.GetBonds() if (b.GetBeginAtom().GetSymbol()=='O' or b.GetEndAtom().GetSymbol()=='O') and b.GetBondTypeAsDouble() == 1])==1]
    fal_idx5 = [atom.GetIdx() for atom in atoms \
        if not atom.GetIsAromatic() and atom.GetSymbol() == 'C' and atom.GetTotalNumHs()==0 \
            and len([n for n in atom.GetNeighbors() if (n.GetSymbol() =='C' or n.GetSymbol() =='S') and not n.IsInRing()]) == 1]
    fal_idx6 = [atom.GetIdx() for atom in atoms \
        if not atom.GetIsAromatic() and atom.GetSymbol() == 'C' and atom.GetTotalNumHs()==0 and atom.IsInRing() \
            and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol() == 'O']) >= 1 ]
    fal_idx7 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.IsInRing() and atom.GetTotalNumHs() == 0 \
            and len([n for n in atom.GetNeighbors()]) == 4]
    fal_idx8 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 0 \
            and len([n for n in atom.GetNeighbors()]) == 3]

    fal_idx = fal_idx1 + fal_idx2 + fal_idx3 + fal_idx4 + fal_idx5 + fal_idx6 + fal_idx7 + fal_idx8
    fal_idx = list(set(fal_idx))
    num_faC, faC_idx = mol_faC(mol)
    fal_idx = [idx for idx in fal_idx if idx not in faC_idx]

    num_fal = len(fal_idx)
    return num_fal, fal_idx

def mol_falH(mol): #CH or CH2
    num_fal = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ringinfo = mol.GetRingInfo()
    aring = [tmp for tmp in ringinfo.AtomRings() if len(tmp) <= 6]
    
    fal_idx0 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 1 \
            and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol()=='O']) == 0 \
            and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol()=='C']) == 1]
    fal_idx1 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 2]
    fal_idx2 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 1 \
            and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol()=='O']) == 0 \
            and len([n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol()=='C']) == 2] 
    fal_idx3 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.IsInRing() and atom.GetTotalNumHs() == 2 \
            and len([r for r in aring if atom.GetIdx() in r]) >= 1 ]
    fal_idx4 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.IsInRing() and atom.GetTotalNumHs() == 1 \
            and len([r for r in aring if atom.GetIdx() in r]) >= 1 ]
    fal_idx5 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and atom.GetTotalNumHs() == 2 \
            and len([r for r in aring if atom.GetIdx() in r]) < 1 ]
    fal_idx6 = [atom.GetIdx() for atom in atoms \
        if atom.GetSymbol() == 'C' and not atom.GetIsAromatic() and not atom.IsInRing() and atom.GetTotalNumHs() == 2 \
            and len([n for n in atom.GetNeighbors()]) == 1]

    fal_idx = fal_idx0 + fal_idx1 + fal_idx2 + fal_idx3 + fal_idx4 + fal_idx5 + fal_idx6
    fal_idx = list(set(fal_idx))
    num_fal = len(fal_idx)
    return num_fal, fal_idx

def mol_falO(mol):
    num_fal = 0
    mol = AllChem.RemoveHs(mol)
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    ringinfo = mol.GetRingInfo()
    aring = ringinfo.AtomRings()

    num_fal, fal_idx = mol_fal(mol)
    num_faC, faC_idx = mol_faC(mol)
    falO_idx1 = [idx for idx in fal_idx \
                if len([n.GetSymbol() for n in atoms[idx].GetNeighbors() if n.GetSymbol()=='O']) >= 1 \
                and idx not in faC_idx ]

    fal_idx = falO_idx1
    fal_idx = list(set(fal_idx))
    num_fal = len(fal_idx)
    return num_fal, fal_idx

####################################

##### Calculate 13C-NMR parameters ######

def carbon_nmr_individual(mol):

    atoms = mol.GetAtoms()
    num_c = len([atom for atom in atoms if atom.GetSymbol()=='C'])
    num_faH, faH_idx = mol_faH(mol)
    num_faP, faP_idx = mol_faP(mol)
    num_faS, faS_idx = mol_faS(mol)
    num_faB, faB_idx = mol_faB(mol)
    num_faC, faC_idx = mol_faC(mol)
    num_fal, fal_idx = mol_fal(mol)
    faH = num_faH
    faP = num_faP
    faS = num_faS
    faB = num_faB
    faC = num_faC
    faN = faP + faS + faB
    fapr = faN + faH

    num_fal, fal_idx = mol_fal(mol)
    num_falas, falas_idx = mol_falas(mol)
    num_falH, falH_idx = mol_falH(mol)
    num_falO, falO_idx = mol_falO(mol)
    fal = num_fal
    falas = num_falas
    falH  = num_falH
    falO  = num_falO

    fal = num_fal
    total_c = num_c
    test = sum([num_faH, num_faP, num_faS, num_faB, num_faC, num_fal])
    result = {'fapr':fapr,
              'faN':faN,
              'faH':faH,
              'faP':faP,
              'faS':faS,
              'faB':faB,
              'faC':faC,
              'fal':fal,
              'fal*':falas,
              'falH':falH,
              'falO':falO,
              'num_C':total_c,
              }

    result_ratio = {'fapr':fapr/total_c,
                    'faN':faN/total_c,
                    'faH':faH/total_c,
                    'faP':faP/total_c,
                    'faS':faS/total_c,
                    'faB':faB/total_c,
                    'faC':faC/total_c,
                    'fal':fal/total_c,
                    'fal*':falas/total_c,
                    'falH':falH/total_c,
                    'falO':falO/total_c,
                    'num_C':total_c,
                    }

    index = {'faH':faH_idx,
            'faP':faP_idx,
            'faS':faS_idx,
            'faB':faB_idx,
            'faC':faC_idx,
            'fal':fal_idx,
            'fal*':falas_idx,
            'falH':falH_idx,
            'falO':falO_idx,
            }
    return result_ratio, index

main()
