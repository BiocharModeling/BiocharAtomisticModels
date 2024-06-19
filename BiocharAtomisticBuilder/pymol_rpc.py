#Initiate Pymol RPC module
import glob,os,sys,subprocess
import xmlrpc.client as xmlrpclib
from collections import Counter
import time
from rdkit.Chem import AllChem
from rdkit import Chem
from Examine_CNMR import *

#!which pip
#!pip install rdkit
#!pip install ipywidgetse
#!pip install ase
#!pip install ipymol
#!pip install networkx

class pymol_jupyter_builder:

    HOST = os.environ.get('PYMOL_RPCHOST', 'localhost')
    PORT = 9123

    def __init__(self, host=HOST, port=PORT):
        self.host = host
        self.port = int(port)
        self._process = None
        self.Grid_mode = False

    def __del__(self):
        self.stop()
    
    #def __getattr__(self,key):
    #    if not self._process_is_running():
    #        self.start()
    #    return getattr(self.server, key)
    
    def _process_is_running(self):
        return self._process is not None and self._process.poll() is None

    def start(self, args=("",), exe="pymol"):
        if self._process_is_running():
            print("A PyMOL RPC server is already running.")
            return

        assert isinstance(args, (list, tuple))
        self._process = subprocess.Popen([exe, "-QR"])
        self.server = xmlrpclib.ServerProxy("http://%s:%d" % (self.host, self.port))
        time.sleep(3)

    def stop(self):
        if self._process_is_running():
            self._process.terminate()

    #Define helper methods

    #Atom Labeling
    def show_label(self):
        self.server.do('label all, index')
        self.server.do('label_size ,20')

    def clear_label(self):
        self.server.do('label all')

    def refresh_model(self):
        #self.server.do('orient')
        self.server.do('set sphere_scale, 0.2, (all)')
        self.server.do('set_bond stick_radius, 0.14, (all), (all)')
        self.server.do('color red,   elem O')
        self.server.do('color cyan,  elem C')
        self.server.do('color white, elem H')
        self.server.do('color blue,  elem N')
        self.server.do('color pink,  elem Si')
        self.server.do('color yellow,  elem S')
        self.server.do('show sticks')
        self.server.do('show spheres')
        return

    #Import files
    def file_import(self,file_paths,label=False):
        self.Grid_mode = False
        object_names = {}
        self.server.do('reinitialize')
        for path in file_paths:
            self.server.do('load '+path+',object_props=*')
            object_names[path.split('/')[-1].split('.')[0]] = 1
        self.server.do('orient')
        self.server.do('set sphere_scale, 0.2, (all)')
        self.server.do('set_bond stick_radius, 0.14, (all), (all)')
        self.server.do('color red,   elem O')
        self.server.do('color cyan,  elem C')
        self.server.do('color white, elem H')
        self.server.do('color blue,  elem N')
        self.server.do('color pink,  elem Si')
        self.server.do('color yellow,  elem S')
        self.server.do('show sticks')
        self.server.do('show spheres')
        if label:
            self.server.do('label all, index')
            self.server.do('label_size ,20')
        return object_names

    #Connect or break bonds in a single object
    def make_bond(self,atom_label1, atom_label2):
        self.server.do('select atom1, name %s' % atom_label1)
        self.server.do('select atom2, name %s' % atom_label2)
        self.server.do('make atom1, atom2')
        return 
    def break_bond(self,atom_label1, atom_label2):
        self.server.do('select atom1, name %s' % atom_label1)
        self.server.do('select atom2, name %s' % atom_label2)
        self.server.do('unbond atom1, atom2')
        return

    def grid_translate(self,object_names,distance=30):
        if self.Grid_mode == True:
            return

            return
        object_names2 = []
        for key, value in object_names.items():
            if value == 1:
                self.server.do('set_name %s, %s%d' % (key,key,1))
                object_names2.append("%s%d" % (key,1))
            else:
                self.server.do('set_name %s, %s%d' % (key,key,1))
                object_names2.append("%s%d" % (key,1))
                for j in range(value-1):
                    self.server.do('copy %s%d, %s%d' % (key,j+2,key,1))
                    object_names2.append("%s%d" % (key,j+2))
        y = 5
        x = int(len(object_names2) / y) + 1   
        pos = []
        
        for x_i in range(x):
            for y_i in range(y):
                x_p = x_i * distance
                y_p = y_i * distance
                pos.append([x_p,y_p,0])

        for i, obj in enumerate(object_names2):
            self.server.do('translate [%d,%d,-1], %s' % (pos[i][0],pos[i][1],obj))
        self.server.do('reset')
        self.Grid_mode = True
        return 

    def Get_Object_Names(self):
        self.server.do('run Get_Object_Names.py')
        f = open('object_name_mw_dbe.txt','r')
        lines = f.readlines()
        lines2 = [line.split()[1] for line in lines]
        object_names = dict(Counter(lines2))
        return object_names

    def Get_Object_Masses(self):
        self.server.do('run Get_Object_Names.py')
        f = open('object_name_mw_dbe.txt','r')
        lines = f.readlines()
        object_masses = {}
        for line in lines:
            line_s = line.split()
            object_name = line_s[1]
            object_mass = line_s[2]
            object_masses[object_name] = float(object_mass)
        return object_masses

    def Get_Atom_Names(self,object_name):
        self.server.do('select "mol1", %s' % object_name)
        self.server.do('run Get_Atom_Names.py')
        f = open('atom_names.txt','r')
        lines = f.readlines()
        atom_names = [line.split()[1] for line in lines]
        self.server.do('delete mol1')
        return atom_names

    def Crosslink_mols1(self,object1,index1,object2,index2):
        self.server.do('select at1, %s & index %s' % (object1, index1))
        self.server.do('select at2, %s & index %s' % (object2, index2))
        self.server.do('edit at1, at2')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % object1)
        self.server.do('rebuild %s' % object2)
        object_masses= self.Get_Object_Masses()
        object_names = list(object_masses.keys())
        # Check if object_names is empty
        while not object_names:
            object_masses = self.Get_Object_Masses()
            object_names = list(object_masses.keys())
        max_num = len(object_names)
        for i in range(max_num):
            new_object_name = 'combined_a%s' % (i+1)  
            if new_object_name not in object_names:
                break
        self.server.do('set_name %s, %s' % (object2,new_object_name))
        self.server.do('zoom %s' % new_object_name)
        self.server.do('clean %s' % new_object_name)
        self.server.do('sculpt %s' % new_object_name)
        self.refresh_model()
        self.clear_label()
        return new_object_name

    def Crosslink_mols2(self,object1,index1,object2,index2):
        self.server.do('select at1, %s & index %s' % (object1, index1))
        self.server.do('select at2, %s & index %s' % (object2, index2))
        self.server.do('edit at1, at2')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % object1)
        self.server.do('rebuild %s' % object2)
        object_masses= self.Get_Object_Masses()
        object_names = list(object_masses.keys())
        # Check if object_names is empty
        while not object_names:
            object_masses = self.Get_Object_Masses()
            object_names = list(object_masses.keys())
        max_num = len(object_names)
        for i in range(max_num):
            new_object_name = 'combined_b%s' % (i+1)  
            if new_object_name not in object_names:
                break
        self.server.do('set_name %s, %s' % (object2,new_object_name))
        self.server.do('zoom %s' % new_object_name)
        self.server.do('clean %s' % new_object_name)
        self.server.do('sculpt %s' % new_object_name)
        self.refresh_model()
        self.clear_label()
        return new_object_name
     
    def Crosslink_mols3(self,object1,index1,object2,index2):
        self.server.do('select at1, %s & index %s' % (object1, index1))
        self.server.do('select at2, %s & index %s' % (object2, index2))
        self.server.do('edit at1, at2')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % object1)
        self.server.do('rebuild %s' % object2)
        object_masses= self.Get_Object_Masses()
        object_names = list(object_masses.keys())
        # Check if object_names is empty
        while not object_names:
            object_masses = self.Get_Object_Masses()
            object_names = list(object_masses.keys())
        max_num = len(object_names)
        for i in range(max_num):
            new_object_name = 'combined_c%s' % (i+1)
            if new_object_name not in object_names:
                break
        self.server.do('set_name %s, %s' % (object2,new_object_name))
        self.server.do('zoom %s' % new_object_name)
        self.server.do('clean %s' % new_object_name)
        self.server.do('sculpt %s' % new_object_name)
        self.refresh_model()
        self.clear_label()
        return new_object_name
    
    def Crosslink_mols4(self,object1,index1,object2,index2):
        self.server.do('select at1, %s & index %s' % (object1, index1))
        self.server.do('select at2, %s & index %s' % (object2, index2))
        self.server.do('edit at1, at2')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % object1)
        self.server.do('rebuild %s' % object2)
        object_masses= self.Get_Object_Masses()
        object_names = list(object_masses.keys())
        # Check if object_names is empty
        while not object_names:
            object_masses = self.Get_Object_Masses()
            object_names = list(object_masses.keys()) 
        max_num = len(object_names)
        for i in range(max_num):
            new_object_name = 'combined_d%s' % (i+1)
            if new_object_name not in object_names:
                break
        self.server.do('set_name %s, %s' % (object2,new_object_name))
        self.server.do('zoom %s' % new_object_name)
        self.server.do('clean %s' % new_object_name)
        self.server.do('sculpt %s' % new_object_name)
        self.refresh_model()
        self.clear_label()
        return new_object_name
    
    def Crosslink_mols5(self,object1,index1,object2,index2):
        self.server.do('select at1, %s & index %s' % (object1, index1))
        self.server.do('select at2, %s & index %s' % (object2, index2))
        self.server.do('edit at1, at2')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % object1)
        self.server.do('rebuild %s' % object2)
        object_masses= self.Get_Object_Masses()
        object_names = list(object_masses.keys())
        # Check if object_names is empty
        while not object_names:
            object_masses = self.Get_Object_Masses()
            object_names = list(object_masses.keys())  
        max_num = len(object_names)
        for i in range(max_num):
            new_object_name = 'combined_e%s' % (i+1)
            if new_object_name not in object_names:
                break
        self.server.do('set_name %s, %s' % (object2,new_object_name))
        self.server.do('zoom %s' % new_object_name)
        self.server.do('clean %s' % new_object_name)
        self.server.do('sculpt %s' % new_object_name)
        self.refresh_model()
        self.clear_label()
        return new_object_name
    
    def Crosslink_mols6(self,object1,index1,object2,index2):
        self.server.do('select at1, %s & index %s' % (object1, index1))
        self.server.do('select at2, %s & index %s' % (object2, index2))
        self.server.do('edit at1, at2')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % object1)
        self.server.do('rebuild %s' % object2)
        object_masses= self.Get_Object_Masses()
        object_names = list(object_masses.keys())
        # Check if object_names is empty
        while not object_names:
            object_masses = self.Get_Object_Masses()
            object_names = list(object_masses.keys())  
        max_num = len(object_names)
        for i in range(max_num):
            new_object_name = 'combined_f%s' % (i+1)
            if new_object_name not in object_names:
                break
        self.server.do('set_name %s, %s' % (object2,new_object_name))
        self.server.do('zoom %s' % new_object_name)
        self.server.do('clean %s' % new_object_name)
        self.server.do('sculpt %s' % new_object_name)
        self.refresh_model()
        self.clear_label()
        return new_object_name

    def Bond_mol(self,object,index1,index2):
        self.server.do('select at1, %s & index %s' % (object, index1))
        self.server.do('select at2, %s & index %s' % (object, index2))
        self.server.do('edit at1, at2')
        self.server.do('bond')
        self.server.do('h_fill')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('set_name %s, %s' % (object,object+'_macrocycle'))
        self.server.do('set_name %s, %s' % (object+'_macrocycle',object))
        self.server.do('rebuild %s' % object)
        self.server.do('zoom %s' % object)
        self.server.do('clean %s' % object)
        self.refresh_model()
        self.show_label()
        return object

    def RemoveH(self,object1,atom_name1):
        self.server.do('select at1, %s & name %s' % (object1, atom_name1))
        self.server.do('remove_picked')
        return

    def Attach_CH3(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/methane.mol2')
        self.server.do('select at2, %s & name %s' % ('Methane', 'H4'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'Methane')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return

    def Attach_CH2CH3(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/ethane.mol2')
        self.server.do('select at2, %s & name %s' % ('ethane', 'H4'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'ethane')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return
    
    def Attach_C2H4CH3(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/propane.mol2')
        self.server.do('select at2, %s & name %s' % ('propane', 'H6'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'propane')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return

    def Attach_OH(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/water.mol2')
        self.server.do('select at2, %s & name %s' % ('Water', 'H2'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'Water')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return
    
    def Attach_SH(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/sulfide.mol2')
        self.server.do('select at2, %s & name %s' % ('Sulfide', 'H2'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'Sulfide')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return

    def Attach_Aldehyde(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/aldehyde.mol2')
        self.server.do('select at2, %s & name %s' % ('Aldehyde', 'H1'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'Aldehyde')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return

    
    def Attach_Ketone(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/Ketone.pdb')
        self.server.do('select at2, %s & name %s' % ('Ketone', 'H1'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'Ketone')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return

    def Attach_Ester(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/cooh.mol2')
        self.server.do('select at2, %s & name %s' % ('cooh', 'H2'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'cooh')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return

    # This would include 3HC-COO-CH2-C group
    def Attach_Ester_2(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/ethyl_acetate.mol2')
        self.server.do('select at2, %s & name %s' % ('ethyl_acetate', 'H1'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'ethyl_acetate')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return
    
    def Attach_Aniline(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/aniline.mol2')
        self.server.do('select at2, %s & name %s' % ('Aniline', 'H1'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'Aniline')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return
    
    def Attach_Carbonyl(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/formaldehyde.mol2')
        self.server.do('select at2, %s & name %s' % ('formaldehyde', 'H2'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'formaldehyde')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return
    def Attach_Carbox(self,object1,h_index):
        self.server.do('select at1, %s & index %s' % (object1, h_index))
        self.server.do('load ./functional_groups/cooh_2.mol2')
        self.server.do('select at2, %s & name %s' % ('cooh_2', 'H1'))
        self.server.do('edit at2, at1')
        self.server.do('fuse')
        self.server.do('unpick')
        self.server.do('delete at1')
        self.server.do('delete at2')
        self.server.do('delete %s' % 'cooh_2')
        self.server.do('rebuild all')
        self.server.do('zoom %s' % object1)
        self.server.do('clean %s' % object1)
        self.refresh_model()
        self.show_label()
        return


    def check_main(self, object, idx):
        self.server.do('save test.mol,%s' % object)
        mol = AllChem.MolFromMolFile('test.mol',sanitize=False,removeHs=False)
        main_atom_idx= [atom.GetIdx()+1 for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']
        main_atom    = [atom for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']
        main_atom_nn = [[n for n in atom.GetNeighbors() if n.GetSymbol() != 'H'] for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']

        index_check = main_atom_idx.index(idx)
        check_atom_symbol = main_atom[index_check].GetSymbol()
        check_atom_nn = main_atom_nn[index_check]
        check_atom_nn_symbol = [n.GetSymbol() for n in main_atom_nn[index_check]]
        check_atom_ar = main_atom[index_check].GetIsAromatic()
        check_atom_inring = main_atom[index_check].IsInRing()
        subprocess.call('rm test.mol',shell=True)
        return check_atom_symbol, check_atom_nn_symbol, check_atom_ar, check_atom_inring

  
    def Change_Element_CtoO(self,object1,index):
        success = False
        self.server.do('zoom %s' % object1)
        s, ns, ar_flag, ring_flag = self.check_main(object1,index)
        if s == 'C' and len(ns) == 2 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('valence 1, at1, elem *')
            self.server.do('replace O, 4, 2')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        elif s == 'C' and len(ns) == 1 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('valence 1, at1, elem *')
            self.server.do('replace O, 4, 2')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        self.refresh_model()
        self.show_label()
        return success
    
    def Change_Element_CtoN(self,object1,index):
        success = False
        self.server.do('zoom %s' % object1)
        s, ns, ar_flag, ring_flag = self.check_main(object1,index)
        if s == 'C' and len(ns) == 3 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('valence 1, at1, elem *')
            self.server.do('replace N, 4, 3')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        elif s == 'C' and len(ns) == 2 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('valence 1, at1, elem *')
            self.server.do('replace N, 4, 3')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        elif s == 'C' and len(ns) == 1 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('valence 1, at1, elem *')
            self.server.do('replace N, 4, 3')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        self.refresh_model()
        self.show_label()
        return success

    def Change_Element_CtoS(self,object1,index):
        success = False
        self.server.do('zoom %s' % object1)
        s, ns, ar_flag, ring_flag = self.check_main(object1,index)
        if s == 'C' and len(ns) == 2 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('valence 1, at1, elem *')
            self.server.do('replace S, 4, 2')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        elif s == 'C' and len(ns) == 1 and ('O' not in ns)  and ('N' not in ns):
            self.server.do('select at1, %s & index %s' % (object1, index))
            self.server.do('edit at1')
            self.server.do('replace S, 4, 2')
            self.server.do('unpick')
            self.server.do('delete at1')
            self.server.do('rebuild all')
            self.server.do('clean %s' % object1)
            success = True
        self.refresh_model()
        self.show_label()
        return success

    def examine_h(self,object):
        self.server.do('zoom %s' % object)
        #self.server.do('rotate y, 0.5, %s' % object)
        self.server.do('save test.mol,%s' % object)
        #self.server.do('rotate y, -0.5, %s' % object)
        self.refresh_model()
 #       mol = AllChem.MolFromMolFile('test.mol',removeHs=False)
        mol = AllChem.MolFromMolFile('test.mol',sanitize=False,removeHs=False)
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
        
        self.show_label()
        #subprocess.call('rm test.mol',shell=True)
        return types

    def examine_main(self,object):
        self.server.do('zoom %s' % object)
        self.server.do('save test.mol,%s' % object)
        self.refresh_model()
     #   mol = AllChem.MolFromMolFile('test.mol',removeHs=False)
        mol = AllChem.MolFromMolFile('test.mol',sanitize=False,removeHs=False)
        main_atom    = [atom for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']
        main_atom_nn = [[n for n in atom.GetNeighbors()] for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']

        types = {'Si_1n':[],'Si_2n':[],'Si_3n':[],'Si_4n':[],
                 'C_1n':[],'C_2n':[],'C_3n':[],'C_4n':[],
                 'N_1n':[],'N_2n':[],'N_3n':[],
                 'O_1n':[],'O_2n':[],
                 'S_1n':[],'S_2n':[],
                }

        for i, matom in enumerate(main_atom):
            nsymbols = [n.GetSymbol() for n in matom.GetNeighbors() if n.GetSymbol() != 'H']
            if matom.GetSymbol() == 'C' and len(nsymbols) == 1:
                types['C_1n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'C' and len(nsymbols) == 2:
                types['C_2n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'C' and len(nsymbols) == 3:
                types['C_3n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'C' and len(nsymbols) == 4:
                types['C_4n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'Si' and len(nsymbols) == 1:
                types['Si_1n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'Si' and len(nsymbols) == 2:
                types['Si_2n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'Si' and len(nsymbols) == 3:
                types['Si_3n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'Si' and len(nsymbols) == 4:
                types['Si_4n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'N' and len(nsymbols) == 1:
                types['N_1n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'N' and len(nsymbols) == 2:
                types['N_2n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'N' and len(nsymbols) == 3:
                types['N_3n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'O' and len(nsymbols) == 1:
                types['O_1n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'O' and len(nsymbols) == 2:
                types['O_2n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'S' and len(nsymbols) == 1:
                types['S_1n'].append(matom.GetIdx()+1)
            elif matom.GetSymbol() == 'S' and len(nsymbols) == 2:
                types['S_2n'].append(matom.GetIdx()+1)

        self.show_label()
        subprocess.call('rm test.mol',shell=True)
        return types

    def examine_cnmr(self,object,color=False):
        self.server.do('save test.mol,%s,0,mol' % object)
        self.refresh_model()
        mol = AllChem.MolFromMolFile('test.mol',sanitize=False,removeHs=False)
        atom_idx = [atom.GetIdx()+1 for atom in mol.GetAtoms() if atom.GetSymbol() != 'H']
        mol = AllChem.RemoveHs(mol,sanitize=False)
        atom_idx2 = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol()]
        cnmr_data, index = carbon_nmr_individual(mol)

        for key, value in index.items():
            new_value = []
            for v in value:
                idx = atom_idx2.index(v)
                new_v = atom_idx[idx]
                new_value.append(new_v)
            index[key] = new_value

        if color: 
            self.server.do('color red,   elem O')
            self.server.do('color cyan,  elem C')
            self.server.do('color white, elem H')
            self.server.do('color blue,  elem N')
            self.server.do('color pink,  elem Si')
            self.server.do('color yellow,  elem S')
            self.server.do('show sticks')
            self.server.do('show spheres')
            for id in index['faH']:
                self.server.do('color skyblue, %s and index %d' % (object,id))
            for id in index['faP']:
                self.server.do('color sand, %s and index %d' % (object,id))
            for id in index['faS']:
                self.server.do('color salmon, %s and index %d' % (object,id))
            for id in index['faB']:
                self.server.do('color violet, %s and index %d' % (object,id))
            for id in index['faC']:
                self.server.do('color raspberry, %s and index %d' % (object,id))
            for id in index['fal']:
                self.server.do('color olive, %s and index %d' % (object,id))
            for id in index['fal*']:
                self.server.do('color lightteal, %s and index %d' % (object,id))
            for id in index['falH']:
                self.server.do('color lightorange, %s and index %d' % (object,id))
            for id in index['falO']:
                self.server.do('color hotpink, %s and index %d' % (object,id))
            self.server.do('set sphere_scale, 0.2, (%s)' % object)
            self.server.do('set_bond stick_radius, 0.14, (%s), (%s)' % (object,object))
        subprocess.call('rm test.mol',shell=True)
        return cnmr_data, index

    def Get_Atom_Index(self,object_name):
        self.server.do('select "mol1", %s' % object_name)
        self.server.do('run Get_Atom_Index.py')
        f = open('atom_index.txt','r')
        lines = f.readlines()
        atom_index = [line.split()[1] for line in lines]
        self.server.do('delete mol1')
        return atom_index

    def Get_Atom_Info(self,object_name):
        self.server.do('select "mol1", %s' % object_name)
        self.server.do('run Get_Atom_Info.py')
        f = open('atom_info.txt','r')
        lines = f.readlines()
        atom_info = [line.split()[1:] for line in lines]
        self.server.do('delete mol1')
        return atom_info
