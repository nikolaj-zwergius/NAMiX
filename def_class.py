
import re #needed for regex split in placing replacement atoms
from copy import deepcopy #only used in non residue class used in non used from_cif_pdb rep_func
import numpy as np
import os
py_path = os.path.dirname(os.path.realpath(__file__))
class Wrong_Base_Error(Exception): # exception for if modification is not for the base is should replace
    def __init__(self,new:str,old:str,name:str ) -> None:
        self.new = new
        self.old = old
        self.name = name
    def __str__(self) -> str:
        s = (f"tried to replace {self.old} with a {self.new} modification: "
        "this will give raise to wrong modifications as only diffence "
        f"from {self.new} will be modified. Try changing the '{self.old} = {self.name}' line to a valid {self.old} modification")

        return s


class modnuc: # overall class for modified nucleotide
    
    def id_atom_fixer(self, atom:str): # code for spitting the given replacement into the parts need to properly use them for pdb file generation
        atom = re.split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",atom)
        return atom


    def __init__(self,name:str,type:str,old_base:str,replace_module,cif_path=None,replacements:list = None, add_sugar:list = None, add_base:list =  None,description="no description") -> None:
      
        self.name = name
        self.old_base = old_base
        self.type = type
        
        
        self.have_additions = False
        self.replacements = []
        self.removal_steps = []
        self.have_replacements = False
        self.rep_module = replace_module
        self.cif_path = cif_path 
        self.description  = description

        self.inverted  = False

        self.calculated_vector = False
        self.mods = []
        if replacements: # convert the list of replacments to the co
            for i in replacements:
                self.add_replacement(i[0],i[1])
        
        if add_sugar:
            self.add_sugar(add_sugar)
        
        if add_base:
            self.add_base(add_base)
  

    def add_replacement(self,start,replacement): # code for adding replacment table
        start = self.id_atom_fixer(start)
        if replacement != None:
            replacement = self.id_atom_fixer(replacement)
        else:
            replacement = [None]
        self.replacements.append((start,replacement) )
        self.have_replacements = True

    def invers_replacement(self):
        for i in self.replacements:
            if i[1][0] != None:
                self.removal_steps.append(i[::-1])
        self.type = "ATOM"
        self.name = self.old_base + "  "
        self.replacements = self.removal_steps
        self.inverted = True

    def add_sugar(self,add_sugar):
        self.core_sugar = add_sugar[0]
        self.additions_sugar = add_sugar[1]
        self.core_coord_sugar = np.empty(shape = (len(add_sugar[0]),3))
        self.additions_coord_sugar = np.empty(shape = (len(add_sugar[1]),3))
        self.add_additions(add_sugar[0],add_sugar[1],self.additions_coord_sugar,self.core_coord_sugar)

    def add_base(self,add_base):
        self.core_base = add_base[0]
        self.additions_base = add_base[1]
        self.core_coord_base = np.empty(shape = (len(add_base[0]),3))
        self.additions_coord_base = np.empty(shape = (len(add_base[1]),3))
        self.add_additions(add_base[0],add_base[1],self.additions_coord_base,self.core_coord_base)


    def add_additions(self,origins,add,add_cord,core_cord):
        pos_dict = self.cif_reader()
        for i,j in enumerate(add):
            add_cord[i] = pos_dict[j]
        for i,j in enumerate(origins):
            #print(self.name,i,j,core_cord)
            core_cord[i] = pos_dict[j]

    def cif_reader (self):
        with open(f"{py_path}{self.cif_path}") as f:
            x_found = False
            y_found = False
            z_found = False
            cif_pos_dict = {}
            for line in f:

                if x_found and y_found and z_found:
                    if line.strip() == "#": break
                    coord_temp = line[40:].split()
                    for coord, i in enumerate(coord_temp):
                        coord_temp[coord] = float(i)
                    cif_pos_dict[line[11:15].strip()] = np.array(coord_temp)
                if line.strip() == "_chem_comp_atom.x": x_found = True
                if line.strip() == "_chem_comp_atom.y": y_found = True
                if line.strip() == "_chem_comp_atom.z": z_found = True
            return cif_pos_dict

    def __str__(self) -> str:
        return f"{self.name} : {self.description}"
    def __repr__(self) -> str:
        return f"{self.name}"
    

class no_mod(modnuc): # supclass for non modified as they are need for base case
    def __init__(self,name,no_replace) -> None:
        super().__init__(name, "ATOM",name,no_replace)
    def gen_atom_id_from_ref_pdb(self):
        return