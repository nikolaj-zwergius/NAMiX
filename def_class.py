
import re #needed for regex split in placing replacment atoms
from copy import deepcopy #only used in non residue class used in non used from_cif_pdb rep_func
import numpy as np
import os
py_path = os.path.dirname(os.path.realpath(__file__))
class Wrong_Base_Error(Exception): # exception for if modication is not for the base is should replace
    def __init__(self,new:str,old:str,name:str ) -> None:
        self.new = new
        self.old = old
        self.name = name
    def __str__(self) -> str:
        s = (f"tried to replace {self.old} with a {self.new} modifcation: "
        "this will give raise to wrong modifications as only diffence "
        f"from {self.new} will be modifed. Try changing the '{self.old} = {self.name}' line to a vaild {self.old} modifcation")

        return s


class modnuc: # overall class for modnucliotide
    
    def id_atom_fixer(self, atom:str): # code for spitting the given replacemet into the parts need to proplery use them for pdb file generation
        atom = re.split("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",atom)
        return atom

    def __init__(self,name:str,type:str,old_base:str,replace_module,cif_path=None,rep:list = None, add:list = None,description="no desribtion") -> None:
        self.name = name
        self.old_base = old_base
        self.type = type
        self.addiations = []
        self.have_addtions = False
        self.replacments = []
        self.removale_steps = []
        self.have_replacments = False
        self.rep_module = replace_module
        self.cif_path = cif_path 
        self.description  = description
        self.inveresed  = False
        self.calculated_vector = False
        self.mods = []
        if rep:
            for i in rep:
                self.add_replacment(i[0],i[1])
        if add:
            for i in add:
                self.add_addtions(i[0],i[1],i[2])
  

    def add_replacment(self,start,replacment): # code for adding replacment table
        start = self.id_atom_fixer(start)
        if replacment != None:
            replacment = self.id_atom_fixer(replacment)
        else:
            replacment = [None]
        self.replacments.append((start,replacment) )
        self.have_replacments = True

    def invers_replacment(self):
        for i in self.replacments:
            if i[1][0] != None:
                self.removale_steps.append(i[::-1])
        self.type = "ATOM"
        self.name = self.old_base
        self.replacments = self.removale_steps
        self.inveresed = True

    def add_addtions(self,origin,compare,add):
        origin = self.id_atom_fixer(origin)
        compare = self.id_atom_fixer(compare)
        add = self.id_atom_fixer(add)
        self.addiations.append([origin,compare,add])
        self.have_addtions = True

    def calculate_vectors_for_addtion(self):
        pos_dict = self.cif_reader()
       
        for i in range(len(self.addiations)):
            
            origin_coor = pos_dict["".join(self.addiations[i][0])]
            compare = pos_dict["".join(self.addiations[i][1])]
            addtion = pos_dict["".join(self.addiations[i][2])]

            compare_vector = compare-origin_coor
            addtion_vector = addtion-origin_coor 
            self.addiations[i].extend([origin_coor,compare_vector,addtion_vector])
        self.calculated_vector = True

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
                    for coor, i in enumerate(coord_temp):
                        coord_temp[coor] = float(i)
                    cif_pos_dict[line[11:15].strip()] = np.array(coord_temp)
                if line.strip() == "_chem_comp_atom.x": x_found = True
                if line.strip() == "_chem_comp_atom.y": y_found = True
                if line.strip() == "_chem_comp_atom.z": z_found = True
            return cif_pos_dict

    def __str__(self) -> str:
        return f"{self.name} : {self.desribstion}"
    def __repr__(self) -> str:
        return f"{self.name}"
    

class no_mod(modnuc): # supclass for non modifed as they are need for base case
    def __init__(self,name,no_replace) -> None:
        super().__init__(name, "ATOM",name,no_replace)
    def gen_atom_id_from_ref_pdb(self):
        return