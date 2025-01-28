# Change the filename.pdb to the target file and the output filename to desired name
# This script works for pdb files that do not have hydrogens included

from def_class import modnuc, no_mod
from rep_func import no_replace,simple_replace_module,addition_atom_generation


    

A_no_mod = no_mod("A",no_replace)
C_no_mod = no_mod("C",no_replace)
U_no_mod = no_mod("U",no_replace)
T_no_mod = no_mod("T",no_replace)
G_no_mod = no_mod("G",no_replace)

UFT = modnuc("UFT","HETATM","U",simple_replace_module,"/mods/UFT.cif", replacements= [("O2'","F2'"),("HO2'",None)],description="2'F Uracil") #  2'F C
CFZ = modnuc("CFZ","HETATM","C",simple_replace_module,"/mods/CFZ.cif", replacements= [("O2'","F2'"),("HO2'",None)],description="2'F Cytosine") #  2'F C
AF2 = modnuc("AF2","HETATM","A",simple_replace_module,None, replacements= [("O2'","F2'"),("HO2'",None)],description="2'F Guanine | Has no cif file yet") #  2'F A
GF2 = modnuc("GF2","HETATM","G",simple_replace_module,None, replacements= [("O2'","F2'"),("HO2'",None)],description="2'F Adenine | Has no cif file yet") #  2'F G

UMX = modnuc("UMX","HETATM","U",addition_atom_generation,None,replacements= [("O2'",None),("HO2'",None)], add=[("C3'","C2'","C6'"),("C3'","C2'","O2'")],description="Locked Uracil | Not implemented yet")# LNA  U NOT DONE
LCC = modnuc("LCC","HETATM","C",addition_atom_generation,None,replacements= [("O2'",None),("HO2'",None)], add=[("C3'","C2'","C6'"),("C3'","C2'","O2'")],description="Locked Cytosine | Not implemented yet")# LNA  C NOT DONE
LCA = modnuc("LCA","HETATM","A",addition_atom_generation,"/mods/LCA.cif",replacements= [("O2'",None),("HO2'",None)], add=[("C3'","C2'","C6'"),("C3'","C2'","O2'")],description="Locked Adenine | cant be reverted by NAMIX")# LNA  A NOT DONE
LCG = modnuc("LCG","HETATM","G",addition_atom_generation,None,replacements= [("O2'",None),("HO2'",None)], add=[("C3'","C2'","C6'"),("C3'","C2'","O2'")],description="Locked Guanine | Not implemented yet")# LNA  G NOT DONE


mod_dict ={UFT.name:UFT,CFZ.name:CFZ,AF2.name:AF2,GF2.name:GF2,
           UMX.name:UMX,LCC.name:LCC,LCA.name:LCA,LCG.name:LCG
           }
