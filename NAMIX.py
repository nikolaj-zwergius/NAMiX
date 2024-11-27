import os
import sys
from io import TextIOWrapper
dir_path = os.path.dirname(os.path.realpath(__file__))
parent_dir_path = os.path.abspath(os.path.join(dir_path, os.pardir))


from def_class import modnuc, Wrong_Base_Error
from define_mods import *
from hydrogen_bond_restrint_gen import restrint_from_road, restrint_from_pb

from shutil import copy

py_path = os.path.dirname(os.path.realpath(__file__))
wd = os.getcwd()

def _base_equality_check(base:modnuc,oldbase:no_mod):
    if base.old_base != oldbase.old_base:
        raise Wrong_Base_Error(base.old_base, oldbase.old_base,base.name)
    return

def _copy_cif(mod:modnuc,dir):
    if mod.cif_path:
        copy(f"{py_path}{mod.cif_path}",f"{wd}/{dir}")
    elif type(mod) == no_mod:
        return
    else:
        print(f"{mod.name} does not have a cif file path define in the define_mods.py for futher use it is recommeded to define this paht")

def run_modules(residue:dict,atom_nr:int,A:modnuc,C:modnuc,G:modnuc,T:modnuc,U:modnuc,out:TextIOWrapper):
    key = list(residue.keys())[0]
    if residue[key][17:20].strip() == "U": atom_nr = U.rep_module(residue,U,atom_nr,out)
    elif residue[key][17:20].strip() == "A":atom_nr =  A.rep_module(residue,A,atom_nr,out)
    elif residue[key][17:20].strip() == "G": atom_nr = G.rep_module(residue,G,atom_nr,out)
    elif residue[key][17:20].strip() == "C": atom_nr =  C.rep_module(residue,C,atom_nr,out)
    elif residue[key][17:20].strip() == "T": atom_nr = T.rep_module(residue,T,atom_nr,out)
    else:
        for atom in residue:
            out.write(residue[atom])
            atom_nr += 1
    return atom_nr


def NAMIX(pdbfile:str,restrin_file:str=None,A:modnuc=A_no_mod, C:modnuc=C_no_mod, G:modnuc=G_no_mod, T:modnuc=T_no_mod, U:modnuc=U_no_mod, Prefix:str ="",overwrite = False,min=False) -> None:
    standard_bases = [A_no_mod,U_no_mod,C_no_mod,G_no_mod,T_no_mod]
    bases = [A,U,C,G,T]
    for i in range(5): # check that modnucs used match the old base they are intented to replace
        _base_equality_check(bases[i],standard_bases[i])
    
    try:
        dir_path = f"{pdbfile}_{Prefix}mod"
        os.mkdir(dir_path)
    except FileExistsError:
        if not overwrite:
            print("\ndictory alread exists set overwrite = True if folder content should be overwritten\n")
            raise
    
    try:
        open(pdbfile + ".pdb")
    except FileNotFoundError:
        if not os.path.getsize(dir_path):
            print("\n\b NO SOURCE PDB FILE \n")
            os.rmdir(dir_path)
        else:
            print("\nNO SOURCE PDB FILE: folder with data for this file found\n ")
        raise

    if restrin_file:
        if restrin_file[0].endswith(".pb"):
            restrint_from_pb(restrin_file[0],dir_path,min)
        else:
            restrint_from_road(restrin_file[0],restrin_file[1],dir_path,min,restrin_file[2])
    if not min:
        for i in range(5):#copy .cif into folder for XNA
            _copy_cif(bases[i],dir_path)


    with open(pdbfile + ".pdb") as f, open(f"{dir_path}/{Prefix}{pdbfile}_mod.pdb", 'w') as out:
            atom_nr = 1
            residue = {}
            current_res = 0
            
            for line in f:
                if line.startswith("HETATM") and line[17] == "R":#QRNA fix
                    line = f"{'ATOM':<6}{atom_nr:>5}{line[11:17]}{line[18]:<4}{line[21:]}"
                
                if line.startswith("SEQRES"):#SEQRES line mod and output NEED fix to work with specifed residues
                    line = line[:16] + line[16:-2].replace("A  ",f"{A.name:<3}") + line[-2:].replace("A",f"{A.name:<3}")
                    line = line[:16] + line[16:-2].replace("U  ",f"{U.name:<3}") + line[-2:].replace("U",f"{U.name:<3}")
                    line = line[:16] + line[16:-2].replace("C  ",f"{C.name:<3}") + line[-2:].replace("C",f"{C.name:<3}")
                    line = line[:16] + line[16:-2].replace("G  ",f"{G.name:<3}") + line[-2:].replace("G",f"{G.name:<3}")
                    line = line[:16] + line[16:-2].replace("T  ",f"{T.name:<3}") + line[-2:].replace("T",f"{T.name:<3}")
                    out.write(line)
                    continue
               
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if current_res == 0: #start case
                        current_res = int(line[22:26].strip())
                    if current_res != int(line[22:26].strip()): # new residue case
                        atom_nr = run_modules(residue,atom_nr,A,C,G,T,U,out)
                        residue = {}
                        current_res = int(line[22:26].strip())
                    if line[12:16].strip() not in residue.keys():# base case
                        residue[line[12:16].strip()] = line

            atom_nr = run_modules(residue,atom_nr,A,C,G,T,U,out)
                         
    if not min:
        copy(f"{pdbfile}.pdb",f"{dir_path}/{pdbfile}.pdb")

def rev_namix(pdbfile:str): # does not work with add mod like LCA rigth now
    mods = mod_dict.keys()
    atom_nr = 1
    with open(pdbfile + ".pdb") as f, open(pdbfile + "_rev.pdb","w") as out:
        for line in f:
            if line[17:20].strip() in mods:
                if not mod_dict[line[17:20].strip()].inveresed:
                    mod_dict[line[17:20].strip()].invers_replacment()
                simple_replace_module({"A":line},mod_dict[line[17:20].strip()],atom_nr,out)
            else:
                if not line.startswith("SEQRES"):
                    out.write(line)
            atom_nr += 1
            
            if line.startswith("SEQRES"):#SEQRES line mod and output
                    for mod in mods:
                        line = line.replace(mod,f"{mod_dict[mod].old_base:<3}")
                    out.write(line)
                    continue

if __name__ == "__main__":
    NAMIX("Relaxed_no_fluorophore",U=UFT,C=CFZ,overwrite=True)



