# Change the filename.pdb to the target file and the output filename to desired name
# This script works for pdb files that do not have hydrogens included
from io import TextIOWrapper
from def_class import modnuc
import numpy as np

def new_pdb_line(org_line:str,type:str,name:str,atom_nr:int,rep:str="") -> str:
    new_line = []
    new_line.append(f"{type:<6}")
    new_line.append(f"{atom_nr:>5}"+" ")
    if rep:
        if len(rep[0]) == 2:
            new_line.append(f"{''.join(rep)<4}")
        else:
            tmp_line = "".join(rep)
            if len(tmp_line) == 4:
                new_line.append(tmp_line + " ")
            elif  len(tmp_line) == 3:
                new_line.append(" " + tmp_line+" ")
            elif  len(tmp_line) == 2:
                new_line.append(" " + tmp_line +"  ")
            else:
                new_line.append(" "+ tmp_line + "   ")
        
        new_line.append(f'{name:>3}'+" ")
        new_line.append(org_line[21:76])
        new_line.append(f"{rep[0]:^4}")
        new_line.append("  \n")
    
    if not rep:
        new_line.append(org_line[12:17])
        new_line.append(f'{name:>3}'+" ")
        new_line.append(org_line[21:])
    #print(new_line)
    return "".join(new_line)

def align_vectors(a, b):
    b = b / np.linalg.norm(b) # normalize a
    a = a / np.linalg.norm(a) # normalize b
    v = np.cross(a, b)
    # s = np.linalg.norm(v)
    c = np.dot(a, b)
    if np.isclose(c, -1.0):
        return -np.eye(3, dtype=np.float64)

    v1, v2, v3 = v
    h = 1 / (1 + c)

    Vmat = np.array([[0, -v3, v2],
                  [v3, 0, -v1],
                  [-v2, v1, 0]])

    R = np.eye(3, dtype=np.float64) + Vmat + (Vmat.dot(Vmat) * h)
    return R

def no_replace(res:dict,mod:modnuc,atom_nr:int,out:TextIOWrapper) -> None:
    for atom in res.keys():
        line = res[atom]
        out.write("".join([line[0:5],f"{atom_nr:>6}",line[11:]]))
        atom_nr += 1
    return atom_nr

def simple_replace_module(res:dict,mod:modnuc,atom_nr:int,out:TextIOWrapper) -> None:
    if (int(res[list(res.keys())[0]][23:26]) not in mod.mods) and (mod.mods != []):
        atom_nr = no_replace(res,mod,atom_nr,out)
        return atom_nr
    for atom in res.keys():
        line = res[atom]
        modifed = False #set defualt at not changed
        if mod.have_replacements: # cheking for replacments in the modnuc
            for rep in mod.replacements: # testing all replacmients in modnuc
                if line[11:16].strip() == "".join(rep[0]): # checking for match in model
                    if rep[1][0] == None:
                        atom_nr -= 1
                    else:
                        out.write(new_pdb_line(line,mod.type,mod.name,atom_nr,rep = rep[1])) #gerneration new line
                    modifed =  True # return to inform change made
            if mod.have_replacements and not modifed:
                out.write(new_pdb_line(line,mod.type,mod.name,atom_nr)) # update type and name to match rest of nuclotide
        atom_nr +=1
    return atom_nr

def addition_atom_generation(res:dict,mod:modnuc,atom_nr:int,out:TextIOWrapper): #placeholder for adddtion logic

    if (int(res[list(res.keys())[0]][23:26]) not in mod.mods) and (mod.mods != []):
        atom_nr = no_replace(res,mod,atom_nr,out)
        return atom_nr

    atom_nr = simple_replace_module(res,mod,atom_nr,out)
    
    if not mod.calculated_vector: mod.calculate_vectors_for_addition()
    
    for add in mod.additions:

        res_compare_coord = res["".join(add[1])][31:54].split()
        res_origin_coord = res["".join(add[0])][31:54].split()

        for i in range(3):
            res_origin_coord[i] = float(res_origin_coord[i])
            res_compare_coord[i] = float(res_compare_coord[i])

        res_vec = np.array(res_compare_coord)-np.array(res_origin_coord)

        R=align_vectors(add[4],res_vec)
        translat = res_origin_coord-np.array(add[3])

        point_coord = add[3]+R.dot(add[5])+translat
        items = [mod.type,atom_nr,"".join(add[2]),mod.name,res["".join(add[0])][21:26],point_coord[0],point_coord[1],point_coord[2],res["".join(add[0])][56:66],add[2][0]]
        spot3 = "".join(add[2])
        out.write(f"{items[0]:<6}{atom_nr:>5}  {spot3:^4}{items[3]:>3} {items[4]}     {items[5]:>7.6g} {items[6]:>7.6g} {items[7]:>7.6g}  {items[8]}           {items[9]}\n")
        atom_nr += 1


    return atom_nr