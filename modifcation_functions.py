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

def umeyama(P, Q):
    assert P.shape == Q.shape
    n, dim = P.shape

    centeredP = P - P.mean(axis=0)
    centeredQ = Q - Q.mean(axis=0)

    C = np.dot(np.transpose(centeredP), centeredQ) / n

    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    R = np.dot(V, W)

    varP = np.var(P, axis=0).sum()
    c = 1/varP * np.sum(S) # scale factor

    t = Q.mean(axis=0) - P.mean(axis=0).dot(c*R)

    return c, R, t

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
    
    
    
    core_coords_sugar = np.empty((len(mod.core_sugar),3))
    for i, add in enumerate(mod.core_sugar):
       core_coords_sugar[i] = res["".join(add)][31:54].split()
    
    for i in mod.additions_sugar:
        res[i] = []
    
    c,R,t = umeyama(mod.core_coord_sugar,core_coords_sugar)
    add_coords_sugar = mod.additions_coord_sugar.dot(c*R)+t
    

    core_coords_base = np.empty((len(mod.core_base),3))
    for i, add in enumerate(mod.core_base):
       core_coords_base[i] = res["".join(add)][31:54].split()
    
    for i in mod.additions_base:
        res[i] = []

    c,R,t = umeyama(mod.core_coord_base,core_coords_base)
    add_coords_base = mod.additions_coord_base.dot(c*R)+t
    mod_index_sugar = 0
    mod_index_base = 0
    for i,j in enumerate(res.keys()):
        
        #print(j)
        if j not in mod.additions_sugar and j not in mod.additions_base:
            #print(j)
            atom = {j:res[j]}
            atom_nr = simple_replace_module(atom,mod,atom_nr,out)
        elif j in mod.additions_sugar:
            items = [mod.type,atom_nr,j,mod.name,res["".join(add)][21:26],add_coords_sugar[mod.addtion_index[j]][0],add_coords_sugar[mod.addtion_index[j]][1],add_coords_sugar[mod.addtion_index[j]][2],"1.00  0.00",mod.additions_sugar[mod.addtion_index[j]][0] ]
            out.write(f"{items[0]:<6}{atom_nr:>5} {j:^5}{items[3]:>3} {items[4]}     {items[5]:>7.6g} {items[6]:>7.6g} {items[7]:>7.6g}  {items[8]}           {items[9]}\n")
            mod_index_sugar += 1
            atom_nr += 1
        else:
            items = [mod.type,atom_nr,j,mod.name,res["".join(add)][21:26],add_coords_base[mod.addtion_index[j]][0],add_coords_base[mod.addtion_index[j]][1],add_coords_base[mod.addtion_index[j]][2],"1.00  0.00",mod.additions_base[mod.addtion_index[j]][0] ]

            out.write(f"{items[0]:<6}{atom_nr:>5} {j:^5}{items[3]:>3} {items[4]}     {items[5]:>7.6g} {items[6]:>7.6g} {items[7]:>7.6g}  {items[8]}           {items[9]}\n")
            mod_index_base += 1
            atom_nr += 1


    return atom_nr

        #res_vec = np.array(res_compare_coord)-np.array(res_origin_coord)

        #R=align_vectors(add[4],res_vec)
        #translat = res_origin_coord-np.array(add[3])

        #point_coord = add[3]+R.dot(add[5])+translat
        #items = [mod.type,atom_nr,"".join(add[2]),mod.name,res["".join(add[0])][21:26],point_coord[0],point_coord[1],point_coord[2],res["".join(add[0])][56:66],add[2][0]]
        #spot3 = "".join(add[2])
        #out.write(f"{items[0]:<6}{atom_nr:>5}  {spot3:^4}{items[3]:>3} {items[4]}     {items[5]:>7.6g} {items[6]:>7.6g} {items[7]:>7.6g}  {items[8]}           {items[9]}\n")
        #atom_nr += 1


    return atom_nr