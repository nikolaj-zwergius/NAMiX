import re
from io import TextIOWrapper
from shutil import copy
from os import remove

def atom_selection_text(out:TextIOWrapper,key1, key2, key3, key4, dist,sigma, slack=0, change = None ,chain1 = None, chain2 = None): # code for formatting the .eff file
    chain1_txt = ""
    chain2_txt = ""
    if chain1:
         chain1_txt = f"chain {chain1} and"
    if chain2:
         chain2_txt = f"chain {chain2} and"
    out.write("    bond {\n")
    if change: out.write("      action = change\n")
    out.write(f"      atom_selection_1 = {chain1_txt} resid {key1} and name {key2}\n")
    out.write(f"      atom_selection_2 = {chain2_txt} resid {key3} and name {key4}\n")
    out.write(f"      distance_ideal = {dist}\n      sigma = {sigma}\n")
    out.write("    }\n")
                    
def parallelity(out:TextIOWrapper,base1,resid1,base2,resid2, angel,sigma,chain1 = None, chain2 = None):
    base_string = {"A":"name C5 or name C6 or name N1 or name C2 or name N3 or name C4",
                   "G":"name C5 or name C6 or name N1 or name C2 or name N3 or name C4",
                   "C":"name N1 or name C2 or name N3 or name C4 or name C5 or name C6",
                   "U":"name N1 or name C2 or name N3 or name C4 or name C5 or name C6",
                   "R":"name C5 or name C6 or name N1 or name C2 or name N3 or name C4"}
    chain1_txt = ""
    chain2_txt = ""
    if chain1:
         chain1_txt = f"chain {chain1} and "
    if chain2:
         chain2_txt = f"chain {chain2} and "
                             
    out.write("    parallelity {\n")
    out.write(f"      atom_selection_1 = {chain1_txt}resid {resid1} and (name C5 or name C6 or name N1 or name C2 or name N3 or name C4)\n")
    out.write(f"      atom_selection_2 = {chain2_txt}resid {resid2} and (name C5 or name C6 or name N1 or name C2 or name N3 or name C4)\n")
    out.write(f"      sigma = {sigma}\n      target_angle_deg = {angel}\n")
    out.write("    }\n")

def stack(dir_path,trace_file,kissing_loops,chain1=None,chain2=None):
    with open(f"{dir_path}/{trace_file[:-4]}_ss.eff","w") as stack_file:
        stack_file.write("pdb_interpretation {\n" +" "*2+"secondary_structure {\n"+ " "*4+"nucleic_acid {\n")
        chain1_txt = ""
        chain2_txt = ""
        if chain1:
            chain1_txt = f"chain {chain1} and"
        if chain2:
            chain2_txt = f"chain {chain2} and"
        for kiss in kissing_loops:
            stack_file.write(" "*6+"stacking_pair {\n")
            stack_file.write(" "*8+f"base1 = {chain1_txt} resid {kiss[0][1]}\n")
            stack_file.write(" "*8+f"base2 = {chain2_txt} resid {kiss[1][1]}\n"+" "*6+"}\n")
            stack_file.write(" "*6+"stacking_pair {\n")
            stack_file.write(" "*8+f"base1 = {chain1_txt} resid {kiss[0][0]+2}\n")
            stack_file.write(" "*8+f"base2 = {chain2_txt} resid {kiss[1][1]-1}\n"+" "*6+"}\n")
            stack_file.write(" "*6+"stacking_pair {\n")
            stack_file.write(" "*8+f"base1 = {chain1_txt} resid {kiss[1][0]+2}\n")
            stack_file.write(" "*8+f"base2 = {chain2_txt} resid {kiss[0][1]-1}\n"+" "*6+"}\n")
        stack_file.write(" "*4+"}\n"+" "*2+"}\n"+"}\n")

def restraints_from_pb(file,dir_path="",minimum = False): #restraints based on chimeraX .pb file

    with open(file) as f, open(f"{dir_path}/{file[:-4]}_pb.eff","w") as out:
        restraints = []
        chains = {}
        plan = {}
        pairs = []
        blocked = ["HO2'","O4'","HO'2","O2'","HO3'","H'O3","O5'","OP1","OP2","O3'","N7"]
        Pu_or_py = {"N6":"R","N1":"R","N2":"R","O6":"R","O2":"Y","N3":"Y","N4":"Y","O4":"Y"}
        # read .pb and split it into the set of donor/accepter
        for line in f:
            
            if line.startswith(";"):
                continue
            if line.startswith("#"):
                re_split = re.split("\s",line)
                rest = [re_split[0].split(":")[0][1:].split("/")[1],re_split[0].split("@")[0][len(re_split[0].split(":")[0])+1:],re_split[0].split("@")[1]]
                rest2 = [re_split[1].split(":")[0][1:].split("/")[1],re_split[1].split("@")[0][len(re_split[0].split(":")[0])+1:],re_split[1].split("@")[1]]
                rest.extend(rest2)
                restraints.append(rest)
            else:
                re_split = re.split("\s",line)
                rest = [re_split[0].split(":")[0][1:],re_split[0].split("@")[0][len(re_split[0].split(":")[0])+1:],re_split[0].split("@")[1]]
                rest2 = [re_split[1].split(":")[0][1:],re_split[1].split("@")[0][len(re_split[0].split(":")[0])+1:],re_split[1].split("@")[1]]
                rest.extend(rest2)
                restraints.append(rest)
            
            if rest[2] not in blocked and rest2[2] not in blocked:
                mini = min(int(rest[1]),int(rest2[1]))
                maxi = max(int(rest[1]),int(rest2[1]))
                pairs.append((mini,maxi))
            #defines length of chains
            if rest[0] in chains:
                chains[rest[0]] = max(chains[rest[0]],int(rest[1]))
            else:
                chains[rest[0]] = int(rest[1])
            if rest[3] in chains:
                chains[rest[3]] = max(chains[rest[3]],int(rest[4]))
            else:
                chains[rest[3]] = int(rest[4])

        # write out file
        out.write("geometry_restraints {\n  edits {\n")
        for res in restraints:
            # remove backbone interactions
            if res[2] in blocked or res[5] in blocked:
                continue
            atom_selection_text(out,res[1],res[2],res[4],res[5],3.4,0.075,0.075,chain1=res[0],chain2=res[3]) # defines all h bonds
            atom_selection_text(out,res[1],res[2],res[1],f"N{res[2][1]}",1,0.01,0.01,chain1=res[0],chain2=res[0]) #fixes H breaking from modnuc if bonded

        for chain in chains.keys():
            for i in range(1,chains[chain]):
                    atom_selection_text(out,i,"O3'",i+1,"P",1.6,0.01,0.01,True,chain1=chain,chain2=chain) # fixes backbone from breaking
        for pair in pairs:
            parallelity(out,"N",pair[0],"N",pair[1],0,0.0335)
        out.write("  }\n}")
    if dir_path and not minimum:copy(file,dir_path)


def restraints_from_road(trace_file,pos=1,dir_path ="",minimum = False,bp_file="",chain=""): #restraints based on ROAD dot bracket file target.txt from trace pattern
    #dict of base-pair interactions
    pos = int(pos)
    basepairs = {"A":{"U":(("N6","O4"),("N1","N3")),"A":(("N1","N6"),("N6","N1"))},
                 "C":{"G":(("N4","O6"),("N3","N1"),("O2","N2"))},
                 "G":{"C":(("O6","N4"),("N1","N3"),("N2","O2")),"U":(("O6","N3"),("N1","O2"))},
                 "U":{"A":(("O4","N6"),("N3","N1")),"G":(("N3","O6"),("O2","N1"))}
                 }
    #dict of with H donors are connected to what Hs
    hydrogen_fix = {"A":[("H61","N6")],
                    "C":[("H41","N4")],
                    "G":[("H1","N1"),("H21","N2")],
                    "U":[("H3","N3")]
                    }
    with open(trace_file) as f, open(f"{dir_path}/{trace_file[:-4]}_bp.eff","w") as out:
        # list for the differet backets
        bracket1 = []
        bracket2 = []
        bracket3 = []
        bracket4 = []

        pairs = [] # list of all basepairs
        
        #read in file
        lines = []
        for line in f:
            lines.append(line)
        dot_bac = lines[1]
        seq  =lines[2]

        # read Dot backet notation
        index = pos
        kissing_loop_id = []
        for char in dot_bac:
            if char == "(":
                bracket1.append(index)
            elif char == ")":
                pairs.append((bracket1.pop(),index))
            elif char == "[":
                bracket2.append(index)
            elif char == "]":
                pairs.append((bracket2.pop(),index))
                if dot_bac[pairs[-1][0]-8] == ".":
                    kissing_loop_id.append(((pairs[-1][0]-7,pairs[-1][0]-6,pairs[-1][0]+1),(pairs[-1][1]-2,pairs[-1][1]-1,pairs[-1][1]+6)))
            elif char == "{":
                bracket3.append(index)
            elif char == "}":
                pairs.append((bracket3.pop(),index))
            elif char == "<":
                bracket4.append(index)
            elif char == ">":
                pairs.append((bracket4.pop(),index))
            index += 1


        # writiing the file
        out.write("geometry_restraints {\n  edits {\n")
        
        for kiss in kissing_loop_id:
            parallelity(out,"A",kiss[0][0],"A",kiss[0][2],0,0.27,chain1=chain,chain2=chain)
            parallelity(out,"A",kiss[1][0],"A",kiss[1][2],0,0.27,chain1=chain,chain2=chain)
            atom_selection_text(out,kiss[0][0],basepairs["A"]["A"][0][0],kiss[0][2],basepairs["A"]["A"][0][1],3,0.02,chain1=chain,chain2=chain)
            atom_selection_text(out,kiss[0][0],basepairs["A"]["A"][1][0],kiss[0][2],basepairs["A"]["A"][1][1],3,0.02,chain1=chain,chain2=chain)
            atom_selection_text(out,kiss[1][0],basepairs["A"]["A"][0][0],kiss[1][2],basepairs["A"]["A"][0][1],3,0.02,chain1=chain,chain2=chain)
            atom_selection_text(out,kiss[1][0],basepairs["A"]["A"][1][0],kiss[1][2],basepairs["A"]["A"][1][1],3,0.02,chain1=chain,chain2=chain)
        
        for pair in pairs:
            base1 = seq[pair[0]-pos]
            base2 = seq[pair[1]-pos]
            parallelity(out,base1,pair[0],base2,pair[1],0,0.027,chain1=chain,chain2=chain)
          
            for interaction in basepairs[base1][base2]: #write out each h-bond
                atom_selection_text(out,pair[0],interaction[0],pair[1],interaction[1],3,0.02,chain1=chain,chain2=chain)
            
        for i in range(pos,pos+len(seq)-2): # fix backbone breaks
            atom_selection_text(out,i,"O3'",i+1,"P",1.61,0.03,change=True,chain1=chain,chain2=chain)
        
        out.write("  }\n}")
    
    if kissing_loop_id:
        stack(dir_path,trace_file,kissing_loop_id,chain,chain)

    if dir_path and not minimum:
        copy(trace_file,dir_path)
        copy(bp_file,dir_path)
        remove(trace_file)

