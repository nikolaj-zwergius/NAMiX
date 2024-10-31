import re
from io import TextIOWrapper
from shutil import copy

def atom_selction_text(out:TextIOWrapper,key1, key2, key3, key4, dist,sigma, slack=0, change = None ,chain1 = None, chain2 = None): # code for formating the .eff file
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

def restrint_from_pb(file,dir_path="",minimum = False): #restrints based on chimira .pb file

    with open(file) as f, open(f"{dir_path}/{file[:-4]}_pb.eff","w") as out:
        restrints = []
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
                restrints.append(rest)
            else:
                re_split = re.split("\s",line)
                rest = [re_split[0].split(":")[0][1:],re_split[0].split("@")[0][len(re_split[0].split(":")[0])+1:],re_split[0].split("@")[1]]
                rest2 = [re_split[1].split(":")[0][1:],re_split[1].split("@")[0][len(re_split[0].split(":")[0])+1:],re_split[1].split("@")[1]]
                rest.extend(rest2)
                restrints.append(rest)
            
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
        for res in restrints:
            # remove backbone interactions
            if res[2] in blocked or res[5] in blocked:
                continue
            atom_selction_text(out,res[1],res[2],res[4],res[5],3.4,0.075,0.075,chain1=res[0],chain2=res[3]) # defines all h bonds
            atom_selction_text(out,res[1],res[2],res[1],f"N{res[2][1]}",1,0.01,0.01,chain1=res[0],chain2=res[0]) #fixes H breakign from modnuc if bonded

        for chain in chains.keys():
            for i in range(1,chains[chain]):
                    atom_selction_text(out,i,"O3'",i+1,"P",1.6,0.01,0.01,True,chain1=chain,chain2=chain) # fixes backbone from breakting
        for pair in pairs:
            parallelity(out,"N",pair[0],"N",pair[1],0,0.0335)
        out.write("  }\n}")
    if dir_path and not minimum:copy(file,dir_path)


def restrint_from_road(ssfile,pos=1,dir_path ="",minimum = False): #restrints based on ROAD dot backet file target.txt from trace patteren
    #dict of basepair interactions
    pos = int(pos)
    basepairs = {"A":{"U":(("N6","O4"),("N1","N3"))},
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
    with open(ssfile) as f, open(f"{dir_path}/{ssfile[:-4]}_db.eff","w") as out:
        # list for the differet backets
        backet1 = []
        backet2 = []
        backet3 = []
        backet4 = []

        pairs = [] # list of all basepairs
        
        #read in file
        lines = []
        for line in f:
            lines.append(line)
        dot_bac = lines[1]
        seq  =lines[2]

        # read Dot backet notation
        index = pos
        for char in dot_bac:
            if char == "(":
                backet1.append(index)
            elif char == ")":
                pairs.append((backet1.pop(),index))
            elif char == "[":
                backet2.append(index)
            elif char == "]":
                pairs.append((backet2.pop(),index))
            elif char == "{":
                backet3.append(index)
            elif char == "}":
                pairs.append((backet3.pop(),index))
            elif char == "<":
                backet4.append(index)
            elif char == ">":
                pairs.append((backet4.pop(),index))
            index += 1
        
        # writiing the file
        out.write("geometry_restraints {\n  edits {\n")
        for pair in pairs:
            base1 = seq[pair[0]-pos]
            base2 = seq[pair[1]-pos]
            parallelity(out,base1,pair[0]+1,base2,pair[1]+1,0,0.027)
          
            for interaction in basepairs[base1][base2]: #write out each hbond
                atom_selction_text(out,pair[0]+1,interaction[0],pair[1]+1,interaction[1],3.4,0.05)
            #
            #for fix in hydrogen_fix[base1]: # fix h breaks for base 1
            #    atom_selction_text(out,pair[0]+1,fix[0],pair[0]+1,fix[1],1,0.01,0.01)

            #for fix in hydrogen_fix[base2]: # fix h breaks for base 2
            #    atom_selction_text(out,pair[1]+1,fix[0],pair[1]+1,fix[1],1,0.01,0.01)
        
        for i in range(pos,pos+len(seq)-1): # fix backbone breaks
            atom_selction_text(out,i,"O3'",i+1,"P",1.6,0.03,change=True)
        
        out.write("  }\n}")
    if dir_path and not minimum:copy(ssfile,dir_path)
