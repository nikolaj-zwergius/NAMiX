import re
from io import TextIOWrapper
from shutil import copy

def atom_selction_text(out:TextIOWrapper,key1, key2, key3, key4, dist,sigma, slack, change = None ,chain1 = None, chain2 = None): # code for formating the .eff file
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
                    out.write(f"      distance_ideal = {dist}\n      sigma = {sigma}\n      slack = {slack}\n")
                    out.write("    }\n")

def restrint_from_pb(file,dir_path="",minimum = False): #restrints based on chimira .pb file

    with open(file) as f, open(f"{dir_path}/{file[:-4]}_pb.eff","w") as out:
        restrints = []
        chains = {}
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
        blocked = ["HO2'","O4'","HO'2","O2'","HO3'","H'O3"]
        for res in restrints:
            # remove backbone interactions
            if res[2] in blocked or res[5] in blocked:
                continue
            atom_selction_text(out,res[1],res[2],res[4],res[5],3.4,0.075,0.075,chain1=res[0],chain2=res[3]) # defines all h bonds
            atom_selction_text(out,res[1],res[2],res[1],f"N{res[2][1]}",1,0.01,0.01,chain1=res[0],chain2=res[0]) #fixes H breakign from modnuc if bonded

        for chain in chains.keys():
            for i in range(1,chains[chain]):
                    atom_selction_text(out,i,"O3'",i+1,"P",1.6,0.01,0.01,True,chain1=chain,chain2=chain) # fixes backbone from breakting
                
        out.write("  }\n}")
    if dir_path and not minimum:copy(file,dir_path)


def restrint_from_ss(ssfile,dir_path ="",minimum = False): #restrints based on ROAD dot backet file target.txt from trace patteren
    #dict of basepair interactions
    basepairs = {"A":{"U":(("H61","O4"),("N1","H3"))},
                 "C":{"G":(("H41","O6"),("N3","H1"),("O2","H21"))},
                 "G":{"C":(("O6","H41"),("H1","N3"),("H21","O2")),"U":(("O6","H3"),("H1","O2"))},
                 "U":{"A":(("O4","H61"),("H3","N1")),"G":(("H3","O6"),("O2","H1"))}
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
        index = 0
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
            base1 = seq[pair[0]]
            base2 = seq[pair[1]]
          
            for interaction in basepairs[base1][base2]: #write out each hbond
                atom_selction_text(out,pair[0]+1,interaction[0],pair[1]+1,interaction[1],3.4,0.075,0.075)
            
            for fix in hydrogen_fix[base1]: # fix h breaks for base 1
                atom_selction_text(out,pair[0]+1,fix[0],pair[0]+1,fix[1],1,0.01,0.01)

            for fix in hydrogen_fix[base2]: # fix h breaks for base 2
                atom_selction_text(out,pair[1]+1,fix[0],pair[1]+1,fix[1],1,0.01,0.01)
        
        for i in range(1,len(seq)-1): # fix backbone breaks
            atom_selction_text(out,i,"O3'",i+1,"P",1.6,0.01,0.01,True)
        
        out.write("  }\n}")
    if dir_path and not minimum:copy(ssfile,dir_path)
