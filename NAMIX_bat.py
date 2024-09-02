from NAMIX import NAMIX,rev_namix
import getopt,sys,os
from define_mods import *
import subprocess

inputs = [None,None,A_no_mod, C_no_mod, G_no_mod, T_no_mod, U_no_mod, "",False,False]
if __name__ == "__main__":
    try:
        opts,args =getopt.getopt(sys.argv[1:], "f:oa:c:u:g:t:r:p:hb:v", ["help""file=","overwrite","amod=","cmod=","umod=","tmod=","restrin=","prefix=","blueprint","min"])

    except getopt.GetoptError:
            print("""
                  
                  Input file can be give as the last argument of the call or as -f or --file
                  
                  -o or --overwrite: for overwrite folder content with same name
                  -r or --restrin: .txt with dot bracket format or .pb from chimira to gennerete basepair restrains for use in phenix generets .eff file
                  -p or --prefix: for give file prefixes and folder suffix
                  -b or --blueprint: make file for restrin based on ROAD blueprint


                  -a or --amod: set modifcation for adenine
                  -c or --cmod: set modifcation for cytosine
                  -g or --gmod: set modifcation for guanine
                  -u or --umod: set modifcation for uracil
                  -t or --tmod: set modifcation for thymine
                  
                   """)
            sys.exit()

    try:
        if not opts:
            if args[0] not in ["help","h","-help"]:
                file = args[0][2:-4]
                inputs[0] = file
            else:
                opts.append(("-h",""))

    except:
        print("\nNo file given\n")
        raise
    


    for i in opts:
        if i[0] == "-h" or "--help" == i[0]:
            print("""
                  
                  Input file can be give as the last argument of the call or as -f or --file
                  
                  -o or --overwrite: for overwrite folder content with same name
                  -r or --restrin: .txt with dot bracket format or .pb from chimira to gennerete basepair restrains for use in phenix generets .eff file
                  -p or --prefix: for give file prefixes and folder suffix
                  -b or --blueprint: make file for restrin based on ROAD blueprint 
                  -v: return mod nuc stucture to RNA

                
                  Mods are in the format "-[base to replace see below] [modifcation_3_lettercode],resid,resid...."
                  if no resid given all bases of the type will be replaced 
                  
                  -a or --amod: set modifcation for adenine
                  -c or --cmod: set modifcation for cytosine
                  -g or --gmod: set modifcation for guanine
                  -u or --umod: set modifcation for uracil
                  -t or --tmod: set modifcation for thymine


                  --min: Make NAMiX output only the modded pdb file and restrint if given
                  
                   """)
            
            print("\t\t  3 letter code for modificationsn\n")
            for key in mod_dict:
                print(f"\t\t  {mod_dict[key]}")
            print("\n")


            sys.exit()
        if i[0] == "-f" or i[0] == "--file":
            file = i[1].strip(".\\")
            file = file.removesuffix(".pdb")
            inputs[0] = file
        if i[0] == "-o" or i[0] == "--overwrite":
            inputs[8] = True
        if i[0] == "-p" or i[0] == "--prefix":
            inputs[7] = i[1]    
        if i[0] == "-r" or i[0] == "--restrin":
            file = i[1].strip(".\\")
            inputs[1] = file
        if i[0] =="-b" or i[0] == "--blueprint":
            pipe = subprocess.run(["perl", f"{os.path.dirname(__file__)}/trace_pattern.pl", f"{os.getcwd()}/{i[1]}"],capture_output=True)
            with open(f"bp_restrint{inputs[0]}.txt","w") as out:
                file = str(pipe.stdout).split("\\n")
                print(file)
                out.write(file[0][2:]+"\n")
                out.write(file[1]+"\n")
                out.write(file[2]+"\n")
            inputs[1] = f"bp_restrint{inputs[0]}.txt"
        
        
        if i[0] == "-u" or i[0] == "--umod":
            split_list = i[1].split(",")
            mod = mod_dict[split_list[0]]
            for base in range(1,len(split_list)):
                mod.mods.append(int(split_list[base]))
            inputs[6] = mod
        if i[0] == "-t" or i[0] == "--tmod":
            split_list = i[1].split(",")
            mod = mod_dict[split_list[0]]
            for base in range(1,len(split_list)):
                mod.mods.append(int(split_list[base]))
            inputs[5] = mod
        if i[0] == "-g" or i[0] == "--gmod":
            split_list = i[1].split(",")
            mod = mod_dict[split_list[0]]
            for base in range(1,len(split_list)):
                mod.mods.append(int(split_list[base]))
            inputs[4] = mod
        if i[0] == "-c" or i[0] == "--cmod":
            split_list = i[1].split(",")
            mod = mod_dict[split_list[0]]
            for base in range(1,len(split_list)):
                mod.mods.append(int(split_list[base]))
            inputs[3] = mod
        if i[0] == "-a" or i[0] == "--amod":
            split_list = i[1].split(",")
            mod = mod_dict[split_list[0]]
            for base in range(1,len(split_list)):
                mod.mods.append(int(split_list[base]))
            inputs[2] = mod
        
        if i[0] == "--min":
            inputs[9] = True

        
        if i[0] == "-v":
            rev_namix(inputs[0])
            sys.exit()
        

    NAMIX(*inputs)
    print(f"assembled folder for {inputs[0]}")
