from NAMIX import NAMIX,rev_namix
import getopt,sys,os
from define_mods import *
import subprocess
from namix_driver import phenix,qrnas
from namix_configer import configer
inputs = [None,None,A_no_mod, C_no_mod, G_no_mod, T_no_mod, U_no_mod, "",False,False]
driver_tags = {}
help_mes =  """
                  
                  Example of use: namix -f [pdb file] -r [restrain file] -a [adenine modification] -c [cytosine modification,resid1,resid2...] -p [prefix]
                  
                  Input file can be give as the last argument of the call or as -f [filename] or --file [filename]
                  
                  
                  -r [file] or --restrain: .txt with dot bracket format(not implemented yet)
                    or .pb from chimera to generate basepair restrains for use in phenix generates .eff file
                  -b [file] or --blueprint: make file for restrains based on ROAD blueprint
                  
                  -p [prefix] or --prefix: for give file prefixes and folder suffix

                  -o or --overwrite: for overwrite folder content with same name 
                  -v: return mod nuc structure to RNA
                  -m or --min: Make NAMiX output only the modded pdb file and restraint if given

                  -q or --qrna [config]  runs QRNAS if possible if -q used default config will be used #not implemented yet
                  -x [map,res] or --phenix [map,res] runs Phenix.real_space_refine if possible #not implemented yet
                
                  --config [filename]: runs config file #not implemented yet
                  --chain [Chain ID]: adds chain id to .eff for multi chain models 

                  Mods are in the format "-[base to replace see below] [modification_3_letter_code],resid,resid...."
                  if no resid given all bases of the type will be replaced 
                  
                  
                  -a [modification,resids] or --amod [modification,resids]: set modification for adenine
                  -c [modification,resids] or --cmod [modification,resids]: set modification for cytosine
                  -g [modification,resids] or --gmod [modification,resids]: set modification for guanine
                  -u [modification,resids] or --umod [modification,resids]: set modification for uracil
                  -t [modification,resids] or --tmod [modification,resids]: set modification for thymine

                  
                   """



if __name__ == "__main__":
    try:
        opts,args =getopt.getopt(sys.argv[1:], "f:a:c:u:g:t:r:p:b:vx:q:hom", ["help""file=","overwrite","amod=","cmod=","umod=","tmod=","restrain=","prefix=","blueprint=","min","phenix=","qrna","config","chain="])

    except getopt.GetoptError:
            print(help_mes)
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
            print(help_mes)
            
            print("\t\t  3 letter codes for modification\n")
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
        if i[0] == "-r" or i[0] == "--restrain":
            file = i[1].strip(".\\")
            inputs[1] = [file]
        if i[0] =="-b" or i[0] == "--blueprint":
            try:
                split_list_bp = i[1].split(",")
                try:
                    id  = split_list_bp[1]
                except:
                    id = 1
                pipe = subprocess.run(["perl", f"{os.path.dirname(__file__)}/trace_pattern.pl", f"{os.getcwd()}/{split_list_bp[0]}"],capture_output=True)
                with open(f"bp_restrain{inputs[0]}.txt","w") as out:
                    file = str(pipe.stdout).split("\\n")
                    out.write(file[0][2:]+"\n")
                    out.write(file[1]+"\n")
                    out.write(file[2]+"\n")
                inputs[1] = (f"bp_restrain{inputs[0]}.txt",id,split_list_bp[0])
            except:
                print("\n\n could not run: perl properly not installed")
        
        if i[0] == "--chain":
            bp_list = list(inputs[1])
            bp_list.append(i[1])
            inputs[1]= tuple(bp_list)
        
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
        
        if i[0] == "--min" or i[0] == "-m":
            inputs[9] = True
        
        if i[0] == "--phenix" or i[0] == "-x":
            maps,res = i[1].split(",")
            driver_tags["phenix"] = (maps,res)
        if i[0] == "-q":
            driver_tags["qrna"] = i[1]
        if i[0] == "--qrna":
            driver_tags["qrna"] = "defualt_qRNA_config"
        if i[0] == "--config":
            configer(i[1])
            sys.exit()

        if i[0] == "-v":
            rev_namix(inputs[0])
            sys.exit()
        

    NAMIX(*inputs)
    print(f"convereted {inputs[0]}")
    
    if "qrna" in driver_tags:
        print("qrna done")
        qrnas(driver_tags["qrna"])
    if "phenix" in driver_tags:
        print("phenix done")
        phenix(driver_tags["phenix"][0],driver_tags["phenix"][0])
        
    
