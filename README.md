# NAMiX
Nuclic acid modelling including xeno-nucleotide is a tool made to make the process of building and refiment of heavly modifed neclic acids stctures into strucule maps.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

Make sure that python is installed

Then copy the git repository
```bash
 git clone https://github.com/nikolaj-zwergius/NAMiX.git
```
After this add the NAMIX folder to your PATH

### Linux:

```bash
 export PATH="path_to_namix_folder:$PATH"
```

### Windows

### Mac

This software has not been tested on Mac OS but show be compatiable, but no guarantee is made 

### If using ROAD style blueprint
If you want to use ROAD style blueprints for the generation of base paring restraints. you need to have perl install in addtion to python, as the code for convereting the blueprint to dot-barcket format, is taken from the [ROAD](https://github.com/esa-lab/ROAD) repository.

## Usage
NAMIX can be called by just writting the following:

```bash
namix -f [filename] [options]
```
The -f [filename] is the only part that is mandatory

For a list of options see below
-r [file]: .txt with dot bracket format(not implemented yet) or .pb from chimira to gennerete basepair restrains for use in phenix generets .eff file. \
-b [file]: make file for restrints based on ROAD blueprint.
                  
-p [prefix]: for give file prefixes and folder suffix

-o: for overwrite folder content with same name.\
-v: return mod nuc stucture to RNA.\
-m: Make NAMiX output only the modded pdb file and restrint if given.\

Mods are in the format "-[base to replace see below] [modifcation_3_lettercode,resid,resid....]"\
if no res ID is given all bases of the type will be replaced                  
                  
 -a [modifcaiton,resids]: set modifcation for adenine\
 -c [modifcaiton,resids]: set modifcation for cytosine\
 -g [modifcaiton,resids]: set modifcation for guanine\
 -u [modifcaiton,resids]: set modifcation for uracil\
 -t [modifcaiton,resids]: set modifcation for thymine

## Contributing
1. Fork the repository.
2. Create a new branch: `git checkout -b feature-name`.
3. Make your changes.
4. Push your branch: `git push origin feature-name`.
5. Create a pull request.

## License
This project is licensed under the [GPL-3.0 license](LICENSE).
