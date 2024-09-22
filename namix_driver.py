import subprocess


def phenix(map,res):
    subprocess.run("phenix.real_space_refine")

namex_config = "path_to_config"
def qrnas(config):
    subprocess.run("QRNA")
    print("done")

