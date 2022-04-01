import subprocess
import sys
import os

def set_parameter(param_path, param_name, value):
    if not os.path.isfile(param_path):
        sys.exit("ERROR: %s is not a file"%param_path)
    if isinstance(value, int):
        param_string = "%s = %d"%(param_name, value)
    elif isinstance(value, str):
        param_string = "%s = %s"%(param_name, value)
    elif isinstance(value, float):
        param_string = "%s = %.8e"%(param_name, value)
    else:
        sys.exit("ERROR: Unrecongnised type %s supplied for %s"%(str(value.type), param_name))
    
    with open(param_path, "r") as f:
        if param_name in f.read():
            subprocess.run(["sed", "-i", "s/^%s .*$/%s/g"%(param_name,param_string), "simulation.par"])
            return
    with open(param_path, "a") as f:
            f.write(param_string + "\n")
