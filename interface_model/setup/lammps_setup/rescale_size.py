"""
LAMMPS molecular system setup and box rescaling utility for aqueous electrolyte simulations.
This file is typically called by other files for its functions.

Typical workflow:
1. Extract box dimensions from XYZ comment lines or input parameters
2. Convert XYZ coordinates to LAMMPS data format via moltemplate
3. Rescale box dimensions based on density/volume corrections from NPT runs
4. Verify structural integrity with single-step LAMMPS calculation
5. Output rescaled XYZ files

Dependencies: LAMMPS executable, moltemplate, subprocess, tempfile
"""

import subprocess
import tempfile
import os
import glob

def create_system_lt(molecule_names, n_molecules, box):
    """ 
        generates the system.lt file
        feeds moltemplate the .lt files for each molecules type and tells it how many molecules we want
    """

    #have to change for the sandwich, look at interface for this
    system_lt = ""
    for molecule_id in range(len(molecule_names)):
        system_lt += "import \"" + molecule_names[molecule_id] + ".lt\" \n"
    for molecule_id in range(len(molecule_names)):
        system_lt += (molecule_names[molecule_id] + "_molecule = new " + molecule_names[molecule_id] + 
                     "["+ str(n_molecules[molecule_id]) + "]\n")
    system_lt += "write_once(\"Data Boundary\") {\n\
      " + str(box[0]) + " " + str(box[1]) + " xlo xhi \n\
      " + str(box[2]) + " " + str(box[3]) + " ylo yhi \n\
      " + str(box[4]) + " " + str(box[5]) + " zlo zhi \n\
    } \n\n"

    print(system_lt)
    file = open("system.lt","w")
    file.write(system_lt)
    file.close()

def run_moltemplate(input_file):
    """ runs moltemplate """
    os.system(f'moltemplate.sh -xyz {input_file} system.lt')

def get_box_from_xyz_comment(xyz_file):
    """Extract box dimensions from XYZ comment line if present"""
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
        
    if len(lines) < 2:
        return None
        
    comment = lines[1].strip()
    
    if "Lattice=" in comment:
        lattice_part = comment.split("Lattice=\"")[1].split()
        box_size = float(lattice_part[0])
        return (0, box_size, 0, box_size, 0, box_size)
    
    return None

def change_box_output_then_verify(lammps_exe, box, input_xyz, output_xyz, x_bounds, y_bounds, z_bounds,scale,ion):
    original_box = box
    molecule_names = ['h2o','cation','nitrate']
    create_system_lt(molecule_names,[64,1*ion,1*ion],original_box)
    run_moltemplate(input_xyz)
    print(original_box)
    script_content = f"""
# Setup
units real
dimension 3
boundary p p p
atom_style full

# Read input XYZ
# region box block {original_box[0]} {original_box[1]} {original_box[2]} {original_box[3]} {original_box[4]} {original_box[5]}
# create_box 1 box
# read_dump {input_xyz} 0 x y z box no format xyz
read_data 'system.data'

# Change the box
change_box all z final ${z_bounds[0]} ${z_bounds[1]}

# Create group for specific atom type
group frozen_type type 6 7

# Scale those atoms back to original positions
displace_atoms frozen_type scale 1.0 1.0 ${scale}

# Output the modified structure to XYZ file
dump output_xyz all xyz 1 {output_xyz}
dump_modify output_xyz element O H Na N O Ti H
run 0
undump output_xyz


# Set up force field for verification step
pair_style lj/cut/coul/long 10.0 10.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  opls 
neigh_modify every 1 delay 0 check no
pair_modify shift yes mix arithmetic
kspace_style pppm 1.0e-5
# import non-bonded (LJ) interaction parameters
include "system.in.settings"

# Set up thermodynamic output for verification
thermo_style custom step atoms vol density pe ke etotal press
thermo 1

print "=== VERIFICATION STEP ==="
timestep 1
run 1

"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False) as f:
        f.write(script_content)
        script_file = f.name
    
    try:
        result = subprocess.run([lammps_exe, '-in', script_file], 
                               capture_output=True, text=True)
        
        if result.returncode == 0:
            print("LAMMPS completed successfully!")
            print("Output:")
            print(result.stdout)
        else:
            print("LAMMPS failed:")
            print(result.stderr)
            
        return result
        
    finally:
        os.unlink(script_file)
