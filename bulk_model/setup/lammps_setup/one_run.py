'''
This script is used to set up a single system for lammps MD.
This is purely used for functions and isn't run itself.
'''

import numpy as np # type: ignore
import os

NA = 6.022e23 # Avogadro's constant
angstrom3_to_liter = 1e-27 # cubic angstrom to liter

def create_ion_xyz(): 
    """ create the cation.xyz file """
    one_atom_xyz = "1\n\n\
    C          0.00000        0.00000        0.00000\n"
    file = open("cation.xyz","w")
    file.write(one_atom_xyz)
    file.close()

def run_packmol(box_size, molecule_names, n_molecules):
    """ 
    runs packmol to pack the electrolyte 
    seed -1 sets a random seed, generating a new starting config each time its run 
    """   

    output_file = "packed_box.xyz"

    packmol_in = "seed -1\n\n\
    tolerance 2.0 \n\n\
    filetype xyz \n\n\
    output " + output_file + "\n\n"

    for molecule_id in range(len(molecule_names)):
        packmol_in += "structure " + molecule_names[molecule_id] + ".xyz\n\
        number " + str(n_molecules[molecule_id]) + "\n\
        inside box 0. 0. 0. " + str(0.9*box_size) + " " + str(0.9*box_size) + " " + str(0.9*box_size) + "\n\
        end structure \n\n"

    print(packmol_in)

    file = open("packmol_in.inp","w")
    file.write(packmol_in)
    file.close()
    
    os.system("packmol < packmol_in.inp") # runs packmol from the generated input script

def create_lt_files(ion_charge,sigma_cation, eps_cation, cation_mass):
    """ 
    creates the .lt file for the cation 
    since we change the cation identity as an independent variable, we generate the file in out script 
        rather than hard coding it up front
    """

    with open("cation-template.lt","rt") as fin:
        with open(f'cation.lt','wt') as fout:
            for line in fin:
                fout.write(line.replace('{CATION_MASS}',str(cation_mass)).replace(
                    '{CATION_EPS}',str(eps_cation)).replace(
                    '{CATION_SIGMA}',str(sigma_cation)).replace(
                    '{ION_CHARGE}',str(ion_charge)))

def create_system_lt(molecule_names, n_molecules, box_size):
    """ 
        generates the system.lt file
        feeds moltemplate the .lt files for each molecules type and tells it how many molecules we want
    """

    system_lt = ""
    for molecule_id in range(len(molecule_names)):
        system_lt += "import \"" + molecule_names[molecule_id] + ".lt\" \n"
    for molecule_id in range(len(molecule_names)):
        system_lt += (molecule_names[molecule_id] + "_molecule = new " + molecule_names[molecule_id] + 
                     "["+ str(n_molecules[molecule_id]) + "]\n")
    system_lt += "write_once(\"Data Boundary\") {\n\
      0 " + str(box_size) + " xlo xhi \n\
      0 " + str(box_size) + " ylo yhi \n\
      0 " + str(box_size) + " zlo zhi \n\
    } \n\n"

    print(system_lt)
    file = open("system.lt","w")
    file.write(system_lt)
    file.close()

def run_moltemplate():
    """ runs moltemplate """
    os.system('moltemplate.sh -xyz packed_box.xyz system.lt')

def create_system_in(coul_cutoff, lj_cutoff, temp, timestep, n_npt_steps):
    with open("system-template-npt.in","rt") as fin:
        with open(f'system.in','wt') as fout:
            for line in fin:
                fout.write(line.replace('{TIMESTEP}',str(timestep)).replace(
                    '{LJ_CUTOFF}',str(lj_cutoff)).replace(
                    '{COUL_CUTOFF}',str(coul_cutoff)).replace(
                    '{TEMP}',str(temp)).replace(
                    '{N_NPT_STEPS}',str(n_npt_steps)))
                
def create_system_nvt_in(coul_cutoff, lj_cutoff, temp, timestep, n_nvt_steps):
    with open("system-template-nvt.in","rt") as fin:
        with open(f'system-nvt.in','wt') as fout:
            for line in fin:
                fout.write(line.replace('{TIMESTEP}',str(timestep)).replace(
                    '{LJ_CUTOFF}',str(lj_cutoff)).replace(
                    '{COUL_CUTOFF}',str(coul_cutoff)).replace(
                    '{TEMP}',str(temp)).replace(
                    '{N_NVT_STEPS}',str(n_nvt_steps)))


def run_setup(box_size, n_molecules, molecule_names, 
    sigma_cation, eps_cation, cation_mass, ion_charge,
    coul_cutoff, lj_cutoff, temp, timestep, n_npt_steps, n_nvt_steps):

    """ series of functions to set up ONE run """
    # 1. create the ion xyz file
    create_ion_xyz()

    # 2. run packmol to create a random packed box --> generates the initial configuration
    run_packmol(box_size, molecule_names, n_molecules)

    # 3. create .lt file for the cation
    create_lt_files(ion_charge,sigma_cation, eps_cation, cation_mass)

    # 4. create system .lt file to be fed to moltemplate 
    create_system_lt(molecule_names, n_molecules, box_size)

    # 5. run moltemplate --> generates the topology file for lammps
    # the topology file defines both the structure AND the interactions
    run_moltemplate()

    # 6. create the system.in file for the NPT simulation  
    create_system_in(coul_cutoff, lj_cutoff, temp, timestep, n_npt_steps)

    # 7. create the system.in file for the NVT simulation 
    create_system_nvt_in(coul_cutoff, lj_cutoff, temp, timestep, n_nvt_steps)