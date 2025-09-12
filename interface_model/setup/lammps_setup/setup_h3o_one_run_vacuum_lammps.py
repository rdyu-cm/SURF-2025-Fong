"""
setup up one run for the piston simulation (h3o vacuum)

-- writes the moltemplate files
-- writes the system.in file
-- copies the necessary files to the run directory

-- adapted to have electrode in the middle sandwiched by electrolyte
"""

import numpy as np
import os
from ase.io import read, write # type: ignore
from ase import Atoms
from ase.build import fcc111, surface, bulk

NA = 6.022e23 # Avogadro's constant
angstrom3_to_liter = 1e-27 # cubic angstrom to liter


def create_slab(slab_size):
    # generate electrode; use lattice constant in Lennard-Jones units
    a = 3.128
    slab_size = int(slab_size)
    bulk_tih2 = Atoms(
        symbols=['Ti', 'Ti', 'Ti', 'Ti', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
        positions=[
            # Ti atom at fcc position
            [1.0*a, 1.0*a, 1.0*a],
            [0.5*a, 0.5*a, 1.0*a],
            [0.5*a, 0.0, 0.5*a],
            [0.0, 0.5*a, 0.5*a],
            # H atoms according to materials project
            [0.25*a, 0.25*a, 0.25*a],
            [0.25*a, 0.25*a, 0.75*a],
            [0.25*a, 0.75*a, 0.75*a],
            [0.75*a, 0.75*a, 0.75*a],
            [0.75*a, 0.25*a, 0.25*a],
            [0.75*a, 0.75*a, 0.25*a],
            [0.25*a, 0.75*a, 0.25*a],
            [0.75*a, 0.25*a, 0.75*a]
        ],
        cell=[a, a, a],
        pbc=True
    )
    slab = surface(bulk_tih2, (1,1,1), layers = 6,vacuum = a, termination=1)
    slab = slab.repeat((slab_size,slab_size,1))
    #slab = hcp0001('Ti', size=(slab_size,slab_size,10), a=2.94, c=4.64, vacuum = 0.0)

    write('slab.xyz', slab)
    return slab

def create_slab2(slab_size):
    # generate electrode; use lattice constant in Lennard-Jones units
    '''
    this one is stoichiometric
    '''

    a = 3.128
    Ti_shift_atoms = Atoms(
    symbols=['Ti', 'Ti', 'Ti', 'Ti', 'H', 'H', 'H', 'H'],
    positions=[
        #Ti
        [1.0*a, 1.0*a, 1.0*a],
        [0.5*a, 0.5*a, 1.0*a],
        [0.5*a, 0.0, 0.5*a],
        [0.0, 0.5*a, 0.5*a],
        # H atoms according to materials project
        [0.25*a, 0.25*a, 0.25*a],
        #[0.25*a, 0.25*a, 0.75*a],
        [0.25*a, 0.75*a, 0.75*a],
        #[0.75*a, 0.75*a, 0.75*a],
        #[0.75*a, 0.25*a, 0.25*a],
        [0.75*a, 0.75*a, 0.25*a],
        #[0.25*a, 0.75*a, 0.25*a],
        [0.75*a, 0.25*a, 0.75*a]
    ],
    cell=[a, a, a],
    pbc=True
    )
    Ti_shift_slab = surface(Ti_shift_atoms, (1,1,1), layers=1, vacuum=0)

    H_shift_atoms = Atoms(
    symbols=['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
    positions=[
        #Ti
        # [1.0*a, 1.0*a, 1.0*a],
        # [0.5*a, 0.5*a, 1.0*a],
        # [0.5*a, 0.0, 0.5*a],
        # [0.0, 0.5*a, 0.5*a],
        # H atoms according to materials project
        [0.25*a, 0.25*a, 0.25*a],
        [0.25*a, 0.25*a, 0.75*a],
        [0.25*a, 0.75*a, 0.75*a],
        [0.75*a, 0.75*a, 0.75*a],
        [0.75*a, 0.25*a, 0.25*a],
        [0.75*a, 0.75*a, 0.25*a],
        [0.25*a, 0.75*a, 0.25*a],
        [0.75*a, 0.25*a, 0.75*a]
    ],
    cell=[a, a, a],
    pbc=True
    )
    H_shift_slab = surface(H_shift_atoms, (1,1,1), layers=1, vacuum=0)



    slab_size = int(slab_size)
    bulk_tih2 = Atoms(
        symbols=['Ti', 'Ti', 'Ti', 'Ti'],#, 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'],
        positions=[
            # Ti atom at fcc position
            [1.0*a, 1.0*a, 1.0*a],
            [0.5*a, 0.5*a, 1.0*a],
            [0.5*a, 0.0, 0.5*a],
            [0.0, 0.5*a, 0.5*a]
            # H atoms according to materials project
            # [0.25*a, 0.25*a, 0.25*a],
            # [0.25*a, 0.25*a, 0.75*a],
            # [0.25*a, 0.75*a, 0.75*a],
            # [0.75*a, 0.75*a, 0.75*a],
            # [0.75*a, 0.25*a, 0.25*a],
            # [0.75*a, 0.75*a, 0.25*a],
            # [0.25*a, 0.75*a, 0.25*a],
            # [0.75*a, 0.25*a, 0.75*a]
        ],
        cell=[a, a, a],
        pbc=True
        )
    
    slab = surface(bulk_tih2, (1,1,1), layers = 6,vacuum = a)
    slab = slab.repeat((slab_size,slab_size,1))
    # positions = slab.get_positions()
    # positions[:, 2] += Ti_shift_slab.cell[2][2]/2  # shift in z-direction
    # slab.set_positions(positions)

    top_atoms = Atoms(
        symbols=['H', 'H', 'H', 'H'],
        positions=[
            # H atoms according to materials project
            [0.25*a, 0.25*a, 0.25*a],
            #[0.25*a, 0.25*a, 0.75*a],
            [0.25*a, 0.75*a, 0.75*a],
            #[0.75*a, 0.75*a, 0.75*a],
            #[0.75*a, 0.25*a, 0.25*a],
            [0.75*a, 0.75*a, 0.25*a],
            #[0.25*a, 0.75*a, 0.25*a],
            [0.75*a, 0.25*a, 0.75*a]
        ],
        cell=[a, a, a],
        pbc=True
    )
    top_H = surface(top_atoms, (1,1,1), layers=6, vacuum=0)
    top_H = top_H.repeat((slab_size,slab_size,1))
    positions = top_H.get_positions()
    positions[:, 2] += a + 1*Ti_shift_slab.cell[2][2]  # shift in z-direction
    top_H.set_positions(positions)

    bottom_atoms = Atoms(
        symbols=['H', 'H', 'H', 'H'],
        positions=[
            # H atoms according to materials project
            [0.25*a, 0.25*a, 0.25*a],
            #[0.25*a, 0.25*a, 0.75*a],
            [0.25*a, 0.75*a, 0.75*a],
            #[0.75*a, 0.75*a, 0.75*a],
            #[0.75*a, 0.25*a, 0.25*a],
            [0.75*a, 0.75*a, 0.25*a],
            #[0.25*a, 0.75*a, 0.25*a],
            [0.75*a, 0.25*a, 0.75*a]
        ],
        cell=[a, a, a],
        pbc=True
    )
    bottom_H = surface(bottom_atoms, (1,1,1), layers=6, vacuum=0)
    bottom_H = bottom_H.repeat((slab_size,slab_size,1))
    positions = bottom_H.get_positions()
    positions[:, 2] += a - Ti_shift_slab.cell[2][2]# shift in z-direction
    bottom_H.set_positions(positions)
    
    slab = slab + top_H + bottom_H
    write('slab.xyz', slab)


    return slab


def run_packmol_elyte(Lx, Ly, Lz, molecule_names, n_molecules, cation_name):   
    output_file = "packed_box_elyte.xyz"
    packmol_in = "seed -1\n\n\
    tolerance 2.0 \n\n\
    filetype xyz \n\n\
    output " + output_file + "\n\n"
    for molecule_id in range(len(molecule_names)):
        packmol_in += "structure " + molecule_names[molecule_id] + ".xyz\n\
        number " + str(int(n_molecules[molecule_id])) + "\n\
        inside box 0. 0. 0. " + str(0.9*Lx) + " " + str(0.9*Ly) + " " + str(0.9*Lz) + "\n\
        end structure \n\n"

    print(packmol_in)

    file = open("packmol_in.inp","w")
    file.write(packmol_in)
    file.close()
    
    os.system("packmol < packmol_in.inp")

def run_combined_packmol(slab_size, elyte_height, slab_height, Lz, Lx, Ly, sigmaM):
    # write packmol file to combine electrodes and electrolyte
    elyte_height_top = elyte_height/2 + slab_height
    #elyte_height_bot = Lz/2 - (elyte_height/2 + slab_height/2)
    with open("packmol_template_vacuum.inp","rt") as fin:
            with open(f'packmol_tower.inp','wt') as fout:
                for line in fin:
                    fout.write(line.replace('{ELYTE_HEIGHT_TOP}',str(elyte_height_top)).replace(
                        '{ELECTRODE_HEIGHT}',str(slab_height/2)).replace(
                        '{SLAB_SIZE}',str(slab_size/2)))

    os.system('packmol < packmol_tower.inp')
    atoms = read('packed_box_tower.xyz')
    atoms.set_cell([Lx,Ly,Lz])
    atoms.set_pbc([True, True, True])

def create_lt_files(ion_charge,sigma_cation, eps_cation, cation_mass, epsM, sigmaM):
    # cation .lt file --> define the cation properties and force field
    with open("cation-template.lt","rt") as fin:
        with open(f'cation.lt','wt') as fout:
            for line in fin:
                fout.write(line.replace('{CATION_MASS}',str(cation_mass)).replace(
                    '{CATION_EPS}',str(eps_cation)).replace(
                    '{CATION_SIGMA}',str(sigma_cation)).replace(
                    '{ION_CHARGE}',str(ion_charge)))

    # with open("anion-template.lt","rt") as fin:
    #     with open(f'anion.lt','wt') as fout:
    #         for line in fin:
    #             fout.write(line.replace('{ANION_MASS}',str(anion_mass)).replace(
    #                 '{ANION_EPS}',str(eps_anion)).replace(
    #                 '{ANION_SIGMA}',str(sigma_anion)).replace(
    #                 '{ION_CHARGE}',str(ion_charge)))

    # metal .lt file --> define the force field of the metal
    with open("metal_forceField-template.lt","rt") as fin:
        with open(f'metal_forceField.lt','wt') as fout:
            for line in fin:
                fout.write(line.replace('{METAL_EPS}',str(epsM)).replace(
                    '{METAL_SIGMA}',str(sigmaM)))

def create_system_lt(molecule_names, n_molecules, Lx, Ly, Lz, slab):
    #modified for the sandwich to ensure atom ordering is correct
    system_lt = ""
    for molecule_id in range(len(molecule_names)):
        system_lt += "import \"" + molecule_names[molecule_id] + ".lt\" \n"
    for molecule_id in range(len(molecule_names)):
        system_lt += (molecule_names[molecule_id] + f"_molecule{molecule_id} = new " + molecule_names[molecule_id] + 
                     "["+ str(n_molecules[molecule_id]) + "]\n")
    system_lt += "write_once(\"Data Boundary\") {\n\
      0 " + str(Lx) + " xlo xhi \n\
      0 " + str(Ly) + " ylo yhi \n\
      0 " + str(Lz) + " zlo zhi \n\
    } \n\n"

    print(system_lt)
    file = open("system.lt","w")
    file.write(system_lt)
    file.close()

def run_moltemplate():
    os.system('moltemplate.sh -xyz packed_box_tower.xyz system.lt ')
    # os.system(f'cp system.in.charges {run_name}')

def create_system_in(cation,coul_cutoff, lj_cutoff, temp, timestep, n_npt_steps):
    with open("system-template-npt-mace.in","rt") as fin:
        with open(f'system.in','wt') as fout:
            for line in fin:
                fout.write(line.replace('{CATION}',str(cation)).replace('{TIMESTEP}',str(timestep)).replace(
                    '{LJ_CUTOFF}',str(lj_cutoff)).replace(
                    '{COUL_CUTOFF}',str(coul_cutoff)).replace(
                    '{TEMP}',str(temp)).replace(
                    '{N_NPT_STEPS}',str(n_npt_steps)))
                
def create_system_nvt_in(cation,coul_cutoff, lj_cutoff, temp, timestep, n_nvt_steps):
    with open("system-template-h3o-nvt-mace.in","rt") as fin:
        with open(f'system-nvt.in','wt') as fout:
            for line in fin:
                fout.write(line.replace('{CATION}',str(cation)).replace('{TIMESTEP}',str(timestep)).replace(
                    '{LJ_CUTOFF}',str(lj_cutoff)).replace(
                    '{COUL_CUTOFF}',str(coul_cutoff)).replace(
                    '{TEMP}',str(temp)).replace(
                    '{N_NVT_STEPS}',str(n_nvt_steps)))

def run_setup(run_name, n_molecules, molecule_names, 
    slab_size, epsM, sigmaM, 
    sigma_cation, eps_cation, cation_mass, cation_name, ion_charge,
    coul_cutoff, lj_cutoff, temp, timestep, n_steps, n_steps_equil):

    os.makedirs(run_name, exist_ok=True)

    slab = create_slab2(slab_size)
    Lx = slab.cell[0][0] 
    Ly = slab.cell[1][1]

    # choose distance between slabs
    #update densities based on cation pair - okay to do this, we're underpacking anyway
    density_guess = 1.04 # g/cm^3
    approx_box_size = ((n_molecules[0])/6.022e23*18.02/density_guess)/Lx/Ly/(1e-8)**3
    elyte_height = approx_box_size # under-pack the box for equilibration

    # create_ion_xyz()
    run_packmol_elyte(Lx, Ly, elyte_height, molecule_names, n_molecules, cation_name)

    #gap_between_electrodes = 25 # Angstrom
    Lz = slab.cell[2][2] + elyte_height + 12 - 3.128 #+ gap_between_electrodes # total z length of simulation box

    run_combined_packmol(slab.cell[0][0], elyte_height, slab.cell[2][2], Lz, Lx, Ly, sigmaM)

    molecule_names = molecule_names + ['metal_Ti','metal_H'] #adding two because of the two electrolytes
    n_molecules = list(n_molecules) + [int(slab.get_number_of_atoms()/3),int(2*slab.get_number_of_atoms()/3)]

    #everything below is for lammps, not ASE
    create_lt_files(ion_charge,sigma_cation, eps_cation, cation_mass, epsM, sigmaM)
    create_system_lt(molecule_names, n_molecules, Lx, Ly, Lz, slab)
    create_system_in(cation_name,coul_cutoff, lj_cutoff, temp, timestep, n_steps_equil)
    create_system_nvt_in(cation_name,coul_cutoff, lj_cutoff, temp, timestep, n_steps)
    run_moltemplate()
    os.system(f'cp system.data {run_name}')
    os.system(f'cp system.in.settings {run_name}')
    os.system(f'cp system.in {run_name}')
    os.system(f'cp system-nvt.in {run_name}')
    os.system(f'cp packed_box_tower.xyz {run_name}/init.xyz')