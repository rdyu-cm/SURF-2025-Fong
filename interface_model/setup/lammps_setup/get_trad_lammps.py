'''
This script is primarily used for the remove_bonds_manual function
It is necessary to remove the bond information for the MLP to run in lammps
'''

import MDAnalysis as mda
import os
import random
import glob
from ase import Atoms
from ase.io import read, write
import numpy as np

def remove_bonds_manual(input_file, output_file):
    """
    Remove bond information by parsing and rewriting the LAMMPS data file
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    skip_section = False
    
    for line in lines:
        line_lower = line.lower().strip()
        original_line = line.strip()
        
        # Skip empty lines when in skip mode
        if skip_section and (not original_line or original_line.isspace()):
            continue
            
        # Header modifications - change counts to 0
        if 'bonds' in line_lower and 'bond types' not in line_lower:
            # Line like "1234 bonds" -> "0 bonds"
            parts = original_line.split()
            if len(parts) >= 2 and parts[1].lower() == 'bonds':
                new_lines.append("0 bonds\n")
                continue
                
        if 'bond types' in line_lower:
            # Line like "2 bond types" -> "0 bond types" 
            parts = original_line.split()
            if len(parts) >= 3 and parts[1].lower() == 'bond' and parts[2].lower() == 'types':
                new_lines.append("0 bond types\n")
                continue
        
        if 'angles' in line_lower and 'angle types' not in line_lower:
            parts = original_line.split()
            if len(parts) >= 2 and parts[1].lower() == 'angles':
                new_lines.append("0 angles\n")
                continue
                
        if 'angle types' in line_lower:
            parts = original_line.split()
            if len(parts) >= 3 and parts[1].lower() == 'angle' and parts[2].lower() == 'types':
                new_lines.append("0 angle types\n")
                continue
        
        if 'dihedrals' in line_lower and 'dihedral types' not in line_lower:
            parts = original_line.split()
            if len(parts) >= 2 and parts[1].lower() == 'dihedrals':
                new_lines.append("0 dihedrals\n")
                continue
                
        if 'dihedral types' in line_lower:
            parts = original_line.split()
            if len(parts) >= 3 and parts[1].lower() == 'dihedral' and parts[2].lower() == 'types':
                new_lines.append("0 dihedral types\n")
                continue
        
        # Detect start of sections to skip
        if (original_line.lower() == 'bonds' or 
            original_line.lower() == 'angles' or 
            original_line.lower() == 'dihedrals' or
            original_line.lower() == 'impropers'):
            skip_section = True
            continue
            
        # Detect start of sections to keep (ends skipping)
        if (original_line.lower() == 'atoms' or 
            original_line.lower() == 'masses' or 
            original_line.lower() == 'pair coeffs' or
            original_line.lower() == 'velocities'):
            skip_section = False
            new_lines.append(line)
            continue
        
        # Skip lines when in a bond/angle/dihedral section
        if skip_section:
            continue
            
        # Keep all other lines
        new_lines.append(line)
    
    # Write the new file
    with open(output_file, 'w') as f:
        f.writelines(new_lines)
    

def mda_to_ase_frame(universe, frame_idx, type_to_element):
    """Convert MDAnalysis frame to ASE Atoms object"""
    universe.trajectory[frame_idx]
    
    # Get positions and elements
    positions = universe.atoms.positions
    elements = [type_to_element.get(int(t), 'X') for t in universe.atoms.types]
    
    # Get cell - handle both orthogonal and triclinic
    box = universe.trajectory.ts.dimensions
    charges = universe.atoms.charges
    resids = universe.atoms.resids
    if np.allclose([box[3], box[4], box[5]], [90.0, 90.0, 90.0]):
        # Orthogonal cell
        cell = [box[0], box[1], box[2]]
    
    # Create ASE Atoms object with PBC
    atoms = Atoms(symbols=elements, positions=positions, cell=cell, pbc=True)
    atoms.set_array('charges', charges)
    atoms.set_array('residues', resids)
    return atoms

def write_positions():
    # Load trajectory
    path = "/home/rdyu/LAMMPS_Tutorial/LAMMPS_Tutorial_runs/"
    os.chdir(path + f'LAMMPS_Tutorial_1_rep0')
    # for file in glob.glob(f"raw_positions"):
    #     os.remove(file)
    u = mda.Universe("end_npt.data", "npt_unwrapped_0.dcd")
    #type_to_element = {1: 'O', 2: 'H', 3: 'Na', 4: 'N', 5: 'O'}

    # Assign elements based on types
    #atom_types = u.atoms.types
    #elements = [type_to_element.get(int(t), 'X') for t in atom_types]

    # Add elements as a new property
    #u.add_TopologyAttr('elements', elements)
    
    # Write specific frame
    frame = random.randint(500,999)
    atoms = mda_to_ase_frame(u, frame)#, type_to_element)
    path = "/home/rdyu/LAMMPS_MACE"
    os.chdir(path)
    write(f"system.data", atoms, atom_style = 'full', format='lammps-data')

def write_positions_ase():
    frame = random.randint(500,999)
    path = "/home/rdyu/LAMMPS_Tutorial/LAMMPS_Tutorial_runs/"
    os.chdir(path + f'LAMMPS_Tutorial_1_rep0')
    atoms = read('npt_unwrapped_0.dcd', index=frame, format='dcd')

    # You'll need the original topology information
    # Option A: If you have the original .data file
    ref_atoms = read('end_npt.data', format='lammps-data', style='full')
    topology_info = {
        'charges': ref_atoms.get_array('charges') if 'charges' in ref_atoms.arrays else None,
        'types': ref_atoms.get_array('numbers'),
        'bonds': ref_atoms.arrays.get('bonds', None),
        'angles': ref_atoms.arrays.get('angles', None),
        'dihedrals': ref_atoms.arrays.get('dihedrals', None)
    }

    # Apply topology to the extracted frame
    atoms.set_array('charges', topology_info['charges'])
    # Add other topology info as needed

    # Write as LAMMPS data file
    write('system.data', atoms, format='lammps-data', atom_style='full')

def write_positions_mdanalysis(repli):
    path = "/home/rdyu/LAMMPS_Tutorial/LAMMPS_Tutorial_runs/"
    os.chdir(path + f'LAMMPS_Tutorial_1_rep0')
    u = mda.Universe("end_npt.data", "npt_unwrapped_0.dcd")
    frame = int(repli)*10+800
    u.trajectory[frame]
    path = "/home/rdyu/LAMMPS_MACE"
    os.chdir(path)
    u.atoms.write('system.data', format='DATA')