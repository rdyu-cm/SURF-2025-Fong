'''
This script prints out energy and force RMSE between dft outputs that were calculated with different numbers of k points

input parameters:
- path: the location of your dft calculations folders
- n_atoms: the number of atoms in your system

output:
- prints: Energy RMSE, Force RMSE
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import ase.io

def access_dft_energies_forces(path,location, length):
    '''
    Reads and outputs forces and energies from a cp2k output file
    '''
    energies = []
    forces = []
    for i in range(*length):
        with open(f'{path}{location}/structure{i}1/ref_calc-1.ener','r') as file:
            lines = file.readlines()
        energies.append((float(lines[1].split()[4])))
        with open(f'{path}{location}/structure{i}1/forces-output.xyz','r') as file:
            lines = file.readlines()
        for j in range(len(lines)-2):
            forces.append(float(lines[j+2].split()[1]))
            forces.append(float(lines[j+2].split()[2]))
            forces.append(float(lines[j+2].split()[3]))
    return energies, forces

def main():
    path = '/anvil/scratch/x-ryu3/LAMMPS_interface/setup/DFT_slab_test/3pbc/'
    location = 'dft_calcs'
    
    #extracting dft energy and forces
    two_energy, two_forces = access_dft_energies_forces(path, location, (2,3))
    three_energy, three_forces = access_dft_energies_forces(path, location, (3,4))
    energy_dif = np.array(three_energy) - np.array(two_energy)
    
    #atomic hartree units to ev
    au_ev = 27.2114079527
    #number of atoms in your system
    n_atoms = 288
    
    #energy conversions and RMSE calculation
    energy_dif = energy_dif * au_ev / n_atoms * 1000
    energy_RMSE = np.sqrt(np.sum(energy_dif**2)/len(energy_dif))
    force_residual = np.array(three_forces) - np.array(two_forces)
    force_residual = force_residual * au_ev / 0.521 * 1000
    force_RMSE = np.sqrt(np.sum(force_residual**2)/len(force_residual))
    
    print(f'energy difference: {energy_RMSE} meV/atom')
    print(f'force RMSE: {force_RMSE} meV/angstrom')
    
    

if __name__ == '__main__':
    main()


