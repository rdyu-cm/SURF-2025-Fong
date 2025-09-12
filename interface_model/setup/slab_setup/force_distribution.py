'''
This script prints out energy and force RMSE between two dft outputs that were calculated with different numbers of k points

input parameters:
- path: the location of your dft calculations folders
- n_atoms: the number of atoms in your system

output:
- saves as png: The forces of the two dft runs together, as well as individually
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

def plot_force_distr_together(parameter2,parameter3):
    '''
    graphs forces for two dft runs together and saves them as a png in a new force_distribution folder
    
    notes:
    - number of forces must be the same
    '''
    x = np.linspace(0,len(parameter2),len(parameter2))
    plt.clf()
    plt.plot(x,parameter2,color = 'purple', label = f'Force Distribution (meV/angstrom), 2x2 k points')
    plt.plot(x,parameter3,color = 'r', label = f'Force Distribution (meV/angstrom), 3x3 k points')
    plt.xlabel('Force number')
    plt.ylabel(f'Force (meV/angstrom)')
    plt.title(f'Forces on Converged Slab')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    os.makedirs('force_distribution', exist_ok = True)
    fig = plt.gcf()
    fig.savefig('force_distribution/force_distribution.png', dpi=300)

def plot_force_distr(parameter, kpoints):
    '''
    graphs forces a single dft run and saves them as a png in a new force_distribution folder
    '''
    x = np.linspace(0,len(parameter),len(parameter))
    plt.clf()
    plt.plot(x,parameter,color = 'purple', label = f'Force Distribution (meV/angstrom), {kpoints}x{kpoints} k points')
    plt.xlabel('Force number')
    plt.ylabel(f'Force (meV/angstrom)')
    plt.title(f'Forces on Converged Slab')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    os.makedirs('force_distribution', exist_ok = True)
    fig = plt.gcf()
    fig.savefig(f'force_distribution/force_distribution{kpoints}.png', dpi=300)

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
    two_forces = np.array(two_forces) * au_ev / 0.521 * 1000
    three_forces = np.array(three_forces) * au_ev / 0.521 * 1000
    
    plot_force_distr_together(two_forces,three_forces)
    plot_force_distr(two_forces, 2)
    plot_force_distr(three_forces,3)
    
    
    
    

if __name__ == '__main__':
    main()


