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
import glob

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

def access_aims_energies_forces(path, location, klength):
    energies = []
    forces = []
    for i in range(*klength):
        filename = glob.glob(f'{path}{location}/structure{i}1/*.own')
        with open(filename[0],'r') as file:
            lines = file.readlines()
        switch = False
        for i, line in enumerate(lines):
            if 'Total atomic forces (unitary forces cleaned) [eV/Ang]:' in line:
                switch = True
            elif 'Start decomposition of the XC Energy' in line:
                switch = False
            elif 'Total energy uncorrected' in line:
                energies.append(float(line.split()[5]))
                energies.append(float(lines[i+1].split()[5]))
                energies.append(float(lines[i+2].split()[5]))
            elif switch:
                if len(line.split()) > 1:
                    forces.append(float(line.split()[2]))
                    forces.append(float(line.split()[3]))
                    forces.append(float(line.split()[4]))
    return energies, forces
    

def plot_force_distr_together(parameters, kpoints):
    '''
    graphs forces for multiple dft runs together and saves them as a png in a new force_distribution folder
    
    notes:
    - number of forces must be the same
    '''
    x = np.linspace(0,len(parameters[0]),len(parameters[0]))
    colors = plt.cm.inferno(np.linspace(0, 0.8, len(parameters)))
    plt.clf()
    for i, parameter in enumerate(parameters):
        plt.plot(x,parameter,color = colors[i], label = f'Force Distribution (meV/angstrom), {kpoints[i]}x{kpoints[i]} k points')
    plt.xlabel('Force number')
    plt.ylabel(f'Force (meV/angstrom)')
    plt.title(f'Forces on Converged Slab')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    os.makedirs('force_distribution', exist_ok = True)
    fig = plt.gcf()
    fig.savefig('force_distribution/force_distribution.png', dpi=300)
    
    #x,y,z directions
    plot_force_distr_together_coordinate(parameters, kpoints, 'x')
    plot_force_distr_together_coordinate(parameters, kpoints, 'y')
    plot_force_distr_together_coordinate(parameters, kpoints, 'z')
    
def plot_force_distr_together_coordinate(parameters,kpoints, direction):
    '''
    graphs forces together in a particular direction
        
    notes:
    - number of forces must be the same
    '''
    direction_conv = np.array(['x','y','z'])
    new_params = []
    for i in range(len(parameters)):
        new_params.append(parameters[i][np.where(direction_conv == direction)[0][0]::3])
    parameters = new_params
    x = np.linspace(0,len(parameters[0]),len(parameters[0]))
    colors = plt.cm.inferno(np.linspace(0, 0.8, len(parameters)))
    plt.clf()
    for i, parameter in enumerate(parameters):
        plt.plot(x,parameter,color = colors[i], label = f'Force Distribution (meV/angstrom), {kpoints[i]}x{kpoints[i]} k points')
    plt.xlabel('Force number')
    plt.ylabel(f'Force (meV/angstrom)')
    plt.title(f'Forces on Converged Slab ({direction})')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    os.makedirs('force_distribution', exist_ok = True)
    fig = plt.gcf()
    fig.savefig(f'force_distribution/force_distribution{direction}.png', dpi=300)
    

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
    
    #x,y,z forces
    plot_force_distr_coordinate(parameter, kpoints, 'x')
    plot_force_distr_coordinate(parameter, kpoints, 'y')
    plot_force_distr_coordinate(parameter, kpoints, 'z')
    
def plot_force_distr_coordinate(parameter, kpoints, direction):
    '''
    graphs forces in a particular direction
    input:
    - direction: x,y, or z
    '''
    direction_conv = np.array(['x','y','z'])
    parameter = parameter[np.where(direction_conv == direction)[0][0]::3]
    x = np.linspace(0,len(parameter),len(parameter))
    plt.clf()
    plt.plot(x,parameter,color = 'purple', label = f'Force Distribution ({direction}) (meV/angstrom), {kpoints}x{kpoints} k points')
    plt.xlabel('Force number')
    plt.ylabel(f'Force (meV/angstrom)')
    plt.title(f'Forces on Converged Slab ({direction})')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    os.makedirs('force_distribution', exist_ok = True)
    fig = plt.gcf()
    fig.savefig(f'force_distribution/force_distribution{kpoints}{direction}.png', dpi=300)

def plot_energy(energies, kpoints, aims):
    '''
    graphs energies
    '''
    plt.clf()
    colors = plt.cm.inferno(np.linspace(0, 0.8, 3))
    if aims:
        energy1 = energies[:,0]
        energy2 = energies[:,1]
        energy3 = energies[:,2]
        x = np.linspace(1,len(energy1),len(energy1))
        plt.plot(kpoints,energy1, color=colors[0], label = f'Energy (meV), total uncorrected energy')
        plt.plot(kpoints,energy1, color=colors[1], label = f'Energy (meV), total corrected energy')
        plt.plot(kpoints,energy1, color=colors[2], label = f'Energy (meV), electronic free energy')
    else:
        x = np.linspace(1,len(energies),len(energies))
        plt.plot(kpoints,energies, color=colors[0], label = f'Energy (meV)')
    plt.xlabel('K points')
    plt.ylabel(f'Energy (meV/angstrom)')
    plt.title(f'Energies On Slab')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    os.makedirs('force_distribution', exist_ok = True)
    fig = plt.gcf()
    fig.savefig(f'force_distribution/energy_distribution.png', dpi=300)

def main():
    path = '/sdcc/u/ryu1/copper_interface/setup/aims_copper/'
    location = 'dft_calcs'
    
    #extracting dft energy and forces
    
    all_energies = []
    all_forces = []
    kpoints = [5,6,7,8,9]
    for kpoint in kpoints:
        energy, force = access_aims_energies_forces(path, location, (kpoint,kpoint+1))
        all_energies.append(energy)
        all_forces.append(force)
    #convert to meV
    all_forces = np.array(all_forces) * 1000
    all_energies = np.array(all_energies) * 1000
    
    # #atomic hartree units to ev
    # au_ev = 27.2114079527
    # #number of atoms in your system
    # n_atoms = 288
    # #energy conversions and RMSE calculation
    # two_forces = np.array(two_forces) * au_ev / 0.521 * 1000
    # three_forces = np.array(three_forces) * au_ev / 0.521 * 1000
    
    plot_force_distr_together(all_forces, kpoints)
    plot_energy(all_energies, kpoints, True)
    #plot_force_distr(two_forces, 2)
    #plot_force_distr(three_forces,3)
    

if __name__ == '__main__':
    main()


