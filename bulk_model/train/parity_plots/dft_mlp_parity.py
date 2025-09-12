'''
This script accesses the previously calculated dft and mlp energy and forces and creates a parity plot

input parameters:
- test.xyz
- MACE_energies.npy
- MACE_forces.npy
- MACE_sizes.npy

output files:
- .png files of energy and force parity plots
- .png files of energy and force residual graphs

notes:
- use access_dft_energies_forces with make_plot_with_train
- use access_test_dft_energies_forces with make_plot_no_train
- for the most unbiased parity plot, use only the test set
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import ase.io

def access_dft_energies_forces(path,location, length):
    '''
    accesses the dft forces from the dft folder directly
    units are in atomic units
    '''
    
    energies = []
    forces = []
    for i in range(*length):
        with open(f'{path}{location}/structure{i}/ref_calc-1.ener','r') as file:
            lines = file.readlines()
        energies.append((float(lines[1].split()[4])))
        with open(f'{path}{location}/structure{i}/forces-output.xyz','r') as file:
            lines = file.readlines()
        for j in range(len(lines)-2):
            forces.append(float(lines[j+2].split()[1]))
            forces.append(float(lines[j+2].split()[2]))
            forces.append(float(lines[j+2].split()[3]))
    return energies, forces

def access_test_dft_energies_forces(path,location):
    '''
    accesses the dft forces through the test.xyz file (and thus has converted units)
    units are in eV, angstrom
    '''
    
    structures = ase.io.read(f'{path}{location}', index=':')
    energies = np.array([atoms.info['dft_energy'] for atoms in structures])
    forces = np.concatenate([atoms.arrays['dft_forces'].flatten() for atoms in structures])
    print(energies.shape)
    print(forces.shape)

    return energies, forces

def access_mlp_energies_forces(path,location,length):
    '''
    used for accessing mlp forces from mlp dump files - usually is not used
    '''
    
    energies = []
    forces = []
    for i in range(*length):
        with open(f'{path}{location}/structure{i}/energy.dat','r') as f:
            lines = f.readlines()
            energies.append(float(lines[2].split()[1]))
        with open(f'{path}{location}/structure{i}/forces.dump','r') as file:
            lines = file.readlines()
        for j in range(len(lines)-9):
            forces.append(float(lines[j+9].split()[5]))
            forces.append(float(lines[j+9].split()[6]))
            forces.append(float(lines[j+9].split()[7]))
    return energies, forces

def make_plot_with_train(path, parameter_dft, parameter_mlp, train_parameter_dft, train_parameter_mlp ,parameter, program):
    '''
    makes an energy and forces parity plot for a specific structure
    unit conversion is necessary for this due to how the dft forces and energies are taken
    '''
    
    #convert parameter units
    au_ev = 27.2114079527
    if parameter == 'Energies':
        #converting from a.u. to meV/atom
        parameter_dft = np.array(parameter_dft) * au_ev / 197 * 1000
        parameter_mlp = np.array(parameter_mlp) / 197 * 1000
        train_parameter_dft = np.array(train_parameter_dft) * au_ev / 197 * 1000
        train_parameter_mlp = np.array(train_parameter_mlp) / 197 * 1000
        units = 'meV/atom'
    if parameter == 'Forces':
        #converting from a.u. to meV/angstrom
        parameter_dft = np.array(parameter_dft) * au_ev / 0.521 * 1000
        train_parameter_dft = np.array(train_parameter_dft) * au_ev / 0.521 * 1000
        parameter_mlp = np.array(parameter_mlp) * 1000
        train_parameter_mlp = np.array(train_parameter_mlp) * 1000
        units = 'meV/A'
    #parity plot
    plt.figure()
    
    plt.scatter(parameter_dft,parameter_mlp, color = 'b', label = f'System {parameter}')
    plt.scatter(train_parameter_dft,train_parameter_mlp, color = 'r', label = f'System {parameter} (Train)')
    
    plt.xlabel(f'DFT {parameter} ({units})')
    plt.ylabel(f'MACE {parameter} ({units})')
    
    min_value = min(min(parameter_dft),min(parameter_mlp),min(train_parameter_dft),min(train_parameter_mlp))
    max_value = max(max(parameter_dft),max(parameter_mlp),max(train_parameter_dft),max(train_parameter_mlp))
    x = np.linspace(min_value,max_value,len(parameter_dft))
    
    plt.plot(x,x,label='y=x')
    plt.legend()
    plt.title('DFT and MLP Parity Plot')
    plt.tight_layout()
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    
    os.makedirs(f'{path}dft_mlp_parity_plots', exist_ok=True)
    fig = plt.gcf()
    fig.savefig(f'{path}dft_mlp_parity_plots/{parameter}_{program}_dft_mlp.png', dpi=300)

    #residual calculation
    residuals = np.array(parameter_dft) - np.array(parameter_mlp)
    residuals_train = np.array(train_parameter_dft) - np.array(train_parameter_mlp)
    RMSE = np.sqrt(np.sum(residuals**2)/len(residuals))
    RMSE_train = np.sqrt(np.sum(residuals_train**2)/len(residuals_train))
    structure_number = np.linspace(1,len(residuals),len(residuals))
    structure_number_train = np.linspace(1+len(residuals),len(residuals_train)+len(residuals),len(residuals_train))
    
    #residual plot
    plt.figure()
    plt.scatter(structure_number,residuals, color='b', label= f'RMSE: {RMSE:.3f} {units}')
    plt.scatter(structure_number_train,residuals_train, color='r', label= f'RMSE: {RMSE_train:.3f} {units} (Train)')
    plt.xlabel('Structure Number (1 ion pair)')
    plt.ylabel(f'Residual of DFT and MLP ({units})')
    plt.title(f'{parameter} Residual graph based on structure')
    plt.tight_layout()
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    fig = plt.gcf()
    fig.savefig(f'{path}dft_mlp_parity_plots/{parameter}_{program}_dft_mlp_residuals.png', dpi=300)

def make_plot_no_train(path, parameter_dft, parameter_mlp, parameter, program, sizes):
    '''
    makes an energy and forces parity plot for a specific structure
    unit conversion for this is unecessary since this has already been done for the test set
    '''
    
    au_ev = 27.2114079527
    if parameter == 'Energies':
        parameter_dft = np.array(parameter_dft) / sizes * 1000
        parameter_mlp = np.array(parameter_mlp) / sizes * 1000
        units = 'meV/atom'
    if parameter == 'Forces':
        parameter_dft = np.array(parameter_dft) * 1000
        parameter_mlp = np.array(parameter_mlp) * 1000
        units = 'meV/A'
    
    #parity plot
    plt.figure()

    #plotting data
    plt.scatter(parameter_dft,parameter_mlp, color = 'b', label = f'System {parameter}')
    plt.xlabel(f'DFT {parameter} ({units})')
    plt.ylabel(f'MACE {parameter} ({units})')
    
    #calculation for the y=x line
    min_value = min(min(parameter_dft),min(parameter_mlp))
    max_value = max(max(parameter_dft),max(parameter_mlp))
    x = np.linspace(min_value,max_value,len(parameter_dft))
    
    #plotting y=x line
    plt.plot(x,x,label='y=x')
    plt.legend()
    plt.title('DFT and MLP Parity Plot')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    
    #saving the figure
    os.makedirs(f'{path}dft_mlp_parity_plots', exist_ok=True)
    fig = plt.gcf()
    fig.savefig(f'{path}dft_mlp_parity_plots/{parameter}_{program}_dft_mlp.png', dpi=300)

    #calculating residuals
    residuals = np.array(parameter_dft) - np.array(parameter_mlp)
    RMSE = np.sqrt(np.sum(residuals**2)/len(residuals))
    structure_number = np.linspace(1,len(residuals),len(residuals))
    
    #residual figure
    plt.figure()
    plt.scatter(structure_number,residuals, color='b', label= f'RMSE: {RMSE:.3f} {units}')
    plt.xlabel('Structure Number (1 ion pair)')
    plt.ylabel(f'Residual of DFT and MLP ({units})')
    plt.title(f'{parameter} Residual graph based on structure')
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(-2,2))
    plt.legend()
    
    #saving the figure
    fig = plt.gcf()
    fig.savefig(f'{path}dft_mlp_parity_plots/{parameter}_{program}_dft_mlp_residuals.png', dpi=300)


def main():
    '''
    This implementation is for the test set only
    '''
    
    from_path = '/anvil/scratch/x-ryu3/Bulk/train/lr_train/'
    to_path = '/anvil/scratch/x-ryu3/Bulk/train/lr_train/train_new_128/verification/'
    mlp_energies = np.load('mlp_ase/MACE_energies.npy')
    mlp_forces = np.load('mlp_ase/MACE_forces.npy', allow_pickle = True)
    sizes = np.load('mlp_ase/MACE_sizes.npy')
    dft_energies, dft_forces = access_test_dft_energies_forces(from_path, 'test.xyz')

    make_plot_no_train(to_path, dft_energies, mlp_energies, 'Energies', 'ase', sizes)
    make_plot_no_train(to_path, dft_forces, mlp_forces, 'Forces', 'ase', sizes)
    

if __name__ == '__main__':
    main()


