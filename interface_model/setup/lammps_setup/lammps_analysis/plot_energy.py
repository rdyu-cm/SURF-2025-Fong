'''
This script plots the energy over time for a NVT run in lammps from the .log file
'''

import matplotlib.pyplot as plt

def plot_energy(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    switch = False
    pot_energy = []
    for line in lines():
        if not line.split():
            continue
        if line.split()[0] == 'Loop':
            switch = False
        if switch:
            if line.split()[0] == 'WARNING:':
                continue
            else:
                pot_energy.append(float(line.split()[6]))
        if line.split()[0] == 'Time':
            switch = True

    plt.plot(pot_energy,label='Energy (eV)')
    plt.xlabel('Time (ns)')
    plt.ylabel('Energy (ev)')
    plt.title('NVT Energy Over Time')
    plt.savefig('pot_energy.png')
    plt.show()

file = 'nvt.log'
plot_energy(file)
    
