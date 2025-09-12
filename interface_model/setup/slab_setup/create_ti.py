'''
This script creates a slab of TiH2

input parameters:
- slab_size: the number of units cells by unit cells for each layer
- layers: the number of layers (z direction) for the slab

output file:
- slab.xyz
'''

from ase import Atoms
from ase.build import fcc111, surface, bulk
from ase.io import read, write

def create_slab2(slab_size, layers):
    '''
    creates a TiH2 slab with the specified slab_size and number of layers
    '''

    a = 3.128 #literature unit cell length
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
    
    slab = surface(bulk_tih2, (1,1,1), layers = layers,vacuum = a)
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
    top_H = surface(top_atoms, (1,1,1), layers=layers, vacuum=0)
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
    bottom_H = surface(bottom_atoms, (1,1,1), layers=layers, vacuum=0)
    bottom_H = bottom_H.repeat((slab_size,slab_size,1))
    positions = bottom_H.get_positions()
    positions[:, 2] += a - Ti_shift_slab.cell[2][2]# shift in z-direction
    bottom_H.set_positions(positions)
    
    slab = slab + top_H + bottom_H

    write('slab.xyz', slab)


    return slab

slab_size = 6
layers = 6
create_slab2(slab_size,layers)