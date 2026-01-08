# Rationalizing Specific Ion Effects in Electrochemical Nitrate Reduction Using Machine Learning Potential Simulations

This repository contains most of the set up and analysis code for this project. Since most .xyz files can be created from the set up scripts, structure files, for the most part, have been omitted. 

## Overview
- This repository contains scripts for a bulk electrolyte system (nitrate, water, cations) as well as a interface system (nitrate, water, cations, TiH2)
- Main methods used: LAMMPS, ASE, MACE, CP2K, FHI-AIMS
- The primary goal of this project is to investigate nitrate adsorption properties in solution with different supporting cations by creating a machine learning potential to be able to simulate bulk electrolyte and interfacial nitrate environments

## Organization
- This repository is organized by system: bulk electrolyte and interface, as well as set up, training, and analysis subdirectories. The relevant scripts are contained in each of those directories
- File System:
    - bulk_model
        - setup
            - dft_setup: Contains scripts to set up a DFT system in CP2K from the (mace_setup) output
            - lammps_setup: Contains scripts to set up an initial system with packmol and moltemplate, as well as run NPT to equilibrate with a traditional forcefield
            - mace_setup: Contains scripts to run NVT in LAMMPS on the systems from (lammps_setup) with a MACE foundation model
        - train: Contains scripts to train a MLP with the MACE architecture (including the LR model) from the dft output in (dft_setup)
            - parity_plots: Contains scripts to generate a parity plot for the test set between dft output and trained MLP output
        - verification: Contains scripts to run the trained MLP in ASE, as well as generate RDFs from the trajectories
    - interface_model
        - setup
            - fhi-aims: Contains a control.in file necessary to run fhi-aims
            - lammps_setup: Contains scripts to generate and run a TiH2 electrolyte interface system in LAMMPS with a MACE foundation model. This includes a sandwich (metal between two layers of electrolyte) set up, as well as a one sided vacuum (slab, then electrolyte, then vacuum) set up. The cations for these systems include most of the alkali metal ions, as well as hydronium. 
                - lammps_analysis: Contains scripts to run and post-process the LAMMPS runs
            - slab_setup: Contains scripts to generate and run a pure metal slab (TiH2) in CP2K and FHI-AIMS. This is used for dft parameter testing, including k points and mixing. 
        - train: Currently empty
        - verification: Currently empty
    - environment_installations: Contains .md files for installing conda environments

## Sample Scripts
- MACE in LAMMPS
    - Relevant Folder: interface_model/setup/lammps_setup
    - Very similar to regular LAMMPS, but with a couple of differences
        - In the system.in file, you have to change the pair_style and pair_potential to point to the model you want to use
        - LAMMPS must be rebuilt to be mace compatible. A build is on the fong repo, but without ml-iap. 
        - Note: .model files have to be converted to .pt files which are compatible with LAMMPS. This can be done with a file on the mace github, mace/cli/create_lammps_file.py
- MACE in ASE
    - Relevant Folder: bulk_model/setup/mace_setup/mace_runs
    - The run_1.py file takes in a .xyz file and runs dynamics with a mace foundation model. A different foundation model can be used with MACECalculator.

