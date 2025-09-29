'''
This script creates rdfs from mlp and aimd runs

notes:
- the first implementation is a manual rdf implementation, which is slower
- the second implementation is from an imported package, which is much faster
- for normal plotting: use run_through_elements (this is a slow implementation)
- for individual atom comparison for mlp vs. aimd: use plot_aimd
- for plotting mlp vs. aimd together in one subplot: use plot_aimd_together
- parameter specification is at the bottom
'''

from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from MDAnalysis.analysis.rdf import InterRDF
import MDAnalysis as mda
from ase.geometry import cell_to_cellpar

#####START OF MLP RDF
def create_structures(path, file):
    '''
    create an ase atoms object, changing the chemical symbols to be correct
    '''
    structures = read(f'{path}{file}', index = ':')
    new_structures = []
    
    for atoms in structures:
        nitrates = False
        new_symbols = []
        for symbol in atoms.get_chemical_symbols():
            if symbol == 'N':
                new_symbols.append('N')
                nitrates = True
            elif nitrates:
                new_symbols.append('O_nit')
            else:
                if symbol == 'O':
                    new_symbols.append('O_wat')
                elif symbol == 'Na':
                    new_symbols.append('NA')
                else:
                    new_symbols.append(symbol)
        new_symbols = np.array(new_symbols)
        atoms.set_array('names', new_symbols)
        new_structures.append(atoms)
    return new_structures



def make_one_rdf(structures, r_max, n_bins, element_1, element_2):
    '''
    Create the bins of an RDF for a single trajectory
    Currently only works for bulk simulations because of the density calculation messing up with 
    an interface vacuum
    '''

    rdf = np.zeros((len(structures),n_bins))
    for i, atoms in enumerate(structures):
        #preprocessing
        element_indices_1 = [i for i, symbol in enumerate(atoms.get_array('names')) 
                        if symbol == element_1]
        element_indices_2 = [i for i, symbol in enumerate(atoms.get_array('names')) 
                        if symbol == element_2]
        element_positions_1 = atoms.get_positions()[element_indices_1]
        element_positions_2 = atoms.get_positions()[element_indices_2]

        cell_x = atoms.get_cell()[0][0]
        cell_y = atoms.get_cell()[1][1]
        cell_z = atoms.get_cell()[2][2]
        box = [cell_x,cell_y,cell_z]

        element_positions_1[:,0] = element_positions_1[:,0] % cell_x
        element_positions_1[:,1] = element_positions_1[:,1] % cell_y
        element_positions_1[:,2] = element_positions_1[:,2] % cell_z

        element_positions_2[:,0] = element_positions_2[:,0] % cell_x
        element_positions_2[:,1] = element_positions_2[:,1] % cell_y
        element_positions_2[:,2] = element_positions_2[:,2] % cell_z

        #calculating the rdf
        dr = r_max/n_bins
        bins = np.arange(0, r_max + dr, dr)
        g_r = []
        for j in range(element_positions_1.shape[0]):
            for k in range(element_positions_2.shape[0]):
                if element_1 == element_2 and j == k:
                    continue
                delta = element_positions_1[j] - element_positions_2[k]
                delta = delta - box * np.round(delta / box)
                r = np.linalg.norm(delta)
                g_r.append(r)
        hist_counts, _ = np.histogram(g_r, bins)
        rdf[i,:] += hist_counts[:n_bins]


        n_atoms = element_positions_2.shape[0]
        n_ref_atoms = element_positions_1.shape[0]
        V = box[0] * box[1] * box[2]
        density = n_atoms / V
        for j in range(n_bins):
            r_i = bins[j] + dr/2
            shell_volume = 4 * np.pi * (r_i**2) * dr
            normalization = density * n_ref_atoms * shell_volume
            rdf[i,j] /= normalization

    rdf = rdf.mean(axis=0)
    bins = bins[:n_bins] + dr/2

    return bins, rdf

def combine_rdf_reps(path, file, runs, r_max, n_bins, element_1, element_2):
    '''
    add up rdfs of different reps for given elements
    '''
    all_rdf = []
    for rep in runs:
        new_path = f'{path}na_runs/na_rep{rep}/'
        structures = create_structures(new_path, file)
        bins, rdf = make_one_rdf(structures, r_max, n_bins, element_1, element_2)
        all_rdf.append(rdf)
    #calculate mean and standard deviation
    all_rdf_mean = np.array(all_rdf).mean(axis=0)
    all_rdf_error = np.array(all_rdf).std(axis=0) #/ np.sqrt(len(runs))
    return bins, all_rdf_mean, all_rdf_error
        
def run_through_elements(path, file, runs, r_max, n_bins, elements, folder):
    '''
    Create RDFs of all of the element combinations
    '''
    for i in range(len(elements)):
        for j in range(len(elements)):
            bins, rdf, rdf_error = combine_rdf_reps(path, file, runs, r_max, n_bins, elements[i], elements[j])
            plot_rdf(path, f'{elements[i]}-{elements[j]}', bins, rdf, rdf_error, folder)

def plot_rdf(path, title, bins, rdf, rdf_error, folder):
    os.makedirs(f'{path}{folder}', exist_ok=True)
    os.chdir(folder)
    plt.figure()
    plt.plot(bins, rdf, 'b-', linewidth=3, label='Mean')
    plt.fill_between(bins, 
                    rdf - rdf_error, 
                    rdf + rdf_error, 
                    alpha=0.3, color='blue', 
                    label=f'±std')
    plt.title(title)
    plt.ylabel('g(r)')
    plt.xlabel('Distance (A)')
    fig = plt.gcf()
    fig.savefig(f'{title} RDF')
    plt.close()
    os.chdir('..')

#######END OF MLP RDF

####TESTING MDANALYSIS vs my own
def plot_mda_comparison(path, file, runs, r_max, n_bins, element_1, element_2, folder):
    os.makedirs(f'{path}{folder}', exist_ok=True)
    os.chdir(folder)
    bins, rdf, rdf_error = combine_rdf_reps(path, file, runs, r_max, n_bins, element_1, element_2)
    mda_bins, mda_rdf, mda_error = combine_mdanalysis_rdf(path, file, runs, r_max, n_bins, element_1, element_2)
    plt.figure()
    plt.plot(bins, rdf, 'b-', linewidth=3, label='own_rdf')
    plt.fill_between(bins, 
                    rdf - rdf_error, 
                    rdf + rdf_error, 
                    alpha=0.3, color='blue', 
                    label=f'±std')
    plt.plot(mda_bins, mda_rdf, 'r-', linewidth=3, label='mda_rdf')
    plt.fill_between(mda_bins, 
                    mda_rdf - mda_error, 
                    mda_rdf + mda_error, 
                    alpha=0.3, color='r', 
                    label=f'±std')
    plt.title(f'{element_1}-{element_2}')
    plt.ylabel('g(r)')
    plt.xlabel('Distance (A)')
    fig = plt.gcf()
    fig.savefig(f'{element_1}-{element_2}-MDA-comp RDF')
    plt.close()
    os.chdir('..')
    


#######START OF MDANALYSIS MLP RDF
def make_mda_universe(new_path, file):
    u = mda.Universe(f'{new_path}{file}')
    u.dimensions = cell_to_cellpar(read(f'{new_path}{file}', index = '0').get_cell())
    return u


def one_rdf_mdanalysis(u, r_max, n_bins, element_1, element_2):
    atom_types = list(set(u.atoms.types))
    #print(atom_types)
    atom_types.remove('O')
    atom_types.append('O_wat')
    atom_types.append('O_nit')
    
    oxygens = u.select_atoms('type O')
    o_nitrate = oxygens[-3:]
    o_water = oxygens[:-3]
    
    atoms = []
    
    for element in [element_1,element_2]:
        if element == 'O_wat':
            atoms.append(o_water)
        elif element == 'O_nit':
            atoms.append(o_nitrate)
        else:
            atoms.append(u.select_atoms(f'type {element}'))

    rdf = InterRDF(atoms[0], atoms[1], nbins = n_bins, range = (0.1, min(u.dimensions[:2]/2)))
    rdf.run()
    return rdf.results.bins, rdf.results.rdf

def combine_mdanalysis_rdf(path, file, runs, r_max, n_bins, element_1, element_2):
    all_rdf = []
    for rep in runs:
        new_path = f'{path}{sub_path}{rep}/'
        #structures = create_structures(new_path, file)
        u = make_mda_universe(new_path, file)
        bins, rdf = one_rdf_mdanalysis(u, r_max, n_bins, element_1, element_2)
        all_rdf.append(rdf)
    #calculate mean and standard deviation
    all_rdf_mean = np.array(all_rdf).mean(axis=0)
    all_rdf_error = np.array(all_rdf).std(axis=0) #/ np.sqrt(len(runs))
    print(bins.shape)
    print(all_rdf_mean.shape)
    return bins, all_rdf_mean, all_rdf_error





########END OF MDANALYSIS MLP RDF

#### START OF AIMD RDF
def load_aimd_data(path):
    with open(path + 'rdfs.pkl', 'rb') as handle:
        rdfs = pickle.load(handle)
    with open(path + 'vdos.pkl', 'rb') as handle:
        vdos = pickle.load(handle)
    return rdfs, vdos

def average_aimd_replicates(paths):
    
#     print(paths)
    rdfs_all = []
    vdos_all = []
    for i, path in enumerate(paths):
#         print(path)
        rdfs, vdos = load_aimd_data(path)
        rdfs_all.append(rdfs)
        vdos_all.append(vdos)
    
    data = [rdfs_all, vdos_all]
    
    rdf_avg = {}
    vdos_avg = {}
    
    rdf_err = {}
    vdos_err = {}
    
    data_avg = [rdf_avg, vdos_avg]
    data_err = [rdf_err, vdos_err]
    
    for i, observable in enumerate(data):
        for key in observable[0].keys():
            values = [d[key] for d in observable]
            avg_values = np.mean(values,axis=0)
    
           
            err_values = np.std([value[1] for value in values],axis=0)
        
        
            data_avg[i][key] = avg_values
            data_err[i][key] = err_values
    
    return tuple(data_avg), tuple(data_err)

### PLOTTING SEPARATELY
def plot_all_f(data, data_err, labels, observable, folder):
    """Helper function to loop over all functions and plot them"""

    for name, d in data[0].items():
        plot_f([dat[name] for dat in data], [dat[name] for dat in data_err],
                          labels,
                          name,
                          observable, folder)

def plot_f(data, data_err, labels, title, observable, folder):
    os.makedirs(f'{path}{folder}', exist_ok=True)
    os.chdir(f'{path}{folder}')
    plt.figure()
    for i, d in enumerate(data):
        x = d[0]
        y = d[1]
        y_err = data_err[i]
        plt.plot(x, y, lw=3, linestyle = '-', color = 'r', label = 'aimd')
        if len(data_err[i]) != 0:
            plt.fill_between(x, y-y_err, y+y_err, alpha = 0.5, color = 'r')

    element_1 = title.split('-')[0]
    element_2 = title.split('-')[1]
    bins, rdf, rdf_error = combine_mdanalysis_rdf(path, file, runs, r_max, n_bins, element_1, element_2)
    
    plt.plot(bins, rdf, 'b', ls = ':', linewidth=3, label='MLP')

    plt.fill_between(bins, 
                    rdf - rdf_error, 
                    rdf + rdf_error, 
                    alpha=0.5, color='blue')
    
    plt.title("Species: " + str(title))
    plt.legend()
    plt.ylabel('g(r)')
    plt.xlabel('Distance (A)')
    fig = plt.gcf()
    fig.savefig(f'{title} RDF')
    plt.close()
    os.chdir(path)

###PLOTTING EVERYTHING ON A SUBPLOT
def plot_together(data, data_err, labels, observable, folder):
    """
    Helper function to loop over all functions and plot them
    This one plots everything on a single subplot
    """

    fig, ax = plt.subplots(3,5)
    i = 0
    for name, d in data[0].items():
        plot_everything([dat[name] for dat in data], [dat[name] for dat in data_err],
                          labels,
                          name,
                          observable, folder, ax[int(i/5),i%5])
        i += 1
    plt.subplots_adjust(
        left=0.08,    # Left padding (8% of figure width)
        right=0.92,   # Right padding (8% from right edge)
        top=0.95,     # Top padding (5% from top)
        bottom=0.08,  # Bottom padding (8% from bottom)
        hspace=0.4,   # Vertical spacing between subplots
        wspace=0.6    # Horizontal spacing between subplots
    )
    
    os.makedirs(f'{path}{folder}', exist_ok=True)
    os.chdir(f'{path}{folder}')
    plt.savefig('all_rdf.png', dpi=300)
    os.chdir(path)

def plot_everything(data, data_err, labels, title, observable, folder, ax):
    '''
    plotting everything on a single subplot
    '''
    for i, d in enumerate(data):
        x = d[0]
        y = d[1]
        y_err = data_err[i]
        ax.plot(x, y, lw=1, linestyle = '-', color = 'r', label = 'aimd')
        if len(data_err[i]) != 0:
            ax.fill_between(x, y-y_err, y+y_err, alpha = 0.5, color = 'r')
    
    element_1 = title.split('-')[0]
    element_2 = title.split('-')[1]
    bins, rdf, rdf_error = combine_mdanalysis_rdf(path, file, runs, r_max, n_bins, element_1, element_2)

    ax.plot(bins, rdf, 'b', ls=':', linewidth=1, label='MLP')

    ax.fill_between(bins, 
                    rdf - rdf_error, 
                    rdf + rdf_error, 
                    alpha=0.5, color='blue')
    
    ax.set_title(f"{title}", fontsize=8)
    #ax.legend(fontsize=8)
        


#### CALLING FUNCTIONS
def plot_aimd(aimd_path, folder, runs):
    '''
    plots aimd vs mlp on separate plots
    '''
    path = aimd_path
    paths = [f'{path}/run_{i}/' for i in runs] 

    data_avg, data_err = average_aimd_replicates(paths)
    rdfs, vdos = data_avg
    rdfs_err, vdos_err = data_err
    plot_all_f([rdfs], [rdfs_err], ['',''], 'RDF', folder)

def plot_aimd_together(aimd_path, folder, runs):
    '''
    plots aimd vs mlp on a single plot
    '''
    path = aimd_path
    paths = [f'{path}/run_{i}/' for i in runs] 

    data_avg, data_err = average_aimd_replicates(paths)
    rdfs, vdos = data_avg
    rdfs_err, vdos_err = data_err
    plot_together([rdfs], [rdfs_err], ['',''], 'RDF', folder)

    
    


path = '/anvil/scratch/x-ryu3/Bulk/verification/NPT/na_aimd_runs/' #keep like this for referencing mlp runs
models = ['lr64_chengucb','lr64_mace','lr128_chengucb','lr128_mace','sr64','sr128']
sub_path = 'Na_1_rep'
aimd_path = '/anvil/scratch/x-ryu3/aimd_nitrate/nano3_64water'
file = 'all_nvt.xyz' #keep like this for referencing mlp runs
#reps = 10
r_max = 6.5
n_bins = 150
elements = ['O_wat','H','NA','N','O_nit']
folder = 'no_run6_mda_rdf'
runs = [0,1,2,3,4,5,7,8,9] #no run 6
#run_through_elements(path, file, reps, r_max, n_bins, elements, folder)
for model in models:
    orig_path = path
    path = f'{path}{model}/'
    plot_aimd(aimd_path, folder, runs)
    plot_aimd_together(aimd_path, folder, runs)
    path = orig_path
#combine_mdanalysis_rdf(path, file, reps, r_max, n_bins, 'NA', 'N')
#plot_mda_comparison(path, file, reps, r_max, n_bins, 'O_wat', 'O_wat', folder)




