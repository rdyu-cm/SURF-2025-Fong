'''
This script converts the numbers output by lammps to atom names
'''

def replace_names(file, cation):
    type_map = {1: 'O',
                2: 'H',
                3: f'{cation}',
                4: 'N',
                5: 'O',
                6: 'Ti',
                7: 'H'}

    with open(file, 'r') as f:
        lines = f.readlines()
    
    modified_lines = []
    
    for line in lines:
        # Skip header lines (number of atoms and comment lines)
        if line.strip().isdigit() or line.startswith('#') or line.strip() == '':
            modified_lines.append(line)
        else:
            # Split the line into components
            parts = line.split()
            if len(parts) >= 4:  # Ensure we have at least atom_type, x, y, z
                try:
                    atom_type = int(parts[0])
                    # Replace the atom type number with atom name
                    if atom_type in type_map:
                        parts[0] = type_map[atom_type]
                    else:
                        print(f"Warning: Unknown atom type {atom_type}")
                    
                    # Reconstruct the line
                    modified_line = ' '.join(parts) + '\n'
                    modified_lines.append(modified_line)
                except ValueError:
                    # If first column isn't a number, keep the line as is
                    modified_lines.append(line)
            else:
                # Keep lines that don't have enough columns
                modified_lines.append(line)
    
    output_file = 'nvt_traj_names.xyz'
    
    with open(output_file, 'w') as f:
        f.writelines(modified_lines)
    
    print(f"Modified file saved as: {output_file}")

file = 'nvt_traj.xyz'
cation = 'Cs'
replace_names(file, cation)