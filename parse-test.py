# Assuming you have added `import numpy as np` at the top of your class code
# Assume the class definition is saved in a file named `fchk_class.py`
from fchk_class import Fchk
import csv


# Create an instance by providing the path to a Gaussian `.fchk` file
# rewrite this so that it takes a prefix and appends it to .fchk and adds it to the path
compound_id = "NiTrans_09-05-23_1"  # Replace with the actual compound name
compound_name = "NiTrans"

path = "/home/user/Documents/wheeler-local/" + compound_id + ".fchk"
fchk_data = Fchk(path, tag="some_tag")
csv_filename = compound_id+ ".csv"
N_orbitals = 5


def main(fchk_data, N_orbitals, compound_name, csv_filename):

    orbitals_info = get_orbitals_around_homo_lumo(fchk_data, N_orbitals)
    for label, info in orbitals_info.items():
        print(f"{label}: Orbital Number = {info['orbital_number']}, Energy = {info['energy']}")

    write_orbitals_to_csv(fchk_data, N_orbitals, compound_name, csv_filename)


def print_attributes(fchk_data):
    # Access the parsed attributes
    print("Title:", fchk_data.title)
    print("Calculation type:", fchk_data.calc)
    print("Basis set:", fchk_data.basis)
    print("Number of atoms:", fchk_data.nat)
    print("Total energy:", fchk_data.total_energy)
    print("Homo:", fchk_data.ahomo)
    print("Lumo:", fchk_data.alumo)
    print("orbital energies:", fchk_data.alpha_energies)


def get_orbitals_around_homo_lumo(fchk_data, N_orbitals):
    result = {}
    
    # HOMO energies and labels
    start_homo = fchk_data.ahomo - N_orbitals + 1
    end_homo = fchk_data.ahomo + 1  # +1 to include HOMO itself
    for i in range(start_homo, end_homo):
        if i < 0:
            continue
        label_diff = fchk_data.ahomo - i
        label = f"HOMO-{label_diff}" if label_diff != 0 else "HOMO"
        result[label] = {
            'orbital_number': i,
            'energy': fchk_data.alpha_energies[i]
        }
    
    # LUMO energies and labels
    start_lumo = fchk_data.alumo
    end_lumo = fchk_data.alumo + N_orbitals  # +N_orbitals to include LUMO itself and N_orbitals above it
    for i in range(start_lumo, end_lumo):
        if i >= len(fchk_data.alpha_energies):
            break
        label_diff = i - fchk_data.alumo
        label = f"LUMO+{label_diff}" if label_diff != 0 else "LUMO"
        result[label] = {
            'orbital_number': i, 
            'energy': fchk_data.alpha_energies[i]
        }
    
    return result


def write_orbitals_to_csv(fchk_data, N_orbitals, compound_name, csv_filename):
    # Get the orbital information
    orbitals_info = get_orbitals_around_homo_lumo(fchk_data, N_orbitals)
    
    # Write to CSV
    with open(csv_filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        # Write the header
        csvwriter.writerow(['compound', 'unround_eV', 'eV', 'symmetry_label', 'orbital_label', 'orbital_num'])
        
        for label, info in orbitals_info.items():
            # For this example, I'm rounding the energy to 2 decimal places
            rounded_energy = round(info['energy'], 2)
            
            # Writing each row. Skipping columns with no data
            csvwriter.writerow([compound_name, info['energy'], rounded_energy, r'" "', label, info['orbital_number']])

main(fchk_data, N_orbitals, compound_name, csv_filename)