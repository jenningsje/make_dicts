import sys
import pandas as pd
import numpy as np
from RadialPsi import *
import sys
from zeff import *
from nucleiterms.current.bifold.matematik import mesh
from nucleiterms.current.bifold.filon import filon

sys.path.insert(0, '/Users/james/Desktop/file_cabinet/work/bioinformatics/github/make_dfs')
print(sys.path)

# import text files

# import items for the Mij matrix Below
# import the acid table
table_file0 = open("./make-tables/txt/schain_prob_table.txt")
lines0 = table_file0.readlines()
table_file0.close()

# import items for the Hij matrix Below
# import the hamiltonian chart
table_file = open("./make-tables/txt/hbond_type_table.txt")
lines = table_file.readlines()
table_file.close()

# import the amino acid list
seq_file = open("./make-tables/txt/list.txt")
seq_string = seq_file.read()
seq_file.close()

# import the periodic table
periodic_file = open("./make-tables/txt/periodic_table.txt")
periodic_table_lines = periodic_file.read()
periodic_file.close()

# import the periodic table where the element symbols are replaced with orbitals
orbital_file = open("./make-tables/txt/orbitals.txt")
orbital_file_lines = orbital_file.read()
orbital_file.close()

# import the atomic radii, atomic numbers and energies
atomic_radii = pd.read_excel("./make-tables/excel/radii.xlsx")
atomic_numbers = pd.read_excel("./make-tables/excel/atomic_numbers.xlsx")
energies = pd.read_excel("./make-tables/excel/energies.xlsx")

# create dictionaries

atom_dict_b = {b'C': 1, b'N': 2, b'O': 3, b'ZN': 4, b'None': 5}
aa_filter_b = {b'ALA', b'ARG', b'ASN', b'ASP', b'CYS', b'GLU', b'GLN', b'GLY', b'HIS', b'ILE', b'LEU', b'LYS', b'MET', b'PHE', b'PRO', b'SER', b'THR', b'TRP', b'TYR', b'VAL', b'HOH'}
aa_dict_b = {b'ALA': 0, b'ARG': 1, b'ASN': 2, b'ASP': 3, b'CYS': 4, b'GLU': 5, b'GLN': 6, b'GLY': 7, b'HIS': 8, b'ILE': 9, b'LEU': 10, b'LYS': 11, b'MET': 12, b'PHE': 13, b'PRO': 14, b'SER': 15, b'THR': 16, b'TRP': 17, b'TYR': 18, b'VAL': 19, b'NAG': 20, b'HOH': 21}
atom_filter = {'C', 'N', 'O', 'ZN'}
aa_dict = {"ALA": 1, "Arg": 2, "ASN": 3, "ASP": 4, "CYS": 5, "GLU": 6, "GLN": 7, "GLY": 8, "HIS": 9, "LIE": 10, "LEU": 11, "LYS": 12, "MET": 13, "PHE": 14, "PRO": 15, "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20}

# coloumb constant
ke = 8.9875517873681764 * 10 ** 9
e = 1.60217663 * 10 * -19
A = []
E = []
Nuclei_V_dict = {}
aa_table = []
h_table = []
all_Zeffs_dict = {}
Zeff_dict = {}
V_nuc = []
periodic_array = []
orbital_array = []

# split the sidechain probability table
acids0 = lines0[0].split()

# split the amino acid list
seq_list = seq_string.split()

# split the periodic table
elements = periodic_table_lines.split('\n')

# split the orbtials
orbitals = orbital_file_lines.split('\n')

# parameters for the Hij matrix below
n = len(seq_list)
seq_index = []

# parameters for the nuclei-nuclei interaction matrix
m = len(atomic_radii)

# construct the periodic table matrix
for periodic_table_line in elements[1:]:
    periodic_table_row = periodic_table_line.split()[1:]
    elements = list(map(str,periodic_table_row))
    periodic_array.append(elements)
    
# construct the periodic thable where the elements are replaced with orbitals
for orbital_file_line in orbitals[1:]:
    orbital_file_row = orbital_file_line.split()[1:]
    orbitals = list(map(str, orbital_file_row))
    orbital_array.append(orbitals)

# construct the matrix for the types of hydrogen bonds
for line in lines[1:]:
    row = line.split()[1:]
    letters = list(map(str,row))
    h_table.append(letters)

# construct the hydrogen bond probability matrix
for line0 in lines0[1:]:
    row0 = line0.split()[1:]
    numbers0 = list(map(float,row0))
    aa_table.append(numbers0)

# construct the matrix for the energy levels
for i in range(n):
    E.append([])
    for j in range(n):
        prob = h_table[i][j]
        if prob == "N":
            prob2 = -1.64013e-22
            E[i].append(prob2)
        elif prob == "O":
            prob2 = -2.09727e-22
            E[i].append(prob2)
        elif prob == "P":
            prob2 = 0.0
            E[i].append(prob2)
        else:
            prob2 = 0.0
            E[i].append(prob2)

# tensor product of acid_table0 and Eij
for i in range(n):
    A.append([])
    for j in range(n):
            prob = aa_table[i][j]*E[i][j]
            A[i].append(prob)

for i in range(m):
    all_Zeffs_dict[atomic_radii["symbol"].iloc[i]] = slater(atomic_radii["symbol"].iloc[i])
    Zeff_dict[atomic_radii["symbol"].iloc[i]] = all_Zeffs_dict[atomic_radii["symbol"].iloc[i]][-1][-1]

r = mesh(500000, 1000000, 5000)
q = mesh(500000, 1000000, 5000)

e_lab = 0 # collission energy which is zero
V_nuc = np.empty(shape=(m, m), dtype='object')
V_nuc_dict = {}

# create nucle-nuclei energy dictionary
for i in range(m):
    Zeff_i = Zeff_dict[atomic_radii["symbol"].iloc[i]]
    radius_i = atomic_radii["atomic radius"].iloc[i]
    a_proj = atomic_radii["atomic number"].iloc[i]
    for j in range(m):
        Zeff_j = Zeff_dict[atomic_radii["symbol"].iloc[j]]
        radius_j = atomic_radii["atomic radius"].iloc[j]
        min_dist = radius_i + radius_j
        if Zeff_i != 0 and Zeff_j != 0:
            rho_p = f_2prm_gaussian(r, ke * (Zeff_i ** 2) / min_dist, 1)
            rho_t = f_2prm_gaussian(r, ke * (Zeff_j ** 2) / min_dist, 1)
            V_nuc[i][j] = u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)
        else:
            V_nuc[i][j] = 0
            
        V_nuc_dict[atomic_radii["symbol"].iloc[i] + atomic_radii["symbol"].iloc[j]] = V_nuc[i][j]
        
V_nuc_df = pd.DataFrame(V_nuc_dict)
        
orbital_to_n = {}
orbital_to_l = {}

for index3, row3 in energies.iterrows():
    orbital_to_n[row3["Orbitals"]] = row3["n"]
    orbital_to_l[row3["Orbitals"]] = row3["l"]

element_to_orbital = {}

for i in range(0, 10):
    for j in range(0, 18):
        element_to_orbital[periodic_array[i][j]] = orbital_array[i][j]

element_to_quantum_numbers = {}

for i in range(len(periodic_array)):
    for j in range(len(periodic_array[i])):
        element = periodic_array[i][j]
        if element != 'X':
            orbital = orbital_array[i][j]
            n = orbital_to_n.get(orbital, 'Unknown')
            print(n)
            l = orbital_to_l.get(orbital, 'Unknown')
            element_to_quantum_numbers[element] = {'n': n, 'l': l}

elements = list(element_to_quantum_numbers.keys())
quantum_numbers_bad = list(element_to_quantum_numbers.values())
quantum_numbers = [item for item in quantum_numbers_bad if not (item == {'n': 'Unknown', 'l': 'Unknown'})]
vector = np.vectorize(np.int_)

