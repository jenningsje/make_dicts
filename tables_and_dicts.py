import pandas as pd
import numpy as np
from RadialPsi import *
from zeff import *
from nucleiterms.current.bifold.matematik import mesh
from nucleiterms.current.bifold.simpson.simpson import *
import json

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
atom_dict = {'H': {'n': 1, 'l': 0}, 'He': {'n': 1, 'l': 0}, 'Li': {'n': 2, 'l': 0}, 'Be': {'n': 2, 'l': 0}, 'B': {'n': 2, 'l': 1}, 'C': {'n': 2, 'l': 1}, 'N': {'n': 2, 'l': 1}, 'O': {'n': 2, 'l': 1}, 'F': {'n': 2, 'l': 1}, 'Ne': {'n': 2, 'l': 1}, 'Na': {'n': 3, 'l': 0}, 'Mg': {'n': 3, 'l': 0}, 'Al': {'n': 3, 'l': 1}, 'Si': {'n': 3, 'l': 1}, 'P': {'n': 3, 'l': 1}, 'S': {'n': 3, 'l': 1}, 'Cl': {'n': 3, 'l': 1}, 'Ar': {'n': 3, 'l': 1}, 'K': {'n': 4, 'l': 0}, 'Ca': {'n': 4, 'l': 0}, 'Sc': {'n': 3, 'l': 2}, 'Ti': {'n': 3, 'l': 2}, 'V': {'n': 3, 'l': 2}, 'Cr': {'n': 3, 'l': 2}, 'Mn': {'n': 3, 'l': 2}, 'Fe': {'n': 3, 'l': 2}, 'Co': {'n': 3, 'l': 2}, 'Ni': {'n': 3, 'l': 2}, 'Cu': {'n': 3, 'l': 2}, 'Zn': {'n': 3, 'l': 2}, 'Ga': {'n': 4, 'l': 1}, 'Ge': {'n': 4, 'l': 1}, 'As': {'n': 4, 'l': 1}, 'Se': {'n': 4, 'l': 1}, 'Br': {'n': 4, 'l': 1}, 'Kr': {'n': 4, 'l': 1}, 'Rb': {'n': 5, 'l': 0}, 'Sr': {'n': 5, 'l': 0}, 'Y': {'n': 4, 'l': 2}, 'Zr': {'n': 4, 'l': 2}, 'Nb': {'n': 4, 'l': 2}, 'Mo': {'n': 4, 'l': 2}, 'Tc': {'n': 4, 'l': 2}, 'Ru': {'n': 4, 'l': 2}, 'Rh': {'n': 4, 'l': 2}, 'Pd': {'n': 4, 'l': 2}, 'Ag': {'n': 4, 'l': 2}, 'Cd': {'n': 4, 'l': 2}, 'In': {'n': 5, 'l': 1}, 'Sn': {'n': 5, 'l': 1}, 'Sb': {'n': 5, 'l': 1}, 'Te': {'n': 5, 'l': 1}, 'I': {'n': 5, 'l': 1}, 'Xe': {'n': 5, 'l': 1}, 'Cs': {'n': 6, 'l': 0}, 'Ba': {'n': 6, 'l': 0}, 'Hf': {'n': 5, 'l': 2}, 'Ta': {'n': 5, 'l': 2}, 'W': {'n': 5, 'l': 2}, 'Re': {'n': 5, 'l': 2}, 'Os': {'n': 5, 'l': 2}, 'Ir': {'n': 5, 'l': 2}, 'Pt': {'n': 5, 'l': 2}, 'Au': {'n': 5, 'l': 2}, 'Hg': {'n': 5, 'l': 2}, 'TL': {'n': 6, 'l': 1}, 'Pb': {'n': 6, 'l': 1}, 'Bi': {'n': 6, 'l': 1}, 'Po': {'n': 6, 'l': 1}, 'At': {'n': 6, 'l': 1}, 'Rn': {'n': 6, 'l': 1}, 'Fr': {'n': 7, 'l': 0}, 'Ra': {'n': 7, 'l': 0}, 'Rf': {'n': 6, 'l': 2}, 'Db': {'n': 6, 'l': 2}, 'Sg': {'n': 6, 'l': 2}, 'Bh': {'n': 6, 'l': 2}, 'Hs': {'n': 6, 'l': 2}, 'Mt': {'n': 6, 'l': 2}, 'Ds': {'n': 6, 'l': 2}, 'Rg': {'n': 6, 'l': 2}, 'Cn': {'n': 6, 'l': 2}, 'Nh': {'n': 7, 'l': 1}, 'FL': {'n': 7, 'l': 1}, 'Mc': {'n': 7, 'l': 1}, 'Lv': {'n': 7, 'l': 1}, 'Ts': {'n': 7, 'l': 1}, 'Og': {'n': 7, 'l': 1}, 'La': {'n': 5, 'l': 2}, 'Ce': {'n': 4, 'l': 3}, 'Pr': {'n': 4, 'l': 3}, 'Nd': {'n': 4, 'l': 3}, 'Pm': {'n': 4, 'l': 3}, 'Sm': {'n': 4, 'l': 3}, 'Eu': {'n': 4, 'l': 3}, 'Gd': {'n': 4, 'l': 3}, 'Tb': {'n': 4, 'l': 3}, 'Dy': {'n': 4, 'l': 3}, 'Ho': {'n': 4, 'l': 3}, 'Er': {'n': 4, 'l': 3}, 'Tm': {'n': 4, 'l': 3}, 'Yb': {'n': 4, 'l': 3}, 'Lu': {'n': 4, 'l': 3}, 'Ac': {'n': 6, 'l': 2}, 'Th': {'n': 6, 'l': 2}, 'Pa': {'n': 5, 'l': 3}, 'U': {'n': 5, 'l': 3}, 'Np': {'n': 5, 'l': 3}, 'Pu': {'n': 5, 'l': 3}, 'Am': {'n': 5, 'l': 3}, 'Cm': {'n': 5, 'l': 3}, 'Bk': {'n': 5, 'l': 3}, 'Cf': {'n': 5, 'l': 3}, 'Es': {'n': 5, 'l': 3}, 'Rm': {'n': 5, 'l': 3}, 'Md': {'n': 5, 'l': 3}, 'No': {'n': 5, 'l': 3}, 'Lr': {'n': 5, 'l': 3}}
atom_list = atomic_radii['symbol'].tolist()

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
q = mesh(500000, 300000, 5000)

e_lab = 141.7 # collission energy which is zero
V_nuc_dict = {}
V_nuc_df = pd.DataFrame({'atom pair': [], 'energy': [], 'min dist': [], 'max dist': []})

# create nucle-nuclei energy dictionary
for i in range(m):
    # retrieve data from nuclei-nuclei repulsion for nuclei i
    Zeff_i = Zeff_dict[atomic_radii["symbol"].iloc[i]]
    radius_i = atomic_radii["atomic radius"].iloc[i]
    a_proj = atomic_radii["atomic number"].iloc[i]
    atom1 = atomic_radii["symbol"]
    for j in range(m):
        # retreive data for nuclei-nuclei repulsion for nuclei j
        Zeff_j = Zeff_dict[atomic_radii["symbol"].iloc[j]]
        radius_j = atomic_radii["atomic radius"].iloc[j]
        min_dist = radius_i + radius_j
        atom2 = atomic_radii["symbol"]

        # obtain nuclei-nuclei repulsion energy for the yukawa potential
        if Zeff_i != 0 and Zeff_j != 0:
            rho_p = f_2prm_gaussian(r, ke * (Zeff_i ** 2) / min_dist, 1)
            rho_t = f_2prm_gaussian(r, ke * (Zeff_j ** 2) / min_dist, 1)
            V_nuc = u_m3y_reid_zr(e_lab, a_proj, rho_p, rho_t, r, q)
        else:
            V_nuc = 0

        # create dictionary containing databases corresponding to each nuclei-nuclei pair
        V_nuc_dict[atomic_radii["symbol"].iloc[i] + atomic_radii["symbol"].iloc[j]] = V_nuc

output_file_path = "V_nuc_dict.pkl"

with open(output_file_path, 'wb') as output:
    json.dump(V_nuc_dict, output)

orbital_to_n = {}
orbital_to_l = {}

for index3, row3 in energies.iterrows():
    orbital_to_n[row3["Orbitals"]] = row3["n"]
    orbital_to_l[row3["Orbitals"]] = row3["l"]

element_to_orbital = {}

for i in range(0, 10):
    for j in range(0, 18):
        element_to_orbital[periodic_array[i][j]] = orbital_array[i][j]

for atom in atom_dict:
    print(str(atom_dict[atom]) + " " + str(atom) + " " + str(atom_dict[atom]['n']))
    n_num = atom_dict[atom]['n']
    l_num = atom_dict[atom]['l']
    R1 = R(n_num, l_num, 5.0)
