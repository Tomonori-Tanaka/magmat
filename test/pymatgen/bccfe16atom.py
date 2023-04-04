from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import symmetry as sym
from pymatgen.symmetry import site_symmetries
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import system

"""
parser = CifParser("./fe_54.cif")
structure = parser.get_structures()[0]
print(structure)
finder = site_symmetries.get_site_symmetries(structure, 0.01)[0][0]["Rot"]
print(finder)
"""

poscar = Poscar.from_file('fe_54.poscar',
                          check_for_POTCAR=False, read_velocities=False)
structure = poscar.structure
#print(structure)
# finder = site_symmetries.get_site_symmetries(structure, 0.01)
# print(type(finder[1][0]))

finder = SpacegroupAnalyzer(structure)
# print(vars(finder)["_space_group_data"]["translations"])

def get_pair_list(pair_mother):
    list = []
    center_atom = pair_mother[0][1]
    for it in pair_mother[1]:
        if len(list) == 0:
            list = [[center_atom, it[1]]]
        else:
            list.append([center_atom, it[1]])
    return list


trans_vec = np.empty(0)
for vec in vars(finder)["_space_group_data"]["translations"]:
    if trans_vec.size == 0:
        trans_vec = np.array([vec])
    else:
        check = True
        for vec2 in trans_vec:
            vec_tmp = np.subtract(vec, vec2)
            if np.linalg.norm(vec_tmp) < 1e-12:
                check = False
                continue
        if check:
            trans_vec = np.append(trans_vec, [vec], axis=0)

# print(trans_vec)

pair_0 = [[0, 0], [[0, 2],
                   [1, 48],
                   [2, 36],
                   [4, 32],
                   [5, 12],
                   [10, 28],
                   [11, 14],
                   [13, 10]]]

pair_1 = [[0, 0], [[0, 5],
           [0, 7],
           [5, 6],
           [10, 8],
           [13, 4]]]

# list = get_pair_list(pair_0)
#print(list)

system = system.System()

a = 2.866
lavec_in = np.array([[2*a, 0, 0], [0, 2*a, 0], [0, 0, 2*a]])
lavec_in = lavec_in.T
# print(lavec_in)
nat_in = 16
kind_in = [1]*16
xf_in =


#for t in vars(finder).items():
#	print(t["_symprec"])