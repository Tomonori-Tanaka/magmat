from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import symmetry as sym
from pymatgen.symmetry import site_symmetries
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np

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


array = np.empty(0)
for vec in vars(finder)["_space_group_data"]["translations"]:
    if array.size == 0:
        array = np.array([vec])
        # array = np.append(array, vec, axis=1)
    else:
        check = True
        for vec2 in array:
            vec_tmp = np.subtract(vec, vec2)
            if np.linalg.norm(vec_tmp) < 1e-12:
                check = False
                continue
        if check:
            array = np.append(array, [vec], axis=0)

print(array)

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

list = get_pair_list(pair_0)
print(list)

#for t in vars(finder).items():
#	print(t["_symprec"])