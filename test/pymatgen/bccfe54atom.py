from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen import symmetry as sym
from pymatgen.symmetry import site_symmetries
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import sys
import system
import libs.mathfunctions
from get_moments import get_moments

"""
parser = CifParser("./fe_54.cif")
structure = parser.get_structures()[0]
print(structure)
finder = site_symmetries.get_site_symmetries(structure, 0.01)[0][0]["Rot"]
print(finder)
"""

poscar = Poscar.from_file('fe_54.vasp',
                          check_for_POTCAR=False, read_velocities=False)
structure = poscar.structure
# print(structure)
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


def is_identical_position(vec1, vec2):
    vec_tmp = np.subtract(vec1, vec2)
    if np.linalg.norm(vec_tmp) < 1e-3:
        return True
    return False


def calc_dot_all_pairs(pair_list, moment_data):
    """

    :param pair_list: list of lists [[0, 1], [0, 2], ...]
    :param moment_data: [[theta, phi, ...]#atom1, [...]#atom2, ...]
    :return: sum of inner products (float)
    """
    in_product = 0.0
    for pair in pair_list:
        at0 = pair[0]
        at1 = pair[1]
        in_product += moment_data[at0][3] * moment_data[at1][3] \
                      + moment_data[at0][4] * moment_data[at1][4] \
                      + moment_data[at0][5] * moment_data[at1][5]
    return in_product


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

pair_0 = [[0, 0], [[0, 1],
                   [1, 47],
                   [2, 35],
                   [4, 31],
                   [5, 11],
                   [10, 27],
                   [11, 13],
                   [13, 9]]]

pair_1 = [[0, 0], [[0, 2],
                   [0, 4],
                   [0, 6],
                   [5, 5],
                   [10, 7],
                   [13, 3]]]

# list1 = get_pair_list(pair_0)
# print(list1)

system = system.System()

a = 2.86600000000
lavec_in = np.array([[3 * a, 0, 0], [0, 3 * a, 0], [0, 0, 3 * a]])
lavec_in = lavec_in.T
# print(lavec_in)
nat_in = 54
kind_in = [1] * 54

xf_in = []
with open("fe_54.positions", mode='r', encoding='utf-8') as f:
    for line in f:
        line = line.split()
        line = list(map(lambda x: float(x), line))
        xf_in.append(line)
# print(xf_in)
system.set_supercell(lavec_in, nat_in, kind_in, xf_in)
system.initialize()
# print(len(system.x_image))
# print(system.x_image[1][47])

pair_0_list = get_pair_list(pair_0)
pair_1_list = get_pair_list(pair_1)

origin_atom_coord = system.x_image[0][0]
# print(origin_atom_coord)
# print(system.supercell.lattice_vector)
# print(system.supercell.x_fractional)

for i, tran_tmp in enumerate(trans_vec):
    trans_vec_cart = libs.mathfunctions.rotvec(tran_tmp, system.supercell.lattice_vector)
    if is_identical_position(trans_vec_cart, np.array([0., 0., 0.])):
        continue
    # print(tran_tmp)
    # print(trans_vec_cart)
    equi_atom_coord = origin_atom_coord + trans_vec_cart
    # print(equi_atom_coord)
    origin_index_tmp = None
    for icell, atom in enumerate(system.x_image):
        for k, dummy in enumerate(atom):
            coord = system.x_image[icell][k]
            # print(coord)
            if is_identical_position(equi_atom_coord, coord):
                origin_index_tmp = k
    if origin_index_tmp == None:
        sys.exit("ERROR: origin_index_tmp is None")
    for pair_0_index in pair_0[1]:
        icell_index = pair_0_index[0]
        atom_index = pair_0_index[1]
        atom_coord = system.x_image[icell_index][atom_index]
        atom_coord_new = atom_coord + trans_vec_cart

        new_index_tmp = None
        for icell, atom in enumerate(system.x_image):
            for k, dummy in enumerate(atom):
                coord = system.x_image[icell][k]
                if is_identical_position(atom_coord_new, coord):
                    new_index_tmp = k
        if new_index_tmp == None:
            sys.exit("ERROR: new_index_tmp is None")
        pair_0_list.append([origin_index_tmp, new_index_tmp])

    for pair_1_index in pair_1[1]:
        icell_index = pair_1_index[0]
        atom_index = pair_1_index[1]
        atom_coord = system.x_image[icell_index][atom_index]
        atom_coord_new = atom_coord + trans_vec_cart

        new_index_tmp = None
        for icell, atom in enumerate(system.x_image):
            for k, dummy in enumerate(atom):
                coord = system.x_image[icell][k]
                if is_identical_position(atom_coord_new, coord):
                    new_index_tmp = k
        if new_index_tmp == None:
            sys.exit("ERROR: new_index_tmp is None")
        pair_1_list.append([origin_index_tmp, new_index_tmp])

# for t in vars(finder).items():
#	print(t["_symprec"])
