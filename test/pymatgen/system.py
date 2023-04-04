'''
system.py
'''
from enum import Enum

import numpy as np
import sys
import libs.mathfunctions as mathfunctions
eps12 = 1e-12


class AtomType:
    def __init__(self):
        self.element = None
        self.magmom = None


class Cell:
    def __init__(self):
        self.lattice_vector = np.zeros((3, 3))  # ndarray((3, 3))
        self.reciprocal_lattice_vector = np.zeros((3, 3))
        # self.set_reciprocal_latt()
        self.volume = None
        self.number_of_atoms = None
        self.number_of_elems = None  # KD
        self.kind = None  # ? element number 1, 2
        self.x_fractional = None  # ? ndarray((number_of_atoms, 3))
        self.x_cartesian = None  # ? ndarray((number_of_atoms, 3))

    """
    def set_reciprocal_latt(self):
        det = np.linalg.det(self.lattice_vector)

        if (abs(det) < eps12):
            sys.exit("set_reciprocal_latt", "Lattice Vector is singular")

        factor = 2.0 * np.pi / det
        self.reciprocal_lattice_vector[0] = np.cross(self.lattice_vector[1], self.lattice_vector[2]) * factor
        self.reciprocal_lattice_vector[1] = np.cross(self.lattice_vector[2], self.lattice_vector[0]) * factor
        self.reciprocal_lattice_vector[2] = np.cross(self.lattice_vector[0], self.lattice_vector[1]) * factor
    """


class Spin:
    def __init__(self):
        lspin = False
        time_reversal = 1
        noncollinear = 0
        magmom = None  # ndarray((number_of_atoms, 3))


class System:
    def __init__(self):

        # Variables for geometric structure
        self.supercell = Cell()
        self.supercell.number_of_atoms = 0
        self.supercell.number_of_elems = 0

        self.kdname = []
        self.is_periodic = [1, 1, 1]  # is_periodic[3]
        self.x_image = None
        self.exist_image = None

        # Variables for spins
        self.spin = Spin()
        self.spin.lspin = False
        self.spin.noncollinear = 0
        self.spin.time_reversal_symm = 1
        self.str_magmom = ""

        self.set_atomtype_group()

        class LatticeType(Enum):
            Direct = 1
            Reciprocal = 2

    def initialize(self):
        """
        This method corresponds to "init" method in alm (system.cpp).
        """
        self.nat = self.supercell.number_of_atoms

        self.set_atomtype_group()

        nneib = 27
        self.x_image = np.zeros((nneib, self.nat, 3))
        self.x_image.fill(np.nan)

        self.exist_image = [None] * nneib

        self.generate_coordinate_of_periodic_images()

    def set_supercell(self, lavec_in, nat_in, kind_in, xf_in):
        """
        lavec_in: 3x3 float
        nat_in: unsigned int
        kind_in: ? list[nat_in] int
        xf_in: nat_in x 3 float
        """
        unique_nums = [0] * nat_in
        wrong_number = False

        nkd = 0
        for i in range(nat_in):
            in_unique_nums = False
            for j in range(nkd):
                if unique_nums[j] == kind_in[i]:
                    in_unique_nums = True
                    break
            if not in_unique_nums:
                unique_nums[nkd] = kind_in[i]
                nkd += 1

        for i in range(nkd):
            if unique_nums[i] > nkd:
                print(" WARNING : integers assigned to atoms are wrong. \n",
                      " The numbers will be resorted.")
                wrong_number = True
                break

        for i in range(3):
            for j in range(3):
                self.supercell.lattice_vector[i][j] = lavec_in[i][j]

        self.set_reciprocal_latt(self.supercell.lattice_vector, self.supercell.reciprocal_lattice_vector)
        self.supercell.volume = self.volume(self.supercell.lattice_vector, "Direct")
        self.supercell.number_of_atoms = nat_in
        self.supercell.number_of_elems = nkd
        self.supercell.kind = [None] * nat_in
        self.supercell.x_fractional = np.zeros((nat_in, 3))
        self.supercell.x_fractional.fill(np.nan)
        self.supercell.x_cartesian = np.zeros((nat_in, 3))
        self.supercell.x_fractional.fill(np.nan)

        if not wrong_number:
            for i in range(nat_in):
                self.supercell.kind.append(kind_in[i])
        else:
            for i in range(nat_in):
                for j in range(nkd):
                    if kind_in[i] == unique_nums[j]:
                        self.supercell.kind.append(j + 1)

        xtmp = [.0] * 3
        for i in range(nat_in):
            for j in range(3):
                xtmp[j] = xf_in[i][j]
            # The fractional coordinate should be in the range of 0<=xf<1
            for j in range(3):
                while xtmp[j] >= 1.0:
                    xtmp[j] -= 1.0
                while xtmp[j] < 0.0:
                    xtmp[j] += 1.0
            self.supercell.x_fractional[i] = xtmp

        xf_tmp = [None] * 3
        xc_tmp = [None] * 3
        for i, xf in enumerate(self.supercell.x_fractional):
            for j in range(3):
                xf_tmp[j] = xf[j]
            xc_tmp = mathfunctions.rotvec(xf_tmp, self.supercell.lattice_vector)
            for j in range(3):
                xtmp[j] = xc_tmp[j]
            self.supercell.x_cartesian[i] = xtmp

        self.spin.magmom = []
        vec = [None] * 3
        for i in range(nat_in):
            for j in range(3):
                vec[j] = 0
            self.spin.magmom.append(vec)

    def get_supercell(self):
        return self.supercell

    def get_x_image(self):
        return self.x_image

    def get_exist_image(self):
        return self.exist_image

    def set_periodicity(self, is_periodic_in):
        """
        is_periodic_in: list(3)
        """
        for i in range(3):
            self.is_periodic[i] = is_periodic_in[i]

    def get_periodicity(self):
        return self.is_periodic

    def set_kdname(self, kdname_in):
        """
        kdname_in: str list 1xnumber_of_elems
        """
        nkd = self.supercell.number_of_elems

        for i in range(nkd):
            self.kdname[i] = kdname_in[i]

    def get_kdname(self):
        return self.kdname

    def set_reciprocal_latt(self, aa, bb):
        """
        Calculate Reciprocal Lattice Vectors

        Here, BB is just the inverse matrix of AA (multiplied by factor 2 Pi)
        BB = 2 Pi AA^{-1},
        = t(b1, b2, b3)

        (b11 b12 b13)
        = (b21 b22 b23)
        (b31 b32 b33),

        b1 = t(b11, b12, b13) etc.

        aa: lattice_vector ndarray 3x3
        bb: reciplocal_lattice_vector ndarray 3x3
        """

        det = np.linalg.det(aa)

        if det < eps12:
            sys.exit("Lattice Vector is singular")

        factor = 2.0 * np.pi / det

        bb[0][0] = (aa[1][1] * aa[2][2] - aa[1][2] * aa[2][1]) * factor
        bb[0][1] = (aa[0][2] * aa[2][1] - aa[0][1] * aa[2][2]) * factor
        bb[0][2] = (aa[0][1] * aa[1][2] - aa[0][2] * aa[1][1]) * factor

        bb[1][0] = (aa[1][2] * aa[2][0] - aa[1][0] * aa[2][2]) * factor
        bb[1][1] = (aa[0][0] * aa[2][2] - aa[0][2] * aa[2][0]) * factor
        bb[1][2] = (aa[0][2] * aa[1][0] - aa[0][0] * aa[1][2]) * factor

        bb[2][0] = (aa[1][0] * aa[2][1] - aa[1][1] * aa[2][0]) * factor
        bb[2][1] = (aa[0][1] * aa[2][0] - aa[0][0] * aa[2][1]) * factor
        bb[2][2] = (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]) * factor

    def frac2cart(self, xf):
        """
        x_cartesian = A x_fractional
        """
        x_tmp = [None] * 3

        for i in range(self.supercell.number_of_atoms):
            x_tmp = mathfunctions.rotvec(xf[i], self.supercell.lattice_vector)

            for j in range(3):
                xf[i][j] = x_tmp[j]

    def volume(self, latt_in, type):
        """
        latt_in: lattice of input. 3x3 list
        type: LatticeType enum
        return: volume of the cell
        """
        mat = np.zeros((3, 3))

        if type == "Direct":
            for i in range(3):
                for j in range(3):
                    mat[i][j] = latt_in[j][i]
        elif type == "Reciprocal":
            for i in range(3):
                for j in range(3):
                    mat[i][j] = latt_in[i][j]
        else:
            sys.exit("Invalid LatticeType is given")

        vol = abs(mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
                  + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2])
                  + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]))

        return vol

    def set_spin_variables(self, nat_in, lspin_in, noncol_in, trev_sym_in, magmom_in):
        self.spin.lspin = lspin_in
        self.spin.noncollinear = noncol_in
        self.spin.time_reversal_symm = trev_sym_in
        self.spin.magmom = np.zeros((nat_in, 3))
        self.spin.magmom = self.spin.magmom.fill(np.nan)

        for i in range(nat_in):
            self.spin.magmom[i] = magmom_in[i]

    def get_spin(self):
        return self.spin

    def set_str_magmom(self, str_magmom_in):
        self.str_magmom = str_magmom_in

    def get_str_magmom(self):
        return self.str_magmom

    def get_atomtype_group(self):
        return self.atomtype_group

    def set_atomtype_group(self):
        type_tmp = AtomType()
        set_type = set()

        for i in range(self.supercell.number_of_atoms):
            type_tmp.element = self.supercell.kind[i]

            if self.spin.noncollinear == 0:
                type_tmp.magmom = self.spin.magmom[i][2]
            else:
                type_tmp.magmom = 0  # maybe this member is not used in noncollinear system.

            set_type.add(type_tmp)

        natomtypes = len(set_type)
        atomtype_group = [None] * natomtypes

        for i in range(self.supercell.number_of_atoms):
            count = 0
            for it in set_type:
                if self.spin.noncollinear:
                    if self.supercell.kind[i] == it.element:
                        if atomtype_group[count] == None:
                            atomtype_group[count] == [i]
                        else:
                            atomtype_group[count].append(i)
                count += 1

    def generate_coordinate_of_periodic_images(self):
        """
        Generate Cartesian coordinates of atoms in the neighboring 27 supercells
        """
        nat = self.supercell.number_of_atoms
        xf_in = self.supercell.x_fractional  # x(means coordinates) of atoms in fractional. "in" means input.

        icell = 0  # icell means index of imaginary (virtual) neighboring cell. 0 means centrally located supercell.
        for i in range(nat):
            for j in range(3):
                self.x_image[0][i][j] = xf_in[i][j]

        # Convert to Cartesian coordinate
        self.frac2cart(self.x_image[0])

        for ia in range(-1, 2):
            for ja in range(-1, 2):
                for ka in range(-1, 2):
                    if ia == 0 and ja == 0 and ka == 0:
                        continue

                    icell += 1
                    for i in range(nat):
                        self.x_image[icell][i][0] = xf_in[i][0] + float(ia)
                        self.x_image[icell][i][1] = xf_in[i][1] + float(ja)
                        self.x_image[icell][i][2] = xf_in[i][2] + float(ka)

                    # Convert to Cartesian coordinate
                    self.frac2cart(self.x_image[icell])

        icell = 0
        self.exist_image[0] = 1

        for ia in range(-1, 2):
            for ja in range(-1, 2):
                for ka in range(-1, 2):
                    if ia == 0 and ja == 0 and ka == 0:
                        continue

                    icell += 1
                    # When periodic flag is zero along an axis,
                    # periodic images along that axis cannot be considered.
                    if ((abs(ia) == 1 and self.is_periodic[0] == 0) or
                            (abs(ja) == 1 and self.is_periodic[1] == 0) or
                            (abs(ia) == 2 and self.is_periodic[2] == 0)):
                        self.exist_image[icell] = 0
                    else:
                        self.exist_image[icell] = 1

    def print_structure_stdout(self, cell):
        print(" SYSTEM")
        print(" ======")

        print("  Lattice Vector")
        print("{:5.10e}".format(cell.lattice_vector[0][0]),
              "{:5.10e}".format(cell.lattice_vector[1][0]),
              "{:5.10e}".format(cell.lattice_vector[2][0]), " : a1")

        print("{:5.10e}".format(cell.lattice_vector[0][1]),
              "{:5.10e}".format(cell.lattice_vector[1][1]),
              "{:5.10e}".format(cell.lattice_vector[2][1]), " : a2")

        print("{:5.10e}".format(cell.lattice_vector[0][2]),
              "{:5.10e}".format(cell.lattice_vector[1][2]),
              "{:5.10e}".format(cell.lattice_vector[2][2]), " : a3\n")

        print("  Cell volume = ", cell.volume, "\n")

        print("  Reciprocal Lattice Vector")
        print("{:5.10e}".format(self.supercell.reciprocal_lattice_vector[0][0]),
              "{:5.10e}".format(self.supercell.reciprocal_lattice_vector[0][1]),
              "{:5.10e}".format(self.supercell.reciprocal_lattice_vector[0][2]), " : b1")

        print("{:5.10e}".format(self.supercell.reciprocal_lattice_vector[1][0]),
              "{:5.10e}".format(self.supercell.reciprocal_lattice_vector[1][1]),
              "{:5.10e}".format(self.supercell.reciprocal_lattice_vector[1][2]), " : b2")

        print("{:5.10e}".format(self.supercell.reciprocal_lattice_vector[2][0]),
              "{:5.10e}".format(self.supercell.reciprocal_lattice_vector[2][1]),
              "{:5.10e}".format(self.supercell.reciprocal_lattice_vector[2][2]), " : b3\n")

        print("  Atomic species:")
        for i in range(cell.number_of_elems):
            print("{:6}".format(i + 1), "{:5}".format(self.kdname[i]))
        print("")

        print("  Atomic positions in fractional basis and atomic species")
        for i in range(cell.number_of_atoms):
            print("{:6}".format(i + 1))
            print("{:15}".format(cell.x_fractional[i][0]),
                  "{:15}".format(cell.x_fractional[i][1]),
                  "{:15}".format(cell.x_fractional[i][2]), cell.kind[i])
        print("\n")

    def print_magmom_stdout(self):
        print("  MAGMOM is given. The magnetic moments of each atom are as follows:")
        for i in range(self.supercell.number_of_atoms):
            print("{:6}".format(i + 1), "{:5}".format(self.spin.magmom[i][0]),
                  "{:5}".format(self.spin.magmom[i][1]),
                  "{:5}".format(self.spin.magmom[i][2]))
        print("")

        if self.spin.noncollinear == 0:
            print("  NONCOLLINEAR = 0: magnetic moments are considered as scalar variables.")
        elif self.spin.noncollinear == 1:
            print("  NONCOLLINEAR = 1: magnetic moments are considered as vector variables.")
            if self.spin.time_reversal_symm:
                print("  TREVSYM = 1: Time-reversal symmetry will be considered for generating magnetic space group")
            else:
                print(
                    "  TREVSYM = 0: Time-reversal symmetry will NOT be considered for generating magnetic space group")
        print("\n")


if __name__ == '__main__':
    latvec = np.array([[0.5, 0.5, 0.], [0., 0.5, 0.5], [0.5, 0., 0.5]])
    cell_test = Cell()
    # print(cell_test.reciprocal_lattice_vector)
    # print(cell_test.lattice_vector)

    system_test = System()
