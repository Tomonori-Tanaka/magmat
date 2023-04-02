import numpy as np
import sys
import warnings
import magmat_main
import files
from input_setter import InputSetter


class InputParser:
    def __init__(self):
        self.input_setter = InputSetter()
        self.input_filename = ""
        self.kdname = []
        self.mode = ""
        self.maxorder = 0
        self.nat = 0
        self.nkd = 0

    def run(self, alm, nargs, arg):
        if nargs == 1:
            self.from_stdin = True
            self.input_filename = arg[0]
        else:
            self.from_stdin = False
            try:
                with open(arg[1]) as f:
                    pass
                self.input_filename = arg[1]
            except:
                sys.exit("No such file or directory: ")

        self.parse_input(alm)

    def get_run_mode(self):
        return self.mode

    def parse_displacement_and_force_files(self, u, f, datfile_in):
        nrequired = None
        if datfile_in.ndata == 0:
            nrequired = -1
        else:
            # Total number of data entries (displacement + force)
            nrequired = 6 * self.nat * datfile_in.ndata

        value_arr = []

        # Open the target file and copy the data to 1D temporary vector
        try:
            with open(datfile_in.filename, encoding='utf-8') as file:
                pass
        except:
            sys.exit("cannot open DFSET file")

        with open(datfile_in.filename, encoding='utf-8') as file:
            # reach_end = False
            for line in file:
                line = line.strip()
                if line[0] == "#":
                    continue
                elif not line:
                    continue
                line = line.split()
                line = map(lambda x: float(x), line)
                value_arr.extend(line)
                if len(value_arr) == nrequired:
                    reach_end = True
                    break

        # Check if the length of the vector is correct.
        # Also, estimate ndata if it is not set.
        n_entries = len(value_arr)
        if nrequired == -1:
            if n_entries % (6 * self.nat) == 0:
                datfile_in.ndata = n_entries / (6 * self.nat)
            else:
                sys.exit("The number of lines in DFSET is too small for the given NDATA = {0}".format(datfile_in.ndata))

        if datfile_in.nstart == 0:
            datfile_in.nstart = 1
        if datfile_in.nend == 0:
            datfile_in.nend = datfile_in.ndata

        # Copy the data into 2D array
        ndata_used = datfile_in.nend - datfile_in.nend - datfile_in.nstart + 1 - datfile_in.skip_e + datfile_in.skip_s

        u = np.zeros((ndata_used, 3 * self.nat))
        u.fill(np.nan)
        f = np.zeros((ndata_used, 3 * self.nat))
        f.fill(np.nan)

        idata = 0
        for i in range(datfile_in.ndata):
            if i < datfile_in.nstart - 1:
                continue
            elif i >= datfile_in.skip_s - 1 and i < datfile_in.skip_e - 1:
                continue
            elif i > datfile_in.nend - 1:
                break

            for j in range(self.nat):
                for k in range(3):
                    u[idata][3 * j + k] = value_arr[6 * self.nat * i + 6 * j + k]
                    f[idata][3 * j + k] = value_arr[6 * self.nat * i + 6 * j + k + 3]

            idata += 1

    def parse_input(self, alm):
        """
        The order of calling methods in this method is important.
        Since following methods rely on variables those already parsed.

        :param alm:
        :return:
        """

        if not self.locate_tag("&general"):
            sys.exit("&general entry not found in the input file")

        # kdname is allocated in this method.
        self.parse_general_vars(alm)

        if not self.locate_tag("&cell"):
            sys.exit("&cell entry not found in the input file")
        self.parse_cell_parameter()

        if not self.locate_tag("&position"):
            sys.exit("&position entry not found in the input file")
        self.parse_atomic_positions()
        self.input_setter.set_geometric_structure(alm)

        if not self.locate_tag("&interaction"):
            sys.exit("&interaction entry not found in the input file")
        self.parse_interaction_vars()

        if not self.locate_tag("&cutoff"):
            sys.exit("&cutoff entry not found in the input file")
        self.parse_cutoff_radii()
        self.input_setter.define(alm)

        if self.mode == "optimize":
            if not self.locate_tag("&optimize"):
                sys.exit("&optimize entry not found in the input file")
            self.parse_optimize_vars(alm)

    def parse_general_vars(self, alm):
        trim_dispsign_for_evenfunc = True
        print_hessian = 0
        print_fcs_alamode = 0
        print_fc2_qefc = 0
        print_fc3_chengbte = 0
        noncollinear = 0
        trevsym = 0
        str_disp_basis = "C"

        kdname_v = []
        periodic_v = []
        magmom_v = []
        str_split = []
        input_list = [
            "PREFIX", "MODE", "NAT", "NKD", "KD", "PERIODIC", "PRINTSYM", "TOLERANCE",
            "DBASIS", "TRIMEVEN", "VERBOSITY",
            "MAGMOM", "NONCOLLINEAR", "TREVSYM", "HESSIAN", "TOL_CONST", "FCSYM_BASIS",
            "NMAXSAVE", "FC3_SHENGBTE", "FC2_QEFC", "FCS_ALAMODE", "FC_ZERO_THR"
        ]
        no_defaults = ["PREFIX", "MODE", "NAT", "NKD", "KD"]
        general_var_dict = {}

        self.get_var_dict(input_list, general_var_dict)

        for it in no_defaults:
            if it in general_var_dict is False:
                sys.exit("The following variable is not found in &general input region: {}".format(it))

        prefix = general_var_dict["PREFIX"]
        mode = general_var_dict["MODE"]
        mode = mode.lower()
        if mode != "suggest" and self.mode != "optimize" and self.mode != "opt":
            sys.exit("Invalid MODE variable")
        if mode == "opt":
            mode = "optimize"
        nat = int(general_var_dict["NAT"])
        nkd = int(general_var_dict["NKD"])

        if not "VERBOSITY" in general_var_dict:
            verbosity = 1
        else:
            verbosity = int(general_var_dict["VERBOSITY"])

        if not "PRINTSYM" in general_var_dict:
            printsymmetry = 0
        else:
            printsymmetry = int(general_var_dict["PRINTSYM"])

        kdname_v = general_var_dict["KD"].split()
        if len(kdname_v) != nkd:
            sys.exit("input_parser.py: The number of entries for KD is inconsistent with NKD")
        else:
            kdname = [None] * nkd
            for i, kd in enumerate(kdname_v):
                kdname[i] = kd

        if not "PERIODIC" in general_var_dict:
            is_periodic = [1, 1, 1]
        else:
            periodic_v = general_var_dict["PERIODIC"].split()
            periodic_v = list(map(lambda x: int(x), periodic_v))
            if len(periodic_v) == 3:
                for i in range(3):
                    try:
                        is_periodic[i] = int(periodic_v[i])
                    except:
                        sys.exit("The PERIODIC tag must be a set of integers.")
            else:
                sys.exit("Invalid number of entries for PERIODIC")

        if not "TOLERANCE" in general_var_dict:
            tolerance = 1.0e-3
        else:
            tolerance = float(general_var_dict["TOLERANCE"])
        if not "TOL_CONST" in general_var_dict:
            tolerance_constraint = 1.0e-6
        else:
            tolerance_constraint = float(general_var_dict["TOL_CONST"])

        if not "FCSYM_BASIS" in general_var_dict:
            basis_force_constant = "Lattice"
        else:
            basis_force_constant = general_var_dict["FCSYM_BASIS"]
            basis_force_constant = basis_force_constant.lower()

            if basis_force_constant[0] != "c" and basis_force_constant[0] != "l":
                sys.exit("Invalid FCSYM_BASIS.")

        if not "NMAXSAVE" in general_var_dict:
            nmaxsave = 5
        else:
            nmaxsave = int(general_var_dict["NMAXSAVE"])

        if not "NONCOLLINEAR" in general_var_dict:
            noncollinear = 0
        else:
            noncollinear = int(general_var_dict["NONCOLLINEAR"])

        if not "TREVSYM" in general_var_dict:
            trevsym = 1
        else:
            trevsym = int(general_var_dict["TREVSYM"])

        if not "HESSIAN" in general_var_dict:
            print_hessian = 0
        else:
            print_hessian = int(general_var_dict["HESSIAN"])

        if not "FCS_ALAMODE" in general_var_dict:
            print_fcs_alamode = 1
        else:
            print_fcs_alamode = int(general_var_dict["FCS_ALAMODE"])

        if not "FC3_SHENGBTE" in general_var_dict:
            print_fc3_shengbte = 0
        else:
            print_fc3_shengbte = int(general_var_dict["FC3_SHENGBTE"])

        if not "FC2_QEFC" in general_var_dict:
            print_fc2_qefc = 0
        else:
            print_fc2_qefc = int(general_var_dict["FC2_QEFC"])

        if not "FC_ZERO_THR" in general_var_dict:
            fc_zero_threshold = 1.0e-12
        else:
            fc_zero_threshold = int(general_var_dict["FC_ZERO_THR"])

        # Convert MAGMOM input to array
        magmom = np.zeros((nat, 3))
        lspin = False

        if "MAGMOM" in general_var_dict:
            lspin = True
            if noncollinear:
                icount = 0
                self.split_str_by_space(general_var_dict["MAGMOM"], magmom_v)
                for it in magmom_v:
                    if "*" in it:
                        sys.exit("Wild card '*' is not supported when NONCOLLINEAR = 1.")
                    else:
                        magmag = float(it)
                        if icount // 3 >= nat:
                            sys.exit("Too many entries for MAGMOM.")
                        magmom[icount // 3][icount % 3] = magmag
                        icount += 1
                if icount != 3 * nat:
                    sys.exit("Number of entries for MAGMOM must be 3*NAT when NONCOLLINEAR = 1.")
            else:
                icount = 0
                magmom_v = general_var_dict["MAGMOM"].split()
                for it in magmom_v:
                    if "*" in it:
                        if it == "*":
                            sys.exit("Please place '*' without space for the MAGMOM-tag.")
                        str_split = it.split('*')
                        if len(str_split) != 2:
                            sys.exit("Invalid format for the MAGMOM-tag.")
                        else:
                            if not str_split[0] or not str_split[1]:
                                sys.exit("Please place '*' without space for the MAGMOM-tag.")
                            ncount = 0
                            try:
                                ncount = int(str_split[0])
                                magmag = float(str_split[1])
                            except:
                                sys.exit("Bad format for MAGMOM.")

                            for i in range(icount, icount + ncount):
                                magmom[i][2] = magmag
                            icount += ncount

                    else:
                        magmag = float(it)
                        if icount == nat:
                            icount = 0
                            break
                        magmom[icount][2] = magmag
                        icount += 1

                if icount != nat:
                    sys.exit("Number of entries for MAGMOM must be NAT.")

        if mode == "suggest":
            if not "DBASIS" in general_var_dict:
                str_disp_basis = "Cart"
            else:
                str_disp_basis = general_var_dict["DBASIS"]
            str_disp_basis = str_disp_basis.upper()
            if str_disp_basis[0] != 'C' and str_disp_basis[1] != 'F':
                sys.exit("Invalid DBASIS")

            if "TRIMEVEN" in general_var_dict:
                trim_dispsign_for_evenfunc = None
                trim_dispsign_for_evenfunc = int(general_var_dict["TRIMEVEN"])

        self.input_setter.set_general_vars(alm,
                                           prefix,
                                           mode,
                                           verbosity,
                                           str_disp_basis,
                                           magmom.flatten(),
                                           nat,
                                           nkd,
                                           printsymmetry,
                                           is_periodic,
                                           trim_dispsign_for_evenfunc,
                                           lspin,
                                           print_hessian,
                                           print_fcs_alamode,
                                           print_fc3_shengbte,
                                           print_fc2_qefc,
                                           noncollinear,
                                           trevsym,
                                           kdname,
                                           magmom,
                                           tolerance,
                                           tolerance_constraint,
                                           basis_force_constant,
                                           nmaxsave,
                                           fc_zero_threshold
                                           )

    def parse_cell_parameter(self):
        a = .0
        lavec_tmp = np.zeros((3, 3))
        line = ""
        line_wo_comment = ""
        line_tmp = ""
        line_vec = []
        line_split = []
        pos_fist_comment_tag = None

        line_cell = None

        # specify &cell entry
        with open(self.input_filename, encoding='utf-8') as f:
            icount = 0
            for line in f:
                if "&cell" in line:
                    line_cell = icount
                    break
                icount += 1

        with open(self.input_filename, encoding='utf-8') as f:
            icount = 0
            for line in f:
                if icount <= line_cell:
                    continue
                line = line.strip()
                if line == "":
                    continue
                elif line[0] == "#":
                    continue
                elif self.is_endof_entry(line):
                    break
                line_vec.append(line)

        if len(line_vec) != 4:
            sys.exit(
                "Too few or too much lines for the &cell field.\n"
                "The number of valid lines for the &cell field should be 4.")

        for i, line in enumerate(line_vec):
            line_split = line.split()
            if i == 0:
                # Lattice factor a
                if len(line_split) == 1:
                    a = float(line_split[0])
                else:
                    sys.exit("Unacceptable format for &cell field.")
            else:
                # Lattice vectors a1, a2, a3
                if len(line_split) == 3:
                    for j in range(3):
                        lavec_tmp[j][i - 1] = float(line_split[j])
                else:
                    sys.exit("Unacceptable format for &cell field.")

        self.input_setter.set_cell_parameter(a, lavec_tmp)

    def parse_interaction_vars(self):
        nbody_included = []
        nbody_v = []
        input_list = ["NORDER", "NBODY"]
        no_defaults = ["NORDER"]
        interaction_var_dict = {}
        maxorder = 0

        self.get_var_dict(input_list, interaction_var_dict)

        for it in no_defaults:
            if it in interaction_var_dict:
                sys.exit("The following variable is not found in &interaction input region: ")

        self.maxorder = int(interaction_var_dict["NORDER"])
        if maxorder < 1:
            sys.exit("maxorder has to be a positive integer")

        nbody_included = [None] * maxorder

        nbody_v = interaction_var_dict["NBODY"].split()

        if not nbody_v[0]:
            for i in range(maxorder):
                nbody_included[i] = i + 2
        elif len(nbody_v) == maxorder:
            for i in range(maxorder):
                try:
                    nbody_included[i] = int(nbody_v[i])
                except:
                    sys.exit("NBODY must be an integer.")
        else:
            sys.exit("The number of entry of NBODY has to be equal to NORDER")

        if nbody_included[0] != 2:
            warnings.warn("Harmonic interaction is always 2 body (except on-site 1 body)")

        self.input_setter.set_interaction_vars(maxorder, nbody_included)

    """
    def parse_optimize_vars(self, alm):
        # This method must not be called before setting NAT by parse_general_vars.
    
        constraint_flag = None
        flag_sparse = 0
        rotation_axis = ""
        
        optcontrol = optimiez
    """

    def parse_atomic_positions(self):
        str_v = []

        # specify &positions entry
        line_positions = None
        with open(self.input_filename, encoding='utf-8') as f:
            icount = 0
            for line in f:
                if "&positions" in line:
                    line_positions = icount
                    break
                icount += 1

        with open(self.input_filename, encoding='utf-8') as f:
            icount = 0
            for line in f:
                if icount <= line_positions:
                    continue
                line = line.strip()
                if line[0] == "#":
                    continue
                if line == "":
                    continue
                if self.is_endof_entry(line):
                    break
                str_v.append(line)

        if len(str_v) != self.nat:
            sys.exit("The number of entries for atomic positions should be NAT")

        xeq = [None] * self.nat
        kd = [None] * self.nat

        for i in range(self.nat):
            pos_line = str_v[i].split()
            if len(pos_line) == 4:
                try:
                    kd[i] = int(pos_line[0])
                except:
                    sys.exit("Invalid entry for the &position field at line {}".format(i + 1))
                for j in range(3):
                    xeq[i][j] = float(pos_line[j + 1])
            else:
                sys.exit("Bad format for &position region")

        self.input_setter.set_atomic_positions(self.nat, kd, xeq)

    def parse_cutoff_radii(self):
        line_cutoff = None
        str_cutoff = []
        with open(self.input_filename, encoding='utf-8') as f:
            icount = 0
            for line in f:
                if "&positions" in line:
                    line_cutoff = icount
                    break
                icount += 1

        with open(self.input_filename, encoding='utf-8') as f:
            icount = 0
            for line in f:
                if icount <= line_cutoff:
                    continue
                line = line.strip()
                if line == "":
                    continue
                elif line[0] == "#":
                    continue
                elif self.is_endof_entry(line):
                    break
                str_cutoff.append(line)

        order = 0
        element_allowed = set()
        str_pair = []
        kd_map = {}

        cutoff_tmp = .0
        cutoff_radii_tmp = .0

        undefined_cutoff = np.zeros((self.maxorder, self.nkd, self.nkd))
        undefined_cutoff.fill(True)

        cutoff_radii_tmp = np.zeros((self.maxorder, self.nkd, self.nkd))

        for i in range(self.nkd):
            element_allowed.add(self.kdname[i])
            kd_map[self.kdname[i]] = i

        element_allowed.add("*")
        kd_map["*"] = -1

        for it in str_cutoff:
            cutoff_line = it.split()
            if len(cutoff_line) < self.maxorder + 1:
                sys.exit("Invalid format for &cutoff entry")

            str_pair = cutoff_line.split("-")

            if len(str_pair) != 2:
                sys.exit("Invalid format for &cutoff entry")

            for i in range(2):
                if not str_pair[i] in element_allowed:
                    sys.exit("Invalid format for &cutoff entry")

            ikd = kd_map[str_pair[0]]
            jkd = kd_map[str_pair[1]]

            for order in range(self.maxorder):
                # Accept any strings starting with 'N' or 'n' as 'None'
                if cutoff_line[order + 1][0] == 'N' or cutoff_line[order + 1][0] == 'n':
                    # Minus value for cutoff radius
                    # This is a flag for neglecting cutoff radius
                    cutoff_tmp = -1.0
                else:
                    cutoff_tmp = float(cutoff_line[order + 1])

                if ikd == -1 and jkd == -1:
                    for i in range(self.nkd):
                        for j in range(self.nkd):
                            cutoff_radii_tmp[order][i][j] = cutoff_tmp
                            undefined_cutoff[order][i][j] = False
                elif ikd == -1:
                    for i in range(self.nkd):
                        cutoff_radii_tmp[order][i][jkd] = cutoff_tmp
                        cutoff_radii_tmp[order][jkd][i] = cutoff_tmp
                        undefined_cutoff[order][i][jkd] = False
                        undefined_cutoff[order][jkd][i] = False
                elif jkd == -1:
                    for j in range(self.nkd):
                        cutoff_radii_tmp[order][j][ikd] = cutoff_tmp
                        cutoff_radii_tmp[order][ikd][j] = cutoff_tmp
                        undefined_cutoff[order][j][ikd] = False
                        undefined_cutoff[order][ikd][j] = False

        for order in range(self.maxorder):
            for j in range(self.nkd):
                for k in range(self.nkd):
                    if undefined_cutoff[order][j][k]:
                        print(" Cutoff radius for ", order + 2, "th-order terms\n",
                              " are not defined between elements ", j + 1, " and ", k + 1)
                        sys.exit("Incomplete cutoff radii")

        cutoff_information_flatten = [None] * (self.maxorder * self.nkd * self.nkd)

        i = 0
        for order in range(self.maxorder):
            for j in range(self.nkd):
                for k in range(self.nkd):
                    cutoff_information_flatten[i] = cutoff_radii_tmp[order][j][k]
                    i += 1

        self.input_setter.set_cutoff_radii(self.maxorder, self.nkd, cutoff_information_flatten)

    def get_var_dict(self, input_list, var_dict):
        keyword_set = set()
        for it in input_list:
            keyword_set.add(it)

        with open(self.input_filename, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                elif line[0] == '#':
                    continue
                elif not '=' in line:
                    continue
                str_entry = line.split(';')

                for it in str_entry:
                    str_tmp = it.strip()
                    if str_tmp:
                        str_varval = str_tmp.split("=")
                        if len(str_varval) != 2:
                            print("Failed to parse")
                            sys.exit("Unacceptable format")
                        key = (str_varval[0].strip()).upper()
                        val = str_varval[1].strip()

                        if key in input_list:
                            var_dict[key] = val

    def is_data_range_consistent(self, datfile_in):
        ndata = datfile_in.ndata
        nstart = datfile_in.nstart
        nend = datfile_in.nend
        skip_s = datfile_in.skip_s
        skip_e = datfile_in.skip_e

        if nstart > 0 and skip_s > 0 and skip_s < nstart:
            return False
        elif nend > 0 and skip_e > 0 and (skip_e - 1) > nend:
            return False
        elif nstart > 0 and nend > 0 and nstart > nend:
            return False
        elif skip_s > 0 and skip_e > 0 and skip_s >= skip_e:
            return False

        if ndata > 0:
            if nstart > 0 and nstart > ndata:
                return False
            if nend > 0 and nend > ndata:
                return False
            if skip_s > 0 and skip_s > ndata:
                return False
            if skip_e > 0 and (skip_e - 1) > ndata:
                return False

        return True

    def locate_tag(self, key):
        ret = 0
        line = ""
        with open(self.input_filename, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line == key:
                    ret = 1
                    break
        return ret

    def is_endof_entry(self, str):
        return str[0] == '/'

    def split_str_by_space(self, str, str_vec):
        str_vec = str.split()

    def assign_val(self, val, key, dict):
        try:
            val = dict[key]
        except:
            sys.exit("Error: assign_val")
