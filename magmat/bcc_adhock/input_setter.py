import numpy as np
import sys

#import memory
import files
#import symmetry
#import optimize
#import constraint
#import patterndisp
import magmat_main
#import error


class InputSetter:
    def __init__(self) -> object:
        self.nat = 0
        self.nkd = 0
        self.kd = None
        self.lavec = np.zeros((3, 3))
        self.xcoords = None
        self.kdname = None
        self.is_periodic = [1, 1, 1]
        self.lspin = False
        self.magmom = None
        self.noncollinear = 0
        self.trevsym = 1
        self.str_magmom = ""
        self.maxorder = None
        self.nbody_include = None
        self.cutoff_radii = None

    def set_cell_parameter(self, a, lavec_in):
        for i in range(3):
            for j in range(3):
                self.lavec[i][j] = a * lavec_in[i][j]

    def set_interaction_vars(self, maxorder_in, nbody_include_in):
        self.maxorder = maxorder_in
        self.nbody_include = [None] * self.maxorder

        for i in range(self.maxorder):
            self.nbody_include[i] = nbody_include_in[i]

    def set_cutoff_radii(self, maxorder_in, nkd_in, cutoff_radii_in):
        if len(cutoff_radii_in) != nkd_in * nkd_in * maxorder_in:
            sys.exit("Incorrect size of the input array cutoff_radii_in")
        self.cutoff_radii = [None] * nkd_in * nkd_in * maxorder_in

        counter = 0
        for i in range(maxorder_in):
            for j in range(nkd_in):
                for k in range(nkd_in):
                    self.cutoff_radii[counter] = cutoff_radii_in[counter]
                    counter += 1

    def set_general_vars(self,
                         alm,
                         prefix,
                         mode,
                         verbosity,
                         str_disp_basis,
                         str_magmom,
                         nat_in,
                         nkd_in,
                         printsymmetery,
                         is_periodic_in,
                         trim_dispsign_for_evenfunc,
                         lspin_in,
                         print_hessian,
                         print_fcs_alamode,
                         print_fc3_shengbte,
                         print_fc2_qefc,
                         noncollinear_in,
                         trevsym_in,
                         kdname_in,
                         magmom_in,
                         tolerance,
                         tolerance_constraint,
                         basis_force_constant,
                         nmaxsave,
                         fc_zero_threshold):

        alm.set_output_filename_prefix(prefix)
        alm.set_verbosity(verbosity)
        self.nat = nat_in
        self.nkd = nkd_in
        alm.set_print_symmetry(printsymmetery)
        alm.set_symmetry_tolerance(tolerance)

        self.kdname = kdname_in
        self.magmom = magmom_in # nat x 3
        self.lspin = lspin_in
        self.noncollinear = noncollinear_in
        self.trevsym = trevsym_in
        self.is_periodic = is_periodic_in

        alm.set_fcs_save_flag("hessian", print_hessian)
        alm.set_fcs_save_flag("alamode", print_fcs_alamode)
        alm.set_fcs_save_flag("shengbte", print_fc3_shengbte)
        alm.set_fcs_save_flag("qefc", print_fc2_qefc)
        alm.set_fc_zero_threshold(fc_zero_threshold)
        alm.set_tolerance_constraint(tolerance_constraint)
        alm.set_forceconstant_basis(basis_force_constant)
        alm.set_nmaxsave(nmaxsave)

        if mode == "suggest":
            alm.set_displacement_basis(str_disp_basis)
            alm.set_displacement_param(trim_dispsign_for_evenfunc)

    def define(self, alm):
        alm.define(self.maxorder, self.nkd, self.nbody_include, self.cutoff_radii)

    def set_optimize_vars(self,
                          alm,
                          u_train_in,
                          f_train_in,
                          u_validation_in,
                          f_validation_in,
                          optcontrol_in):
        alm.set_u_train(u_train_in)
        alm.set_f_train(f_train_in)
        alm.set_validation_data(u_validation_in, f_validation_in)
        alm.set_optimizer_control(optcontrol_in)

    def set_file_vars(self, alm, datfile_train_in, datfile_validation_in):
        alm.set_datfile_train(datfile_train_in)
        alm.set_datfile_validation(datfile_validation_in)

    def set_constraint_vars(self,
                            alm,
                            constraint_flag,
                            rotation_axis,
                            fc2_file,
                            fc3_file,
                            fix_harmonic,
                            fix_cubic):
        alm.set_constraint_mode(constraint_flag)
        alm.set_rotation_axis(rotation_axis)
        alm.set_fc_file(2, fc2_file)
        alm.set_fc_fix(2, fix_harmonic)
        alm.set_fc_file(3, fc3_file)
        alm.set_fc_fix(3, fix_cubic)
        use_algebraic_constraint = constraint_flag / 10
        alm.set_algebrainc_constraint(use_algebraic_constraint)

    def set_atomic_positions(self, nat_in, kd_in, xcoord_in):
        self.kd = [None] * nat_in
        self.xcoord = np.zeros((nat_in, 3))
        self.xcoord.fill(np.nan)

        for i in range(nat_in):
            self.kd[i] = kd_in[i]
            for j in range(3):
                self.xcoord[i][j] = xcoord_in[i][j]

    def set_geometric_structure(self, alm):
        alm.set_cell(self.nat, self.lavec, self.xcoord, self.kd, self.kdname)
        alm.set_periodicity(self.is_periodic)
        alm.set_magnetic_params(self.nat, self.magmom, self.lspin, self.noncollinear, self.trevsym, self.str_magmom)


