"""
generate_conf.py
Generate various (random) spin configurations
"""
import argparse
import numpy as np
from interface.OpenMX import OpenmxParser

def check_code_options(args):
    conditions = [args.OpenMX is None]

    if conditions.count(True) == len(conditions):
        raise RuntimeError("Error : Either --VASP or --OpenMX option must be given.")
    elif len(conditions) - conditions.count(True) > 1:
        raise RuntimeError("Error : --VASP and "
                           "--OpenMX cannot be given simultaneously.")

    elif args.OpenMX:
        code = "OpenMX"
        struct_format = "OpenMX dat format"
        str_outfiles = "%s{counter}.dat" % args.prefix
        file_original = args.OpenMX

    return code, file_original, struct_format, str_outfiles

def get_code_object(code):
    if code == "VASP":
        pass
        # return VaspParser()

    if code == "OpenMX":
        return OpenmxParser()

def generate_spin_conf_openmx(codeobj, args):
    theta_min = args.theta_min
    theta_max = args.theta_max
    phi_min = args.phi_min
    phi_max = args.phi_max
    nat = codeobj.nat

    moment_directions = []
    for i in range(nat):
        # theta = np.arccos(-2*np.random.rand() + 1)
        # theta = np.rad2deg(theta)
        while True:
            theta = np.arccos(-2*np.random.rand() + 1)
            theta = np.rad2deg(theta)
            if theta >= theta_min and theta <= theta_max:
                break

        while True:
            phi = 2*np.pi*np.random.rand()
            phi = np.rad2deg(phi)
            if phi >= phi_min and phi <= phi_max:
                break
        moment_directions.append([theta, phi, theta, phi])

    return moment_directions
    # codeobj.moment_direction = moment_directions




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--prefix',
                        type=str, default="pattern",
                        help="Prefix of the files to be created. (default: pattern)")

    parser.add_argument('--OpenMX',
                        metavar='supercell.dat', type=str,
                        help="dat file")

    parser.add_argument('--theta_min',
                        type=float, default=.0,
                        help='minimum value of theta')

    parser.add_argument('--theta_max',
                        type=float, default=180.,
                        help='maxmum value of theta')

    parser.add_argument('--phi_min',
                        type=float, default=0.,
                        help='minimum value of phi')

    parser.add_argument('--phi_max',
                        type=float, default=360.,
                        help='maxmum value of phi')

    parser.add_argument('-nc',
                        type=int, default=1,
                        help='Specify the number of spin configuration patterns')

    parser.add_argument('-sn',
                       type=int, default=0,
                       help='Specify the start number of patterns')

    args = parser.parse_args()

    code, file_original, struct_format, str_outfiles = check_code_options(args)
    codeobj = get_code_object(code)
    codeobj.load_initial_structure(file_original)

    header_list = [i for i in range(args.sn, args.sn + args.nc)]
    for i in range(args.nc):
        moment_directions = generate_spin_conf_openmx(codeobj, args)
        # print(moment_directions)
        codeobj.generate_input(args.prefix, header_list[i], moment_directions)
