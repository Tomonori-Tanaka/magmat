"""
generate_conf_vasp.py
Generate various (random) spin configurations for vasp INCAR file.
"""

import argparse
import numpy as np
import sys
import re



def generate_spin_conf_vasp(args):
    theta_min = args.theta_min
    theta_max = args.theta_max
    phi_min = args.phi_min
    phi_max = args.phi_max
    nat = args.nat
    moment = args.moment

    moment_directions = []
    for i in range(nat):
        # theta = np.arccos(-2*np.random.rand() + 1)
        # theta = np.rad2deg(theta)
        while True:
            theta = np.arccos(-2*np.random.rand() + 1)
            theta_deg = np.rad2deg(theta)
            if theta_deg >= theta_min and theta_deg <= theta_max:
                break

        while True:
            phi = 2*np.pi*np.random.rand()
            phi_deg = np.rad2deg(phi)
            if phi_deg >= phi_min and phi_deg <= phi_max:
                break

        moment_x = moment * np.sin(theta) * np.cos(phi)
        moment_y = moment * np.sin(theta) * np.sin(phi)
        moment_z = moment * np.cos(theta)

        # convert to str
        moment_x = '{:.3}'.format(moment_x)
        moment_y = '{:.3}'.format(moment_y)
        moment_z = '{:.3}'.format(moment_z)


        moment_directions.append(moment_x)
        moment_directions.append(moment_y)
        moment_directions.append(moment_z)

    # check
    if len(moment_directions) != 3*nat:
        sys.exit("ERROR: moment_directions is not 3*nat. len(moment_directions) = {}".format(len(moment_directions)))

    return moment_directions

def generate_input(infile, prefix, header, moment_directions):
    pattern_magmom = re.compile(r'^MAGMOM')
    pattern_mconst = re.compile(r'^M_CONSTR')

    outfile = prefix + header
    with open(outfile, mode='w', encoding='utf-8') as f_out:
        with open(infile, encoding='utf-8') as f_in:
            for line in f_in:
                line = line.strip()
                if bool(pattern_magmom.search(line)):
                    line = "MAGMOM = " + ' '.join(moment_directions) + "\n"
                if bool(pattern_mconst.search(line)):
                    line = "M_CONSTR = " + ' '.join(moment_directions) + "\n"

                f_out.write(line + '\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('nat',
                        type=int,
                        help='the number of atoms')

    parser.add_argument('--prefix',
                        type=str, default="pattern",
                        help="Prefix of the files to be created. (default: pattern)")

    parser.add_argument('--Vasp',
                        metavar='incar', type=str,
                        help="INCAR file")

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

    parser.add_argument('--moment',
                        type=float, default=1.0,
                        help='magnitude of magnetic moment')


    parser.add_argument('-nc',
                        type=int, default=1,
                        help='Specify the number of spin configuration patterns')

    parser.add_argument('-sn',
                       type=int, default=0,
                       help='Specify the start number of patterns')

    args = parser.parse_args()

    pattern_magmom = re.compile(r'^MAGMOM')
    pattern_mconst = re.compile(r'^M_CONSTR')

    header_list = [str(i) for i in range(args.sn, args.sn + args.nc)]

    for i in range(args.nc):
        moment_directions = generate_spin_conf_vasp(args)
        # print(moment_directions)
        generate_input(args.Vasp, args.prefix, header_list[i], moment_directions)


