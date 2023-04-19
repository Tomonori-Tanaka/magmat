import argparse
import re
import numpy as np
import sys
from itertools import islice
from itertools import tee

parser = argparse.ArgumentParser()

parser.add_argument('target_file', metavar='file_to_parse', type=str, nargs='+',
                    help="Output file of DFT codes, e.g., OUTCAR.")

parser.add_argument('--oszicar', type=str, nargs='+',
                    help='OSZICAR file.')

args = parser.parse_args()


def my_sign(x):
    value = (x > 0) - (x < 0)
    if value == 0:
        value = 1
    return value


def run_parse(file_results, oszicar):
    penalty_e_list = []
    e_p = None
    for file in oszicar:
        with open(file, encoding='utf-8') as f:
            for line in f:
                if "E_p" in line:
                    line = line.split()
                    e_p = float(line[2])
        if e_p == None:
            sys.exit("ERROR: penalty energy cannot be found.")
        else:
            penalty_e_list.append(e_p)


    line_num_x = 0
    line_num_y = 0
    line_num_z = 0
    toten = .0
    toten_list = []
    e_wo_entropy = .0
    e_wo_entropy_list = []
    e_sigma0 = .0
    e_sigma0_list = []
    nat = 0
    momx_listoflist = []
    momy_listoflist = []
    momz_listoflist = []

    for i, file in enumerate(file_results):
        # print("# {}".format(i))
        with open(file, encoding='utf-8') as f:
            for i, line in enumerate(f):
                if "magnetization (x)" in line:
                    line_num_x = i + 4
                elif "magnetization (y)" in line:
                    line_num_y = i + 4
                elif "magnetization (z)" in line:
                    line_num_z = i + 4
                elif "free  energy   TOTEN" in line:
                    line = line.split()
                    toten = float(line[4])
                elif "energy  without entropy=" in line:
                    line = line.split()
                    e_wo_entropy = float(line[3])
                    e_sigma0 = float(line[6])
                elif "number of ions     NIONS =" in line:
                    line = line.split()
                    nat = int(line[11])
        toten_list.append(toten)
        e_wo_entropy_list.append(e_wo_entropy)
        e_sigma0_list.append(e_sigma0)


        if line_num_x == line_num_y == line_num_z == 0:
            sys.exit("ERROR: line_num_{xyz} is zero.")

        with open(file, encoding='utf-8') as f:
            f1, f2, f3 = tee(f, 3)
            mom_x_list = []
            mom_y_list = []
            mom_z_list = []
            for line in islice(f1, line_num_x, line_num_x + nat):
                line = line.split()
                mom_x_list.append(float(line[4]))
            for line in islice(f2, line_num_y, line_num_y + nat):
                line = line.split()
                mom_y_list.append(float(line[4]))
            for line in islice(f3, line_num_z, line_num_z + nat):
                line = line.split()
                mom_z_list.append(float(line[4]))
        momx_listoflist.append(mom_x_list)
        momy_listoflist.append(mom_y_list)
        momz_listoflist.append(mom_z_list)

    for i in range(len(file_results)):
        print("# {}".format(i))
        print(e_sigma0_list[i]-penalty_e_list[i])

        for mom_x, mom_y, mom_z in zip(momx_listoflist[i], momy_listoflist[i], momz_listoflist[i]):
            magnitude = np.sqrt(mom_x ** 2 + mom_y ** 2 + mom_z ** 2)
            theta = np.rad2deg(np.arccos(mom_z / magnitude))
            if np.sqrt(mom_x ** 2 + mom_y ** 2) < 1e-10:
                phi = .0
            else:
                phi = np.rad2deg(my_sign(mom_y) * np.arccos(mom_x / np.sqrt(mom_x ** 2 + mom_y ** 2)))

            l_direction_cos = mom_x / magnitude
            m_direction_cos = mom_y / magnitude
            n_direction_cos = mom_z / magnitude

            print('{:15.9f}'.format(magnitude),
                  '{:20.9f}'.format(theta),
                  '{:15.9f}'.format(phi),
                  '{:20.9f}'.format(l_direction_cos),
                  '{:15.9f}'.format(m_direction_cos),
                  '{:15.9f}'.format(n_direction_cos))


if __name__ == "__main__":
    if len(args.target_file) != len(args.oszicar):
        print(len(args.target_file), "  ", len(args.oszicar))
        sys.exit("ERROR: the number of OUTCAR files and OSZICAR files.")

    run_parse(args.target_file, args.oszicar)
