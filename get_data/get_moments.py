import numpy as np


def get_moments(moments_file, num_atom):
    '''
    get the information of magnetic moment from a dataset on magentic moments (MOMENTS file).
    :param moments_file: dataset of magnetic moments
    :param num_atom: number of atoms in the unit cell
    :return: np.array [energy#1, energy#2, ...] and
                      [[[magnetic moment, theta, phi, l, m, n]#atom1, [magnetic moment, ...]#atom2, ...]#1, ...]
            where l, m, and n mean direction cosines, #1, #2, ... means data set number.
    '''
    energies = np.empty(0)
    moments = np.empty(0)
    moments_list = []
    moments_on_a_pattern = np.empty(0)
    with open(moments_file, encoding='utf-8') as file:
        line_counter = 0
        for line in file:
            line = line.strip()

            # skip unnecessary line
            if not line:
                continue
            if line[0] == '#':
                continue

            if line_counter % (num_atom + 1) == 0:
                energy = float(line)
                energies = np.append(energies, energy)
            else:
                line = line.split()
                moment_magnitude = line[0]
                theta = line[1]
                phi = line[2]
                cosine_l = line[3]
                cosine_m = line[4]
                cosine_n = line[5]

                moment_on_an_atom = np.array([moment_magnitude,
                                              theta,
                                              phi,
                                              cosine_l,
                                              cosine_m,
                                              cosine_n])
                if moments_on_a_pattern.size == 0:
                    moments_on_a_pattern = np.append(moments_on_a_pattern, moment_on_an_atom)
                else:
                    moments_on_a_pattern = np.stack([moments_on_a_pattern, moment_on_an_atom])

                if line_counter % (num_atom + 1) == num_atom:
                    """
                    if moments.size == 0:
                        moments = np.array([moments_on_a_pattern])
                        print(moments)
                    else:
                        moments = np.stack([moments, np.array([moments_on_a_pattern])])
                    """
                    moments_list.append(moments_on_a_pattern)

                    # initialize moments_on_a_pattern
                    moments_on_a_pattern = np.empty(0)
            line_counter += 1

    # convert list to ndarray
    moments = np.array(moments_list)
    return energies, moments


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('moments', metavar='MOMENTS file',
                        help="Dataset file on magnetic moments, e.g., MOMENTS.")
    parser.add_argument('num_atom', metavar='number of atoms', type=int,
                        help="Number of atoms in the unit cell.")

    args = parser.parse_args()

    energies, moments = get_moments(args.moments, args.num_atom)
    # print(energies)
    print(moments)
