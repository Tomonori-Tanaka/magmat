import argparse
import re
import numpy as np
import sys

parser = argparse.ArgumentParser()

parser.add_argument('target_file', metavar='file_to_parse', type=str, nargs='+',
                    help="Output file of DFT codes, e.g., *.out.")

def run_parse(file_results):
    # pattern = '\s*\d+\s+\w+\s+(\d*\.?\d+(?:(?<!\.\d+)\.)?\s+){3}'
    # pattern = '^\s*\d+\s+\w+\s+(\d+\.?\d{5}s+){3}[+-]?(?:\d\.?\d{5})\s+\d+\.?\d{5}|\.\d+\s+[+-]?(?:\d+\.?\d*|\.\d+)'
    # pattern = '\s+\d+\s+\S+\s+(\d+\.\d{5}\s+){3}[+-]?\d+\.\d{5}\s+\d+\.\d{5}\s+[+-]?\d+\.\d{5}\\\\n'
    pattern_spin = '^\d+\s+\S+\s+(\d+\.\d{5}\s+){3}[+-]?\d+\.\d{5}\s+\d+\.\d{5}\s+[+-]?\d+\.\d{5}$'
    pattern_spin = re.compile(pattern_spin)
    pattern_etot = '^Utot\.\s+[+-]?\d+\.\d+'
    pattern_etot = re.compile(pattern_etot)
    pattern_ecs = re.compile('^Ucs\.\s+[+-]?\d+\.\d+')

    for i, file in enumerate(file_results):
        with open(file, encoding='utf-8') as f:
            print("# {}".format(i))
            for line in f:
                line = line.strip()
                # print(line.strip())

                if pattern_ecs.search(line):
                    line_ecs = line.split()
                    ecs = float(line_ecs[1])

                if pattern_etot.search(line):
                    line_etot = line.split()
                    etot = float(line_etot[1])
                    try:
                        etot = etot - ecs
                        print(etot)
                        del(ecs)
                    except:
                        sys.exit("Probably, program cannot catch Ucs. line.")

                if pattern_spin.search(line):
                    line_spin = line.split()
                    sdiff = float(line_spin[5])
                    theta = float(line_spin[6])
                    phi = float(line_spin[7])
                    l_direction_cos = np.sin(np.deg2rad(theta)) * np.cos(np.deg2rad(phi))
                    # l_direction_cos = round(l_direction_cos, 8)
                    m_direction_cos = np.sin(np.deg2rad(theta)) * np.sin(np.deg2rad(phi))
                    # m_direction_cos = round(m_direction_cos, 8)
                    n_direction_cos = np.cos(np.deg2rad(theta))
                    # n_direction_cos = round(n_direction_cos, 8)

                    print('{:15.9f}'.format(sdiff),
                          '{:20.9f}'.format(theta),
                          '{:15.9f}'.format(phi),
                          '{:20.9f}'.format(l_direction_cos),
                          '{:15.9f}'.format(m_direction_cos),
                          '{:15.9f}'.format(n_direction_cos))


if __name__ == '__main__':

    args = parser.parse_args()
    file_results = args.target_file
    run_parse(file_results)