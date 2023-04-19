from matplotlib import pyplot as plt
import sys
sys.path.append("/Users/tomorin/PycharmProjects/magmat/magmat")
from get_moments import get_moments
import numpy as np

args = sys.argv
input_file = args[1]
nat_in = int(args[2])

energies, moments = get_moments(input_file, nat_in)
energies = energies - energies.mean()

index = []
for i, energy in enumerate(energies):
    index.append(i)
plt.bar(index, energies)
plt.show()

