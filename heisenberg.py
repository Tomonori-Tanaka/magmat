import numpy as np
import pymc as pm
import matplotlib.pyplot as plt
import re
from get_data.get_moments import get_moments

# https://qiita.com/takubb/items/a2e02fe61aad80f8f2f3

num_atoms = 2
moments_file = "MOMENTS_small"
energies, moments = get_moments(moments_file, num_atoms)


energies = energies - np.mean(energies)
x = np.array([])
for pattern in moments:
    vec_a = pattern[0][3:6]
    vec_b = pattern[1][3:6]
    x = np.append(x, np.dot(vec_a, vec_b))

energies = energies*27.2114*1000/2


with pm.Model() as model:
    jij = pm.Normal('jij', mu=0, sigma=1000)
    # noise = pm.Normal('noise', mu=0, sigma=1)
    noise = pm.HalfFlat('noise')
    b = pm.Normal('b', mu=0, sigma=np.std(energies))
    y_pred = pm.Normal('y_pred', mu=-jij*(1-x)+b, sigma=noise, observed=energies)
    trace = pm.sample(draws=5000, chains=4) 

#trace_n = trace[1000:]
pm.plot_trace(trace)
pm.summary(trace)
plt.figure()


x_deg = np.degrees(np.arccos(x))
x_deg_sort = np.sort(x_deg)

jij_m = trace.posterior.jij.values.mean()
b_m = trace.posterior.b.values.mean()
print(jij_m)
print(b_m)
print(np.std(energies))
plt.plot(x_deg, energies, '.')
x_mesh = np.linspace(0, 180, 180)
plt.plot(x_mesh, -jij_m*(1-np.cos(np.radians(x_mesh)))+b_m, 'red', label='Estimated')
plt.figure()
plt.show()
