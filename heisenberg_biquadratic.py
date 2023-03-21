import numpy as np
import pymc as pm
import matplotlib.pyplot as plt
import re
np.random.seed(123)

# https://qiita.com/takubb/items/a2e02fe61aad80f8f2f3

x = np.array([])
y = np.array([])
with open("jij_perpendicular.txt", encoding="utf-8") as f:
  pattern = re.compile(r'^\#')
  for line in f:
    if pattern.match(line):
      continue
    line = line.split()
    x = np.append(x, float(line[0]))
    y = np.append(y, float(line[2]))

x = np.cos(np.deg2rad(x))
x2 = x**2
y = y*27.2114*1000/2


with pm.Model() as model:
    jij = pm.Normal('jij', mu=200, sigma=100)
    lambda_eff = pm.Normal('lambda_eff', mu = 0, sigma=1000)
    # noise = pm.Normal('noise', mu=0, sigma=1)
    noise = pm.HalfFlat('noise')
    b = pm.Normal('b', mu=0, sigma=500)
    # y_pred = pm.Normal('y_pred', mu=-jij*x, sigma=noise, observed=y)
    y_pred = pm.Normal('y_pred', mu=-jij*x+lambda_eff*x2+b, sigma=noise, observed=y)
    trace = pm.sample(draws=5000, chains=4) 

#trace_n = trace[1000:]
pm.plot_trace(trace)
pm.summary(trace)
plt.show()
print(trace.posterior.jij.values.mean())
print(trace.posterior.lambda_eff.values.mean())

x_deg = np.degrees(np.arccos(x))
jij_m = trace.posterior.jij.values.mean()
lambda_eff_m = trace.posterior.lambda_eff.values.mean()
b_m = trace.posterior.b.values.mean()
plt.plot(x_deg, y, '.')
plt.plot(x_deg, -jij_m*np.cos(np.radians(x_deg))+lambda_eff_m*np.cos(np.radians(x_deg))**2+b_m, 'red', label='Estimated')
plt.show()
