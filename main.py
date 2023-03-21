import numpy as np
import pymc as pm
import matplotlib.pyplot as plt
np.random.seed(123)

# https://qiita.com/takubb/items/a2e02fe61aad80f8f2f3
N = 100
alpha_real = 2.5
beta_real = 0.9
eps_real = np.random.normal(0, 0.5, size=N)

x = np.random.normal(10, 1, N)
y_real = alpha_real + beta_real *x
y = y_real + eps_real

with pm.Model() as model:
    a = pm.Normal('alpha', mu=0, sigma=10)
    b = pm.Normal('beta', mu=0, sigma=1) 
    noise = pm.Normal('noise', mu=0, sigma=1)
    y_pred = pm.Normal('y_pred', mu=a*x+b, sigma=noise, observed=y)
    trace = pm.sample(draws=5000, chains=2) 

#trace_n = trace[1000:]
# pm.traceplot(trace_n)
pm.plot_trace(trace)
#pm.summary(trace_n)
pm.summary(trace)
plt.show()
