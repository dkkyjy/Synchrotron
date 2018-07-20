import time
import numpy as np
import matplotlib.pyplot as plt
import Synchrotron

delta = 10
dL = 1e27
R = 1e16
B = 1
z = 0

K = 1
n = 2
gamma_min = 1e2
gamma_max = 1e5
SpectrumType = 1
SpectrumPars = (K, n, gamma_min, gamma_max)

nu = np.logspace(10, 20, 100)
for B in [1, 10, 100]:
    print(B)
    j = Synchrotron.J(nu, B, SpectrumType, SpectrumPars)
    #k = Synchrotron.K(nu, B, SpectrumType, SpectrumPars)
    #I = Synchrotron.I(list(nu), B, R, SpectrumType, SpectrumPars)
    #L = Synchrotron.L(list(nu), B, R, SpectrumType, SpectrumPars)
    #F = Synchrotron.F(list(nu), B, R, delta, dL, z, SpectrumType, SpectrumPars)

    #tau = 2 * R * np.array(k)
    #I_python = np.array(j)/np.array(k) * (1 - 2/tau**2 * (1 - np.e**(-tau) * (tau + 1)))
    #plt.plot(I - I_python)

    #nu_obs = Synchrotron.Nu_obs(list(nu), delta, z)
    #nu_obs = np.array(nu_obs)
    plt.loglog(nu, nu*j, label=f'B={B}T')
    #plt.loglog(nu, nu*k, label=f'B={B}T')
    #plt.loglog(nu_obs, nu_obs*F, label=f'B={B}T')

'''
nu, j = np.load('../SEDmodel/j_1.npy')
plt.loglog(nu, nu*j, '--')
nu, j = np.load('../SEDmodel/j_10.npy')
plt.loglog(nu, nu*j, '--')
nu, j = np.load('../SEDmodel/j_100.npy')
plt.loglog(nu, nu*j, '--')
nu, k = np.load('../SEDmodel/alpha_1.npy')
plt.loglog(nu, nu*k, '--')
nu, k = np.load('../SEDmodel/alpha_10.npy')
plt.loglog(nu, nu*k, '--')
nu, k = np.load('../SEDmodel/alpha_100.npy')
plt.loglog(nu, nu*k, '--')
nu, F = np.load('../SEDmodel/flux_1.npy')
plt.loglog(nu, nu*F, '--')
nu, F = np.load('../SEDmodel/flux_10.npy')
plt.loglog(nu, nu*F, '--')
nu, F = np.load('../SEDmodel/flux_100.npy')
plt.loglog(nu, nu*F, '--')
'''

plt.legend()
plt.xlim(1e10, 1e20)
plt.ylim(1e-17, 1e-7)
plt.xlabel(r'$\nu$ [Hz]')
plt.ylabel(r'$\nu$ j$_{\nu}$ [erg cm$^{-3}$ s$^{-1}$]')
plt.show()
assert False

tmin = time.time()
nu_obs = Synchrotron.Nu_obs(list(nu), delta, z)
tmax = time.time()
print(f'C use: {tmax-tmin}s')

tmin = time.time()
nu_obs_python = nu * delta / (1+z)
tmax = time.time()
print(f'Python use: {tmax-tmin}s')
