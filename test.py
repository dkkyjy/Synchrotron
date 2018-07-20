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
    j = Synchrotron.J(nu, B, SpectrumType, SpectrumPars)
    k = Synchrotron.K(nu, B, SpectrumType, SpectrumPars)
    #I = Synchrotron.I(nu, B, R, SpectrumType, SpectrumPars)
    #L = Synchrotron.L(nu, B, R, SpectrumType, SpectrumPars)
    #F = Synchrotron.F(nu, B, R, delta, dL, z, SpectrumType, SpectrumPars)

    tau = 2 * R * np.array(k)
    plt.loglog(nu, tau, label=f'B={B}T')
    I = np.array(j)/np.array(k) * (1 - 2/tau**2 * (1 - np.e**(-tau) * (tau + 1)))
    L = 4 * np.pi**2 * R**2 * I
    F = L / (4*np.pi*dL**2)

    nu_obs = nu * delta / (1+z)
    #plt.loglog(nu, nu*j, label=f'B={B}T')
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
