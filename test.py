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

def tau(R, k):
    return 2 * R * k

def I(j, k, tau):
    #return j / k * (1 - np.exp(-tau))
    #return j / k * (1 - 2/tau**2 * (1 - np.exp(-tau) * (tau+1)))
    I = np.zeros_like(tau)
    I[tau < 1e-5] = j[tau < 1e-5] * (tau[tau <= 1e-5] / k[tau <= 1e-5])
    I[tau >= 1e-5] = j[tau >= 1e-5] / k[tau >= 1e-5] * (1 - np.exp(-tau[tau >= 1e-5]))
    return I


def L(I, R):
    return 4 * np.pi**2 * R**2 * I

def F(L, rL):
    return L / (4*np.pi * dL**2)

def F_obs(F, delta, z):
    return F * delta**3 * (1+z)

def Nu_obs(nu, delta, z):
    return nu * delta / (1+z)

nu = np.logspace(5, 25, 1000)
for B in [1, 10, 100]:
    j = Synchrotron.J(nu, B, SpectrumType, SpectrumPars)
    k = Synchrotron.K(nu, B, SpectrumType, SpectrumPars)

    tau = 2 * R * k
    I_source = I(j, k, tau)
    L_source = L(I_source, R)
    Flux = F(L_source, dL)
    Flux_obs = F_obs(Flux, delta, z)
    nu_obs = Nu_obs(nu, delta, z)
    #Flux_obs = np.pi * R**2 * j / k / dL**2 * (1 - 2/tau**2 * (1 - np.exp(-tau) * (tau+1)))
    #plt.loglog(nu, nu*j, '-', label=f'B={B}T')
    #plt.loglog(nu, nu*k, '-', label=f'B={B}T')
    #plt.loglog(nu, nu*I_source, '-', label=f'B={B}T')
    #plt.loglog(nu, nu*L_source, label=f'B={B}T')
    plt.loglog(nu_obs, nu_obs*Flux_obs, label=f'B={B}T')
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

nu, I = np.load('../SEDmodel/I_1.npy')
plt.loglog(nu, nu*I, '--')
nu, I = np.load('../SEDmodel/I_10.npy')
plt.loglog(nu, nu*I, '--')
nu, I = np.load('../SEDmodel/I_100.npy')
plt.loglog(nu, nu*I, '--')

nu, L = np.load('../SEDmodel/flux_1.npy')
plt.loglog(nu, nu*L, '--')
nu, L = np.load('../SEDmodel/flux_10.npy')
plt.loglog(nu, nu*L, '--')
nu, L = np.load('../SEDmodel/flux_100.npy')
plt.loglog(nu, nu*L, '--')
'''
nu, F = np.load('../SEDmodel/flux_1.npy')
plt.loglog(nu, nu*F, '--')
nu, F = np.load('../SEDmodel/flux_10.npy')
plt.loglog(nu, nu*F, '--')
nu, F = np.load('../SEDmodel/flux_100.npy')
plt.loglog(nu, nu*F, '--')


plt.legend()
plt.xlim(1e5, 1e25)
plt.ylim(1e-17, 1e-7)
plt.xlabel(r'$\nu$ [Hz]')
plt.ylabel(r'$\nu$ j$_{\nu}$ [erg cm$^{-3}$ s$^{-1}$]')
plt.show()
