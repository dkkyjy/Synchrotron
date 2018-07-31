import time
import numpy as np
import matplotlib.pyplot as plt
import Synchrotron

'''
  Follow the paper A&A 367, 809-825 (2001)
  The multifrequency emission of Mrk 501 from radio to TeV gamma-rays
  The emission and absorption coefficient calculation using the approximation (Appendix A)
  The equation (3) has some difficult in high energy, donot use it.

x = np.logspace(-4, 1, 100)
F_x = Synchrotron.F(x)
R_x = Synchrotron.R(x)
plt.figure()
plt.loglog(x, F_x, 'r-')
plt.xlim(1e-4, 1e1)
plt.ylim(1e-4, 1e0)
plt.xlabel('x')
plt.ylabel('F(x)')
plt.show()
'''

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

BList = np.logspace(-3, 3, 10)
KList = np.logspace(0, 5, 10)
nList = np.linspace(1.5, 4, 10)
#SpectrumParsList = [(K, n, gamma_min, gamma_max) for K in KList]
#SpectrumParsList = [(K, n, gamma_min, gamma_max) for n in nList]
#SpectrumParsList = [(K, n, gamma_min, gamma_max) for K in KList for n in nList]

def tau(k, R):
    return 2 * R * k

def I(j, k, R):
    #return j / k * (1 - np.exp(-tau))
    #return j / k * (1 - 2/tau**2 * (1 - np.exp(-tau) * (tau+1)))
    tau = 2 * R * k
    I = np.zeros_like(tau)
    I[tau < 1e-5] = j[tau < 1e-5] * 2 * R
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
np.save('nu.npy', nu)
np.savetxt('nu.dat', nu)
for B in [1, 10, 100]:
#for B in BList:
#for SpectrumPars in SpectrumParsList:
    tmin = time.time()
    j = Synchrotron.J(nu, B, SpectrumType, SpectrumPars)
    k = Synchrotron.K(nu, B, SpectrumType, SpectrumPars)
    tmax = time.time()
    print(f'Synchrotron emission use {tmax - tmin}s')

    np.save(f'j_{B}.npy', j)
    np.save(f'k_{B}.npy', k)
    np.savetxt(f'j_{B}.dat', j)
    np.savetxt(f'k_{B}.dat', k)

    I_source = I(j, k, R)
    L_source = L(I_source, R)
    Flux = F(L_source, dL)
    Flux_obs = F_obs(Flux, delta, z)
    nu_obs = Nu_obs(nu, delta, z)

    plt.loglog(nu, nu*j, '-', label=f'B={B}T')
    plt.loglog(nu, nu*k, '-', label=f'B={B}T')
    #plt.loglog(nu, nu*I_source, '-', label=f'B={B}T')
    #plt.loglog(nu, nu*L_source, label=f'B={B}T')
    #plt.loglog(nu_obs, nu_obs*Flux_obs, label=f'B={B:.2f}T')


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
'''
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

nu, F = np.load('../SEDmodel/flux_1.npy')
plt.loglog(nu, nu*F, '--')
nu, F = np.load('../SEDmodel/flux_10.npy')
plt.loglog(nu, nu*F, '--')
nu, F = np.load('../SEDmodel/flux_100.npy')
plt.loglog(nu, nu*F, '--')
'''


plt.legend()
plt.xlim(1e5, 1e25)
plt.ylim(1e-25, 1e-5)
plt.xlabel(r'$\nu$ [Hz]')
plt.ylabel(r'$\nu$ j$_{\nu}$ [erg cm$^{-3}$ s$^{-1}$]')
plt.show()
