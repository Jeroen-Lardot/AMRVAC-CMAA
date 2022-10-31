import numpy as np
import matplotlib.pyplot as plt

c = 2.99e10
arad = 7.5646e-15

#file = '16Msol_deJager_profile20.data'
#file = '16Msol_Leuven_profile18.data'
#file = '20Msol_deJager_profile19.data'
#file = '20Msol_Leuven_profile13.data'


file = '16Msol_deJager_profile20.data'

DATA = np.loadtxt(file,skiprows=6)

I = np.transpose(DATA)[0]
M = np.transpose(DATA)[1]
logR = np.transpose(DATA)[2]
logT = np.transpose(DATA)[3]
logRho =  np.transpose(DATA)[4]
pgas = np.transpose(DATA)[11]
prad = np.transpose(DATA)[10]
kappa = np.transpose(DATA)[16]
Lum = np.transpose(DATA)[17]
Lrad = np.transpose(DATA)[18]

X = np.transpose(DATA)[6]
Y = np.transpose(DATA)[7]
Z = np.transpose(DATA)[8]


Diffcoef = c/(kappa*10**logRho)

plt.figure()
plt.plot(10**logR, M)
plt.xlim([300,700])
plt.xlabel('R')
plt.ylabel('mass')

plt.figure()
plt.plot(10**logR, 10**logRho)
plt.xlim([300,700])
plt.ylim([1.e-12,1.e-7])
plt.xlabel('R')
plt.ylabel('Rho')

plt.figure()
plt.plot(10**logR, X, label='X')
plt.plot(10**logR, Y, label='Y')
plt.plot(10**logR, Z, label='Z')
plt.xlim([300,700])
# plt.ylim([1.e-12,1.e-7])
plt.xlabel('R')
plt.ylabel('Fraction')
plt.legend()

plt.figure()
plt.plot(10**logR, 10**logT)
plt.plot(10**logR, (3*prad/arad)**0.25)
plt.ylim([100,100000])
plt.xlim([300,700])
plt.xlabel('R')
plt.ylabel('T')


plt.figure()
plt.plot(10**logR, kappa)
plt.xlim([300,700])
plt.xlabel('R')
plt.ylabel('kappa')

plt.figure()
plt.plot(10**logR, 10**Lum)
plt.plot(10**logR, 10**Lrad)
plt.xlabel('R')
plt.ylabel('L')
plt.xlim([0,600])
plt.ylim([0,10**6])

plt.plot(10**logR, Diffcoef)
np.savetxt('first_profile', np.transpose([10**logR[::-1],10**logRho[::-1],pgas[::-1],prad[::-1]*3,Diffcoef[::-1]]))

plt.show()
