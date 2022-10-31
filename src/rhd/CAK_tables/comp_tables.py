import numpy as np
import matplotlib.pyplot as plt

def Get_table(path):
    D = np.loadtxt(path,usecols=[0])[1:]
    T = np.loadtxt(path)[0][1:]

    table= np.loadtxt(path,skiprows=1)[:,1:]
    return D,T,table.T

mins = {'al':0.4,'Ke':0.0,'Q0':100,'Qb':300}
maxs = {'al':0.99,'Ke':0.4,'Q0':4000,'Qb':3000}

fig,axs = plt.subplots(4,3)

for ii,param in enumerate(['al','Ke','Q0','Qb']):
    for jj,met in enumerate(['Y09800/','Y09900/','Y09960/']):
        path = met + param + '_TD'
        D,T,table = Get_table(path)
        pc = axs[ii,jj].pcolor(D,T,table,vmin=mins[param],vmax=maxs[param])
        cb = plt.colorbar(pc,ax=axs[ii,jj])

        axs[0,jj].set_title(met)
    axs[ii,0].set_ylabel(param)


plt.show()
