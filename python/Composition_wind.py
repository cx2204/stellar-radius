import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
import mesaPlot as mp
import matplotlib
Msun = 1.989 * 10**33 # solar mass in kg.
plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}
fontsize = {'fontsize':30}

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

def MS(h):
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    L = L[index:]
    R = h.log_R[index:]
    age = h.star_age[index:]
    mu_ave = h.mu_int[index:]
    abs_mdot = h.log_abs_mdot[index:]
    # return h.center_h1[index:]
    return R, L

################ FIGURE 5 ################ ################ ################
h7 = mr.MesaData('data/Z0.04_hardcode_wind0.0/LOGS/history.data')
R7, L7 = MS(h7)
h6 = mr.MesaData('data/Z0.04_wind0.02_use_min/LOGS/history.data')
R6, L6 = MS(h6)
# h1 = mr.MesaData('data/Z0.02_wind0.02_use_min/LOGS/history.data')
# R1, L1 = MS(h1)
h5 = mr.MesaData('data/Z0.04_wind0.04_use_min/LOGS/history.data')
R5, L5 = MS(h5)
# h4 = mr.MesaData('data/Z2e-2/LOGS/history.data')
# R4, L4 = MS(h4)

fig, ax1 = plt.subplots(ncols=1,nrows=1,sharex=True,sharey=True,figsize=(6.5,4.5))

width_fiducial = 2.5
width_mp = 2.

# color = plt.cm.Reds(np.linspace(0,1,12))

ax1.plot(R5, L5, c='green',linestyle='--',label=r'$Z=Z_{\rm wind}=0.04$')
ax1.plot(R6, L6, c='k',linewidth=2.,label=r'$Z=0.04$, $Z_{\rm wind}=0.02$')
ax1.plot(R7, L7, c='brown',linestyle='-.',label=r'$Z=0.04$, $Z_{\rm wind}=0.0$')

ax1.set_ylim(4.1,5.2)
ax1.set_xlim(0.5,3.2)

ax1.legend(fontsize=10.7, ncol=1, framealpha=0.5,loc=2, frameon=False)

ax1.set_yticks([4.2,4.4,4.6,4.8,5.0,5.2])
ax1.set_xticks([0.5,1.,1.5,2.,2.5,3.])

ax1.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax1.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)

plt.tight_layout()
# plt.savefig('figures/HR_wind.pdf')
plt.show()

################ FIGURE 6 ################ ################ ################

h0 = mr.MesaData('data/Z2e-2/LOGS/history.data')
R0, L0 = MS(h0)
h1 = mr.MesaData('data/Z_I/LOGS/history.data')
R1, L1 = MS(h1)
h2 = mr.MesaData('data/Z_II/LOGS/history.data')
R2, L2 = MS(h2)
h3 = mr.MesaData('data/Z0.04_wind0.02_use_min/LOGS/history.data')
R3, L3 = MS(h3)

fig, ax1 = plt.subplots(ncols=1,nrows=1,sharex=True,sharey=True,figsize=(6.5,4.5))

# width_fiducial = 2.5
# width_mp = 2.

color = plt.cm.Reds(np.linspace(0,1,12))

ax1.plot(R0,L0,c='k',label='X=0.75,Y=0.23,Z=0.02',linewidth=width_fiducial)
ax1.plot(R1,L1,c=color[11],label='X=0.64,Y=0.32,Z=0.04 (I)',linewidth=width_mp,linestyle=':')
ax1.plot(R2,L2,c=color[9],label='X=0.73,Y=0.23,Z=0.04 (II)',linewidth=width_mp,linestyle='--')
ax1.plot(R3,L3,c=color[4],label='X=0.75,Y=0.21,Z=0.04 (III)',linewidth=width_mp,linestyle='-.')

ax1.set_ylim(4.1,5.2)
ax1.set_xlim(0.5,3.2)

ax1.legend(fontsize=10.7, ncol=1, framealpha=0.5,loc=2, frameon=False)

ax1.set_yticks([4.2,4.4,4.6,4.8,5.0,5.2])
ax1.set_xticks([0.5,1.,1.5,2.,2.5,3.])

ax1.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax1.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)

plt.tight_layout()
# plt.savefig('figures/HR_Y.pdf')
plt.show()