import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.interpolate import interp1d
import os
from Epochs_LR import max_drdt

plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}
fontsize = {'fontsize':30}

matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)

Msun = 1.989 * 10**33 # solar mass in kg.

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def choseMS(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    Teff = Teff[index:]
    L = L[index:]
    R = h.log_R[index:]
    t = h.star_age[index:]
    return t, R, L

def EndOfMS(h):
    # compute radius, luminosity and effective temperature of star ONLY on the MS
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    r = h.log_R[index:]
    L_ms = L[index:]
    t = h.star_age[index:]
    center_h1 = h.center_h1[index:]
    ind2 = np.where(center_h1 > (1e-3))[0]
    # return r[ind2][-1], L_ms[ind2][-1], Teff[ind2][-1]
    ind_max = np.where(r[ind2] == max(r[ind2]))[0]
    return t[ind2][ind_max],r[ind2][ind_max],L_ms[ind2][ind_max]

h1 = mr.MesaData('data/high_Z_2Z/Z2e-2/LOGS/history.data')
h2 = mr.MesaData('data/low_Z_2Z/Z1e-3/LOGS/history.data')

dir1 = [x[1] for x in os.walk('data/high_Z_2Z')][0]
dir2 = [x[1] for x in os.walk('data/low_Z_2Z')][0]

t1, R1, L1 = choseMS(h1)
t2, R2, L2 = choseMS(h2)

R_hg_1, L_hg_1 = max_drdt(h1)
R_hg_2, L_hg_2 = max_drdt(h2)
ind_hg_1, ind_hg_2 = find_nearest(R1,R_hg_1), find_nearest(R2,R_hg_2)
t_hg_1, t_hg_2 = t1[ind_hg_1], t2[ind_hg_2]

t_ms_1, R_ms_1, L_ms_1 = EndOfMS(h1)
t_ms_2, R_ms_2, L_ms_2 = EndOfMS(h2)

R_begin_1, R_end_1 = 1.416, 2.544
ind_begin_1, ind_end_1 = find_nearest(R1,R_begin_1), find_nearest(R1,R_end_1)

R_begin_2, R_end_2 = 1.181, 2.713
ind_begin_2, ind_end_2 = find_nearest(R2,R_begin_2), find_nearest(R2,R_end_2)

fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(6.5,4.5))

for d1 in dir1:
    h = mr.MesaData('data/high_Z_2Z/{}/LOGS/history.data'.format(d1))
    t, R, L = choseMS(h)
    ax.plot(t, R, c='k', linewidth=0.3)

for d2 in dir2:
    h = mr.MesaData('data/low_Z_2Z/{}/LOGS/history.data'.format(d2))
    t, R, L = choseMS(h)
    ax.plot(t, R, c='mediumblue', linewidth=0.3)

ax.plot(t1, R1, label='$Z=0.02$', c='k',linewidth=2)
ax.plot(t2, R2, label='$Z=10^{-3}$',c='mediumblue',linewidth=2)

# ax.plot(t1[ind_end_1:], R1[ind_end_1:], c='k',linewidth=2,linestyle='--')
# ax.plot(t2[ind_end_2:], R2[ind_end_2:], c='mediumblue',linewidth=2,linestyle='--')

ax.plot(t1[ind_begin_1:ind_end_1], R1[ind_begin_1:ind_end_1], c='lime',alpha=0.3,linewidth=10,zorder=0,label='Hertzsprung Gap')
ax.plot(t2[ind_begin_2:ind_end_2], R2[ind_begin_2:ind_end_2], c='lime',alpha=0.3,linewidth=10,zorder=0)

ax.scatter([t_hg_1,t_hg_2],[R_hg_1,R_hg_2],c='lime',edgecolors='k',marker='o',s=70,zorder=5,label='Maximum dlogR/dt')
ax.scatter([t_ms_1,t_ms_2],[R_ms_1,R_ms_2],c='cyan',edgecolors='k',marker='*',s=200,zorder=5,label='MS')


# ax.set_yscale('log')
ax.set_xlabel('age [yr]',fontsize=15)
ax.set_ylabel('log$_{10}$(R/R$_{\odot}$)',fontsize=15)
ax.legend(fontsize=15, ncol=1, framealpha=0.5,loc=2, frameon=True)
plt.tight_layout()
# plt.savefig('figures/age_radius_HG.pdf')
plt.show()