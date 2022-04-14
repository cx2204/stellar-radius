import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
import os
import matplotlib

Msun = 1.989 * 10**33 # solar mass in kg.

plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}
fontsize = {'fontsize':30}

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

metallicity = ['Z1e-4','Z1e-3','Z1e-2','Z2e-2','Z3e-2','Z4e-2']
color = plt.cm.Greys_r(np.linspace(0,1,20))
# c = color[[2,4,6,8,10,12,14,16,18]]

def choseMS(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    # print(np.where(np.round(L, 3) == np.round(Lnuc, 3))[0])
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    R = h.log_R
    age = h.star_age
    cxhe = h.center_he4
    return Teff[index:], R[index:], L[index:], age[index:], cxhe[index:]

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def He_ignition(h):
    # compute radius, luminosity and effective temperature at He ignition
    center_h1 = h.center_h1
    ind_tams = np.where(center_h1 <= 1e-3)[0] # this is after TAMS
    center_he = h.center_he4[ind_tams]
    r = h.log_R[ind_tams]
    age = h.star_age[ind_tams]
    lum = h.log_L[ind_tams]
    ind_he = np.where(center_he <= 0.90)[0] # he ignition
    # print(center_he[ind_max], age[ind_max])
    return r[ind_he][0], lum[ind_he][0], age[ind_he][0]

fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(6.5,4.5))

r_he = []
t_he = []
l_he = []
for z in metallicity:
    # print(float(z))
    h = mr.MesaData('data/{}/LOGS/history.data'.format(z))

    # print(float(z))
    cyclers = cycler(color=color)
    plt.rc('axes', prop_cycle=cyclers)
    Teff, R, L, age, cxhe = choseMS(h)

    rad_he,lum_he,age_he = He_ignition(h)
    r_he.append(rad_he)
    t_he.append(age_he)
    l_he.append(lum_he)

    sc = ax.scatter(R,L,c=cxhe, cmap='viridis_r',s=5)
    # plt.plot(R, L, label=r'Z = {}'.format(float(z)), zorder=5)

fig.colorbar(sc)
# ax.set_clim(0,1)
ax.text(3.2,5.25,'X$_c(^4He)$',fontsize=15)
ax.scatter(r_he,l_he,marker='D',s=40,color='orangered',label=r'X$_c(^4He)=90\%$',zorder=10)
ax.legend(fontsize=14, ncol=2, framealpha=0.5,loc=2, frameon=False)
# plt.title('Z=[0.0001,0.0002,0.001,0.01,0.02,0.03,0.04]\n(top to bottom)')
# plt.xlabel('Age [yr]')
ax.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)
ax.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
plt.tight_layout()
# plt.savefig('figures/HeIgn-HR-colored.pdf')
plt.show()

# z_norm = 0.0




# if z_norm == 0.02:
#     import os
#
#     root2 = 'data/microphysics_0.02_0.04/'
#     dirlist2 = [item for item in os.listdir(root2) if os.path.isdir(os.path.join(root2, item))]
#
#     color = plt.cm.Greys_r(np.linspace(0, 1, 20))
#     c = color[[2, 4, 6, 8, 10, 12, 14, 16, 18]]
#
#     r_he2 = []
#     t_he2 = []
#     l_he2 = []
#     for i in range(len(dirlist2)):
#         h2 = mr.MesaData('data/{}/LOGS/history.data'.format(dirlist2[i]))
#         Teff2, R2, L2, age2, cxhe2 = choseMS(h2)
#
#         cyclers = cycler(color=c)
#         plt.rc('axes', prop_cycle=cyclers)
#
#         rad_he2,lum_he2,age_he2 = He_ignition(h2)
#         r_he2.append(rad_he2)
#         t_he2.append(age_he2)
#         l_he2.append(lum_he2)
#         #
#         plt.scatter(R2,L2,c=cxhe2, cmap='jet',s=1)
#         plt.plot(R2, L2, zorder=5)
#
#     plt.colorbar()
#     plt.clim(0,1)
#     plt.plot(r_he2,l_he2,'o',color='silver',label=r'CX$_{\rm He}$ = 90%',zorder=10)
#     plt.legend()
#     # plt.title('Z=[0.0001,0.0002,0.001,0.01,0.02,0.03,0.04]\n(top to bottom)')
#     # plt.xlabel('Age [yr]')
#     plt.xlabel('log(R/R$_{\odot}$)')
#     plt.ylabel('log(L/L$_{\odot}$)')
#
#     # plt.savefig('figures/HeIgn-HR-0.02.pdf')
#     plt.show()