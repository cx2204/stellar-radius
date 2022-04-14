import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from Epochs_LR import choseMS, EndOfMS, max_drdt, He_ignition, ConvectiveFraction, Hayashi
import os

plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}
fontsize = {'fontsize':30}

matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)

dir1 = [x[1] for x in os.walk('data/High_Z_null')][0]
dir2 = [x[1] for x in os.walk('data/Low_Z_null')][0]
dir3 = [x[1] for x in os.walk('data/Resolution')][0]
dir4 = [x[1] for x in os.walk('data/Resolution_highZ')][0]

l_hay_h = 4.98
l_hay_l = 5.02

fig, ((ax1,ax2),(ax4,ax3)) = plt.subplots(ncols=2,nrows=2,figsize=(11,7.5),sharey=True,sharex=True)

lw = 0.5
c_main = 'k'
c_epoch = 'cyan', 'lime', 'yellow', 'orchid', 'pink'
marker_epoch = '*', 'o', 'P', 'X', 'D'
size_epoch = 130, 50, 80, 50, 40
edge_c = 'k'

# print(dir1)

r1 = []
for d1 in dir1:

    h = mr.MesaData('data/High_Z_null/{}/LOGS/history.data'.format(d1))
    R, L = choseMS(h)
    ax1.plot(R, L, c=c_main, linewidth=lw)

    r_ems, l_ems = EndOfMS(h)
    r_hg, l_hg = max_drdt(h)
    r_he, l_he = He_ignition(h)
    r_conv, l_conv = ConvectiveFraction(h,50)
    r_hay, l_hay = Hayashi(h,l_hay_h)
    r_epoch = [r_ems[0],r_hg,r_he,r_conv,r_hay]
    # print(r_epoch)
    l_epoch = [l_ems, l_hg, l_he, l_conv, l_hay]
    r1.append(r_epoch)
    for i in range(len(r_epoch)):
        ax1.scatter(r_epoch[i],l_epoch[i],c=c_epoch[i],edgecolors=edge_c,marker=marker_epoch[i],s=size_epoch[i],zorder=5)

r1_T = np.array(r1).T
# print(r1_T)
for i in range(len(r1_T)):
    r_h = max(r1_T[i])
    r_l = min(r1_T[i])
    ax1.axvspan(r_l,r_h,color=c_epoch[i],alpha=0.3,zorder=0)

r2 = []
for d2 in dir2:

    h = mr.MesaData('data/Low_Z_null/{}/LOGS/history.data'.format(d2))
    R, L = choseMS(h)
    ax2.plot(R, L, c=c_main, linewidth=lw)

    r_ems, l_ems = EndOfMS(h)
    r_hg, l_hg = max_drdt(h)
    r_he, l_he = He_ignition(h)
    r_conv, l_conv = ConvectiveFraction(h,50)
    r_hay, l_hay = Hayashi(h,l_hay_l)
    r_epoch = [r_ems[0],r_hg,r_he,r_conv,r_hay]
    l_epoch = [l_ems, l_hg, l_he, l_conv, l_hay]
    r2.append(r_epoch)
    for i in range(len(r_epoch)):
        ax2.scatter(r_epoch[i],l_epoch[i],c=c_epoch[i],edgecolors=edge_c,marker=marker_epoch[i],s=size_epoch[i],zorder=5)

r2_T = np.array(r2).T
# print(r1_T)
for i in range(len(r2_T)):
    r_h = max(r2_T[i])
    r_l = min(r2_T[i])
    ax2.axvspan(r_l,r_h,color=c_epoch[i],alpha=0.3,zorder=0)

r3 = []
for d3 in dir3:
    # print(d3)
    h = mr.MesaData('data/Resolution/{}/LOGS/history.data'.format(d3))
    R, L = choseMS(h)
    ax3.plot(R, L, c=c_main, linewidth=lw)

    r_ems, l_ems = EndOfMS(h)
    r_hg, l_hg = max_drdt(h)
    r_he, l_he = He_ignition(h)
    r_conv, l_conv = ConvectiveFraction(h,50)
    r_hay, l_hay = Hayashi(h,l_hay_l)
    r_epoch = [r_ems[0],r_hg,r_he,r_conv,r_hay]
    l_epoch = [l_ems, l_hg, l_he, l_conv, l_hay]
    r3.append(r_epoch)
    for i in range(len(r_epoch)):
        ax3.scatter(r_epoch[i],l_epoch[i],c=c_epoch[i],edgecolors=edge_c,marker=marker_epoch[i],s=size_epoch[i],zorder=5)

r3_T = np.array(r3).T
# print(r1_T)
for i in range(len(r3_T)):
    r_h = max(r3_T[i])
    r_l = min(r3_T[i])
    ax3.axvspan(r_l,r_h,color=c_epoch[i],alpha=0.3,zorder=0)

r4 = []
for d4 in dir4:
    # print(d4)
    h = mr.MesaData('data/Resolution_highZ/{}/LOGS/history.data'.format(d4))
    R, L = choseMS(h)
    ax4.plot(R, L, c=c_main, linewidth=lw)

    r_ems, l_ems = EndOfMS(h)
    r_hg, l_hg = max_drdt(h)
    r_he, l_he = He_ignition(h)
    r_conv, l_conv = ConvectiveFraction(h,50)
    r_hay, l_hay = Hayashi(h,l_hay_l)
    r_epoch = [r_ems[0],r_hg,r_he,r_conv,r_hay]
    l_epoch = [l_ems, l_hg, l_he, l_conv, l_hay]
    r4.append(r_epoch)
    for i in range(len(r_epoch)):
        ax4.scatter(r_epoch[i],l_epoch[i],c=c_epoch[i],edgecolors=edge_c,marker=marker_epoch[i],s=size_epoch[i],zorder=5)

r4_T = np.array(r4).T
# print(r1_T)
for i in range(len(r4_T)):
    r_h = max(r4_T[i])
    r_l = min(r4_T[i])
    ax4.axvspan(r_l,r_h,color=c_epoch[i],alpha=0.3,zorder=0)

ax1.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
# ax2.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=15)
ax4.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax3.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)
ax4.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)

ax1.set_title('High Z',fontsize=20)
ax2.set_title('Low Z',fontsize=20)

ax1.text(.6,5.2,'Implementation of microphysics',fontsize=14)
# ax1.text(0.7,5.,r'$Z=Z_{\kappa,\mu,\epsilon}=Z_{\rm wind}=0.02$',fontsize=12)
ax2.text(.6,5.2,'Implementation of microphysics',fontsize=14)
# ax2.text(1.7,4.6,r'$Z=Z_{\kappa,\mu,\epsilon}=Z_{\rm wind}=10^{-3}$',fontsize=12)
ax3.text(.6,5.2,'Resolution',fontsize=14)
# ax3.text(1.4,4.6,'$d_x=0.5-1.0$,$vct=10^{-4}-10^{-3}$',fontsize=12)
# ax3.text(1.4,4.5,r'N$_{\rm max}$=1932-4136,dt$_{\rm max}$=10$^{4.3-5.6}$ yr',fontsize=12)
ax4.text(.6,5.2,'Resolution',fontsize=14)

ax1.set_xlim(0.5,3.1)
ax1.set_ylim(4.1,5.3)
ax1.set_xticks([0.5,1.,1.5,2.,2.5,3.])
ax1.set_yticks([4.2,4.4,4.6,4.8,5.,5.2])
# fig.text(0.2,0.85,'$Z=10^{-3}$',fontsize=15)
# fig.text(0.2,0.85,'$Z=10^{-3}$',fontsize=15)

plt.tight_layout()
# plt.savefig('figures/appendix.pdf')
plt.show()