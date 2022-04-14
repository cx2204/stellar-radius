import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import os

Msun = 1.989 * 10**33 # solar mass in kg.

plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}
fontsize = {'fontsize':30}

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

def choseMS(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    Teff = Teff[index:]
    L = L[index:]
    R = h.log_R[index:]
    return Teff, R, L

h0l = mr.MesaData('data/Z1e-3/LOGS/history.data')
Teff0l, R0l, L0l = choseMS(h0l)

fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(7.5,5.5))

dir1 = [x[1] for x in os.walk('data/low_Z_2Z')][0]

for d1 in dir1:
    h = mr.MesaData('data/low_Z_2Z/{}/LOGS/history.data'.format(d1))
    T, R, L = choseMS(h)
    ax.plot(R, L, c='gray', linewidth=0.3,zorder=0)

ax.plot(R0l, L0l, color='k', zorder=1., linewidth=2.5,label='$Z=10^{-3}$')

p_mu, p_kap, p_eps = 100*0.02/8.8, 100*3.53/8.8, 100*5/8.8
un_m, un_k = 100*0.03/8.8, 100*0.0/8.8
err = [un_m, un_k, 0]
ll = [1,2,3]
label = ['μ', 'κ', 'ε']
sizes1 = [p_mu, p_kap, p_eps]
bar_ax1 = ax.inset_axes([0.01,0.65,0.2,0.2])
bar_ax1.errorbar(ll[:-1],sizes1[:-1],yerr=err[:-1],fmt='o',capsize=5,markersize=1,ecolor='k',c='k')
mybars = bar_ax1.bar(ll,sizes1,color='blue',alpha=0.6,tick_label=label) # add ' autopct='%1.2f%%' ' to show percentage in pie-chart.
bar_ax1.axhline(0,color='k')
bar_ax1.set_title('MS',c='brown')
bar_ax1.set_ylim(0,70)

p4_mu, p4_kap, p4_eps = 100*-0.7/26.9, 100*-33.5/26.9, 100*16.3/26.9
un4_m, un4_k = 100*6.4/26.9, 100*42/26.9
err4 = [un4_m, un4_k, 0]
ll4 = [1,2,3]
label4 = ['μ', 'κ', 'ε']
sizes4 = [p4_mu, p4_kap, p4_eps]
bar_ax4 = ax.inset_axes([0.25,0.3,0.2,0.2])
bar_ax4.errorbar(ll4[:-1],sizes4[:-1],yerr=err4[:-1],fmt='o',capsize=5,markersize=1,ecolor='k',c='k')
mybars4 = bar_ax4.bar(ll4,sizes4,color='blue',alpha=0.5,tick_label=label4) # add ' autopct='%1.2f%%' ' to show percentage in pie-chart.
bar_ax4.axhline(0,color='k')
bar_ax4.set_title('Hertzsprung Gap',c='brown')
bar_ax4.set_ylim(-350,200)

p3_mu, p3_kap, p3_eps = 100*-10.2/51.3, 100*-7.1/51.3, 100*25.6/51.3
un3_m, un3_k = 100*13.4/51.3, 100*24.6/51.3
err3 = [un3_m, un3_k, 0]
ll3 = [1,2,3]
label3 = ['μ', 'κ', 'ε']
sizes3 = [p3_mu, p3_kap, p3_eps]
bar_ax3 = ax.inset_axes([0.45,0.7,0.2,0.2])
bar_ax3.errorbar(ll3[:-1],sizes3[:-1],yerr=err3[:-1],fmt='o',capsize=5,markersize=1,ecolor='k',c='k')
mybars3 = bar_ax3.bar(ll3,sizes3,color='blue',alpha=0.5,tick_label=label3) # add ' autopct='%1.2f%%' ' to show percentage in pie-chart.
bar_ax3.axhline(0,color='k')
bar_ax3.set_title('He burning',c='brown')
bar_ax3.set_ylim(-100,100)

p5_mu, p5_kap, p5_eps = 100*-0.3/1.1, 100*0.06/1.1, 100*0.6/1.1
un5_m, un5_k = 100*0.4/1.1, 100*1.2/1.1
err5 = [un5_m, un5_k, 0]
ll5 = [1,2,3]
label5 = ['μ', 'κ', 'ε']
sizes5 = [p5_mu, p5_kap, p5_eps]
bar_ax5 = ax.inset_axes([0.65,0.3,0.2,0.2])
bar_ax5.errorbar(ll5[:-1],sizes5[:-1],yerr=err5[:-1],fmt='o',capsize=5,markersize=1,ecolor='k',c='k')
mybars5 = bar_ax5.bar(ll5,sizes5,color='blue',alpha=0.5,tick_label=label5) # add ' autopct='%1.2f%%' ' to show percentage in pie-chart.
bar_ax5.axhline(0,color='k')
bar_ax5.set_title('Early Hayashi',c='brown')
bar_ax5.set_ylim(-150,200)

p2_mu, p2_kap, p2_eps = 100*-0.11/6.8, 100*6.2/6.8, 100*-0.8/6.8
un2_m, un2_k = 100*0.06/6.8, 100*0.9/6.8
err2 = [un2_m, un2_k, 0]
ll2 = [1,2,3]
label2 = ['μ', 'κ', 'ε']
sizes2 = [p2_mu, p2_kap, p2_eps]
bar_ax2 = ax.inset_axes([0.8,0.75,0.2,0.2])
bar_ax2.errorbar(ll2[:-1],sizes2[:-1],yerr=err2[:-1],fmt='o',capsize=5,markersize=1,ecolor='k',c='k')
mybars2 = bar_ax2.bar(ll2,sizes2,color='blue',alpha=0.6,tick_label=label2) # add ' autopct='%1.2f%%' ' to show percentage in pie-chart.
bar_ax2.axhline(0,color='k')
bar_ax2.set_title('Late Hayashi',c='brown')
bar_ax2.set_ylim(-15,150)

ax.set_xticks([1.,1.5,2.,2.5,3.])
ax.set_xlim([0.5,3.3])
ax.set_ylim([3.9,5.5])
ax.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)
ax.text(0.6,4.0,'$Z=10^{-3}$',fontsize=20)
# ax.legend(fontsize=12, ncol=1, framealpha=0.5,loc=3)
# ax.text(1.6,4.8,'... ...',fontsize=50)

# get rid of the frame
for spine in bar_ax1.spines.values():
    spine.set_visible(False)

# bar_ax1.set_xlabel([1,2,3],labelpad=10)
bar_ax1.tick_params(which='both',top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
bar_ax1.xaxis.set_ticks_position('none')
bar_ax1.yaxis.set_ticks_position('none')
bar_ax1.get_yaxis().set_visible(False)
bar_ax1.xaxis.labelpad = -5
bar_ax1.patch.set_alpha(0.)

for i,bari in zip(range(3),mybars):
    height = bari.get_height()
    bar_ax1.text(bari.get_x() + bari.get_width()/2, bari.get_height()+err[i]+5, str('%1d%s' % (height,'%')),
                 ha='center', color='blue', fontsize=12)

# get rid of the frame
for spine in bar_ax2.spines.values():
    spine.set_visible(False)

# bar_ax1.set_xlabel([1,2,3],labelpad=10)
bar_ax2.tick_params(which='both',top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
bar_ax2.xaxis.set_ticks_position('none')
bar_ax2.yaxis.set_ticks_position('none')
bar_ax2.get_yaxis().set_visible(False)
bar_ax2.xaxis.labelpad = -5
bar_ax2.patch.set_alpha(0.)
# bar_ax2.text(2.8,20,'?',fontsize=20)

for i,bari in zip(range(3),mybars2):
    height = bari.get_height()
    bar_ax2.text(bari.get_x() + bari.get_width()/2, bari.get_height()+err2[i]+20, str('%1d%s' % (height,'%')),
                 ha='center', color='blue', fontsize=12)

# get rid of the frame
for spine in bar_ax3.spines.values():
    spine.set_visible(False)

# bar_ax1.set_xlabel([1,2,3],labelpad=10)
bar_ax3.tick_params(which='both',top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
bar_ax3.xaxis.set_ticks_position('none')
bar_ax3.yaxis.set_ticks_position('none')
bar_ax3.get_yaxis().set_visible(False)
bar_ax3.xaxis.labelpad = -5
bar_ax3.patch.set_alpha(0.)
# bar_ax2.text(2.8,20,'?',fontsize=20)

for i,bari in zip(range(3),mybars3):
    height = bari.get_height()
    bar_ax3.text(bari.get_x() + bari.get_width()/2, bari.get_height()+err3[i]+10, str('%1d%s' % (height,'%')),
                 ha='center', color='blue', fontsize=12)

for spine in bar_ax4.spines.values():
    spine.set_visible(False)

# bar_ax1.set_xlabel([1,2,3],labelpad=10)
bar_ax4.tick_params(which='both',top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
bar_ax4.xaxis.set_ticks_position('none')
bar_ax4.yaxis.set_ticks_position('none')
bar_ax4.get_yaxis().set_visible(False)
bar_ax4.xaxis.labelpad = -5
bar_ax4.patch.set_alpha(0.)
# bar_ax2.text(2.8,20,'?',fontsize=20)

for i,bari in zip(range(3),mybars4):
    height = bari.get_height()
    bar_ax4.text(bari.get_x() + bari.get_width()/2, bari.get_height()+err4[i]+50, str('%1d%s' % (height,'%')),
                 ha='center', color='blue', fontsize=12)

for spine in bar_ax5.spines.values():
    spine.set_visible(False)

# bar_ax1.set_xlabel([1,2,3],labelpad=10)
bar_ax5.tick_params(which='both',top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
bar_ax5.xaxis.set_ticks_position('none')
bar_ax5.yaxis.set_ticks_position('none')
bar_ax5.get_yaxis().set_visible(False)
bar_ax5.xaxis.labelpad = -5
bar_ax5.patch.set_alpha(0.)
# bar_ax2.text(2.8,20,'?',fontsize=20)

for i,bari in zip(range(3),mybars5):
    height = bari.get_height()
    bar_ax5.text(bari.get_x() + bari.get_width()/2, bari.get_height()+err5[i]+20, str('%1d%s' % (height,'%')),
                 ha='center', color='blue', fontsize=12)

plt.tight_layout()
plt.savefig('figures/bar_low_z_2z.pdf')
plt.show()