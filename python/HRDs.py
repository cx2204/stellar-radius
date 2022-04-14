import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.interpolate import interp1d

Msun = 1.989 * 10**33 # solar mass in kg.

plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}
fontsize = {'fontsize':30}

matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)

def EndOfMS(h):

    # compute radius, luminosity and effective temperature of star ONLY on the MS
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    r = h.log_R[index:]
    L_ms = L[index:]
    Teff = h.log_Teff[index:]
    center_h1 = h.center_h1[index:]
    ind2 = np.where(center_h1 > (1e-3))[0]
    # return r[ind2][-1], L_ms[ind2][-1], Teff[ind2][-1]
    ind_max = np.where(r[ind2] == max(r[ind2]))[0]
    return r[ind2][ind_max],L_ms[ind2][ind_max],Teff[ind2][ind_max]

def He_ignition(h):

    # compute radius, luminosity and effective temperature at He ignition
    center_h1 = h.center_h1
    ind_tams = np.where(center_h1 <= (1e-3))[0] # this is after TAMS
    center_he = h.center_he4[ind_tams]
    r = h.log_R[ind_tams]
    lum = h.log_L[ind_tams]
    he4, ind_he = find_nearest(center_he, 0.9)
    r = r[:ind_he]
    lum = lum[:ind_he]
    center_he = center_he[:ind_he]
    ind_max = np.where(r == max(r))[0]
    # print(center_he[ind_max], age[ind_max])
    return r[ind_max], lum[ind_max]

def Hayashi(h,L_const):
    # compute radius, luminosity and effective temperature on the hayashi track
    radius = h.log_R
    L = h.log_L
    R_int = np.interp(L_const,L,radius)
    # print(R_int)
    return R_int

def Hayashi_end(h):
    radius = h.log_R
    ind_m = np.where(radius == max(radius))[0]
    L = h.log_L
    return max(radius), L[ind_m]

def choseMS(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    Teff = Teff[index:]
    L = L[index:]
    R = h.log_R[index:]
    # center_h1 = h.center_h1[index:]
    # ind2 = np.where(center_h1 < (1e-4 + (Z_eps - Z_fid)))[0]
    # return Teff[ind2], R[ind2], L[ind2]
    return Teff, R, L

def choseMS2(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    Teff = Teff[index:]
    L = L[index:]
    R = h.log_R[index:]
    CXhe = h.center_he4[index:]
    CXh = h.center_h1[index:]
    age = h.star_age[index:]
    logLhe = h.log_LHe[index:]
    logLh = h.log_LH[index:]
    return CXhe, CXh, age, logLhe,logLh
    # return Teff, R, L

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def choseMS3(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    Teff = Teff[index:]
    L = L[index:]
    R = h.log_R[index:]
    CXHe = h.center_he4[index:]
    age = h.star_age[index:]
    m_conv_tot = h.m_conv_tot[index:]
    m_conv_core = h.mass_conv_core[index:]
    m_star = h.star_mass[index:]
    m_h1 = h.total_mass_h1[index:]
    m_tot, m_he_core = h.star_mass[index:], h.he_core_mass[index:]
    return age, CXHe, m_conv_tot, m_conv_core, R, L, m_h1/m_star, m_tot, m_he_core

def ConvectiveFraction(h, p):
    # p = fraction of convective envelope mass to total envelope mass in %
    age, CXHe, m_conv_tot, m_conv_core, R, L, xh, m_tot, m_he_core = choseMS3(h)
    # m_envel = []
    # for l in range(len(R)):
    #     m_c_tot = m_conv_tot[l] / Msun
    #     m_core = m_conv_core[l]
    #     if m_c_tot >= m_core:
    #         m_envel.append(m_c_tot - m_core)
    #     else:
    #         m_envel.append(0)
    #
    # conv_frac = np.array(m_envel) / (m_tot - m_he_core)
    conv_frac = (m_conv_tot / Msun) / (m_tot - m_he_core)
    f_50, ind_50 = find_nearest(conv_frac, p/100)
    return R[ind_50], L[ind_50]

def max_drdt(h):

    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    L = L[index:]
    R = h.log_R[index:]
    t = h.star_age[index:]
    center_h1 = h.center_h1[index:]
    ind2 = np.where(center_h1 >= (1e-2))[0][-1]
    L, R, t = L[ind2:], R[ind2:], t[ind2:]
    # he_c = h.center_he4[index:][ind2:]
    ind3 = np.where(R < 2.5)[0]
    L, R, t = L[ind3], R[ind3], t[ind3]
    dr_dt = abs(np.diff(R)/np.diff(t))
    ind = np.where(dr_dt == max(dr_dt))[0]
    # print(ind)
    return R[ind][0], L[ind][0]

# plt.style.use('dark_background')

# ############## FIGURE 1 ##############
l_hay_h = 4.98
l_hay_l = 5.02

h0 = mr.MesaData('data/Z0.001_overshoot_any/LOGS/history.data')
h1 = mr.MesaData('data/Z2e-3_wind1e-3_use_min/LOGS/history.data')
h2 = mr.MesaData('data/Z2e-2/LOGS/history.data')
h3 = mr.MesaData('data/Z0.04_wind0.02_use_min/LOGS/history.data')
# HRDs
T0, R0, L0 = choseMS(h0)
T1, R1, L1 = choseMS(h1)
T2, R2, L2 = choseMS(h2)
T3, R3, L3 = choseMS(h3)
# End of MS
r0_ms, l0_ms, t0_ms = EndOfMS(h0)
r1_ms, l1_ms, t1_ms = EndOfMS(h1)
r2_ms, l2_ms, t2_ms = EndOfMS(h2)
r3_ms, l3_ms, t3_ms = EndOfMS(h3)
R_ms = [r0_ms, r1_ms, r2_ms, r3_ms]
L_ms = [l0_ms, l1_ms, l2_ms, l3_ms]
# He ignition
r0_he, l0_he = He_ignition(h0)
r1_he, l1_he = He_ignition(h1)
r2_he, l2_he = He_ignition(h2)
r3_he, l3_he = He_ignition(h3)
R_he = [r0_he, r1_he, r2_he, r3_he]
L_he = [l0_he, l1_he, l2_he, l3_he]
# Max dlogR/dt
r0_max, l0_max = max_drdt(h0)
r1_max, l1_max = max_drdt(h1)
r2_max, l2_max = max_drdt(h2)
r3_max, l3_max = max_drdt(h3)
R_max = [r0_max, r1_max, r2_max, r3_max]
L_max = [l0_max, l1_max, l2_max, l3_max]
# Hayashi -- constant L
r0_cons = Hayashi(h0,l_hay_l)
r1_cons = Hayashi(h1,l_hay_l)
r2_cons = Hayashi(h2,l_hay_h)
r3_cons = Hayashi(h3,l_hay_h)
R_cons = [r0_cons, r1_cons, r2_cons, r3_cons]
L_cons = [l_hay_l, l_hay_l, l_hay_h, l_hay_h]
# Hayashi -- end of life
r0_end, l0_end = Hayashi_end(h0)
r1_end, l1_end = Hayashi_end(h1)
r2_end, l2_end = Hayashi_end(h2)
r3_end, l3_end = Hayashi_end(h3)
R_end = [r0_end, r1_end, r2_end, r3_end]
L_end = [l0_end, l1_end, l2_end, l3_end]
# Convective fraction = 50 %
p = 50
r0_conv, l0_conv = ConvectiveFraction(h0, p)
r1_conv, l1_conv = ConvectiveFraction(h1, p)
r2_conv, l2_conv = ConvectiveFraction(h2, p)
r3_conv, l3_conv = ConvectiveFraction(h3, p)
R_conv = [r0_conv, r1_conv, r2_conv, r3_conv]
L_conv = [l0_conv, l1_conv, l2_conv, l3_conv]
#
color = plt.cm.Greys(np.linspace(0,1,20))
fig, ax = plt.subplots(ncols=1,nrows=1,sharex=True,sharey=True,figsize=(6.5,4.5))
ax.plot(R0, L0, c=color[10])
ax.plot(R1, L1, c=color[13])
ax.plot(R2, L2, c=color[15])
ax.plot(R3, L3, c=color[19])

ax.scatter(R_ms, L_ms, edgecolors='k', c='cyan',marker='*',s=200, label='End of MS', zorder=5)
ax.scatter(R_he, L_he, edgecolors='k', c='yellow',marker='P',s=150, label='Helium Ignition' '\n' r'$X_C( ^4{\rm He})=90\%$', zorder=5)
ax.scatter(R_max, L_max, edgecolors='k', c='lime',marker='o',s=100, label='Maximum dlogR/dt', zorder=5)
ax.scatter(R_conv, L_conv, edgecolors='k', c='orchid',marker='X',s=100, label=r'$f_{\rm conv,env}=0.5$', zorder=5)
ax.scatter(R_cons, L_cons, edgecolors='k', c='pink',marker='D',s=70, label='Hayashi Track $(L=10^5 L_{\odot})$', zorder=5)

plt.text(1.814,5.04,r'$Z=10^{-3}$',c=color[10],fontsize=15)
plt.text(1.814,4.86,r'$Z=2\times 10^{-3}$',c=color[13],fontsize=15)
plt.text(2.0,4.73,'$Z=0.02$',c=color[15],fontsize=15)
plt.text(1.65,4.46,'$Z=0.04$',c=color[19],fontsize=15)

ax.set_ylim(3.8,5.2)
ax.set_xlim(0.4,3.2)
ax.legend(fontsize=12, ncol=2, framealpha=0.5,loc=8, frameon=False)
# ax.set_xticks(fontsize=20)
# ax.set_yticks(fontsize=20)
ax.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)
plt.tight_layout()
# plt.savefig('figures/HRD_physics.pdf')
plt.show()


#################### ############## ############## ############## ##############

############## FIGURE 2 #################
z_norm = 0.02
fig, ((ax1,ax3),(ax2,ax4)) = plt.subplots(ncols=2,nrows=2,sharex=True,sharey=True,figsize=(11,7.5))
#
# # for z_norm in [0.02,0.0001]:
fig.subplots_adjust(wspace=0, hspace=0)
##### 4-panel HRD with epochs

width_fiducial = 2.5
width_mp = 2.

l_hay_high_z = 4.98
l_hay_low_z = 5.02

p=50

h0 = mr.MesaData('data/Z2e-2/LOGS/history.data')
Teff0, R0, L0 = choseMS(h0)
ax1.plot(R0, L0, color='k', zorder=0., label='$Z=0.02$ (model 1)', linewidth=width_fiducial)
ax2.plot(R0, L0, color='k', zorder=0., linewidth=width_fiducial)

h_mu = mr.MesaData('data/mu_0.04/LOGS/history.data')
Teff_mu, R_mu, L_mu = choseMS(h_mu)
ax1.plot(R_mu, L_mu, color='limegreen', zorder=1, linestyle='--',linewidth=width_mp, label=r'$Z_{\mu}=0.04$ (model 2)')

h_kap = mr.MesaData('data/kap_0.04/LOGS/history.data')
Teff_k, R_k, L_k = choseMS(h_kap)
ax1.plot(R_k, L_k, color='mediumblue', linestyle='--',zorder=1, linewidth=width_mp, label=r'$Z_{\kappa}=0.04$ (model 3)')

h_eps_0 = mr.MesaData('data/eps_f_burn/LOGS/history.data')
c_h1 = h_eps_0.center_h1

Teff_t, R_t, L_t = choseMS(h_eps_0)
ax1.plot(R_t, L_t, color='palevioletred', zorder=1, linestyle='--',linewidth=width_mp, label=r'$Z_{\epsilon}=0.04$ (model 4)')

h_km = mr.MesaData('data/kap_mu_0.04/LOGS/history.data')
Teff_km, R_km, L_km = choseMS(h_km)
ax2.plot(R_km, L_km, color='darkorange', zorder=1, linestyle='--',linewidth=width_mp, label=r'$Z_{\mu,\kappa}=0.04$ (model 5)')

h_em = mr.MesaData('data/eps_mu_f_burn/LOGS/history.data')
Teff_em, R_em, L_em = choseMS(h_em)
ax2.plot(R_em, L_em, color='magenta', zorder=1, linestyle='--',linewidth=width_mp, label=r'$Z_{\mu,\epsilon}=0.04$ (model 6)')

h_ke = mr.MesaData('data/eps_kap_f_burn/LOGS/history.data')
# print(min(h_ke.center_h1))
Teff_ke, R_ke, L_ke = choseMS(h_ke)
ax2.plot(R_ke, L_ke, color='dodgerblue', zorder=1, linestyle='--',linewidth=width_mp, label=r'$Z_{\kappa,\epsilon}=0.04$ (model 7)')

h_3 = mr.MesaData('data/eps_kap_mu_f_burn/LOGS/history.data')
Teff_3, R_3, L_3 = choseMS(h_3)
ax2.plot(R_3, L_3, color='dimgray', zorder=1, linestyle='--',linewidth=width_mp, label=r'$Z_{\mu,\kappa,\epsilon}=0.04$ (model 8)')

h1 = mr.MesaData('data/Z0.04_wind0.02_use_min/LOGS/history.data')
Teff1, R1, L1 = choseMS(h1)
ax1.plot(R1, L1, color='maroon', zorder=0, label='$Z=0.04$ (model 9)', linewidth=width_fiducial)
ax2.plot(R1, L1, color='maroon', zorder=0, linewidth=width_fiducial)

epoch_l = ['EMS','He ignition','Hertzsprung','Convective Env Fraction','Hayashi const L','C depletion']
for epoch in epoch_l:
    if epoch == 'EMS':
        r_0_ms, l_0_ms, t_0_ms = EndOfMS(h0)
        r_m_ms, l_m_ms, t_m_ms = EndOfMS(h_mu)
        r_ke_ms, l_ke_ms, t_ke_ms = EndOfMS(h_ke)
        r_k_ms, l_k_ms, t_k_ms = EndOfMS(h_kap)
        r_t_ms, l_t_ms, t_t_ms = EndOfMS(h_eps_0)
        r_km_ms, l_km_ms, t_km_ms = EndOfMS(h_km)
        r_em_ms, l_em_ms, t_em_ms = EndOfMS(h_em)
        r_3_ms, l_3_ms, t_3_ms = EndOfMS(h_3)
        r_1_ms, l_1_ms, t_1_ms = EndOfMS(h1)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$R_{\rm max, EMS}$'
        ax1.scatter(r_ms_1, l_ms_1, edgecolors='k', c='cyan', marker='*', s=250, zorder=5)
        ax2.scatter(r_ms_2, l_ms_2, edgecolors='k', c='cyan', marker='*', s=250, zorder=5)
    elif epoch == 'He ignition':
        r_0_ms, l_0_ms = He_ignition(h0)
        r_m_ms, l_m_ms = He_ignition(h_mu)
        r_k_ms, l_k_ms = He_ignition(h_kap)
        r_t_ms, l_t_ms = He_ignition(h_eps_0)
        r_ke_ms, l_ke_ms = He_ignition(h_ke)
        r_km_ms, l_km_ms = He_ignition(h_km)
        r_em_ms, l_em_ms = He_ignition(h_em)
        r_3_ms, l_3_ms = He_ignition(h_3)
        r_1_ms, l_1_ms = He_ignition(h1)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$X_c(^4{\rm He})=90\%$'
        ax1.scatter(r_ms_1, l_ms_1, edgecolors='k', c='yellow',marker='P',s=100, zorder=5)
        ax2.scatter(r_ms_2, l_ms_2, edgecolors='k', c='yellow',marker='P',s=100, zorder=5)
    elif epoch == 'Hertzsprung':
        r_0_ms, l_0_ms = max_drdt(h0)
        r_m_ms, l_m_ms = max_drdt(h_mu)
        r_ke_ms, l_ke_ms = max_drdt(h_ke)
        r_k_ms, l_k_ms = max_drdt(h_kap)
        r_t_ms, l_t_ms = max_drdt(h_eps_0)
        r_km_ms, l_km_ms = max_drdt(h_km)
        r_em_ms, l_em_ms = max_drdt(h_em)
        r_3_ms, l_3_ms = max_drdt(h_3)
        r_1_ms, l_1_ms = max_drdt(h1)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = 'Maximum $dlogR/dt$'
        ax1.scatter(r_ms_1,l_ms_1, edgecolors='k', c='lime',marker='o',s=100,zorder=5)
        ax2.scatter(r_ms_2,l_ms_2, edgecolors='k', c='lime',marker='o',s=100,zorder=5)
    elif epoch == 'Convective Env Fraction':
        r_0_ms, l_0_ms = ConvectiveFraction(h0,p)
        r_m_ms, l_m_ms = ConvectiveFraction(h_mu,p)
        r_ke_ms, l_ke_ms = ConvectiveFraction(h_ke,p)
        r_k_ms, l_k_ms = ConvectiveFraction(h_kap,p)
        r_t_ms, l_t_ms = ConvectiveFraction(h_eps_0,p)
        r_km_ms, l_km_ms = ConvectiveFraction(h_km,p)
        r_em_ms, l_em_ms = ConvectiveFraction(h_em,p)
        r_3_ms, l_3_ms = ConvectiveFraction(h_3,p)
        r_1_ms, l_1_ms = ConvectiveFraction(h1,p)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$f_{\rm conv, env}=0.5$'
        ax1.scatter(r_ms_1,l_ms_1, edgecolors='k', c='orchid',marker='X',s=100,zorder=5)
        ax2.scatter(r_ms_2,l_ms_2, edgecolors='k', c='orchid',marker='X',s=100,zorder=5)
    elif epoch == 'Hayashi const L':
        r_0_ms = Hayashi(h0,l_hay_high_z)
        r_m_ms = Hayashi(h_mu,l_hay_high_z)
        r_ke_ms = Hayashi(h_ke,l_hay_high_z)
        r_k_ms = Hayashi(h_kap,l_hay_high_z)
        r_t_ms = Hayashi(h_eps_0,l_hay_high_z)
        r_km_ms = Hayashi(h_km,l_hay_high_z)
        r_em_ms = Hayashi(h_em,l_hay_high_z)
        r_3_ms = Hayashi(h_3,l_hay_high_z)
        r_1_ms = Hayashi(h1,l_hay_high_z)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_hay_high_z, l_hay_high_z, l_hay_high_z, l_hay_high_z, l_hay_high_z])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_hay_high_z, l_hay_high_z, l_hay_high_z, l_hay_high_z, l_hay_high_z, l_hay_high_z])
        labell = '$L=10^{5}L_{\odot}$'
        ax1.scatter(r_ms_1,l_ms_1, edgecolors='k', c='pink',marker='D',s=80,zorder=5)
        ax2.scatter(r_ms_2,l_ms_2, edgecolors='k', c='pink',marker='D',s=80,zorder=5)
    elif epoch == 'C depletion':
        r_0_ms, l_0_ms = Hayashi_end(h0)
        r_m_ms, l_m_ms = Hayashi_end(h_mu)
        r_ke_ms, l_ke_ms = Hayashi_end(h_ke)
        r_k_ms, l_k_ms = Hayashi_end(h_kap)
        r_t_ms, l_t_ms = Hayashi_end(h_eps_0)
        r_km_ms, l_km_ms = Hayashi_end(h_km)
        r_em_ms, l_em_ms = Hayashi_end(h_em)
        r_3_ms, l_3_ms = Hayashi_end(h_3)
        r_1_ms, l_1_ms = Hayashi_end(h1)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$X_{c}(^{12}C)<10^{-8}$'


ax1.set_ylim(4.1,5.3)
ax1.set_xlim(0.5,3.2)
ax2.set_ylim(4.1,5.3)
ax2.set_xlim(0.5,3.2)

ax1.legend(fontsize=12, ncol=2, framealpha=0.5,loc=2, frameon=False)
ax2.legend(fontsize=12, ncol=1, framealpha=0.5,loc=2, frameon=False)

ax1.set_yticks([4.2,4.4,4.6,4.8,5.0,5.2])
ax2.set_yticks([4.2,4.4,4.6,4.8,5.0,5.2])
ax1.set_xticks([0.5,1.,1.5,2.,2.5,3.])
ax2.set_xticks([0.5,1.,1.5,2.,2.5,3.])


ax1.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax1.set_title('High Z',fontsize=20)
ax2.set_ylabel('log$_{10}$(L/L$_{\odot}$)',fontsize=20)
ax2.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)


h0l = mr.MesaData('data/Z0.001_overshoot_any/LOGS/history.data')
Teff0l, R0l, L0l = choseMS(h0l)
ax3.plot(R0l, L0l, color='k', zorder=0, label='$Z=10^{-3}$ (model 10)', linewidth=width_fiducial)
ax4.plot(R0l, L0l, color='k', zorder=0., linewidth=width_fiducial)

hl_mu = mr.MesaData('data/Z1e-3_Zmu_2e-3/LOGS/history.data')
Teff_mul, R_mul, L_mul = choseMS(hl_mu)
ax3.plot(R_mul, L_mul, color='limegreen', zorder=1, linestyle='--', linewidth=width_mp, label=r'$Z_{\mu}=2\times10^{-3}$ (model 11)')

hl_kap = mr.MesaData('data/Z1e-3_Zkap_2e-3/LOGS/history.data')
Teff_kl, R_kl, L_kl = choseMS(hl_kap)
ax3.plot(R_kl, L_kl, color='mediumblue', zorder=1, linestyle='--', linewidth=width_mp, label=r'$Z_{\kappa}=2\times10^{-3}$ (model 12)')

hl_eps_0 = mr.MesaData('data/Z1e-3_Zeps_2e-3/LOGS/history.data')
Teff_tl, R_tl, L_tl = choseMS(hl_eps_0)
ax3.plot(R_tl, L_tl, color='palevioletred', zorder=1,linestyle='--',  linewidth=width_mp, label=r'$Z_{\epsilon}=2\times10^{-3}$ (model 13)')

hl_km = mr.MesaData('data/Z1e-3_Zkap_mu_2e-3/LOGS/history.data')
Teff_kml, R_kml, L_kml = choseMS(hl_km)
ax4.plot(R_kml, L_kml, color='darkorange', zorder=1, linestyle='--', linewidth=width_mp, label=r'$Z_{\mu,\kappa}=2\times10^{-3}$ (model 14)')

hl_em = mr.MesaData('data/Z1e-3_Zmu_eps_2e-3/LOGS/history.data')
Teff_eml, R_eml, L_eml = choseMS(hl_em)
ax4.plot(R_eml, L_eml, color='magenta', zorder=1, linestyle='--', linewidth=width_mp, label=r'$Z_{\mu,\epsilon}=2\times10^{-3}$ (model 15)')

hl_ke = mr.MesaData('data/Z1e-3_Zkap_eps_2e-3/LOGS/history.data')
# print(len(hl_ke.center_c12))
# print(hl_ke.center_c12)
Teff_kel, R_kel, L_kel = choseMS(hl_ke)
ax4.plot(R_kel, L_kel, color='dodgerblue', zorder=1, linestyle='--', linewidth=width_mp, label=r'$Z_{\kappa,\epsilon}=2\times10^{-3}$ (model 16)')

hl_3 = mr.MesaData('data/Z1e-3_Zkap_mu_eps_2e-3/LOGS/history.data')
Teff_3l, R_3l, L_3l = choseMS(hl_3)
ax4.plot(R_3l, L_3l, color='dimgray', zorder=1, linestyle='--', linewidth=width_mp, label=r'$Z_{\mu,\kappa,\epsilon}=2\times10^{-3}$ (model 17)')

h1l = mr.MesaData('data/Z2e-3_wind1e-3_use_min/LOGS/history.data')

Teff1l, R1l, L1l = choseMS(h1l)
ax3.plot(R1l, L1l, color='maroon', zorder=0, label=r'$Z=2\times10^{-3}$ (model 18)', linewidth=width_fiducial)
ax4.plot(R1l, L1l, color='maroon', zorder=0, linewidth=width_fiducial)

for epoch in epoch_l:
    if epoch == 'EMS':
        r_0_ms, l_0_ms, t_0_ms = EndOfMS(h0l)
        r_m_ms, l_m_ms, t_m_ms = EndOfMS(hl_mu)
        r_ke_ms, l_ke_ms, t_ke_ms = EndOfMS(hl_ke)
        r_k_ms, l_k_ms, t_k_ms = EndOfMS(hl_kap)
        r_t_ms, l_t_ms, t_t_ms = EndOfMS(hl_eps_0)
        r_km_ms, l_km_ms, t_km_ms = EndOfMS(hl_km)
        r_em_ms, l_em_ms, t_em_ms = EndOfMS(hl_em)
        r_3_ms, l_3_ms, t_3_ms = EndOfMS(hl_3)
        r_1_ms, l_1_ms, t_1_ms = EndOfMS(h1l)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$R_{\rm max, EMS}$'
        ax3.scatter(r_ms_1, l_ms_1, edgecolors='k', c='cyan', marker='*', s=250, zorder=5)
        ax4.scatter(r_ms_2, l_ms_2, edgecolors='k', c='cyan', marker='*', s=250, zorder=5)
    elif epoch == 'He ignition':
        r_0_ms, l_0_ms = He_ignition(h0l)
        r_m_ms, l_m_ms = He_ignition(hl_mu)
        r_ke_ms, l_ke_ms = He_ignition(hl_ke)
        r_k_ms, l_k_ms = He_ignition(hl_kap)
        r_t_ms, l_t_ms = He_ignition(hl_eps_0)
        r_km_ms, l_km_ms = He_ignition(hl_km)
        r_em_ms, l_em_ms = He_ignition(hl_em)
        r_3_ms, l_3_ms = He_ignition(hl_3)
        r_1_ms, l_1_ms = He_ignition(h1l)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$X_c(^4{\rm He})=90\%$'
        ax3.scatter(r_ms_1, l_ms_1, edgecolors='k', c='yellow',marker='P',s=100, zorder=5)
        ax4.scatter(r_ms_2, l_ms_2, edgecolors='k', c='yellow',marker='P',s=100, zorder=5)
    elif epoch == 'Hertzsprung':
        r_0_ms, l_0_ms = max_drdt(h0l)
        r_m_ms, l_m_ms = max_drdt(hl_mu)
        r_ke_ms, l_ke_ms = max_drdt(hl_ke)
        r_k_ms, l_k_ms = max_drdt(hl_kap)
        r_t_ms, l_t_ms = max_drdt(hl_eps_0)
        r_km_ms, l_km_ms = max_drdt(hl_km)
        r_em_ms, l_em_ms = max_drdt(hl_em)
        r_3_ms, l_3_ms = max_drdt(hl_3)
        r_1_ms, l_1_ms = max_drdt(h1l)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = 'Maximum $dlogR/dt$'
        ax3.scatter(r_ms_1,l_ms_1, edgecolors='k', c='lime',marker='o',s=100,zorder=5)
        ax4.scatter(r_ms_2,l_ms_2, edgecolors='k', c='lime',marker='o',s=100,zorder=5)
    elif epoch == 'Convective Env Fraction':
        r_0_ms, l_0_ms = ConvectiveFraction(h0l,p)
        r_m_ms, l_m_ms = ConvectiveFraction(hl_mu,p)
        r_ke_ms, l_ke_ms = ConvectiveFraction(hl_ke,p)
        r_k_ms, l_k_ms = ConvectiveFraction(hl_kap,p)
        r_t_ms, l_t_ms = ConvectiveFraction(hl_eps_0,p)
        r_km_ms, l_km_ms = ConvectiveFraction(hl_km,p)
        r_em_ms, l_em_ms = ConvectiveFraction(hl_em,p)
        r_3_ms, l_3_ms = ConvectiveFraction(hl_3,p)
        r_1_ms, l_1_ms = ConvectiveFraction(h1l,p)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$f_{\rm conv, env}=0.5$'
        ax3.scatter(r_ms_1,l_ms_1, edgecolors='k', c='orchid',marker='X',s=100,zorder=5)
        ax4.scatter(r_ms_2,l_ms_2, edgecolors='k', c='orchid',marker='X',s=100,zorder=5)
    elif epoch == 'Hayashi const L':
        r_0_ms = Hayashi(h0l,l_hay_low_z)
        r_m_ms = Hayashi(hl_mu,l_hay_low_z)
        r_ke_ms = Hayashi(hl_ke,l_hay_low_z)
        r_k_ms = Hayashi(hl_kap,l_hay_low_z)
        r_t_ms = Hayashi(hl_eps_0,l_hay_low_z)
        r_km_ms = Hayashi(hl_km,l_hay_low_z)
        r_em_ms = Hayashi(hl_em,l_hay_low_z)
        r_3_ms = Hayashi(hl_3,l_hay_low_z)
        r_1_ms = Hayashi(h1l,l_hay_low_z)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_hay_low_z, l_hay_low_z, l_hay_low_z, l_hay_low_z, l_hay_low_z])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_hay_low_z, l_hay_low_z, l_hay_low_z, l_hay_low_z, l_hay_low_z, l_hay_low_z])
        labell = '$L=10^{5}L_{\odot}$'
        ax3.scatter(r_ms_1,l_ms_1, edgecolors='k', c='pink',marker='D',s=80,zorder=5)
        ax4.scatter(r_ms_2,l_ms_2, edgecolors='k', c='pink',marker='D',s=80,zorder=5)
    elif epoch == 'C depletion':
        r_0_ms, l_0_ms = Hayashi_end(h0l)
        r_m_ms, l_m_ms = Hayashi_end(hl_mu)
        r_ke_ms, l_ke_ms = Hayashi_end(hl_ke)
        r_k_ms, l_k_ms = Hayashi_end(hl_kap)
        r_t_ms, l_t_ms = Hayashi_end(hl_eps_0)
        r_km_ms, l_km_ms = Hayashi_end(hl_km)
        r_em_ms, l_em_ms = Hayashi_end(hl_em)
        r_3_ms, l_3_ms = Hayashi_end(hl_3)
        r_1_ms, l_1_ms = Hayashi_end(h1l)
        r_ms_1 = np.array([r_0_ms, r_1_ms, r_m_ms, r_k_ms, r_t_ms])
        l_ms_1 = np.array([l_0_ms, l_1_ms, l_m_ms, l_k_ms, l_t_ms])
        r_ms_2 = np.array([r_0_ms, r_1_ms, r_km_ms, r_em_ms, r_ke_ms, r_3_ms])
        l_ms_2 = np.array([l_0_ms, l_1_ms, l_km_ms, l_em_ms, l_ke_ms, l_3_ms])
        labell = r'$X_{c}(^{12}C)<10^{-8}$'


ax3.legend(fontsize=12, ncol=1, framealpha=0.5,loc=4, frameon=False)
ax4.legend(fontsize=12, ncol=1, framealpha=0.5,loc=4, frameon=False)

ax3.set_xticks([0.5,1.,1.5,2.,2.5,3.])
ax3.set_title('Low Z',fontsize=20)
ax4.set_xticks([0.5,1.,1.5,2.,2.5,3.])

ax4.set_xlabel('log$_{10}$(R/R$_{\odot}$)',fontsize=20)

plt.tight_layout()
# plt.savefig('figures/HR-4panel.pdf')
plt.show()