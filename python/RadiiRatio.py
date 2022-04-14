import matplotlib.pyplot as plt
import numpy as np
import mesa_reader as mr
import matplotlib

plt.rcParams["font.family"] = "Times New Roman"
font = {'fontname':'Times New Roman'}

Msun = 1.989 * 10**33 # solar mass in kg.

def EndOfMS(h,Z_eps,Z_fid):
    # compute radius, luminosity and effective temperature of star ONLY on the MS
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    r = h.log_R[index:]
    center_h1 = h.center_h1[index:]
    ind2 = np.where(center_h1 > (1e-3 + (Z_eps - Z_fid)))[0]
    return 10 ** np.max(r[ind2])

def He_ignition(h,Z_eps,Z_fid):
    # compute radius, luminosity and effective temperature at He ignition
    center_h1 = h.center_h1
    ind_tams = np.where(center_h1 <= (1e-3 + (Z_eps - Z_fid)))[0] # this is after TAMS
    center_he = h.center_he4[ind_tams]
    r = h.log_R[ind_tams]
    age = h.star_age[ind_tams]
    # he4, ind_he = find_nearest(center_he, 0.90)
    # r = r[:ind_he]
    # ind_max = np.where(r == max(r))[0]
    ind_he = np.where(center_he <= 0.90)[0] # he ignition
    # return 10 ** r[ind_he][0], age[ind_he][0]
    return 10 ** r[ind_he][0]

# h = mr.MesaData('data/Z2e-2/LOGS/history.data')
# r, age = He_ignition(h,0.02,0.02)
# print(r)

def Hayashi(h,L_const):
    # compute radius, luminosity and effective temperature on the hayashi track
    radius = h.log_R
    L = h.log_L
    Teff = h.log_Teff
    c1_top, c1_bot = h.conv_mx1_top_r, h.conv_mx1_bot_r
    c2_top, c2_bot = h.conv_mx2_top_r, h.conv_mx2_bot_r
    R_int = np.interp(L_const,L,radius)
    T_int = np.interp(L_const,L,Teff)
    # print(T_int)
    c1_top_int, c1_bot_int = np.interp(L_const,L,c1_top), np.interp(L_const,L,c1_bot)
    c2_top_int, c2_bot_int = np.interp(L_const,L,c2_top), np.interp(L_const,L,c2_bot)
    R_conv = abs(c1_top_int-c1_bot_int)+abs(c2_top_int-c2_bot_int)
    # print(R_conv/(10 ** R_int))

    # print(T_int)
    return 10 ** R_int

# def Hayashi_end(h):
#     radius = h.log_R
#     return 10 ** max(radius)

def Hayashi_end(h):
    L = h.log_L
    Lnuc = h.log_Lnuc
    ind_ms = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    c12_core, o16_core = h.center_c12[ind_ms:], h.center_o16[ind_ms:]
    R, L = h.log_R[ind_ms:], h.log_L[ind_ms:]
    ind_c_dep = np.where((c12_core<0.001) & (o16_core>0.4))[0]
    # print(ind_c_dep)
    R_cd = R[ind_c_dep][0]
    # print(R_cd)
    return 10 ** R_cd


def max_drdt(h,Z_eps,Z_fid):
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    L = L[index:]
    R = h.log_R[index:]
    t = h.star_age[index:]
    center_h1 = h.center_h1[index:]
    ind2 = np.where(center_h1 < (1e-3 + (Z_eps - Z_fid)))[0]
    L, R, t = L[ind2], R[ind2], t[ind2]
    ind3 = np.where(R < 2.5)[0]
    L, R, t = L[ind3], R[ind3], t[ind3]
    dr_dt = abs(np.diff(R)/np.diff(t))
    ind = np.where(dr_dt == max(dr_dt))[0]
    # print(ind)
    return 10 ** R[ind][0]

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
    m_envel = []
    for l in range(len(R)):
        m_c_tot = m_conv_tot[l] / Msun
        m_core = m_conv_core[l]
        if m_c_tot >= m_core:
            m_envel.append(m_c_tot - m_core)
        else:
            m_envel.append(0)

    conv_frac = np.array(m_envel) / (m_tot - m_he_core)
    f_50, ind_50 = find_nearest(conv_frac, p/100)
    return R[ind_50]


# z_kap = [0.0001,0.0002,0.001,0.01,0.02,0.03,0.04,0.05]
# h = mr.MesaData('data/15m_ze-4_null_eos/LOGS/history.data')
# r_k = EndOfMS(h)
# r_he = He_ignition(h)
# L_const = 4.97
# r_hay = Hayashi(h, L_const)
#
# r_0_ms, r_kap_ms, r_mu_ms = [1.], [1.], [1.]
# r_0_he, r_kap_he, r_mu_he = [1.], [1.], [1.]
# r_0_hay, r_kap_hay, r_mu_hay = [1.], [1.], [1.]
#
# for z in z_kap[1:]:
#     h_0 = mr.MesaData('data/15m_z{}_x0.75/LOGS/history.data'.format(float(z)))
#     h_kap = mr.MesaData('data/kap_changeY_{}/LOGS/history.data'.format(float(z)))
#     h_mu = mr.MesaData('data/eos_{}/LOGS/history.data'.format(float(z)))
#
#     r_0_ms.append(10 ** EndOfMS(h_0)/10 ** r_k)
#     r_kap_ms.append(10 ** EndOfMS(h_kap)/10 ** r_k)
#     r_mu_ms.append(10 ** EndOfMS(h_mu)/10 ** r_k)
#
#     r_0_he.append(10 ** He_ignition(h_0)/10 ** r_he)
#     r_kap_he.append(10 ** He_ignition(h_kap) / 10 ** r_he)
#     r_mu_he.append(10 ** He_ignition(h_mu) / 10 ** r_he)
#
#     r_0_hay.append(10 ** Hayashi(h_0, L_const)/10 ** r_hay)
#     r_kap_hay.append(10 ** Hayashi(h_kap, L_const) / 10 ** r_hay)
#     r_mu_hay.append(10 ** Hayashi(h_mu, L_const) / 10 ** r_hay)
#
# # h_eps = mr.MesaData('data/eps_2e-4_dt1/LOGS/history.data')
# # r_eps_ms = [1., 10 ** EndOfMS(h_eps) / 10 ** r_k]
# # z_eps = [0.0001,0.0002]
# # plt.plot(z_kap,r_kap_ms,'-bo',label='r$_{\kappa}$/r$_0$')
# # plt.plot(z_kap,r_mu_ms,'-go',label='r$_{\mu}$/r$_0$')
# # plt.plot(z_eps,r_eps_ms,'-ko',label='r$_{\epsilon}$/r$_0$ + extrapolation')
# # plt.plot(z_kap,r_0_ms,'-ro',label='r$_z$/r$_0$')
# # plt.ylim(1,50)
# # plt.title('End of Main Sequence')
#
# # plt.plot(z_kap,r_kap_he,'-bo',label='r$_{\kappa}$/r$_0$')
# # plt.plot(z_kap,r_mu_he,'-go',label='r$_{\mu}$/r$_0$')
# # plt.plot(z_kap[:-1],r_0_he[:-1],'-ro',label='r$_z$/r$_0$')
# # plt.title('He Ignition')
#
# pht_data = np.loadtxt('R_phot.txt')
# z_analy, r_analy = pht_data[:,0], pht_data[:,1]
# plt.plot(z_kap,r_kap_hay,'-bo',label='r$_{\kappa}$/r$_0$')
# plt.plot(z_kap,r_mu_hay,'-go',label='r$_{\mu}$/r$_0$')
# plt.plot(z_kap[:-1],r_0_hay[:-1],'-ro',label='r$_z$/r$_0$')
# plt.plot(z_analy,r_analy / r_analy[0],c='k',label='analytical soln.')
# plt.title('Hayashi Track')
#
# plt.xscale('log')
# plt.xlabel('Z',fontsize=13)
# plt.ylabel('Ratios of Radii',fontsize=13,labelpad=10)
# plt.xticks([0.0001,0.0002,0.001,0.01,0.03,0.05],labels=['1e-4','2e-4','1e-3','0.01','0.03','0.05'])
# plt.legend(fontsize=15)
#
# # plt.savefig('figures/r_ratio_hayashi.pdf')
# plt.show()
#
