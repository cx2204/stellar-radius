import mesa_reader as mr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from scipy.interpolate import interp1d

Msun = 1.989 * 10**33 # solar mass in kg.

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
    return r[ind2][ind_max],L_ms[ind2][ind_max]

def He_ignition(h):
    # compute radius, luminosity and effective temperature at He ignition
    center_h1 = h.center_h1
    ind_tams = np.where(center_h1 <= 1e-3)[0] # this is after TAMS
    center_he = h.center_he4[ind_tams]
    r = h.log_R[ind_tams]
    age = h.star_age[ind_tams]
    L = h.log_L[ind_tams]
    # he4, ind_he = find_nearest(center_he, 0.90)
    # r = r[:ind_he]
    # ind_max = np.where(r == max(r))[0]
    ind_he = np.where(center_he <= 0.90)[0] # he ignition
    # return 10 ** r[ind_he][0], age[ind_he][0]
    return r[ind_he][0], L[ind_he][0]

# def He_ignition(h):
#
#     # compute radius, luminosity and effective temperature at He ignition
#     center_h1 = h.center_h1
#     ind_tams = np.where(center_h1 <= (1e-3))[0] # this is after TAMS
#     center_he = h.center_he4[ind_tams]
#     r = h.log_R[ind_tams]
#     lum = h.log_L[ind_tams]
#     he4, ind_he = find_nearest(center_he, 0.9)
#     r = r[:ind_he]
#     lum = lum[:ind_he]
#     center_he = center_he[:ind_he]
#     ind_max = np.where(r == max(r))[0]
#     # print(center_he[ind_max], age[ind_max])
#     return r[ind_max], lum[ind_max]

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
    return R_int, L_const

def choseMS(h):
    Teff = h.log_Teff
    L = h.log_L
    Lnuc = h.log_Lnuc
    index = np.where(np.round(L, 3) == np.round(Lnuc, 3))[0][0]
    Teff = Teff[index:]
    L = L[index:]
    R = h.log_R[index:]
    t = h.star_age[index:]
    # center_h1 = h.center_h1[index:]
    # ind2 = np.where(center_h1 < (1e-4 + (Z_eps - Z_fid)))[0]
    # return Teff[ind2], R[ind2], L[ind2]
    return R, L

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
    ind2 = np.where(center_h1 < 1e-3)[0]
    L, R, t = L[ind2], R[ind2], t[ind2]
    ind3 = np.where(R < 2.5)[0]
    L, R, t = L[ind3], R[ind3], t[ind3]
    dr_dt = abs(np.diff(R)/np.diff(t))
    ind = np.where(dr_dt == max(dr_dt))[0]
    # print(ind)
    return R[ind][0], L[ind][0]

# h1 = mr.MesaData('data/high_Z_2Z/Z2e-2/LOGS/history.data')
# print(h1.bulk_names)
# h2 = mr.MesaData('data/low_Z_2Z/Z1e-3/LOGS/history.data')
# R1, L1 = choseMS(h1)
# R2, L2 = choseMS(h2)
#
# plt.plot(R1, L1, label='high Z')
# plt.plot(R2, L2, label='low Z')
# # plt.xscale('log')
# plt.legend()
# plt.show()