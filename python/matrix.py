import matplotlib.pyplot as plt
import numpy as np
import mesa_reader as mr
from RadiiRatio import EndOfMS, He_ignition, Hayashi, Hayashi_end, max_drdt, ConvectiveFraction
import pandas as pd

z = 0.02 # Z=0.001 for low-Z matrices

def entry(r_h_til,r_h,r_l):
    return ((r_h_til - r_l) / r_l) * 100
    # return (r_h_til - r_l) / (r_h - r_l)

if z == 0.001:
    L_c = 5.02
    h0 = mr.MesaData('data/Z1e-3/LOGS/history.data') # load low Z model
    #
    # ---------- Z=1e-3 (null matrices)
    # h1 = mr.MesaData('data/Low_Z_null/Z1e-3_wind1e-3_use_min/LOGS/history.data') # load high Z model
    # h_ee = mr.MesaData('data/Low_Z_null/Z1e-3_eps/LOGS/history.data')
    # h_kk = mr.MesaData('data/Low_Z_null/Z1e-3_kap/LOGS/history.data')
    # h_mm = mr.MesaData('data/Low_Z_null/Z1e-3_mu/LOGS/history.data')
    # h_mk = mr.MesaData('data/Low_Z_null/Z1e-3_kap_mu/LOGS/history.data')
    # h_em = mr.MesaData('data/Low_Z_null/Z1e-3_mu_eps/LOGS/history.data')
    # h_ek = mr.MesaData('data/Low_Z_null/Z1e-3_kap_eps/LOGS/history.data')
    # #
    # h_ekm = mr.MesaData('data/Low_Z_null/Z1e-3_kap_mu_eps/LOGS/history.data')

    # ---------- Z=2e-3 (result matrices)
    h1 = mr.MesaData('data/low_Z_2Z/Z2e-3_wind1e-3_use_min/LOGS/history.data') # load high Z model
    # h1 = mr.MesaData('data/Z_V/LOGS/history.data') # change h1 file to produce numbers in Table 2 and 3
    h_ee = mr.MesaData('data/low_Z_2Z/Z1e-3_Zeps_2e-3/LOGS/history.data')
    h_kk = mr.MesaData('data/low_Z_2Z/Z1e-3_Zkap_2e-3/LOGS/history.data')
    h_mm = mr.MesaData('data/low_Z_2Z/Z1e-3_Zmu_2e-3/LOGS/history.data')
    h_mk = mr.MesaData('data/low_Z_2Z/Z1e-3_Zkap_mu_2e-3/LOGS/history.data')
    h_em = mr.MesaData('data/low_Z_2Z/Z1e-3_Zmu_eps_2e-3/LOGS/history.data')
    h_ek = mr.MesaData('data/low_Z_2Z/Z1e-3_Zkap_eps_2e-3/LOGS/history.data')
    #
    h_ekm = mr.MesaData('data/low_Z_2Z/Z1e-3_Zkap_mu_eps_2e-3/LOGS/history.data')
    
    ############## uncomment the epoch to produce the corresponding matrix

    # ######### End of MS matrix
    # epoch = 'MS'
    # R_l, R_h = EndOfMS(h0,1e-3,1e-3), EndOfMS(h1,1e-3,1e-3) # scalars
    # #matrix entries
    # r_mm, r_kk, r_ee = EndOfMS(h_mm,1e-3,1e-3), EndOfMS(h_kk,1e-3,1e-3), EndOfMS(h_ee,1e-3,1e-3)
    # r_mk, r_em, r_ek = EndOfMS(h_mk,1e-3,1e-3), EndOfMS(h_em,1e-3,1e-3), EndOfMS(h_ek,1e-3,1e-3)
    # r_ekm = EndOfMS(h_ekm,1e-3,1e-3)

    # ######### max dlogR/dt matrix (Hertzsprung Gap)
    # epoch = 'Hertzsprung Gap'
    # R_l, R_h = max_drdt(h0,1e-3,1e-3), max_drdt(h1,1e-3,1e-3) # scalars
    # r_mm, r_kk = max_drdt(h_mm,1e-3,1e-3), max_drdt(h_kk,1e-3,1e-3)
    # r_ee = max_drdt(h_ee,1e-3,1e-3)
    # r_mk = max_drdt(h_mk,1e-3,1e-3)
    # r_em, r_ek = max_drdt(h_em,1e-3,1e-3), max_drdt(h_ek,1e-3,1e-3)
    # r_ekm = max_drdt(h_ekm,1e-3,1e-3)

    # ######### He ignition matrix
    # epoch = 'He ignition'
    # R_l, R_h = He_ignition(h0,1e-3,1e-3), He_ignition(h1,1e-3,1e-3) # scalars
    # r_mm, r_kk = He_ignition(h_mm,1e-3,1e-3), He_ignition(h_kk,1e-3,1e-3)
    # r_ee = He_ignition(h_ee,1e-3,1e-3)
    # r_mk = He_ignition(h_mk,1e-3,1e-3)
    # r_em, r_ek = He_ignition(h_em,1e-3,1e-3), He_ignition(h_ek,1e-3,1e-3)
    # r_ekm = He_ignition(h_ekm,1e-3,1e-3)

    # ######### convective envelope matrix (Early Hayashi)
    # epoch = 'Early Hayashi'
    # p = 50
    # R_l, R_h = ConvectiveFraction(h0,p), ConvectiveFraction(h1,p) # scalars
    # r_mm, r_kk = ConvectiveFraction(h_mm,p), ConvectiveFraction(h_kk,p)
    # r_ee = ConvectiveFraction(h_ee,p)
    # r_mk = ConvectiveFraction(h_mk,p)
    # r_em, r_ek = ConvectiveFraction(h_em,p), ConvectiveFraction(h_ek,p)
    # r_ekm = ConvectiveFraction(h_ekm,p)

    # ######### Hayashi matrix -- constant L
    epoch = 'Hayashi track L=1e5Lsun'
    R_l, R_h = Hayashi(h0,L_c), Hayashi(h1,L_c) # scalars
    r_mm, r_kk, r_ee = Hayashi(h_mm,L_c), Hayashi(h_kk,L_c), Hayashi(h_ee,L_c)
    r_mk, r_em, r_ek = Hayashi(h_mk,L_c), Hayashi(h_em,L_c), Hayashi(h_ek,L_c)
    r_ekm = Hayashi(h_ekm,L_c)

    # construct matrix
    T_r = [(entry(r_mm, R_h, R_l), entry(r_mk, R_h, R_l), entry(r_em, R_h, R_l)),
           ('-', entry(r_kk, R_h, R_l), entry(r_ek, R_h, R_l)),
           ('-', '-', entry(r_ee, R_h, R_l))]

    dfObj = pd.DataFrame(T_r, columns=['μ', 'κ', 'ε']
                         , index=['μ', 'κ', 'ε'])

    print("Z=2d-3; {}".format(str(epoch)), dfObj, sep='\n')
    print("Δ_{μ,κ,ε} =",entry(r_ekm, R_h, R_l))
    print("Δ_μ+Δ_κ+Δ_ε =",entry(r_mm, R_h, R_l) + entry(r_kk, R_h, R_l) + entry(r_ee, R_h, R_l))
    print("Δ_2Z =",100 * (R_h - R_l) / R_l)

if z == 0.02:
    L_c = 4.98
    h0 = mr.MesaData('data/Z2e-2/LOGS/history.data') # load low Z model

    # # ---------- Z=0.04 (result matrices)

    # h1 = mr.MesaData('data/Z0.04_Zwind0.04/LOGS/history.data')  # load high Z model
    # h1 = mr.MesaData('data/Z0.04_Zwind0.0/LOGS/history.data')  # load high Z model
    # h1 = mr.MesaData('data/Z_II/LOGS/history.data')  # load high Z model
    # h1 = mr.MesaData('data/Zwind_0.04/LOGS/history.data')  # load high Z model
    h1 = mr.MesaData('data/high_Z_2Z/Z0.04_wind0.02_use_min/LOGS/history.data')  # load high Z model
    # h1 = mr.MesaData('data/high_Z_2Z/Z0.04_hardcode_wind0.0/LOGS/history.data')  # load high Z model
    h_ee = mr.MesaData('data/high_Z_2Z/Z2e-2_Zeps4e-2/LOGS/history.data') # finished at C-depletion
    h_kk = mr.MesaData('data/high_Z_2Z/Z2e-2_Zkap4e-2/LOGS/history.data') # X(C12) at the end was <1e-6
    h_mm = mr.MesaData('data/high_Z_2Z/Z2e-2_Zmu4e-2/LOGS/history.data') # X(C12) at the end was 2e-6
    h_mk = mr.MesaData('data/high_Z_2Z/Z2e-2_Zkap_mu4e-2/LOGS/history.data') # finished at C-depletion
    h_em = mr.MesaData('data/high_Z_2Z/Z2e-2_Zmu_eps4e-2/LOGS/history.data') # finished at C-depletion
    h_ek = mr.MesaData('data/high_Z_2Z/Z2e-2_Zkap_eps4e-2/LOGS/history.data') # finished at C-depletion
    h_ekm = mr.MesaData('data/high_Z_2Z/Z2e-2_Zkap_mu_eps4e-2/LOGS/history.data') # finished at C-depletion

    # ---------- Z=0.02 (null matrices)
    # h1 = mr.MesaData('data/High_Z_null/Z0.02_wind0.02_use_min/LOGS/history.data')  # load high Z model
    # h_ee = mr.MesaData('data/High_Z_null/Z2e-2_Zeps2e-2/LOGS/history.data')
    # h_kk = mr.MesaData('data/High_Z_null/Z2e-2_Zkap2e-2/LOGS/history.data')
    # h_mm = mr.MesaData('data/High_Z_null/Z2e-2_Zmu2e-2/LOGS/history.data')
    #
    # h_mk = mr.MesaData('data/High_Z_null/Z2e-2_Zkap_mu2e-2/LOGS/history.data')
    # h_em = mr.MesaData('data/High_Z_null/Z2e-2_Zmu_eps2e-2/LOGS/history.data')
    # h_ek = mr.MesaData('data/High_Z_null/Z2e-2_Zkap_eps2e-2/LOGS/history.data')
    # #
    # h_ekm = mr.MesaData('data/High_Z_null/Z2e-2_Zkap_mu_eps2e-2/LOGS/history.data')
    
    ############## uncomment the epoch to produce the corresponding matrix

    # ######### End of MS matrix
    # epoch = 'MS'
    # R_l, R_h = EndOfMS(h0,0.02,0.02), EndOfMS(h1,0.02,0.02) # scalars
    # #matrix entries
    # r_mm, r_kk, r_ee = EndOfMS(h_mm,0.02,0.02), EndOfMS(h_kk,0.02,0.02), EndOfMS(h_ee,0.02,0.02)
    # r_mk, r_em, r_ek = EndOfMS(h_mk,0.02,0.02), EndOfMS(h_em,0.02,0.02), EndOfMS(h_ek,0.02,0.02)
    # r_ekm = EndOfMS(h_ekm,0.02,0.02)

    # ######### max dlogR/dt matrix (Hertzsprung Gap)
    epoch = 'Hertzsprung Gap'
    R_l, R_h = max_drdt(h0,0.02,0.02), max_drdt(h1,0.04,0.04) # scalars
    r_mm, r_kk = max_drdt(h_mm,0.02,0.02), max_drdt(h_kk,0.02,0.02)
    r_ee = max_drdt(h_ee,0.02,0.02)
    r_mk = max_drdt(h_mk,0.02,0.02)
    r_em, r_ek = max_drdt(h_em,0.02,0.02), max_drdt(h_ek,0.02,0.02)
    r_ekm = max_drdt(h_ekm,0.02,0.02)

    # ######### He ignition matrix
    # epoch = 'He ignition'
    # R_l, R_h = He_ignition(h0,0.02,0.02), He_ignition(h1,0.02,0.02) # scalars
    # r_mm, r_kk = He_ignition(h_mm,0.02,0.02), He_ignition(h_kk,0.02,0.02)
    # r_ee = He_ignition(h_ee,0.02,0.02)
    # r_mk = He_ignition(h_mk,0.02,0.02)
    # r_em, r_ek = He_ignition(h_em,0.02,0.02), He_ignition(h_ek,0.02,0.02)
    # r_ekm = He_ignition(h_ekm,0.02,0.02)

    # ######### convective envelope matrix (Early Hayashi track)
    # epoch = 'Early Hayashi'
    # p = 50
    # R_l, R_h = ConvectiveFraction(h0,p), ConvectiveFraction(h1,p) # scalars
    # r_mm, r_kk = ConvectiveFraction(h_mm,p), ConvectiveFraction(h_kk,p)
    # r_ee = ConvectiveFraction(h_ee,p)
    # r_mk = ConvectiveFraction(h_mk,p)
    # r_em, r_ek = ConvectiveFraction(h_em,p), ConvectiveFraction(h_ek,p)
    # r_ekm = ConvectiveFraction(h_ekm,p)

    # ######### Hayashi matrix -- constant L
    # epoch = 'Hayashi Track L=1e5Lsun'
    # R_l, R_h = Hayashi(h0,L_c), Hayashi(h1,L_c) # scalars
    # r_mm, r_kk, r_ee = Hayashi(h_mm,L_c), Hayashi(h_kk,L_c), Hayashi(h_ee,L_c)
    # r_mk, r_em, r_ek = Hayashi(h_mk,L_c), Hayashi(h_em,L_c), Hayashi(h_ek,L_c)
    # r_ekm = Hayashi(h_ekm,L_c)

    # ######### C-depletion -- end of evolution
    # epoch = 'C-depletion'
    # R_l, R_h = Hayashi_end(h0), Hayashi_end(h1) # scalars
    # r_mm, r_kk, r_ee = Hayashi_end(h_mm), Hayashi_end(h_kk), Hayashi_end(h_ee)
    # r_mk, r_em, r_ek = Hayashi_end(h_mk), Hayashi_end(h_em), Hayashi_end(h_ek)
    # r_ekm = Hayashi_end(h_ekm)

    # construct matrix
    T_r = [(entry(r_mm,R_h,R_l),entry(r_mk,R_h,R_l),entry(r_em,R_h,R_l)),
           ('-',entry(r_kk,R_h,R_l),entry(r_ek,R_h,R_l)),
           ('-','-',entry(r_ee,R_h,R_l))]

    dfObj = pd.DataFrame(T_r, columns=['μ', 'κ', 'ε']
                         , index=['μ', 'κ', 'ε'])

    print("Z=0.02; {}".format(str(epoch)), dfObj, sep='\n')
    print("Δ_{μ,κ,ε} =",entry(r_ekm, R_h, R_l))
    print("Δ_μ+Δ_κ+Δ_ε =",entry(r_mm, R_h, R_l) + entry(r_kk, R_h, R_l) + entry(r_ee, R_h, R_l))
    print("Δ_2Z =",100 * (R_h - R_l) / R_l)
