!  ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      ! use eos_lib
      ! use chem_def
      use star_data_def

      implicit none

      ! these routines are called by the standard run_star check_model
      contains
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         ! use eos_lib, only: eosDT_get_Rho

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items


         ! edit the extras_controls routine to set the procedure pointers
         ! e.g.,
         s% other_eosDT_get => my_eosDT_get
         s% other_eosDT_get_Rho => my_eosDT_get_Rho
         s% other_eosPT_get => my_eosPT_get
         s% other_eosDT_get_T => my_eosDT_get_T
         s% other_kap_get => artificially_high_Z_kap_get
         s% other_net_get => my_net_get

         ! hack the wind routine metallicity
         s% other_wind => Dutch_wind_artificially_high_Z


end subroutine extras_controls


! wind subroutine
      subroutine Dutch_wind_artificially_high_Z(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
         use star_def
         integer, intent(in) :: id
         real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
         real(dp) :: Zwind, T_high, T_low, alfa, w1, w2
         real(dp) :: L1, M1, R1, T1
         real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula
         ! NOTE: surface is outermost cell. not necessarily at photosphere.
         ! NOTE: don't assume that vars are set at this point.
         ! so if you want values other than those given as args,
         ! you should use values from s% xh(:,:) and s% xa(:,:) only.
         ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
         real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         w = 0
         ierr = 0

         L1 = Lsurf
         M1 = Msurf
         R1 = Rsurf
         T1 = Tsurf

         Zwind = min(1.d-3, Z)
         print *, "Hardcoded Z_wind =", Zwind, Z

         T_high = 11000
         T_low = 10000
         if (s% Dutch_scaling_factor == 0) then
            w = 0
         else if (T1 <= T_low) then
            call eval_lowT_Dutch(w)
         else if (T1 >= T_high) then
            call eval_highT_Dutch(w)
         else ! transition
            call eval_lowT_Dutch(w1)
            call eval_highT_Dutch(w2)
            alfa = (T1 - T_low)/(T_high - T_low)
            w = (1-alfa)*w1 + alfa*w2
         end if
         w = s% Dutch_scaling_factor * w


       contains

         subroutine eval_lowT_Dutch(w)
           real(dp), intent(out) :: w
           include 'formats'
           if (s% Dutch_wind_lowT_scheme == 'de Jager') then
              call eval_de_Jager_wind(w)
           else if (s% Dutch_wind_lowT_scheme == 'van Loon') then
              call eval_van_Loon_wind(w)
           else if (s% Dutch_wind_lowT_scheme == 'Nieuwenhuijzen') then
              call eval_Nieuwenhuijzen_wind(w)
           else
              write(*,*) 'unknown value for Dutch_wind_lowT_scheme ' // &
                   trim(s% Dutch_wind_lowT_scheme)
              w = 0
           end if
         end subroutine eval_lowT_Dutch

         subroutine eval_de_Jager_wind(w)
            ! de Jager, C., Nieuwenhuijzen, H., & van der Hucht, K. A. 1988, A&AS, 72, 259.
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
            w = exp10(log10w)
         end subroutine eval_de_Jager_wind


         subroutine eval_van_Loon_wind(w)
            ! van Loon et al. 2005, A&A, 438, 273
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -5.65d0 + 1.05d0*log10(L1/(1d4*Lsun)) - 6.3d0*log10(T1/35d2)
            w = exp10(log10w)
         end subroutine eval_van_Loon_wind


         subroutine eval_Nieuwenhuijzen_wind(w)
            ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
            real(dp), intent(out) :: w
            real(dp) :: log10w
            include 'formats'
            log10w = -14.02d0 + &
                     1.24d0*log10(L1/Lsun) + &
                     0.16d0*log10(M1/Msun) + &
                     0.81d0*log10(R1/Rsun)
            w = exp10(log10w)
         end subroutine eval_Nieuwenhuijzen_wind


         subroutine eval_highT_Dutch(w)
            real(dp), intent(out) :: w
            include 'formats'
            if (s% surface_h1 < 0.4d0) then ! helium rich Wolf-Rayet star: Nugis & Lamers
               w = 1d-11 * pow(L1/Lsun,1.29d0) * pow(Y,1.7d0) * sqrt(Zwind)
            else
               call eval_Vink_wind(w)
            end if
         end subroutine eval_highT_Dutch



         subroutine eval_Vink_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc

            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
               Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Zwind/Zsolar)))
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = (T1 - (Teff_jump - dT)) / (2*dT)
               end if
            end if

            if (alfa > 0) then ! eval hot side wind (eqn 24)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Zwind/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Zwind/Zsolar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if

            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Zwind/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.688d0 &
                  + 2.210d0*log10(L1/Lsun/1d5) &
                  - 1.339d0*log10(M1/Msun/30) &
                  - 1.601d0*log10(vinf_div_vesc/2d0) &
                  + 1.07d0*log10(T1/2d4) &
                  + 0.85d0*log10(Zwind/Zsolar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if

            w = alfa*w1 + (1 - alfa)*w2
         end subroutine eval_Vink_wind

      end subroutine Dutch_wind_artificially_high_Z


! kap subroutine.
         subroutine artificially_high_Z_kap_get( &
               id, k, handle, zbar, X, Z, Zbase, XC, XN, XO, XNe, &
               log10_rho, log10_T, species, chem_id, net_iso, xa, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)
 ! #######################
            use const_def, only: dp
            ! the line below was added so that we can call the kap_get function
            use kap_lib, only: kap_get
! ########################
            ! INPUT
            integer, intent(in) :: id ! star id if available; 0 otherwise
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell
            integer, intent(in) :: handle ! from alloc_kap_handle
            real(dp), intent(in) :: zbar ! average ion charge
            real(dp), intent(in) :: X, Z, Zbase, XC, XN, XO, XNe ! composition
            real(dp), intent(in) :: log10_rho ! density
            real(dp), intent(in) :: log10_T ! temperature
            real(dp), intent(in) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
               ! free_e := total combined number per nucleon of free electrons and positrons

            ! define new variables
            real(dp) :: Z_mod, Zbase_mod, XC_mod, XN_mod, XO_mod, XNe_mod

            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
               ! index from 1 to species
               ! value is between 1 and num_chem_isos
            integer, pointer :: net_iso(:) ! maps chem id to species number
               ! index from 1 to num_chem_isos (defined in chem_def)
               ! value is 0 if the iso is not in the current net
               ! else is value between 1 and number of species in current net
            real(dp), intent(in) :: xa(:) ! mass fractions

            ! OUTPUT
            real(dp), intent(out) :: frac_Type2
            real(dp), intent(out) :: kap ! opacity
            real(dp), intent(out) :: dln_kap_dlnRho ! partial derivative at constant T
            real(dp), intent(out) :: dln_kap_dlnT   ! partial derivative at constant Rho
            integer, intent(out) :: ierr ! 0 means AOK.

!  frac_Type2 = 0; kap = 0; dln_kap_dlnRho = 0; dln_kap_dlnT = 0

!  write(*,*) 'no implementation for other_kap_get'
!  ierr = -1
! #######################
            ! try hardcode these from the initial values of the high Z models
            ! Z, Zbase, XC, XN, XO, XNe
            
            ! Z_mod = max(0.00019999999999997797d0,Z)
            ! Zbase_mod = 0.00019999999999997797d0
            ! XC_mod = 3.4416073257056d-5
            ! XN_mod = 1.0081585205318563d-5
            ! XO_mod = 9.360449279773093d-5
            ! XNe_mod = 2.0994604783338605d-5

            ! Z_mod = max(0.040000000000000036d0,Z)
            ! Zbase_mod = 0.040000000000000036d0
            ! XC_mod = 0.006883214651411344d0
            ! XN_mod = 0.002016317041064625d0
            ! XO_mod = 0.01872089855954843d0
            ! XNe_mod = 0.004198920956668185d0
            
            ! Z_mod = max(1.d-4,Z)
            ! Zbase_mod = 1.d-4
            ! XC_mod = 1.7208036628528303d-5
            ! XN_mod = 5.040792602659346d-6
            ! XO_mod = 4.680224639886629d-5
            ! XNe_mod = 1.0497302391669425d-5    
            
            ! Z_mod = max(0.020000000000000018,Z)
            ! Zbase_mod =  0.020000000000000018
            ! XC_mod = 0.003441607325705952
            ! XN_mod = 0.001008158520531976
            ! XO_mod = 0.009360449279774226
            ! XNe_mod = 0.002099460478334096

            Z_mod = max(0.0010000000000000009,Z)
            Zbase_mod = 0.0010000000000000009
            XC_mod = 0.00017208036628530164
            XN_mod = 5.040792602659834d-5
            XO_mod = 0.00046802246398871075
            XNe_mod = 0.00010497302391670443
            
            ! this you want to evolve
            
            ! X

            !try not touch *lnfree_e* at first

            ! This is to do exactly what MESA would normally do
            call kap_get( &
                 id, zbar, X, Z_mod, Zbase_mod, XC_mod, XN_mod, XO_mod, XNe_mod, &
                 log10_Rho, log10_T, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
                 frac_Type2, kap, dln_kap_dlnRho, dln_kap_dlnT, ierr)

         end subroutine artificially_high_Z_kap_get
! kap routine end.

! epsilon nuc routine begins

      subroutine my_net_get(  &
            id, k, net_handle, just_dxdt, n, num_isos, num_reactions,  &
            x, temp, log10temp, rho, log10rho,  &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode, theta_e_for_graboske_et_al,  &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, ierr)

         ! use const_def, only: dp
         use net_lib, only: net_get
         use net_def, only: Net_Info
         ! logical :: exist
         ! use chem_def

         integer, intent(in) :: id ! id for star
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell
         integer, intent(in) :: net_handle
         logical, intent(in) :: just_dxdt
         type (Net_Info), pointer:: n
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in)  :: x(:) ! (num_isos)
         real(dp), intent(in)  :: temp, log10temp ! log10 of temp
         real(dp), intent(in)  :: rho, log10rho ! log10 of rho
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! electron degeneracy from eos.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         logical, intent(in) :: reuse_rate_raw, reuse_rate_screened ! if true. use given rate_screened

         real(dp), intent(out) :: eps_nuc ! ergs/g/s from burning after subtract reaction neutrinos
         real(dp), intent(out) :: d_eps_nuc_dT
         real(dp), intent(out) :: d_eps_nuc_dRho
         real(dp), intent(inout) :: d_eps_nuc_dx(:) ! (num_isos)
         real(dp), intent(inout) :: dxdt(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dRho(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dT(:) ! (num_isos)
         real(dp), intent(inout) :: d_dxdt_dx(:,:) ! (num_isos, num_isos)
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: eps_neu_total ! ergs/g/s neutrinos from weak reactions
         integer, intent(in) :: screening_mode
         real(dp), intent(in)  :: theta_e_for_graboske_et_al
         integer, intent(in) :: lwork ! size of work >= result from calling net_work_size
         real(dp), pointer :: work(:) ! (lwork)
         integer, intent(out) :: ierr ! 0 means okay
         real(dp) :: abar_new, zbar_new, z2bar_new
         real(dp) :: r_h
         ! real(dp), allocatable :: d_eps_nuc_dx_prev(:), dxdt_prev(:), d_dxdt_dRho_prev(:), &
         !         d_dxdt_dT_prev(:), d_dxdt_dx_prev(:,:), eps_nuc_categories_prev(:)
         real(dp) :: metals, z_high, r_nonmetal, r_metal

         type (star_info), pointer :: s

         ! define new arguments
         integer :: LENGTH

         real(dp) ,allocatable :: x_new(:), num_nucleon(:), num_charge(:), num_charge_sq(:) 
         real(dp), allocatable :: ion_abun(:), num_density(:)

         ierr = 0
         call star_ptr(id, s, ierr)        

         if (s% net_name == 'basic.net') then
                 LENGTH = 8
         
         else if (s% net_name == 'co_burn.net') then 
                 LENGTH = 9

         else if (s% net_name == 'approx21.net') then
                 LENGTH = 21

         endif 
         
         allocate (x_new(LENGTH))
         allocate (num_nucleon(LENGTH))
         allocate (num_charge(LENGTH))
         allocate (num_charge_sq(LENGTH))
         allocate (ion_abun(LENGTH))
         allocate (num_density(LENGTH))

         ! allocate (d_eps_nuc_dx_prev(SIZE(d_eps_nuc_dx)))
         ! allocate (dxdt_prev(SIZE(dxdt)))
         ! allocate (d_dxdt_dRho_prev(SIZE(d_dxdt_dRho))) 
         ! allocate (d_dxdt_dT_prev(SIZE(d_dxdt_dT)))
         ! allocate (d_dxdt_dx_prev(SIZE(d_dxdt_dx(:,:))))
         ! allocate (eps_nuc_categories_prev(SIZE(eps_nuc_categories_prev)))

          if (LENGTH == 8) then
                  metals = SUM(x(4:))
                  z_high = max(0.04d0,metals)
                  r_metal = z_high / metals
                  num_nucleon = (/1,3,4,12,14,16,20,24/) ! A(i)
                  num_charge = (/1,2,2,6,7,8,10,12/) ! Z(i)
                  num_charge_sq = (/1,4,4,36,49,64,100,144/) ! Z(i)^2
                  x_new = [(1 - z_high - x(2) - x(3)), x(2), x(3), x(4:) * r_metal]
          
          else if (LENGTH == 9) then
                  ! STOP "change net"
                  ! metals = SUM(x(4:))
                  ! z_high = max(0.04d0,metals)
                  ! r_metal = z_high / metals
                  ! num_nucleon = (/1,3,4,12,14,16,20,24,28/) ! A(i)
                  ! num_charge = (/1,2,2,6,7,8,10,12,14/) ! Z(i)
                  ! num_charge_sq = (/1,4,4,36,49,64,100,144,196/) ! Z(i)^2
                  ! x_new = [(1 - z_high - x(2) - x(3)), x(2), x(3), x(4:) * r_metal]
                  x_new = x
          else if (LENGTH == 21) then
                  ! metals = SUM(x(6:))
                  ! z_high = max(0.04d0,metals)
                  ! r_metal = z_high / metals
                  ! num_nucleon = (/1,1,1,3,4,12,14,16,20,24,28,32,36,40,44,48,56,52,54,56,56/) ! A(i)
                  ! num_charge = (/0,1,1,2,2,6,7,8,10,12,14,16,18,20,22,24,24,26,26,26,28/) ! Z(i)
                  ! num_charge_sq = (/0,1,1,4,4,36,49,64,100,144,194,256,324,400,484,576,576,676,676,676,784/) ! Z(i)^2
                  ! x_new = [x(1), (1 - z_high - x(1) - x(3) - x(4) - x(5)), x(3), x(4), x(5), x(6:) * r_metal]
                  x_new = x
          endif
          
          ! definitions from chem_lib.f90, get_composition_info subroutine.

          ion_abun = x_new / num_nucleon ! Y(i)
          num_density = ion_abun * avo * rho ! n(i)
          
          abar_new = SUM(num_density * num_nucleon) / SUM(num_density)
          zbar_new = SUM(num_density * num_charge) / SUM(num_density)
          z2bar_new = SUM(num_density * num_charge ** 2) / SUM(num_density)
          
          ! print *, x
          ! print *, x_new

          ! if (s% model_number /= 1) then
          !      ! print *, s% time_step
          !      d_eps_nuc_dx_prev = s% xtra_old(1)
          !      dxdt_prev = s% xtra_old(2)
          !      d_dxdt_dRho_prev = s% xtra_old(3)
          !      d_dxdt_dT_prev = s% xtra_old(4)
          !      d_dxdt_dx_prev = s% xtra_old(5)
          !      eps_nuc_categories_prev = s% xtra_old(6)


          !      ! STOP "check"
          ! else if (s% model_number == 1) then
          !      d_eps_nuc_dx_prev = d_eps_nuc_dx
          !      dxdt_prev = dxdt
          !      d_dxdt_dRho_prev = d_dxdt_dRho
          !      d_dxdt_dT_prev = d_dxdt_dT
          !      d_dxdt_dx_prev = d_dxdt_dx
          !      eps_nuc_categories_prev = eps_nuc_categories
          ! endif

         call net_get( &
              net_handle, just_dxdt, n, num_isos, num_reactions,  &
              x_new, temp, log10temp, rho, log10rho,  &
              abar_new, zbar_new, z2bar_new, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
              rate_factors, weak_rate_factor, &
              reaction_Qs, reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
              eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
              dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
              screening_mode, theta_e_for_graboske_et_al,  &
              eps_nuc_categories, eps_neu_total, &
              lwork, work, ierr)


         deallocate (x_new)
         deallocate (num_nucleon)
         deallocate (num_charge)
         deallocate (num_charge_sq)
         deallocate (ion_abun)
         deallocate (num_density)


         ! s% x_ctrl(1) = d_eps_nuc_dx
         ! s% x_ctrl(2) = dxdt
         ! s% x_ctrl(3) = d_dxdt_dRho
         ! s% x_ctrl(4) = d_dxdt_dT
         ! s% x_ctrl(5) = d_dxdt_dx
         ! s% x_ctrl(6) = eps_nuc_categories

         ! deallocate (d_eps_nuc_dx_prev)
         ! deallocate (dxdt_prev)
         ! deallocate (d_dxdt_dRho_prev)
         ! deallocate (d_dxdt_dT_prev)
         ! deallocate (d_dxdt_dx_prev)
         ! deallocate (eps_nuc_categories_prev)

         ! STOP "check"

   end subroutine my_net_get


! epsilon nuc routine ends

! eos routine begins.

      subroutine my_eosDT_get( &
              id, k, handle, Z, X, abar, zbar, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, &
              res, d_dlnRho_const_T, d_dlnT_const_Rho, &
              d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
         ! use const_def, only: dp

         ! INPUT
         use chem_def, only: num_chem_isos
         use eos_lib
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell
         integer, intent(in) :: handle ! eos handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction

         real(dp), intent(in) :: abar
            ! mean atomic number (nucleons per nucleus; grams per mole)
         real(dp), intent(in) :: zbar ! mean charge per nucleus

         real(dp), intent(in) :: xa(:), Rho, log10Rho, T, log10T

         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number

         real(dp) :: Z_mod

         ! OUTPUT

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results)

         integer, intent(out) :: ierr ! 0 means AOK.

         ! print *, s% gamma_law_hydro, ">0?"
         ! print *, s% use_eosDT_ideal_gas, "true?"
         ! print *, s% use_eosDT_HELMEOS, "true?"

         ! print *, z
         ! print *, res

         Z_mod = max(1.0d-3,Z)
         ! print *, Z_mod

       call eosDT_get( &
            handle, Z_mod, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Rho, log10Rho, T, log10T, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
        ! print *, res
        ! STOP

      end subroutine my_eosDT_get

         subroutine my_eosDT_get_Rho( &
                  id, k, handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  logT, which_other, other_value, &
                  logRho_tol, other_tol, max_iter, logRho_guess,  &
                  logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
                  logRho_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                  d_dabar_const_TRho, d_dzbar_const_TRho, &
                  eos_calls, ierr)

            ! finds log10 Rho given values for temperature and 'other', and initial guess for density.
            ! does up to max_iter attempts until logRho changes by less than tol.

            ! 'other' can be any of the basic result variables for the eos
            ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)

            use eos_lib, only: eosDT_get_Rho
            use const_def, only: dp

            integer, intent(in) :: id ! star id if available; 0 otherwise
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell
            integer, intent(in) :: handle

            real(dp), intent(in) :: Z ! the metals mass fraction
            real(dp), intent(in) :: X ! the hydrogen mass fraction

            real(dp), intent(in) :: abar
               ! mean atomic number (nucleons per nucleus; grams per mole)
            real(dp), intent(in) :: zbar ! mean charge per nucleus

            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
               ! index from 1 to species
               ! value is between 1 and num_chem_isos
            integer, pointer :: net_iso(:) ! maps chem id to species number
               ! index from 1 to num_chem_isos (defined in chem_def)
               ! value is 0 if the iso is not in the current net
               ! else is value between 1 and number of species in current net
            real(dp), intent(in) :: xa(:) ! mass fractions

            real(dp), intent(in) :: logT ! log10 of temperature

            integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
            real(dp), intent(in) :: other_value ! desired value for the other variable
            real(dp), intent(in) :: other_tol

            real(dp), intent(in) :: logRho_tol

            integer, intent(in) :: max_iter ! max number of Newton iterations

            real(dp), intent(in) :: logRho_guess ! log10 of density
            real(dp), intent(in) :: logRho_bnd1, logRho_bnd2 ! bounds for logRho
               ! if don't know bounds, just set to arg_not_provided (defined in const_def)
            real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
               ! if don't know these values, just set to arg_not_provided (defined in const_def)

            real(dp), intent(out) :: logRho_result ! log10 of density

            real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results)

            integer, intent(out) :: eos_calls
            integer, intent(out) :: ierr ! 0 means AOK.

            ! print *, res

            call eosDT_get_Rho( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               !eos_calls, ierr)
               d_dabar_const_TRho, d_dzbar_const_TRho, eos_calls, ierr)
            ! STOP
         end subroutine my_eosDT_get_Rho

         subroutine my_eosPT_get(&
                  id, k, handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa,&
                  Pgas, log10Pgas, T, log10T, &
                  Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
                  res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                  d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
            use const_def, only: dp
            use eos_lib, only: eosPT_get, eos_ptr
            ! INPUT

            integer, intent(in) :: id ! star id if available; 0 otherwise
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell
            integer, intent(in) :: handle

            real(dp), intent(in) :: Z ! the metals mass fraction
            real(dp), intent(in) :: X ! the hydrogen mass fraction

            real(dp), intent(in) :: abar
               ! mean atomic number (nucleons per nucleus; grams per mole)
            real(dp), intent(in) :: zbar ! mean charge per nucleus

            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
               ! index from 1 to species
               ! value is between 1 and num_chem_isos
            integer, pointer :: net_iso(:) ! maps chem id to species number
               ! index from 1 to num_chem_isos (defined in chem_def)
               ! value is 0 if the iso is not in the current net
               ! else is value between 1 and number of species in current net
            real(dp), intent(in) :: xa(:) ! mass fractions

            real(dp), intent(in) :: Pgas, log10Pgas ! the gas pressure
               ! provide both if you have them.  else pass one and set the other to arg_not_provided
               ! "arg_not_provided" is defined in mesa const_def

            real(dp), intent(in) :: T, log10T ! the temperature
               ! provide both if you have them.  else pass one and set the other to arg_not_provided

            type (EoS_General_Info), pointer :: rq

            ! OUTPUT

            real(dp), intent(out) :: Rho, log10Rho ! density
            real(dp), intent(out) :: dlnRho_dlnPgas_const_T
            real(dp), intent(out) :: dlnRho_dlnT_const_Pgas
            real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
            ! partial derivatives of the basic results wrt lnd and lnT
            real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
            ! d_dlnRho_const_T(i) = d(res(i))/dlnd|T
            real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
            ! d_dlnT_const_Rho(i) = d(res(i))/dlnT|Rho
            real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results)


            integer, intent(out) :: ierr ! 0 means AOK.
          call eos_ptr(handle,rq,ierr)
            ! print *, res

          call eosPT_get( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Pgas, log10Pgas, T, log10T, &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               !ierr)
               d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
           ! STOP
         end subroutine my_eosPT_get

         subroutine my_eosDT_get_T( &
                  id, k, handle, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  logRho, which_other, other_value, &
                  logT_tol, other_tol, max_iter, logT_guess, &
                  logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
                  logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                  d_dabar_const_TRho, d_dzbar_const_TRho, &
                  eos_calls, ierr)
            use const_def, only: dp
            use eos_lib, only: eosDT_get_T

            integer, intent(in) :: id ! star id if available; 0 otherwise
            integer, intent(in) :: k ! cell number or 0 if not for a particular cell
            integer, intent(in) :: handle

            real(dp), intent(in) :: Z ! the metals mass fraction
            real(dp), intent(in) :: X ! the hydrogen mass fraction

            real(dp), intent(in) :: abar
            real(dp), intent(in) :: zbar ! mean charge per nucleus

            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
            integer, pointer :: net_iso(:) ! maps chem id to species number
            real(dp), intent(in) :: xa(:) ! mass fractions

            real(dp), intent(in) :: logRho ! log10 of density
            integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
            real(dp), intent(in) :: other_value ! desired value for the other variable
            real(dp), intent(in) :: other_tol

            real(dp), intent(in) :: logT_tol
            integer, intent(in) :: max_iter ! max number of iterations

            real(dp), intent(in) :: logT_guess ! log10 of temperature
            real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds

            real(dp), intent(out) :: logT_result ! log10 of temperature
            real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results)

            integer, intent(out) :: eos_calls
            integer, intent(out) :: ierr ! 0 means AOK.

            call eosDT_get_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               !eos_calls, ierr)
               d_dabar_const_TRho, d_dzbar_const_TRho, eos_calls, ierr)

         end subroutine my_eosDT_get_T

! eos routine ends.


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         ! character (len=60) :: file_name = "/Users/cxin/Documents/Stars/data/eps_nuc_z0.02_ms.txt"
         ! real :: eps_nuc_data(924)
         ! integer :: unit

         ! inquire(file=file_name,exist=exist)

         ! i = 0
         ! if (exist) then
         ! OPEN(unit=11, file=file_name,status='old',action='read',iostat=ierr)
         ! do i=1, size(eps_nuc_data)
               ! i = i + 1
         !      READ(11,'es12.6',iostat=ierr) eps_nuc_data(i)
         !      print *, 'a', eps_nuc_data(i)
         !      s% xtra1_array(i) = eps_nuc_data(i)
         ! enddo
         ! close(11)
         ! end if

      end subroutine extras_startup


      integer function extras_start_step(id)

         use chem_def, only: ih1

         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         integer :: i

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
         ! print *, s% model_number

         !if (s% model_number >= 1005) then
         !       s% use_other_net_get = .true.
                ! print *, s% model_number
         !else
                ! print *, "nope"
         !endif

         !print *, s% model_number

         ! do i = 1, s% nz, 1
         !        if (s% xa(s% net_iso(ih1), i) < 1d-4) then
         !                s% max_years_for_timestep = 3.d2
         !        endif
         ! enddo

      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 4
      end function how_many_extra_history_columns

         ! I started editing here.
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)

           use math_lib, only: safe_log10
           use chem_def, only: ih1

           integer, intent(in) :: id, n
           character (len=maxlen_history_column_name) :: names(n)
           real(dp) :: vals(n)
           integer, intent(out) :: ierr
           type (star_info), pointer :: s

           integer :: i
           real(dp) :: mu_ave ! mass average mean molecular weight
           real(dp) :: kap_ave ! flux averaged opacity in radiative region
           real(dp) :: f_rad, f_tot, m_conv_tot, aux_var

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           mu_ave = 0
           do i = s% nz, 1, -1
              mu_ave = mu_ave + s% mu(i) * s% dm(i)
           end do

           kap_ave = 0
           f_tot = 0
           do i = s% nz, -1, 1
              if (s% mixing_type(i) == 0) then
                 f_rad = s% L(i) / (4 * pi * (s% r(i) ** 2))
                 kap_ave = kap_ave + f_rad * s% opacity(i)
                 f_tot = f_tot + f_rad
              endif
           end do

           ! r_conv_tot = 0

           aux_var = 0
           m_conv_tot = 0
           ! if (s% xa(s% net_iso(ih1), s% nz) > 0.6) then
           do i = 1, s% nz, 1
                if (s% xa(s% net_iso(ih1), i) > 0.6) then
                        if (s% mixing_type(i) == 1) then
                                m_conv_tot = m_conv_tot + s% dm(i)
                        aux_var = 0
                        else
                                aux_var = aux_var + s% dm(i)
                        endif
                else
                        exit
                endif

                if (aux_var >= 3 * msol) then
                        exit
                endif
           end do


           ! endif
           ! kap_ave = kap_ave / f_tot

           ! print *, mu_ave
           names(1) = "mu_int"
           ! vals(1) = mu_ave / SUM(s% dm,DIM=1)
           vals(1) = mu_ave / (s% star_mass * msol)

           names(2) = "mu_tot"
           vals(2) = mu_ave

           names(3) = "kap_ave"
           vals(3) = kap_ave / f_tot

           names(4) = "m_conv_tot"
           vals(4) = m_conv_tot

           ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve



      end module run_star_extras
