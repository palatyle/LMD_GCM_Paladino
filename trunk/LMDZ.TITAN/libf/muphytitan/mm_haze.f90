! Copyright 2013-2015,2017 Université de Reims Champagne-Ardenne 
! Contributor: J. Burgalat (GSMA, URCA)
! email of the author : jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute
! microphysics processes using a two-moments scheme.
! 
! This library is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.

!! file: mm_haze.f90
!! summary: Haze microphysics module.
!! author: J. Burgalat
!! date: 2013-2015,2017
MODULE MM_HAZE
  !! Haze microphysics module.
  !!
  !! This module contains all definitions of the microphysics processes related to aerosols:
  !!
  !! - [coagulation](page/haze.html#coagulation)
  !! - [sedimentation](page/haze.html#sedimentation)
  !! - [production](page/haze.html#production)
  !!
  !! @note
  !! The production function is specific to Titan, where aerosols are created above the detached 
  !! haze layer. No other source is taken into account.  This process is controled by two parameters, 
  !! the pressure level of production and the production rate. Then both M0 and M3 of the aerosols 
  !! distribution are updated in the production zone by addition of the production rate along a
  !! gaussian shape.
  !!
  !! @note
  !! The interface methods always uses the global variables defined in [[mm_globals(module)]] when 
  !! values (any kind, temperature, pressure, moments...) over the vertical grid are required.
  !!
  !! @warning
  !! The tendencies returned by the method are always defined over the vertical grid from __TOP__ 
  !! to __GROUND__.
  !!
  !! @todo 
  !! Modify tests on tendencies vectors to get sure that allocation is done:
  !! Currently, we assume the compiler handles automatic allocation of arrays.
  USE MM_MPREC
  USE MM_GLOBALS
  USE MM_INTERFACES
  USE MM_METHODS
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_haze_microphysics, mm_haze_coagulation, mm_haze_sedimentation, &
            mm_haze_production

  CONTAINS

  !============================================================================
  ! HAZE MICROPHYSICS INTERFACE SUBROUTINE
  !============================================================================

  SUBROUTINE mm_haze_microphysics(dm0a_s,dm3a_s,dm0a_f,dm3a_f)
    !! Get the evolution of moments tracers through haze microphysics processes.
    !!
    !! The subroutine is a wrapper to the haze microphysics methods. It computes the tendencies 
    !! of moments tracers for coagulation, sedimentation and production processes for the 
    !! atmospheric column.
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0a_s
      !! Tendency of the 0th order moment of the spherical mode distribution (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3a_s
      !! Tendency of the 3rd order moment of the spherical mode distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0a_f
      !! Tendency of the 0th order moment of the fractal mode distribution (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3a_f
      !! Tendency of the 3rd order moment of the fractal mode distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: zdm0as
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: zdm3as
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: zdm0af
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: zdm3af

    dm0a_s = 0._mm_wp ; dm3a_s = 0._mm_wp ; dm0a_f = 0._mm_wp ; dm3a_f = 0._mm_wp

    ALLOCATE(zdm0as(mm_nla),zdm3as(mm_nla),zdm0af(mm_nla),zdm3af(mm_nla))
    zdm0as(1:mm_nla) = 0._mm_wp
    zdm3as(1:mm_nla) = 0._mm_wp
    zdm0af(1:mm_nla) = 0._mm_wp
    zdm3af(1:mm_nla) = 0._mm_wp

    IF (mm_w_haze_coag) THEN 
      ! Calls coagulation
      call mm_haze_coagulation(dm0a_s,dm3a_s,dm0a_f,dm3a_f)
    ENDIF

    IF (mm_w_haze_sed) THEN
      ! Calls sedimentation
      call mm_haze_sedimentation(zdm0as,zdm3as,zdm0af,zdm3af)

      ! Computes precipitations
      mm_aer_prec = SUM(zdm3as*mm_dzlev) + SUM(zdm3af*mm_dzlev)

      ! Updates tendencies
      dm0a_s=dm0a_s+zdm0as ; dm3a_s=dm3a_s+zdm3as 
      dm0a_f=dm0a_f+zdm0af ; dm3a_f=dm3a_f+zdm3af
    ENDIF

    IF (mm_w_haze_prod) THEN
      call mm_haze_production(zdm0as,zdm3as)
      ! We only produce spherical aerosols
      dm0a_s=dm0a_s+zdm0as ; dm3a_s=dm3a_s+zdm3as 
    ENDIF

    RETURN
  END SUBROUTINE mm_haze_microphysics


  !============================================================================
  ! COAGULATION PROCESS RELATED METHODS
  !============================================================================
 
  SUBROUTINE mm_haze_coagulation(dM0s,dM3s,dM0f,dM3f)
    !! Get the evolution of the aerosols moments vertical column due to coagulation process.
    !! 
    !! This is main method of the coagulation process:
    !!
    !! 1. Computes gamma pre-factor for each parts of the coagulation equation(s)
    !! 2. Applies the electic correction on the gamma pre-factor
    !! 3. Computes the specific flow regime "kernels"
    !! 4. Computes the harmonic mean of the kernels
    !! 5. Finally computes the tendencies of the moments.
    !!
    !! All arguments are assumed vectors of __N__ elements where __N__ is the total number of 
    !! vertical __layers__.
    !!
    !! @note
    !! The method uses directly the global variables related to the vertical atmospheric structure 
    !! stored in [[mm_globals(module)]]. Consequently they must be updated before calling the subroutine. 
    !!
    !! @bug
    !! If the transfert probabilities are set to 1 for the two flow regimes (pco and pfm), 
    !! a floating point exception occured (i.e. a NaN) as we perform a division by zero
    !!
    !! @todo
    !! Get rid of the fu\*\*\*\* STOP statement...
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dM0s
      !! Tendency of the 0th order moment of the spherical size-distribution over a time step (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dM3s
      !! Tendency of the 3rd order moment of the spherical size-distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dM0f
      !! Tendency of the 0th order moment of the fractal size-distribution over a time step (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dM3f
      !! Tendency of the 3rd order moment of the fractal size-distribution over a time step (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: c_kco,c_kfm,c_slf,tmp, &
                                                   kco,kfm,pco,pfm,mq
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: a_ss,a_sf,b_ss,b_ff,c_ss,c_sf
    INTEGER                                     :: i

    IF (mm_coag_choice < 0 .OR. mm_coag_choice > mm_coag_ss+mm_coag_sf+mm_coag_ff) &
      STOP "invalid choice for coagulation mode interaction activation"

    ! Alloctes local arrays
    ALLOCATE(kco(mm_nla),kfm(mm_nla),c_slf(mm_nla), &
             c_kco(mm_nla),c_kfm(mm_nla),mq(mm_nla), &
             pco(mm_nla),pfm(mm_nla))
    ALLOCATE(a_ss(mm_nla),a_sf(mm_nla), &
             b_ss(mm_nla),b_ff(mm_nla), &
             c_ss(mm_nla),c_sf(mm_nla))
             
    a_ss(:) = 0._mm_wp ; a_sf(:) = 0._mm_wp
    b_ss(:) = 0._mm_wp ; b_ff(:) = 0._mm_wp 
    c_ss(:) = 0._mm_wp ; c_sf(:) = 0._mm_wp

    ! gets kco, kfm pre-factors
    c_kco(:) = mm_get_kco(mm_temp) ; c_kfm(:) = mm_get_kfm(mm_temp)
    ! get slf (slip-flow factor)
    c_slf(:) = mm_akn * mm_lambda_g(mm_temp,mm_play) 

    DO i=1,mm_nla
      ! SS interactions
      IF (mm_rcs(i) > mm_rc_min .AND. IAND(mm_coag_choice,mm_coag_ss) /= 0) THEN
        ! compute probability for M0/CO and M0/FM (resp.)
        pco(i) = mm_ps2s(mm_rcs(i),0,0,mm_temp(i),mm_play(i))
        pfm(i) = mm_ps2s(mm_rcs(i),0,1,mm_temp(i),mm_play(i))
        ! (A_SS_CO x A_SS_FM) / (A_SS_CO + A_SS_FM)
        kco(i) = g0ssco(mm_rcs(i),c_slf(i),c_kco(i)) 
        kfm(i) = g0ssfm(mm_rcs(i),c_kfm(i)) 
        IF (kco(i)*(pco(i)-2._mm_wp)+kfm(i)*(pfm(i)-2._mm_wp) /=0) THEN
          a_ss(i) = (kco(i)*(pco(i)-2._mm_wp)*kfm(i)*(pfm(i)-2._mm_wp))/(kco(i)*(pco(i)-2._mm_wp)+kfm(i)*(pfm(i)-2._mm_wp))
        ENDIF
        ! (B_SS_CO x B_SS_FM) / (B_SS_CO + B_SS_FM)
        kco(i) = kco(i) * (1._mm_wp-pco(i)) ; kfm(i) = kfm(i) * (1._mm_wp-pfm(i))
        IF (kco(i) + kfm(i) /= 0._mm_wp) THEN
          b_ss(i) = kco(i)*kfm(i)/(kco(i)+kfm(i))
        ENDIF
        ! compute and apply eletric charge correction for M0/SS interactions
        mq(i) = mm_qmean(mm_rcs(i),mm_rcs(i),0,'SS',mm_temp(i),mm_play(i))

        a_ss(i) = a_ss(i) * mq(i)
        b_ss(i) = b_ss(i) * mq(i)
        kco(i) = 0._mm_wp ; kfm(i) = 0._mm_wp ; mq(i) = 1._mm_wp
        ! compute probability for M3/CO and M3/FM (resp.)
        pco(i) = mm_ps2s(mm_rcs(i),3,0,mm_temp(i),mm_play(i))
        pfm(i) = mm_ps2s(mm_rcs(i),3,1,mm_temp(i),mm_play(i))
        ! (C_SS_CO x C_SS_FM) / (C_SS_CO + C_SS_FM)
        kco(i) = g3ssco(mm_rcs(i),c_slf(i),c_kco(i))*(pco(i)-1._mm_wp)
        kfm(i) = g3ssfm(mm_rcs(i),c_kfm(i))*(pfm(i)-1._mm_wp)
        IF (kco(i) + kfm(i) /= 0._mm_wp) THEN
          c_ss(i) = (kco(i)*kfm(i))/(kco(i)+kfm(i))
        ENDIF
        IF (b_ss(i) <= 0._mm_wp) c_ss(i) = 0._mm_wp
        ! compute and apply eletric charge correction for M3/SS interactions
        mq(i) = mm_qmean(mm_rcs(i),mm_rcs(i),3,'SS',mm_temp(i),mm_play(i))
        c_ss(i) = c_ss(i) * mq(i)
      ENDIF
      kco(i) = 0._mm_wp ; kfm(i) = 0._mm_wp ; mq(i) = 1._mm_wp

      ! SF interactions
      IF (mm_rcs(i) > mm_rc_min .AND. mm_rcf(i) > mm_rc_min .AND. IAND(mm_coag_choice,mm_coag_sf) /= 0) THEN
        ! (A_SF_CO x A_SF_FM) / (A_SF_CO + A_SF_FM)
        kco(i) = g0sfco(mm_rcs(i),mm_rcf(i),c_slf(i),c_kco(i))
        kfm(i) = g0sffm(mm_rcs(i),mm_rcf(i),c_kfm(i))
        IF(kco(i)+kfm(i) /= 0._mm_wp) THEN
          a_sf(i) = (kco(i)*kfm(i))/(kco(i)+kfm(i))
        ENDIF
        ! compute and apply eletric charge correction for M0/SF interactions
        mq(i) = mm_qmean(mm_rcs(i),mm_rcf(i),0,'SF',mm_temp(i),mm_play(i))
        a_sf(i) = a_sf(i) * mq(i)
        ! (C_SF_CO x C_SF_FM) / (C_SF_CO + C_SF_FM)
        kco(i) = g3sfco(mm_rcs(i),mm_rcf(i),c_slf(i),c_kco(i))
        kfm(i) = g3sffm(mm_rcs(i),mm_rcf(i),c_kfm(i))
        IF (kco(i)+kfm(i) /= 0._mm_wp) THEN
          c_sf(i) = (kco(i)*kfm(i))/(kco(i)+kfm(i))
        ENDIF
        ! compute and apply eletric charge correction for M3/SF interactions
        mq(i) = mm_qmean(mm_rcs(i),mm_rcf(i),3,'SF',mm_temp(i),mm_play(i))
        c_sf(i) = c_sf(i) * mq(i)
      ENDIF
      kco(i) = 0._mm_wp ; kfm(i) = 0._mm_wp ; mq(i) = 1._mm_wp
      ! FF interactions
      IF(mm_rcf(i) > mm_rc_min .AND. IAND(mm_coag_choice,mm_coag_sf) /= 0) THEN
        ! (B_FF_CO x B_FF_FM) / (B_FF_CO + B_FF_FM)
        kco(i) = g0ffco(mm_rcf(i),c_slf(i),c_kco(i))
        kfm(i) = g0fffm(mm_rcf(i),c_kfm(i))
        b_ff(i) = (kco(i)*kfm(i))/(kco(i)+kfm(i))
        ! compute and apply eletric charge correction for M0/FF interactions
        mq(i) = mm_qmean(mm_rcf(i),mm_rcf(i),0,'FF',mm_temp(i),mm_play(i))
        b_ff(i) = b_ff(i) * mq(i)
      ENDIF
    ENDDO
   
    DEALLOCATE(kco,kfm,c_kco,c_kfm,pco,pfm,c_slf)

    ! Now we will use the kharm two by two to compute :
    ! dm_0_S/mm_dt = kharm(1) * m_0_S^2 - kharm(2) * m_0_S * m_0_F
    ! dm_0_F/mm_dt = kharm(3) * m_0_S^2 - kharm(4) * m_0_F^2
    ! dm_3_S/mm_dt = kharm(5) * m_3_S^2 - kharm(6) * m_3_S * m_3_F
    ! ... and finally :
    ! dm_3_F/mm_dt = - dm_3_S/mm_dt
    !
    ! We use a (semi) implicit scheme : when X appears as square we set one X
    ! at t+1, the other a t
    ALLOCATE(tmp(mm_nla))
    ! --- dm0s
    tmp(:) = mm_dt*(a_ss*mm_m0aer_s - a_sf*mm_m0aer_f)
    dm0s(:) =  mm_m0aer_s * (tmp/(1._mm_wp - tmp))
    ! --- dm0f
    tmp(:) = b_ff*mm_dt*mm_m0aer_f
    dm0f(:) = (b_ss*mm_dt*mm_m0aer_s**2 - tmp*mm_m0aer_f)/(1._mm_wp + tmp)
    ! --- dm3s
    tmp(:) = mm_dt*(c_ss*mm_m3aer_s - c_sf*mm_m3aer_f)
    dm3s(:) =  mm_m3aer_s * (tmp/(1._mm_wp - tmp))
    ! --- dmm3f
    dm3f(:) = -dm3s

    ! Deallocates memory explicitly ... another obsolete statement :)
    DEALLOCATE(a_ss,a_sf,b_ss,b_ff,c_ss,c_sf,tmp)

    ! Time to do something else !
    RETURN
  END SUBROUTINE mm_haze_coagulation

  ELEMENTAL FUNCTION g0ssco(rcs,c_slf,c_kco) RESULT(res)
    !! Get &gamma; pre-factor for the 0th order moment with SS interactions in the continuous flow regime.
    !!
    !! @note 
    !! If __rcs__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_slf !! Slip-Flow correction pre-factor.
    REAL(kind=mm_wp), INTENT(in) :: c_kco !! Thermodynamic continuous flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    a1=mm_alpha_s(1._mm_wp) ; a2=mm_alpha_s(-1._mm_wp) ; a3=mm_alpha_s(-2._mm_wp)
    ! Computes gamma pre-factor
    res = (1._mm_wp + a1*a2 + c_slf/rcs *(a2+a1*a3))*c_kco
    RETURN
  END FUNCTION g0ssco

  ELEMENTAL FUNCTION g0sfco(rcs,rcf,c_slf,c_kco) RESULT(res)
    !! Get &gamma; pre-factor for the 0th order moment with SF interactions in the continuous flow regime.
    !!
    !! @note 
    !! If __rcs__ or __rcf__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: rcf   !! Characteristic radius of the fractal size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_slf !! Slip-Flow correction pre-factor.
    REAL(kind=mm_wp), INTENT(in) :: c_kco !! Thermodynamic continuous flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, a6, e, rcff
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp .OR. rcf <= 0._mm_wp) RETURN
    e=3._mm_wp/mm_df ; rcff = rcf**e*mm_rb2ra 
    ! computes mm_alpha coefficients
    a1=mm_alpha_s(1._mm_wp)   ; a2=mm_alpha_f(-e)   ; a3=mm_alpha_f(e) 
    a4=mm_alpha_s(-1._mm_wp) ; a5=mm_alpha_s(-2._mm_wp) ; a6=mm_alpha_f(-2._mm_wp*e)
    ! Computes gamma pre-factor
    res = c_kco*( 2._mm_wp + a1*a2*rcs/rcff + a4*a3*rcff/rcs + c_slf*( a4/rcs + &
                         a2/rcff + a5*a3*rcff/rcs**2 + a1*a6*rcs/rcff**2))
    RETURN
  END FUNCTION g0sfco

  ELEMENTAL FUNCTION g0ffco(rcf,c_slf,c_kco) RESULT(res)
    !! Get &gamma; pre-factor for the 0th order moment with FF interactions in the continuous flow regime.
    !!
    !! @note 
    !! If __rcf__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcf   !! Characteristic radius of the fractal size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_slf !! Slip-Flow correction pre-factor.
    REAL(kind=mm_wp), INTENT(in) :: c_kco !! Thermodynamic continuous flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, e, rcff
    res = 0._mm_wp ; IF (rcf <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    e = 3._mm_wp/mm_df ; rcff = rcf**e*mm_rb2ra
    a1=mm_alpha_f(e) ; a2=mm_alpha_f(-e) ; a3=mm_alpha_s(-2._mm_wp*e)
    ! Computes gamma pre-factor
    res = (1._mm_wp + a1*a2 + c_slf/rcff *(a2+a1*a3))*c_kco
    RETURN
  END FUNCTION g0ffco

  ELEMENTAL FUNCTION g3ssco(rcs, c_slf, c_kco) RESULT(res)
    !! Get &gamma; pre-factor for the 3rd order moment with SS interactions in the continuous flow regime.
    !!
    !! @note 
    !! If __rcs__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_slf !! Slip-Flow correction pre-factor.
    REAL(kind=mm_wp), INTENT(in) :: c_kco !! Thermodynamic continuous flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, a6
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    a1=mm_alpha_s(3._mm_wp) ; a2=mm_alpha_s(2._mm_wp)  ; a3=mm_alpha_s(1._mm_wp) 
    a4=mm_alpha_s(4._mm_wp) ; a5=mm_alpha_s(-1._mm_wp) ; a6=mm_alpha_s(-2._mm_wp)

    ! Computes gamma pre-factor
    res = (2._mm_wp*a1 + a2*a3 + a4*a5 + c_slf/rcs*(a3**2 + a4*a6 + a1*a5 + a2))* & 
          c_kco/(a1**2*rcs**3)
    RETURN
  END FUNCTION g3ssco

  ELEMENTAL FUNCTION g3sfco(rcs, rcf, c_slf, c_kco) RESULT(res)
    !! Get &gamma; pre-factor for the 3rd order moment with SF interactions in the continuous flow regime.
    !!
    !! @note 
    !! If __rcs__ or __rcf__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: rcf   !! Characteristic radius of the fractal size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_slf !! Slip-Flow correction pre-factor.
    REAL(kind=mm_wp), INTENT(in) :: c_kco !! Thermodynamic continuous flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, a6, a7, a8, e, rcff
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp .OR. rcf <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    e=3._mm_wp/mm_df ; rcff = rcf**e*mm_rb2ra 
    a1=mm_alpha_s(3._mm_wp)    ; a2=mm_alpha_s(4._mm_wp) ; a3=mm_alpha_f(-e) 
    a4=mm_alpha_s(2._mm_wp)    ; a5=mm_alpha_f(e)    ; a6=mm_alpha_s(1._mm_wp)
    a7=mm_alpha_f(-2._mm_wp*e) ; a8=mm_alpha_f(3._mm_wp)
    ! Computes gamma pre-factor
    res = (2._mm_wp*a1*rcs**3 + a2*rcs**4*a3/rcff + a4*rcs**2*a5*rcff + &
                c_slf *( a4*rcs**2 + a1*rcs**3*a3/rcff + a6*rcs*a5*rcff + &
                a2*rcs**4*a7/rcff**2))* c_kco/(a1*a8*(rcs*rcf)**3)
    RETURN
  END FUNCTION g3sfco

  ELEMENTAL FUNCTION g0ssfm(rcs, c_kfm) RESULT(res)
    !! Get &gamma; pre-factor for the 0th order moment with SS interactions in the Free Molecular flow regime.
    !!
    !! @note 
    !! If __rcs__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_kfm !! Thermodynamic free molecular flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    a1=mm_alpha_s(0.5_mm_wp) ; a2=mm_alpha_s(1._mm_wp) ; a3=mm_alpha_s(-0.5_mm_wp)
    a4=mm_alpha_s(2._mm_wp)   ; a5=mm_alpha_s(-1.5_mm_wp)
    ! Computes gamma pre-factor
    res = (a1 + 2._mm_wp*a2*a3 + a4*a5)*rcs**0.5_mm_wp*mm_get_btk(1,0)*c_kfm
    RETURN
  END FUNCTION g0ssfm

  ELEMENTAL FUNCTION g0sffm(rcs, rcf, c_kfm) RESULT(res)
    !> Get &gamma; pre-factor for the 0th order moment with SF interactions in the Free Molecular flow regime. 
    !!
    !! @note 
    !! If __rcs__ or __rcf__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: rcf   !! Characteristic radius of the fractal size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_kfm !! Thermodynamic free molecular flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
    REAL(kind=mm_wp) :: e1, e2, e3, e4
    REAL(kind=mm_wp) :: rcff1, rcff2, rcff3, rcff4
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp .OR. rcf <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    e1 = 3._mm_wp/mm_df 
    e2 = 6._mm_wp/mm_df 
    e3 = (6._mm_wp-3._mm_wp*mm_df)/(2._mm_wp*mm_df)
    e4 = (12._mm_wp-3._mm_wp*mm_df)/(2._mm_wp*mm_df)

    rcff1 = mm_rb2ra * rcf**e1 
    rcff2 = rcff1**2
    rcff3 = mm_rb2ra * rcf**e3 
    rcff4 = mm_rb2ra**2 * rcf**e4

    a1=mm_alpha_s(0.5_mm_wp)  ; a2=mm_alpha_s(-0.5_mm_wp) ; a3=mm_alpha_f(e1) 
    a4=mm_alpha_s(-1.5_mm_wp) ; a5=mm_alpha_f(e2)      ; a6=mm_alpha_s(2._mm_wp)
    a7=mm_alpha_f(-1.5_mm_wp) ; a8=mm_alpha_s(1._mm_wp)   ; a9=mm_alpha_f(e3) 
    a10=mm_alpha_f(e4)

    ! Computes gamma pre-factor
    res = (a1*rcs**0.5_mm_wp + 2._mm_wp*rcff1*a2*a3/rcs**0.5_mm_wp + a4*a5*rcff2/rcs**1.5_mm_wp + &
           a6*a7*rcs**2/rcf**1.5_mm_wp + 2._mm_wp*a8*a9*rcs*rcff3 + a10*rcff4 &
          )*mm_get_btk(4,0)*c_kfm
    RETURN
  END FUNCTION g0sffm

  ELEMENTAL FUNCTION g0fffm(rcf, c_kfm) RESULT(res)
    !! Get &gamma; pre-factor for the 0th order moment with FF interactions in the Free Molecular flow regime. 
    !!
    !! @note 
    !! If __rcf__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcf   !! Characteristic radius of the fractal size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_kfm !! Thermodynamic free molecular flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor. 
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, e1, e2, e3, rcff
    res = 0._mm_wp ; IF (rcf <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    e1=3._mm_wp/mm_df ; e2=(6_mm_wp-3._mm_wp*mm_df)/(2._mm_wp*mm_df) 
    e3=(12._mm_wp-3._mm_wp*mm_df)/(2._mm_wp*mm_df)
    rcff=mm_rb2ra**2*rcf**e3
    a1=mm_alpha_f(e3)      ; a2=mm_alpha_f(e1) ; a3=mm_alpha_f(e2)
    a4=mm_alpha_f(-1.5_mm_wp) ; a5=mm_alpha_f(2._mm_wp*e1)
    ! Computes gamma pre-factor
    res = (a1 + 2._mm_wp*a2*a3 + a4*a5)*rcff*mm_get_btk(3,0)*c_kfm
    RETURN
  END FUNCTION g0fffm

  ELEMENTAL FUNCTION g3ssfm(rcs, c_kfm) RESULT(res)
    !! Get &gamma; pre-factor for the 3rd order moment with SS interactions in the Free Molecular flow regime.
    !!
    !! @note
    !! If __rcs__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_kfm !! Thermodynamic free molecular flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    a1=mm_alpha_s(3.5_mm_wp)  ; a2=mm_alpha_s(1._mm_wp)  ; a3=mm_alpha_s(2.5_mm_wp) 
    a4=mm_alpha_s(2._mm_wp)   ; a5=mm_alpha_s(1.5_mm_wp) ; a6=mm_alpha_s(5._mm_wp)
    a7=mm_alpha_s(-1.5_mm_wp) ; a8=mm_alpha_s(4._mm_wp)  ; a9=mm_alpha_s(-0.5_mm_wp) 
    a10=mm_alpha_s(3._mm_wp) ; a11=mm_alpha_s(0.5_mm_wp)
    ! Computes gamma pre-factor
    res = (a1 + 2._mm_wp*a2*a3 + a4*a5 + a6*a7 + 2._mm_wp*a8*a9 + a10*a11) &
          *mm_get_btk(1,3)*c_kfm/(a10**2*rcs**2.5_mm_wp)
    RETURN
  END FUNCTION g3ssfm

  ELEMENTAL FUNCTION g3sffm(rcs, rcf, c_kfm) RESULT(res)
    !! Get &gamma; pre-factor for the 3rd order moment with SF interactions in the Free Molecular flow regime.
    !!
    !! @note 
    !! If __rcs__ or __rcf__ is 0, the function returns 0.
    REAL(kind=mm_wp), INTENT(in) :: rcs   !! Characteristic radius of the spherical size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: rcf   !! Characteristic radius of the fractal size-distribution.
    REAL(kind=mm_wp), INTENT(in) :: c_kfm !! Thermodynamic free molecular flow regime pre-factor.
    REAL(kind=mm_wp) :: res               !! &gamma; coagulation kernel pre-factor.
    REAL(kind=mm_wp) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12
    REAL(kind=mm_wp) :: e1, e2, e3, rcff1, rcff2, rcff3
    res = 0._mm_wp ; IF (rcs <= 0._mm_wp .OR. rcf <= 0._mm_wp) RETURN
    ! computes mm_alpha coefficients
    e1=3._mm_wp/mm_df 
    e2=(6._mm_wp-3._mm_wp*mm_df)/(2._mm_wp*mm_df) 
    e3=(12._mm_wp-3._mm_wp*mm_df)/(2._mm_wp*mm_df)
    rcff1=mm_rb2ra*rcf**e1 ; rcff2=mm_rb2ra*rcf**e2 ; rcff3=mm_rb2ra**2*rcf**e3
    a1=mm_alpha_s(3.5_mm_wp)  ; a2=mm_alpha_s(2.5_mm_wp)   ; a3=mm_alpha_f(e1) 
    a4=mm_alpha_s(1.5_mm_wp)  ; a5=mm_alpha_f(2._mm_wp*e1) ; a6=mm_alpha_s(5._mm_wp) 
    a7=mm_alpha_f(-1.5_mm_wp) ; a8=mm_alpha_s(4._mm_wp)    ; a9=mm_alpha_f(e2) 
    a10=mm_alpha_s(3._mm_wp)  ; a11=mm_alpha_f(e3)         ; a12=mm_alpha_f(3._mm_wp)
    ! Computes gamma pre-factor
    res = (a1*rcs**3.5_mm_wp + 2._mm_wp*a2*a3*rcs**2.5_mm_wp*rcff1 + a4*a5*rcs**1.5_mm_wp*rcff1**2 + &
          a6*a7*rcs**5/rcf**1.5_mm_wp + 2._mm_wp*a8*a9*rcs**4*rcff2 + & 
          a10*a11*rcs**3*rcff3)*mm_get_btk(4,3)*c_kfm/(a10*a12*(rcs*rcf)**3)
    RETURN
  END FUNCTION g3sffm

  !============================================================================
  ! SEDIMENTATION PROCESS RELATED METHODS
  !============================================================================
  
  SUBROUTINE mm_haze_sedimentation(dm0s,dm3s,dm0f,dm3f)
    !! Interface to sedimentation algorithm.
    !!
    !! The subroutine computes the evolution of each moment of the aerosols tracers
    !! through sedimentation process and returns their tendencies for a timestep. 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0s
      !! Tendency of the 0th order moment of the spherical mode (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3s
      !! Tendency of the 3rd order moment of the spherical mode (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0f
      !! Tendency of the 0th order moment of the fractal mode (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3f
      !! Tendency of the 3rd order moment of the fractal mode (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: ft,fdcor,wth
    REAL(kind=mm_wp)                            :: m,n,p 
    REAL(kind=mm_wp), PARAMETER                 :: fac = 4._mm_wp/3._mm_wp * mm_pi 

    ALLOCATE(ft(mm_nle),wth(mm_nle),fdcor(mm_nle))

    !mm_aer_s_flux(:) = 0._mm_wp ; mm_aer_f_flux(:) = 0._mm_wp
    IF (mm_wsed_m0) THEN
      ! Spherical particles
      ! M0
      call get_weff(mm_m0aer_s,0._mm_wp,3._mm_wp,mm_rcs,mm_dt,mm_alpha_s,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m0as_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m0aer_s,-ft,mm_dt,dm0s)
      ! M3
      mm_m3as_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m3aer_s,-ft,mm_dt,dm3s)
      ! get mass flux
      mm_aer_s_flux(:) = fac * mm_rhoaer * ft(2:) * mm_m3aer_s
      ! Fractal particles
      ! M0
      call get_weff(mm_m0aer_f,0._mm_wp,mm_df,mm_rcf,mm_dt,mm_alpha_f,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m0af_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m0aer_f,-ft,mm_dt,dm0f)
      ! M3
      mm_m3af_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m3aer_f,-ft,mm_dt,dm3f)
      ! get mass flux
      mm_aer_f_flux(:) = fac * mm_rhoaer * ft(2:) * mm_m3aer_f
    ELSEIF (mm_wsed_m3) THEN
      ! Spherical particles
      ! M3
      call get_weff(mm_m3aer_s,3._mm_wp,3._mm_wp,mm_rcs,mm_dt,mm_alpha_s,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m3as_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m3aer_s,-ft,mm_dt,dm3s)
      ! get mass flux
      mm_aer_s_flux(:) = fac * mm_rhoaer * ft(2:) * mm_m3aer_s
      ! M0
      mm_m0as_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m0aer_s,-ft,mm_dt,dm0s)
      ! Fractal particles
      ! M3
      call get_weff(mm_m3aer_f,3._mm_wp,mm_df,mm_rcf,mm_dt,mm_alpha_f,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m3af_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m3aer_f,-ft,mm_dt,dm3f)
      ! get mass flux
      mm_aer_f_flux(:) = fac * mm_rhoaer * ft(2:) * mm_m3aer_f
      ! M0
      mm_m0af_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m0aer_f,-ft,mm_dt,dm0f)
    ELSE
      ! Spherical particles
      ! M0
      call get_weff(mm_m0aer_s,0._mm_wp,3._mm_wp,mm_rcs,mm_dt,mm_alpha_s,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m0as_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m0aer_s,-ft,mm_dt,dm0s)
      ! M3
      call get_weff(mm_m3aer_s,3._mm_wp,3._mm_wp,mm_rcs,mm_dt,mm_alpha_s,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m3as_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m3aer_s,-ft,mm_dt,dm3s)
      ! get mass flux
      mm_aer_s_flux(:) = fac * mm_rhoaer * ft(2:) * mm_m3aer_s
      ! Fractal particles
      ! M0
      call get_weff(mm_m0aer_f,0._mm_wp,mm_df,mm_rcf,mm_dt,mm_alpha_f,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m0af_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m0aer_f,-ft,mm_dt,dm0f)
      ! M3
      call get_weff(mm_m3aer_f,3._mm_wp,mm_df,mm_rcf,mm_dt,mm_alpha_f,wth,fdcor)
      ft(:) = wth(:) * fdcor(:)  ; mm_m3af_vsed(:) = ft(1:mm_nla)
      call let_me_fall_in_peace(mm_m3aer_f,-ft,mm_dt,dm3f)
      ! get mass flux
      mm_aer_f_flux(:) = fac * mm_rhoaer * ft(2:) * mm_m3aer_f
    ENDIF
    DEALLOCATE(ft,wth,fdcor)
    RETURN
  END SUBROUTINE mm_haze_sedimentation

  SUBROUTINE let_me_fall_in_peace(mk,ft,dt,dmk)
    !! Compute the tendency of the moment through sedimentation process.
    !!
    !! 
    !! The method computes the time evolution of the \(k^{th}\) order moment through sedimentation:
    !!
    !! $$ \dfrac{dM_{k}}{dt} = \dfrac{\Phi_{k}}{dz} $$
    !!
    !! The equation is resolved using a [Crank-Nicolson algorithm](http://en.wikipedia.org/wiki/Crank-Nicolson_method).
    !! 
    !! Sedimentation algorithm is quite messy. It appeals to the dark side of the Force and uses evil black magic spells 
    !! from ancient times. It is based on \cite{toon1988b,fiadeiro1977,turco1979} and is an update of the algorithm
    !! originally implemented in the LMDZ-Titan 2D GCM.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: mk  !! \(k^{th}\) order moment to sediment (in \(m^{k}\)).
    REAL(kind=mm_wp), INTENT(in)                :: dt  !! Time step (s).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: ft  !! Downward sedimentation flux  (effective velocity of the moment).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dmk !! Tendency of \(k^{th}\) order moment (in \(m^{k}.m^{-3}\)).
    INTEGER                                 :: i 
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: as,bs,cs,mko
    ALLOCATE(as(mm_nla), bs(mm_nla), cs(mm_nla), mko(mm_nla))
    mko(1:mm_nla) = 0._mm_wp
    cs(1:mm_nla) = ft(2:mm_nla+1) - mm_dzlay(1:mm_nla)/dt
    IF (ANY(cs > 0._mm_wp)) THEN
      ! implicit case
      as(1:mm_nla) = ft(1:mm_nla)
      bs(1:mm_nla) = -(ft(2:mm_nle)+mm_dzlay(1:mm_nla)/dt)
      cs(1:mm_nla) = -mm_dzlay(1:mm_nla)/dt*mk(1:mm_nla)
      ! (Tri)diagonal matrix inversion 
      mko(1) = cs(1)/bs(1)
      DO i=2,mm_nla ; mko(i) = (cs(i)-mko(i-1)*as(i))/bs(i) ; ENDDO
    ELSE
      ! explicit case
      as(1:mm_nla)=-mm_dzlay(1:mm_nla)/dt
      bs(1:mm_nla)=-ft(1:mm_nla)
      ! boundaries
      mko(1)=cs(1)*mk(1)/as(1)
      mko(mm_nla)=(bs(mm_nla)*mk(mm_nla-1)+cs(mm_nla)*mk(mm_nla))/as(mm_nla)
      ! interior
      mko(2:mm_nla-1)=(bs(2:mm_nla-1)*mk(1:mm_nla-2) + &
                       cs(2:mm_nla-1)*mk(2:mm_nla-1)   &
                      )/as(2:mm_nla-1)
    ENDIF
    dmk = mko - mk
    DEALLOCATE(as,bs,cs,mko)
    RETURN
  END SUBROUTINE let_me_fall_in_peace

  SUBROUTINE get_weff(mk,k,df,rc,dt,afun,wth,corf)
    !! Get the effective settling velocity for aerosols moments.
    !! 
    !! This method computes the effective settling velocity of the \(k^{th}\) order moment of aerosol 
    !! tracers. The basic settling velocity (\(v^{eff}_{M_{k}}\)) is computed using the following 
    !! equation:
    !! 
    !! $$ 
    !! \begin{eqnarray*}
    !! \Phi^{sed}_{M_{k}} &=& \int_{0}^{\infty} n(r) r^{k} \times w(r) dr 
    !!                     == M_{k}  \times v^{eff}_{M_{k}} \\
    !!     v^{eff}_{M_{k} &=& \dfrac{2 \rho g r_{c}^{\dfrac{3D_{f}-3}{D_{f}}}}
    !!     {r_{m}^{D_{f}-3}/D_{f}} \times \alpha(k)} \times \left( \alpha \right(
    !!     \frac{D_{f}(k+3)-3}{D_{f}}\left) + \dfrac{A_{kn}\lambda_{g}}{r_{c}^{
    !!     3/D_{f}}} \alpha \right( \frac{D_{f}(k+3)-6}{D_{f}}\left)\left)
    !! \end{eqnarray*}
    !! $$
    !!
    !! \(v^{eff}_{M_{k}\) is then corrected to reduce numerical diffusion of the sedimentation algorithm 
    !! as defined in \cite{toon1988b}.
    !!
    !! @warning
    !! Both __df__, __rc__ and __afun__ must be consistent with each other otherwise wrong values will 
    !! be computed.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(mm_nla)            :: mk
      !! Moment of order __k__ (\(m^{k}.m^{-3}\)) at each layer.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(mm_nla)            :: rc
      !! Characteristic radius associated to the moment at each layer.
    REAL(kind=mm_wp), INTENT(in)                               :: k
      !! The order of the moment.
    REAL(kind=mm_wp), INTENT(in)                               :: df
      !! Fractal dimension of the aersols.
    REAL(kind=mm_wp), INTENT(in)                               :: dt
      !! Time step (s).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(mm_nle)           :: wth
      !! Theoretical Settling velocity at each vertical __levels__ (\( wth \times corf = weff\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(mm_nle), OPTIONAL :: corf 
      !! _Fiadero_ correction factor applied to the theoretical settling velocity at each vertical __levels__.
    INTERFACE
      FUNCTION afun(order)
        !! Inter-moment relation function (see [[mm_interfaces(module):mm_alpha_s(interface)]]).
        IMPORT mm_wp
        REAL(kind=mm_wp), INTENT(in) :: order !! Order of the moment.
        REAL(kind=mm_wp) :: afun              !! Alpha value.
      END FUNCTION afun
    END INTERFACE 
    INTEGER                             :: i
    REAL(kind=mm_wp)                    :: af1,af2,ar1,ar2
    REAL(kind=mm_wp)                    :: csto,cslf,ratio,wdt,dzb
    REAL(kind=mm_wp)                    :: rb2ra
    REAL(kind=mm_wp), DIMENSION(mm_nle) :: zcorf
    ! ------------------
    
    wth(:) = 0._mm_wp ; zcorf(:) = 1._mm_wp

    ar1 = (3._mm_wp*df -3._mm_wp)/df    ; ar2 = (3._mm_wp*df -6._mm_wp)/df
    af1 = (df*(k+3._mm_wp)-3._mm_wp)/df ; af2 = (df*(k+3._mm_wp)-6._mm_wp)/df 
    rb2ra = mm_rm**((df-3._mm_wp)/df)
    DO i=2,mm_nla
      IF (rc(i-1) <= 0._mm_wp) CYCLE 
      dzb = (mm_dzlay(i)+mm_dzlay(i-1))/2._mm_wp
      csto = 2._mm_wp*mm_rhoaer*mm_effg(mm_zlev(i))/(9._mm_wp*mm_eta_g(mm_btemp(i)))
      cslf = mm_akn * mm_lambda_g(mm_btemp(i),mm_plev(i))
      wth(i) = - csto/(rb2ra*afun(k)) * (rc(i-1)**ar1 * afun(af1) + cslf/rb2ra * rc(i-1)**ar2 * afun(af2)) 
      ! now correct velocity to reduce numerical diffusion
      IF (.NOT.mm_no_fiadero_w) THEN
        IF (mk(i) <= 0._mm_wp) THEN
          ratio = mm_fiadero_max 
        ELSE
          ratio = MAX(MIN(mk(i-1)/mk(i),mm_fiadero_max),mm_fiadero_min) 
        ENDIF
        ! apply correction
        IF ((ratio <= 0.9_mm_wp .OR. ratio >= 1.1_mm_wp) .AND. wth(i) /= 0._mm_wp) THEN
          wdt = wth(i)*dt
          ! bugfix: max exponential arg to 30)
          zcorf(i) = dzb/wdt * (exp(MIN(30._mm_wp,-wdt*log(ratio)/dzb))-1._mm_wp) / (1._mm_wp-ratio)
          !zcorf(i) = dzb/wdt * (exp(-wdt*log(ratio)/dzb)-1._mm_wp) / (1._mm_wp-ratio)
        ENDIF
      ENDIF
    ENDDO
    ! last value (ground) set to first layer value: arbitrary :)
    wth(i) = wth(i-1)
    IF (PRESENT(corf)) corf(:) = zcorf(:)
    RETURN
  END SUBROUTINE get_weff

  !============================================================================
  ! PRODUCTION PROCESS RELATED METHOD
  !============================================================================

  SUBROUTINE mm_haze_production(dm0s,dm3s)
    !! Compute the production of aerosols moments.
    !! 
    !! The method computes the tendencies of M0 and M3 for the spherical mode through production process. 
    !! Production values are distributed along a normal law in altitude, centered  at 
    !! [[mm_globals(module):mm_p_prod(variable)]] pressure level with a fixed sigma of 20km.
    !!
    !! First M3 tendency is computed and M0 is retrieved using the inter-moments relation a spherical 
    !! characteristic radius set to [[mm_globals(module):mm_rc_prod(variable)]].
    !!
    !! If [[mm_globals(module):mm_var_prod(variable)]] is set to .true., the method computes time-dependent
    !! tendencies using a sine wave of anuglar frequency [[mm_globals(module):mm_w_prod(variable)]]
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0s !! Tendency of M0 (\(m^{-3}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3s !! Tendency of M3 (\(m^{3}.m^{-3}\)).
    INTEGER                     :: i
    REAL(kind=mm_wp)            :: zprod,cprod,timefact 
    REAL(kind=mm_wp), PARAMETER :: sigz  = 20e3_mm_wp, &
                                   fnorm = 1._mm_wp/(dsqrt(2._mm_wp*mm_pi)*sigz), &
                                   znorm = dsqrt(2._mm_wp)*sigz
    REAL(kind=mm_wp), SAVE      :: ctime = 0._mm_wp
    !$OMP THREADPRIVATE (ctime)
    zprod = -1._mm_wp
    ! locate production altitude
    DO i=1, mm_nla-1
      IF (mm_plev(i) < mm_p_prod.AND.mm_plev(i+1) >= mm_p_prod) THEN
        zprod = mm_zlay(i) ; EXIT
      ENDIF
    ENDDO
    IF (zprod < 0._mm_wp) THEN
      WRITE(*,'(a)') "cannot find aerosols production altitude"
      call EXIT(11)
    ENDIF

    dm3s(:)= mm_tx_prod *0.75_mm_wp/mm_pi *mm_dt / mm_rhoaer / 2._mm_wp / mm_dzlev(1:mm_nla) * &
             (erf((mm_zlev(1:mm_nla)-zprod)/znorm) - &
             erf((mm_zlev(2:mm_nla+1)-zprod)/znorm)) 
    dm0s(:) = dm3s(:)/(mm_rc_prod**3*mm_alpha_s(3._mm_wp))

    IF (mm_var_prod) THEN
      timefact = mm_d_prod*sin(mm_w_prod*ctime)+1._mm_wp
      dm3s(:) = timefact*dm3s(:)
      dm0s(:) = timefact*dm0s(:)
      ctime = ctime + mm_dt
    ENDIF


  END SUBROUTINE mm_haze_production

END MODULE MM_HAZE
