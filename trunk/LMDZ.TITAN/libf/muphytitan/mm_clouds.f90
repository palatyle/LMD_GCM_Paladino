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

!! file: mm_clouds.f90
!! summary: Clouds microphysics module
!! author: J. Burgalat
!! date: 2013-2015,2017

MODULE MM_CLOUDS
  !! Clouds microphysics module.
  !!
  !! The module contains all definitions of the microphysics processes related to clouds: 
  !!
  !! - [nucleation](page/clouds.html#nucleation)
  !! - [condensation](page/clouds.html#condensation)
  !! - [sedimentation](page/clouds.html#sedimentation)
  !!
  !! 
  !! The interface methods always use the global variables defined in [[mm_globals(module)]] when values 
  !! (any kind, temperature, pressure, moments...) over the vertical grid are required.
  !! Consequently, all these functions only deals with output argument which are most of the time the 
  !! tendencies of relevant variables on the atmospheric column.
  !!
  !! @note 
  !! Tendencies returned by public methods are always defined from  __TOP__ of the atmosphere to the 
  !! __GROUND__.
  USE MM_MPREC
  USE MM_GLOBALS
  USE MM_METHODS
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_cloud_microphysics, mm_cloud_sedimentation, mm_cloud_nucond 

  CONTAINS

  !============================================================================
  ! CLOUDS MICROPHYSICS INTERFACE SUBROUTINE
  !============================================================================

  SUBROUTINE mm_cloud_microphysics(dm0a,dm3a,dm0n,dm3n,dm3i,dgazs)
    !! Get the evolution of moments tracers through clouds microphysics processes.
    !!
    !! The subroutine is a wrapper to the clouds microphysics methods. It computes the tendencies of moments 
    !! tracers for nucleation, condensation and sedimentation processes for the atmospheric column.
    !!
    !! @note 
    !! Both __dm3i__ and __dgazs__ are 2D-array with the vertical layers in first dimension and the number 
    !! of ice components in the second. 
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm0a
      !! Tendency of the 0th order moment of the aerosols (fractal mode) (\(m^{-3}\)). 
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm3a
      !! Tendency of the 3rd order moment of the aerosols distribution (fractal mode) (\(m^{3}.m^{-3}\)) .
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm0n
      !! Tendency of the 0th order moment of the aerosols distribution (fractal mode) (\(m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm3n
      !! Tendency of the 3rd order moment of the ccn distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(out) :: dm3i
      !! Tendencies of the 3rd order moments of each ice components (\(m^{3}.m^{-3}\)). 
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(out) :: dgazs
      !! Tendencies of each condensible gaz species (\(mol.mol^{-1}\)). 
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE   :: zdm0n,zdm3n 
    REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE :: zdm3i
    INTEGER                                    :: i
    dm0a = 0._mm_wp ; dm3a = 0._mm_wp 
    dm0n = 0._mm_wp ; dm3n = 0._mm_wp 
    dm3i = 0._mm_wp ; dgazs = 0._mm_wp

    IF (mm_w_cloud_nucond) THEN
      ! Calls condensation/nucleation (and update saturation ratio diagnostic)
      call mm_cloud_nucond(dm0a,dm3a,dm0n,dm3n,dm3i,dgazs,mm_gazs_sat)
    ENDIF

    IF (mm_w_cloud_sed) THEN
      ! Calls sedimentation
      ALLOCATE(zdm0n(mm_nla),zdm3n(mm_nla),zdm3i(mm_nla,mm_nesp))
      call mm_cloud_sedimentation(zdm0n,zdm3n,zdm3i)

      ! computes precipitations / ice fluxes
      mm_ccn_prec = SUM(zdm3n*mm_dzlev)
      mm_ccn_flux(:) = get_mass_flux(mm_rhoaer,mm_m3ccn(:)) 

      DO i=1, mm_nesp 
        mm_ice_prec(i) = SUM(zdm3i(:,i)*mm_dzlev)
        mm_ice_fluxes(:,i) = get_mass_flux(mm_xESPS(i)%rho,mm_m3ice(:,i)) 
      ENDDO 
      ! updates tendencies
      dm0n = dm0n + zdm0n
      dm3n = dm3n + zdm3n
      dm3i = dm3i + zdm3i
    ENDIF

  END SUBROUTINE mm_cloud_microphysics

  !-----------------------------------------------------------------------------
  ! NUCLEATION/CONDENSATION PROCESS RELATED METHODS
  !-----------------------------------------------------------------------------

  SUBROUTINE mm_cloud_nucond(dm0a,dm3a,dm0n,dm3n,dm3i,dgazs,gazsat)
    !! Get moments tendencies through nucleation/condensation/evaporation.
    !!
    !! The method is a wrapper of [[mm_clouds(module):nc_esp(subroutine)]] which computes the 
    !! tendencies of tracers for all the condensible species given in the vector __xESPS__.
    !!
    !! @warning
    !! __xESPS__, __m3i__ and __gazes__ must share the same indexing. For example if __xESPS(IDX)__ 
    !! corresponds to \(CH_{4}\) properties then __m3i(IDX)__ must be the total volume of solid 
    !! \(CH_{4}\) (ice) and  __gazs(IDX)__ its vapor mole fraction.
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm0a
      !! Tendency of the 0th order moment of the aerosols (fractal mode) (\(m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm3a
      !! Tendency of the 3rd order moment of the aerosols distribution (fractal mode) (\(m^{3}.m^{-3}\)). 
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm0n
      !! Tendency of the 0th order moment of the aerosols distribution (fractal mode) (\(m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(out)   :: dm3n
      !! Tendency of the 3rd order moment of the ccn distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(out) :: dm3i
      !! Tendencies of the 3rd order moments of each ice components (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(out) :: dgazs
      !! Tendencies of each condensible gaz species (\(mol.mol^{-1}\)) .
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(out) :: gazsat
      !! Saturation ratio of each condensible specie.
    INTEGER                                   :: i,idx,ng
    TYPE(mm_esp)                                 :: xESP
    REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE :: zdm0a,zdm3a,zdm0n,zdm3n
    ALLOCATE(zdm0a(mm_nla,mm_nesp),zdm3a(mm_nla,mm_nesp), &
             zdm0n(mm_nla,mm_nesp),zdm3n(mm_nla,mm_nesp))
    zdm0a(:,:) = 0._mm_wp ; zdm3a(:,:) = 0._mm_wp
    zdm0n(:,:) = 0._mm_wp ; zdm3n(:,:) = 0._mm_wp
    DO i = 1, mm_nesp
      call nc_esp(mm_xESPS(i),mm_gazs(:,i),mm_m3ice(:,i),dgazs(:,i),dm3i(:,i), &
                  zdm0a(:,i),zdm3a(:,i),zdm0n(:,i),zdm3n(:,i),gazsat(:,i))
    ENDDO

    ! Computes balance :
    ! Each ice components has been treated independently from the others, and
    ! their tendencies are returned as is (just converted in input units).
    !
    ! Aerosols and CCN distribution evolution depends on the ice components:
    !   - For nucleation only creation of CCN can occur.
    !   - For condensation only loss of CCN can occur.
    ! We use the simple following rule :
    !   The global variation of CCN (and thus aerosols) is determined from the
    !   most intense activity of the ice components.
    !   that is the maximum value of the CCN tendencies regardless of its sign.
    DO i=1, mm_nla
      idx = MAXLOC(zdm0n(i,:),DIM=1) ! WARNING this is not the definition above (should be in ABS() func)
      dm0n(i) = zdm0n(i,idx)
      dm3n(i) = zdm3n(i,idx)
      dm0a(i) = zdm0a(i,idx)
      dm3a(i) = zdm3a(i,idx)
      ! all ice are returned but we must convert their units
      dm3i(i,:) = dm3i(i,:)
    ENDDO
  END SUBROUTINE mm_cloud_nucond

  SUBROUTINE nc_esp(xESP,vapX,m3iX,dvapX,dm3iX,dm0aer,dm3aer,dm0ccn,dm3ccn,Xsat)
    !! Get moments tendencies through nucleation/condensation/evaporation of a given condensible specie.
    !!
    !! The method computes the global tendencies of the aerosols, ccn and "ice" moments through cloud 
    !! microphysics processes (nucleation & condensation).
    !!
    !! @warning
    !! Input quantities __m3iX__,__m3iO__, __m0aer__,__m3aer__, __m0ccn__,__m3ccn__ are assumed to be in 
    !! \(X.kg^{-1}\) (where X is the unit of the moment that is, a number for M0 and a volume - \(m^3\) 
    !! for M3) ; __vapX__ must be expressed in term of molar fraction.
    TYPE(mm_esp), INTENT(in)                   :: xESP
      !! Condensate specie properties.
    REAL(kind=mm_wp),INTENT(in), DIMENSION(:)  :: vapX
      !! Gas specie molar fraction on the vertical grid from __TOP__ to __GROUND__ (\(mol.mol^{-1}\)).
    REAL(kind=mm_wp),INTENT(in), DIMENSION(:)  :: m3iX
      !! 3rd order moment of the ice component (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: dvapX
      !! Tendency of gas specie (\(mol.mol^{-1}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: dm3iX
      !! Tendency of the 3rd order moment of the ice component (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: dm0aer
      !! Tendency of the 0th order moment of the fractal mode distribution (\(m^{-3}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: dm3aer
      !! Tendency of the 3rd order moment of the fractal mode distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: dm0ccn
      !! Tendency of the 0th order moment of the ccn distribution (\(m^{-3}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: dm3ccn
      !! Tendency of the 3rd order moment of the ccn distribution (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp),INTENT(out), DIMENSION(:) :: Xsat
      !! Saturation ratio values on the vertical grid (--).
    INTEGER                                 :: i
    REAL(kind=mm_wp)                        :: bef,aft
    REAL(kind=mm_wp), DIMENSION(SIZE(vapX)) :: sm0a,sm3a,sm0n,sm3n,sm3iX
    REAL(kind=mm_wp), DIMENSION(SIZE(vapX)) :: zm0a,zm3a,zm0n,zm3n,zm3iX,zvapX
    REAL(kind=mm_wp), DIMENSION(SIZE(vapX)) :: pX,sig,qsat,lheat,seq,up,down, &
                                               ctot,newvap,nucr,grate,cm0,cm3
    ! Initialization : 
    ! Copy input argument and convert units X.m-3 -> X.kg-1) 
    ! sXXX is the initial converted value saved
    sm3iX = m3iX/mm_rhoair 
    sm0a = mm_m0aer_f/mm_rhoair ; sm3a = mm_m3aer_f/mm_rhoair
    sm0n = mm_m0ccn/mm_rhoair ; sm3n = mm_m3ccn/mm_rhoair
    ! zXXX is our working copy
    zm3ix = sm3ix ; zm0a = sm0a ; zm3a = sm3a ; zm0n = sm0n ; zm3n = sm3n

    ! Molar fraction of X specie is set in mass mixing ratio 
    zvapX  = vapX  * xESP%fmol2fmas
    ! Surface tension
    sig = mm_sigX(mm_temp,xESP)
    ! X specie mass mixing ratio at saturation
    qsat = mm_qsatX(mm_temp,mm_play,xESP)
    ! partial pressure of X specie
    pX = vapX * mm_play
    ! Saturation ratio
    Xsat = zvapX / qsat
    ! Equilibrium saturation near the drop 
    seq = dexp(2._mm_wp*sig*xESP%masmol/(xESP%rho*mm_rgas*mm_temp*mm_drad))
    ! Latent heat released
    lheat = mm_lheatX(mm_temp,xESP)
    ! Gets nucleation rate (ccn radius is the monomer !)
    call nuc_rate((/(mm_rm, i=1,mm_nla)/),mm_temp,xESP,pX,Xsat,nucr)
    ! Gets growth rate
    call growth_rate(mm_temp,mm_play,pX/Xsat,xESP,seq,mm_drad,grate)
    ctot = zvapx + xESP%rho * m3iX
    up = vapx + mm_dt * grate * 4._mm_wp * mm_pi * xESP%rho * mm_drad * seq * zm0n
    down = 1._mm_wp + mm_dt * grate * 4._mm_wp * mm_pi * xESP%rho * mm_drad / qsat * zm0n
    ! gets new vapor X specie mass mixing ratio : cannot be greater than the
    ! total gas + ice and lower than nothing :)
    newvap = max(min(up/down,ctot),0._mm_wp)
    ! gets "true" growth rate
    grate = grate * (newvap/qsat - seq)
    
    ! computes tendencies through condensation
    ! 1) check for the specific case : NO ICE and SUBLIMATION
    WHERE (zm3iX <= 0._mm_wp .AND. grate <= 0._mm_wp) 
      ! no ice and sublimation : reset ice to 0
      zm3iX = 0._mm_wp
    ELSEWHERE
      ! update ice volume ...
      zm3iX = zm3iX + mm_dt*grate*4._mm_wp*mm_pi*mm_drad*zm0n
      ! ... and check if there ice left in the ccn
      WHERE (zm3ix <= 0._mm_wp)
        zm3ix = 0._mm_wp
        zm0a = zm0a + zm0n ; zm0n = 0._mm_wp
        zm3a = zm3a + zm3n ; zm3n = 0._mm_wp
      ENDWHERE
    ENDWHERE
    
    ! computes tendencies
    ! all of these tendencies are in X.kg-1
    dm0aer = zm0a - sm0a
    dm3aer = zm3a - sm3a
    dm0ccn = zm0n - sm0n
    dm3ccn = zm3n - sm3n
    dm3ix  = zm3ix - sm3ix
    ! and this one is in mol.mol-1
    dvapx  = -xESP%rho * dm3ix / xESP%fmol2fmas

    ! reset temporary arrays to initial values
    zm3ix = sm3ix ; zm0a = sm0a ; zm3a = sm3a ; zm0n = sm0n ; zm3n = sm3n
  
    ! Computes global tendencies (nucleation && condensation)
    !
    ! for nucleation we have the following equations:
    !   dMaer_(k)/dt = - dMccn_(k)/dt (conservation of aerosols+ccn)       (1)
    !   dMaer_(k)/dt = - 4*PI*nucr/rm * Maer_(k+3)                         (2)
    !                = - 4*PI*nucr/rm * alpha(k+3)/alpha(k) * Maer_(k)
    !   where alpha(k+3) = Maer_(k+3)/Maer_(0) /  rc**(k+3)
    ! We solve (implicit scheme) :
    !   CONST_M(k) = 4*PI*nucr/rm * alpha(k+3)/alpha(k)*rc**3 * dt
    !   Maer_(k)[t+dt] = 1/(1+CONST_M(k)) * Maer_(k)[t]                    (3)
    ! Then, from eq. 2:
    ! Mccn_(k)[t+dt] = Mccn_(k)[t] + CONST_M(k)/(1+CONST_M(k))*Maer_(k)[t] (4)
    ! 
    cm0 = 4._mm_wp*mm_pi*nucr/mm_rm*mm_alpha_f(3._mm_wp)*mm_rcf**3*mm_dt
    cm3 = 4._mm_wp*mm_pi*nucr/mm_rm*mm_alpha_f(6._mm_wp)/mm_alpha_f(3._mm_wp)*mm_rcf**3*mm_dt
    zm0a = 1._mm_wp/(1._mm_wp+cm0) * zm0a
    zm3a = 1._mm_wp/(1._mm_wp+cm0) * zm3a
    WHERE (zm0a <= 0._mm_wp .OR. zm3a <= 0._mm_wp) 
      zm0a=0._mm_wp
      zm3a=0._mm_wp
      zm0n = zm0n + sm0a
      zm3n = zm3n + sm3a
    ELSEWHERE
      zm0n = zm0n + cm0/(1.+cm0)*zm0a
      zm3n = zm3n + cm3/(1.+cm3)*zm3a
    ENDWHERE
  
    ! Adds condensation tendencies
    zm0a = zm0a + dm0aer
    zm3a = zm3a + dm3aer
    zm0n = zm0n + dm0ccn
    zm3n = zm3n + dm3ccn

    ! Computes balance
    IF (mm_debug) THEN
      WRITE(*,'(a)') "Condensation/nucleation balance :"
      DO i=1,mm_nla
        bef = sm0a(i) + sm0n(i)
        aft = zm0a(i)  + zm0n(i)
        IF (ABS(bef-aft)/bef > 1e-10_mm_wp) WRITE(*,'((a),I2.2,(a))') &
        "[WARNING] nc_esp speaking: Total number not conserved (z=",i,")"
        bef = sm3a(i) + sm3n(i)
        aft = zm3a(i)  + zm3n(i)
        IF (ABS(bef-aft)/bef > 1e-10_mm_wp) WRITE(*,'((a),I2.2,(a))') &
        "[WARNING] nc_esp speaking: Total volume not conserved (z=",i,")"
      ENDDO
    ENDIF

    ! Now updates tendencies
    dm0aer = (zm0a-sm0a)*mm_rhoair
    dm3aer = (zm3a-sm3a)*mm_rhoair
    dm0ccn = (zm0n-sm0n)*mm_rhoair
    dm3ccn = (zm3n-sm3n)*mm_rhoair

  END SUBROUTINE nc_esp
  
  SUBROUTINE nuc_rate(rccn,temp,xESP,pvp,sat,rate)
    !! Get nucleation rate.
    !!
    !! The method computes the heterogeneous nucleation rate for the given specie on a fractal particle 
    !! of size __rccn__.
    !! Except __xESP__, all arguments are vectors of the same size (vertical grid).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: rccn !! Radius of the cloud condensation nuclei (m).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: temp !! Temperature (K).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: pvp  !! Partial vapor pressure of X specie (Pa).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: sat  !! Saturation ratio of given specie (--).
    TYPE(mm_esp), INTENT(in)                    :: xESP !! X specie properties (--).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: rate !! The nucleation rate (\(m^-{2}.s^{-1}\)).
    INTEGER                                 :: nv
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: nX,rstar,gstar,x,zeldov,deltaf, &
                                               fsh,fstar,sig
    nv = SIZE(rccn)
    ALLOCATE(nX(nv), rstar(nv), gstar(nv), x(nv), zeldov(nv), &
             deltaf(nv), fsh(nv), fstar(nv))
    ! Activation condition
    WHERE (sat > 1._mm_wp)
      sig = mm_sigX(temp,xESP)
      nX    = pvp/mm_kboltz/temp
      rstar = 2._mm_wp*sig*xESP%vol/(mm_kboltz*temp*dlog(sat))
      ! curvature radius
      x     = rccn/rstar
      fsh   = mm_fshape(xESP%mteta,x)
      fstar = (4._mm_wp/3._mm_wp*mm_pi)*sig*(rstar**2.)*fsh
      deltaf=MIN(MAX((2.*mm_fdes-mm_fdif-fstar)/(mm_kboltz*temp),-100._mm_wp),100._mm_wp)
      WHERE (deltaf == -100._mm_wp) 
        rate = 0._mm_wp
      ELSEWHERE
        gstar  = 4._mm_wp*mm_pi*(rstar**3)/(3._mm_wp*xESP%vol)
        zeldov = dsqrt(fstar/(3._mm_wp*mm_pi*mm_kboltz*temp*(gstar**2)))
        rate   = zeldov*mm_kboltz*temp*(nX*rstar)**2._mm_wp*dexp(deltaf)/ &
                  (fsh*mm_nus*xESP%mas)
      ENDWHERE
    ELSEWHERE
      rate = 0._mm_wp
    ENDWHERE
    DEALLOCATE(nX,rstar,gstar,x,zeldov,deltaf,fsh,fstar)
    RETURN
  END SUBROUTINE nuc_rate

  SUBROUTINE growth_rate(temp,pres,pXsat,xESP,seq,drad,rate)
    !! Get growth rate through condensation/evaporation process.
    !!
    !! The method computes the growth rate a drop through condensation/evaporation processes:
    !! 
    !! $$ r \times \frac{dr}{dt} = g_{rate} \times (S - S_{eq}) $$
    !!
    !! Except __xESP__ which is a scalar, all arguments are vectors of the same size (vertical grid).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: temp  !! Temperature (K).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: pres  !! Pressure level (Pa).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: seq   !! Saturation vapor pressure of specie (Pa).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: drad  !! Specie properties.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)  :: pXsat !! Equilibrium saturation near the drop.
    TYPE(mm_esp), INTENT(in)                    :: xESP  !! Drop radius (m).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: rate  !! Growth rate (\(m^{2}.s^{-1}\)).
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: k,knu,slf,rkc,rdc,l,dv
    INTEGER                                     :: n
    n = SIZE(temp)
    ALLOCATE(k(n),knu(n),slf(n),rkc(n),rdc(n),l(n),dv(n))

    ! N2 (air) Thermal conductivity (where does it come from ?)
    k(:) = ( 2.857e-2_mm_wp*temp-0.5428_mm_wp)*4.184e-3_mm_wp
    ! Gas mean free path
    l(:) = mm_lambda_g(temp,pres)
    ! Diffusion coefficient of X gas
    Dv(:) = 1._mm_wp/3._mm_wp*dsqrt(8._mm_wp*mm_rgas*temp(:)/(mm_pi*xESP%masmol))*mm_kboltz*temp(:) / &
         (mm_pi*pres(:)*(mm_air_rad+xESP%ray)**2* dsqrt(1._mm_wp+xESP%fmol2fmas))
    knu(:) = l(:)/drad(:)                                         ! The knudsen number of the drop
    slf(:) = (1.333_mm_wp+0.71_mm_wp/knu(:))/(1._mm_wp+1._mm_wp/knu(:)) ! Slip flow correction
    Dv(:)  = Dv(:)/(1._mm_wp+slf(:)*knu(:))
    ! latent heat resistance coefficient
    rkc(:) = mm_lheatX(temp(:),xESP)**2 * xESP%rho * xESP%masmol / (k(:)*mm_rgas*temp(:)**2)
    ! Diffusion resistance coefficient
    rdc(:) = mm_rgas * temp(:) * xESP%rho / (Dv(:)*pXsat(:)*xESP%masmol)
    ! Growth rate: rdr/dt = rate * (S-Seq) ; rate is returned
    rate(:) = 1._mm_wp / (seq(:)*rkc(:)+rdc(:))
    RETURN
  END SUBROUTINE growth_rate
 

  !-----------------------------------------------------------------------------
  ! SEDIMENTATION PROCESS RELATED METHODS
  !-----------------------------------------------------------------------------

  SUBROUTINE mm_cloud_sedimentation(dm0n,dm3n,dm3i)
    !! Compute the tendency of _clouds_ related moments through sedimentation process.
    !!
    !! The method computes the tendencies of moments related to cloud microphysics through 
    !! sedimentation process. The algorithm used here differs from 
    !! [[mm_haze(module):mm_haze_sedimentation(subroutine)]] as all moments settle with the same 
    !! terminal velocity which is computed with the average drop radius of the size distribution. 
    !! We simply compute an _exchange matrix_ that stores the new positions of each cells through 
    !! sedimentation process and then computes the matrix
    !! product with input moments values to get final tendencies.
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm0n
      !! Tendency of the 0th order moment of the ccn distribution (\(m^{-3}\)). 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm3n
      !! Tendency of the 3rd order moment of the ccn distribution (\(m^{3}.m^{-3}\)). 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:,:) :: dm3i
      !! Tendencies of the 3rd order moment of each ice component of the cloud (\(m^{3}m^{-3}\)). 
    INTEGER                                   :: im,nm
    REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE :: moms, momsf,chg_matrix
    nm = 2 + mm_nesp
    ALLOCATE(moms(mm_nla,nm),momsf(mm_nla,nm),chg_matrix(mm_nla,mm_nla))
    ! Initializes moms 
    moms(:,1)  = mm_m0ccn * mm_dzlev
    moms(:,2)  = mm_m3ccn * mm_dzlev
    DO im=1,mm_nesp
      moms(:,2+im) = mm_m3ice(:,im) * mm_dzlev
    ENDDO
    ! Computes exchange matrix
    CALL exchange(mm_drad,mm_drho,mm_dt,chg_matrix)
    ! Computes final moments values
    DO im=1,nm ; momsf(:,im) = MATMUL(chg_matrix,moms(:,im)) ; ENDDO
    ! Computes tendencies (converted in X.m-3)
    dm0n = (momsf(:,1)-moms(:,1))/mm_dzlev
    dm3n = (momsf(:,2)-moms(:,2))/mm_dzlev
    DO im=1,mm_nesp
      dm3i(:,im) = (momsf(:,2+im)-moms(:,2+im))/mm_dzlev
    ENDDO
    RETURN
  END SUBROUTINE mm_cloud_sedimentation

  SUBROUTINE exchange(rad,rhog,dt,matrix)
    !! Compute the exchange matrix.
    !!
    !! The subroutine computes the matrix exchange used by 
    !! [[mm_clouds(module):mm_cloud_sedimentation(subroutine)]] to compute moments tendencies 
    !! through sedimentation process. Both __rad__ and __rhog__ must be vector with relevant 
    !! values over the atmospheric vertical structure. __matrix__ is square 2D-array with same
    !! dimension size than __rad__.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: rad
      !! Cloud drop radius over the atmospheric vertical structure (m).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: rhog
      !! Cloud drop density over the atmospheric vertical structure (\(kg.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(in)               :: dt
      !! Timestep (s).
    REAL(kind=mm_wp), INTENT(out)              :: matrix(:,:)
      !! The output _exchange matrix_.
    INTEGER                                 :: nz,i,j,jj,jinf,jsup
    REAL(kind=mm_wp)                            :: zni,znip1,xf,xft,xcnt
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: puit
    REAL(kind=mm_wp)                            :: cpte,cpte2
    INTEGER, PARAMETER                      :: ichx  = 1
    matrix = 0._mm_wp ; nz = SIZE(rad) ; ALLOCATE(puit(nz))
    ! compute exchange matrix 
    DO i=1,nz
      puit(i) = 0._mm_wp
      xcnt = 0._mm_wp
      ! computes drop move (i.e. its new positions) 
      CALL getnzs(ichx,i,rad,rhog,dt,zni,znip1) 

      ! Peculiar case : Ground level precipitation [znip1<0 && (zni<0 || zni>0)]
      ! Complete precipitation [ znip1 <= 0 && zni <= 0 ] :
      IF(zni<=0._mm_wp.and.znip1<=0._mm_wp) THEN 
        xft=0._mm_wp
        xf=1._mm_wp
        xcnt=xcnt+xf
        puit(i)=puit(i)+xf
      ENDIF
      ! partial precipitation [ znip1 <= 0 && zni > 0 ] :
      IF (zni>0._mm_wp .and. znip1 <= 0._mm_wp) THEN
        xft=zni/(zni-znip1)
        xf=(1.-xft)
        xcnt=xcnt+xf
        puit(i)=puit(i)+xf
      ENDIF
      ! General case: no ground precipitation [ znip1 > 0 && zni > 0 ]
      IF (zni>0._mm_wp.and.znip1>0._mm_wp) THEN
        xft = 1._mm_wp       ! on a la totalite de la case
        xf  = 0._mm_wp
        xcnt=xcnt+xf
        puit(i)=puit(i)+xf
      ENDIF
      ! Fix minimum level to the ground 
      znip1 = MAX(znip1,0.)
      zni   = MAX(zni,0.)
      ! Locate new "drop" position in the verical grid
      jsup=nz+1
      jinf=nz+1
      DO j=1,nz
        IF (zni<=mm_zlev(j).and.zni>=mm_zlev(j+1)) jsup=j
        IF (znip1<=mm_zlev(j).and.znip1>=mm_zlev(j+1)) jinf=j
      ENDDO
      ! Volume is out of range: (all drops have touched the ground!)
      ! Note: can happen here,it has been treated previously :)
      IF (jsup>=nz+1.and.jinf==jsup) THEN 
        WRITE(*,'(a)') "[EXCHANGE] speaking: The impossible happened !"
        call EXIT(666)
      ENDIF
      ! Volume inside a single level
      IF (jsup==jinf.and.jsup<=nz) THEN 
        xf=1._mm_wp
        xcnt=xcnt+xft*xf
        matrix(jinf,i)=matrix(jinf,i)+xft*xf
      ENDIF 
  
      ! Volume over 2 levels 
      IF (jinf==jsup+1) THEN 
        xf=(zni-mm_zlev(jinf))/(zni-znip1)
        xcnt=xcnt+xf*xft
        IF(jsup<=nz) THEN
          matrix(jsup,i)=matrix(jsup,i)+xft*xf
        ENDIF
        xf=(mm_zlev(jinf)-znip1)/(zni-znip1)
        xcnt=xcnt+xf*xft
        IF (jinf<=nz) THEN
          matrix(jinf,i)=matrix(jinf,i)+xft*xf
        ENDIF
      ENDIF
  
      ! Volume over 3 or more levels
      IF (jinf > jsup+1) THEN
        ! first cell
        xf=(zni-mm_zlev(jsup+1))/(zni-znip1)
        xcnt=xcnt+xf*xft
        matrix(jsup,i)=matrix(jsup,i)+xft*xf
        ! last cell
        xf=(mm_zlev(jinf)-znip1)/(zni-znip1)
        xcnt=xcnt+xf*xft
        matrix(jinf,i)=matrix(jinf,i)+xft*xf
        ! other :)
        DO jj=jsup+1,jinf-1
          xf=(mm_zlev(jj)-mm_zlev(jj+1))/(zni-znip1)
          xcnt=xcnt+xf*xft
          matrix(jj,i)=matrix(jj,i)+xft*xf
        ENDDO
      ENDIF 
    ENDDO 
    ! checking if everything is alright if debug enabled... 
    IF (mm_debug) THEN
      cpte=0._mm_wp ; cpte2=0._mm_wp
      DO j=1,nz
        DO jj=1,nz
          cpte=cpte+matrix(jj,j)
        ENDDO
        cpte2=cpte+puit(j)
      ENDDO
      IF (abs(cpte2-nz)>1.e-4_mm_wp) THEN
        WRITE(*,'(a)')"[EXCHANGE] speaking :"
        WRITE(*,'("tx expl (/nz):",2(2X,ES10.3))') cpte,cpte2
      ENDIF
    ENDIF
    RETURN 
  END SUBROUTINE exchange

  SUBROUTINE getnzs(ichx,idx,rad,rho,dt,zni,zns)
    !! Compute displacement of a cell under sedimentation process.
    !!
    !! The method computes the new position of a _drop cell_ through sedimentation process as 
    !! descibed in the following scheme:
    !!
    !! ![Cloud sedimentation scheme](|media|/cloud_sed_scheme.svg)
    !!
    !! New positions are returned in __zni__ and __zns__ ouptut arguments.
    !!
    !! @note
    !! The method uses directly [[mm_globals(module):mm_play(variable)]], [[mm_globals(module):mm_plev(variable)]],
    !! [[mm_globals(module):mm_temp(variable)]],[[mm_globals(module):mm_btemp(variable)]], 
    !! [[mm_globals(module):mm_zlay(variable)]] and [[mm_globals(module):mm_zlev(variable)]] and uses __idx__ to
    !! get the relevant value to use on the vertical grid.
    INTEGER, INTENT(in)                        :: ichx
      !! Velocity extrapolation control flag (0 for linear, 1 for exponential -preferred -).
    INTEGER, INTENT(in)                        :: idx
      !! Initial position of the drop (subscript of vertical layers vectors).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: rad
      !! Cloud drop radius over the atmospheric vertical structure (m).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: rho
      !! Cloud drop density over the atmospheric vertical structure (\(kg.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(in)               :: dt
      !! Timestep (s).
    REAL(kind=mm_wp), INTENT(out)              :: zni
      !! Final layer center position (m).
    REAL(kind=mm_wp), INTENT(out)              :: zns
      !! Final layer upper boundary position (m).
    REAL(kind=mm_wp)            :: ws,wi,w,zi,zs 
    REAL(kind=mm_wp)            :: alpha,argexp,v0,arg1,arg2
    INTEGER                 :: i,nz
    REAL(kind=mm_wp), PARAMETER :: es = 30._mm_wp
    nz = SIZE(rad)
    ! Linear extrapolation of velocity
    IF (ichx==0) THEN
      ! velocity upper interface
      ws = wsettle(mm_plev(idx),mm_btemp(idx),mm_zlev(idx),rho(idx),rad(idx))
      IF (idx==nz) THEN
        ! veloctity center layer
        wi = wsettle(mm_play(idx),mm_temp(idx),mm_zlay(idx),rho(idx),rad(idx))
      ELSEIF(idx<nz) THEN
        ! velocity lower interface
        wi = wsettle(mm_plev(idx+1),mm_btemp(idx+1),mm_zlev(idx+1), &
                     rho(idx+1),rad(idx+1))
      ELSE 
        WRITE(*,'(a)') "[getnzs] speaking:" 
        WRITE(*,'(a)') "This is the fatal error..."
        WRITE(*,'(a)') "index is higher than number of levels"
        call EXIT(111)
      ENDIF
      w   = (ws+wi)/2._mm_wp
      zni = mm_zlev(idx)-w*dt
      zns = mm_zlev(idx)-mm_dzlev(idx)-w*dt
      RETURN
    ! Exponential extrapolation of velocity
    ELSEIF(ichx==1) THEN
      ws = wsettle(mm_plev(idx),mm_btemp(idx),mm_zlev(idx),rho(idx),rad(idx))
      zs = mm_zlev(idx)
      wi = wsettle(mm_play(idx),mm_temp(idx),mm_zlay(idx),rho(idx),rad(idx))
      zi=mm_zlay(idx)
      ! ws & wi must be different !
      IF(dabs(wi-ws)/wi <= 1.e-3_mm_wp)  wi=ws/1.001_mm_wp
      IF (wi /= 0._mm_wp) alpha = dlog(ws/wi)/(zs-zi)  ! alpha < 0 if wi > ws
      !   -es < argexp < es 
      argexp=MAX(MIN(alpha*zs,es),-es) 
      v0 = ws/dexp(argexp)
      arg1=1._mm_wp+v0*alpha*dexp(argexp)*dt
      argexp=MAX(MIN(alpha*(mm_zlev(idx)-mm_dzlev(idx)),es),-es)
      arg2=1._mm_wp+v0*alpha*dexp(argexp)*dt
      IF (arg1<=0._mm_wp.OR.arg2<=0._mm_wp) THEN 
        ! correct velocity
        ! divides the velocity argument in arg1 and arg2 :
        ! argX=1+alpha*v0*exp(alpha*z)*dt <==> argX-1=alpha*v0*exp(alpha*z)*dt
        ! ==> argX' = 1 + (argX-1)/2  <==> argX' = (1+argX)/2.
        DO i=1,25
          IF (arg1<=0._mm_wp.OR.arg2<=0._mm_wp) THEN
            IF (mm_debug) & 
            WRITE(*,'((a),I2.2,(a))') "[getnzs] must adjust velocity (",i,"/25)"
            arg1=(arg1+1._mm_wp)/2._mm_wp ; arg2=(arg2+1._mm_wp)/2._mm_wp
          ELSE
            EXIT 
          ENDIF
        ENDDO
        ! sh** we have to stop
        IF (i>25) THEN
          WRITE(*,'(a)')"[getnzs] speaking:"
          WRITE(*,'(a)') "Cannot adjust velocities"
          call EXIT(111)
        ENDIF
      ENDIF
      zni = mm_zlev(idx)-dlog(arg1)/alpha
      zns = mm_zlev(idx)-mm_dzlev(idx)-dlog(arg2)/alpha
      RETURN
    ENDIF  
  END SUBROUTINE getnzs

  ELEMENTAL FUNCTION wsettle(p,t,z,rho,rad) RESULT(w)
    !! Compute the settling velocity of a spherical particle. 
    !!
    !! The method computes the effective settling velocity of spherical particle of 
    !! radius __rad__. It accounts for the slip-flow transition (no approximation).
    REAL(kind=mm_wp), INTENT(in) :: p   !! The pressure level (Pa).
    REAL(kind=mm_wp), INTENT(in) :: t   !! The temperature (K).
    REAL(kind=mm_wp), INTENT(in) :: z   !! The altitude level (m).
    REAL(kind=mm_wp), INTENT(in) :: rho !! Density of the particle (\(kg.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(in) :: rad !! Radius of the particle (m).
    REAL(kind=mm_wp) :: w               !! Settling velocity (\(m.s^{-1}\)).
    REAL(kind=mm_wp)            :: g,a,kn,nu 
    REAL(kind=mm_wp), PARAMETER :: ra = 1.75e-10_mm_wp, nu0 = 1.74e-4_mm_wp, c = 109._mm_wp
    ! Computes corrected gravity
    g = mm_effg(z)
    ! Knudsen number
    kn = mm_kboltz*t/(p*4._mm_wp*sqrt(2._mm_wp)*mm_pi*ra**2)/rad
    ! Air viscosity
    nu=nu0*sqrt(t/293._mm_wp)*(1._mm_wp+c/293._mm_wp)/(1._mm_wp+c/t)
    ! Computes settling velocity
    w = 2._mm_wp/9._mm_wp * rad**2*g*rho/nu
    ! apply slip-flow correction
    w = w*(1._mm_wp+1.2517_mm_wp*kn+0.4_mm_wp*kn*dexp(-1.1_mm_wp/kn)) 
  END FUNCTION wsettle

  FUNCTION get_mass_flux(rho,m3) RESULT(flx)
    !> Get the mass flux of (clouds related) moment through sedimention.
    !!
    !! @warning
    !! The method is __only__ valid for cloud moments (i.e. ice or ccn). It calls
    !! [[mm_clouds(module):wsettle(function)]] that compute the _mean_ settling velocity of a
    !! cloud drop.
    !!
    !! @note
    !! The computed flux is always positive. 
    REAL(kind=mm_wp), INTENT(in)               :: rho
      !! Tracer density (\(kg.m^{-3}\)).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: m3
      !! Vertical profile of the total volume of tracer (i.e. M3) from __TOP__ to __GROUND__ (\(m^{3}.m^{-3}\)).
    REAL(kind=mm_wp), DIMENSION(SIZE(m3)) :: flx
      !! Mass sedimentation fluxes at each layer from __TOP__ to __GROUND__ (\(kg.m^{-2}.s^{-1}\)).
    REAL(kind=mm_wp), SAVE :: fac = 4._mm_wp/3._mm_wp * mm_pi
    flx = fac * rho * m3 * wsettle(mm_play,mm_temp,mm_zlay,mm_drho,mm_drad)
    RETURN
  END FUNCTION get_mass_flux

END MODULE MM_CLOUDS
