SUBROUTINE calchim(ngrid,qy_c,declin,dtchim,            &
     ctemp,cpphi,cpphis,cplay,cplev,czlay,czlev,dqyc)

  !---------------------------------------------------------------------------------------------------------
  !  
  ! Purpose : Interface subroutine to photochemical C model for Titan GCM.
  ! -------
  !           The subroutine computes the chemical processes for a single vertical column.
  !
  ! - Only tendencies are returned.
  ! - With moyzon_ch=.true. and input vectors zonally averaged
  !   the calculation is done only once per lat. band
  !
  ! Authors: + S. Lebonnois : 01/2000 | 09/2003 
  ! -------                  adaptation for Titan 3D : 02/2009
  !                          adaptation for // : 04/2013
  !                          extension chemistry up to 1300km : 10/2013
  !      
  !          + J. Vatant d'Ollone
  !               + 02/17 - adaptation for the new generic-forked physics
  !               + 01/18 - 03/18 - Major transformations :
  !                   - Upper chemistry fields are now stored in startfi 
  !                     and defined on a pressure grid from Vervack profile
  !                   - These modifs enables to run chemistry with others resolution than 32x48x55 !
  !                   - Only the actinic fluxes are still read in a 49-lat input but interp. on lat grid
  !                   - Chemistry can still be done in 2D 
  !                      -> Calcul. once per band lat and put same tendency in all longi. 
  !                         Check for negs in physiq_mod.
  !                      -> If procs sharing a lat band, no problem, the calcul will just be done twice.
  !                      -> Will not work with Dynamico, where the chemistry will have to be done in 3D.
  !                          ( and there'll be work to do to get rid of averaged fields )
  !
  !               + 02/19 : To always have correct photodissociations rates, altitudes sent here by physiq are always 
  !                         calculated with effective g - and with reference to the body not the local surface -
  !                          even if in physiq we keep altitudes coherent with dynamics !
  !
  ! + STILL TO DO : + Replug the interaction with haze (cf titan.old) -> to see with JB.
  !                 + Use iso_c_binding for the fortran-C exchanges.
  !---------------------------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------
  ! Structure :
  ! -----------
  !   0.  Declarations
  !   I.  Init and firstcall
  !         1. Read and store Vervack profile
  !         2. Compute planetar averaged atm. properties
  !         3. Init compound caracteristics
  !         4. Init photodissociations rates from actinic fluxes
  !         5. Init chemical reactions
  !         6. Init eddy diffusion coeff
  !   II. Loop on latitudes/grid-points
  !         0. Check on 2D chemistry
  !           1. Compute atm. properties at grid points
  !           2. Interpolate photodissociation rates at lat,alt,dec
  !           3. Read composition
  !           4. Call main solver gptitan C routine
  !           5. Calculate output tendencies on advected tracers
  !           6. Update upper chemistry fields
  !         0bis. If 2D chemsitry, don't recalculate if needed
  ! -----------------------------------------------------------------

  USE, INTRINSIC :: iso_c_binding
  USE comchem_h
  USE dimphy
  USE datafile_mod, ONLY: datadir
  USE comcstfi_mod, ONLY: g, rad, pi, r, kbol
  USE geometry_mod, ONLY: latitude
  USE logic_mod, ONLY: moyzon_ch
  USE moyzon_mod, ONLY: tmoy, playmoy

  IMPLICIT NONE

! ------------------------------------------
! *********** 0. Declarations *************
! ------------------------------------------

  ! Arguments
  ! ---------

  INTEGER, INTENT(IN)                              :: ngrid       ! Number of atmospheric columns.
  REAL*8, DIMENSION(ngrid,klev,nkim), INTENT(IN)   :: qy_c        ! Chemical species on GCM layers after adv.+diss. (mol/mol).
  REAL*8, INTENT(IN)                               :: declin      ! Solar declination (rad).
  REAL*8, INTENT(IN)                               :: dtchim      ! Chemistry timsetep (s).
  REAL*8, DIMENSION(ngrid,klev),      INTENT(IN)   :: ctemp       ! Mid-layer temperature (K).
  REAL*8, DIMENSION(ngrid,klev),      INTENT(IN)   :: cpphi       ! Mid-layer geopotential (m2.s-2).
  REAL*8, DIMENSION(ngrid),           INTENT(IN)   :: cpphis      ! Surface geopotential (m2.s-2).
  REAL*8, DIMENSION(ngrid,klev),      INTENT(IN)   :: cplay       ! Mid-layer pressure (Pa).
  REAL*8, DIMENSION(ngrid,klev+1),    INTENT(IN)   :: cplev       ! Inter-layer pressure (Pa).
  REAL*8, DIMENSION(ngrid,klev),      INTENT(IN)   :: czlay       ! Mid-layer effective altitude (m) : ref = geoid.
  REAL*8, DIMENSION(ngrid,klev+1),    INTENT(IN)   :: czlev       ! Inter-layer effective altitude (m) ref = geoid.

  REAL*8, DIMENSION(ngrid,klev,nkim), INTENT(OUT)  :: dqyc        ! Chemical species tendencies on GCM layers (mol/mol/s).

  ! Local variables :
  ! -----------------

  INTEGER :: i , l, ic, ig, igm1

  INTEGER :: dec, idec, ipres, ialt, klat

  REAL*8  :: declin_c  ! Declination (deg).
  REAL*8  :: factp, factalt, factdec, factlat, krpddec, krpddecp1, krpddecm1
  REAL*8  :: temp1, logp

  ! Variables sent into chemistry module (must be in double precision)
  ! ------------------------------------------------------------------

  REAL*8, DIMENSION(nlaykim_tot) :: temp_c  ! Temperature (K).
  REAL*8, DIMENSION(nlaykim_tot) :: press_c ! Pressure (Pa).
  REAL*8, DIMENSION(nlaykim_tot) :: phi_c   ! Geopotential (m2.s-2) - actually not sent in chem. module but used to compute alts.
  REAL*8, DIMENSION(nlaykim_tot) :: nb      ! Density (cm-3).

  REAL*8, DIMENSION(nlaykim_tot,nkim) :: cqy   ! Chemical species in whole column (mol/mol) sent to chem. module.
  REAL*8, DIMENSION(nlaykim_tot,nkim) :: cqy0  !     "      "     "   "      "        "     before modifs.

  REAL*8  :: surfhaze(nlaykim_tot)
  REAL*8  :: cprodaer(nlaykim_tot,4), cmaer(nlaykim_tot,4)
  REAL*8  :: ccsn(nlaykim_tot,4), ccsh(nlaykim_tot,4)

  REAL*8, DIMENSION(nlaykim_tot) :: rmil   ! Mid-layer distance (km) to planetographic center.
  REAL*8, DIMENSION(nlaykim_tot) :: rinter ! Inter-layer distance (km) to planetographic center (RA grid in chem. module).
  ! NB : rinter is on nlaykim_tot too, we don't care of the uppermost layer upper boundary altitude.

  ! Saved variables initialized at firstcall
  ! ----------------------------------------

  LOGICAL, SAVE :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall)

  REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: kedd ! Eddy mixing coefficient for upper chemistry (cm^2.s-1)
!$OMP THREADPRIVATE(kedd)

  REAL*8, DIMENSION(:,:), ALLOCATABLE, SAVE :: md   ! Mean molecular diffusion coefficients (cm^2.s-1)
  REAL*8, DIMENSION(:),   ALLOCATABLE, SAVE :: mass ! Molar mass of the compounds (g.mol-1)
!$OMP THREADPRIVATE(mass,md)

  REAL*8, DIMENSION(:), ALLOCATABLE, SAVE :: r1d, ct1d, p1d, t1d ! Vervack profile
  ! JVO 18 : No threadprivate for those as they'll be read in tcp.ver by master 

  REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: krpd  ! Photodissociations rate table
  REAL*8, DIMENSION(:,:)    , ALLOCATABLE, SAVE :: krate ! Reactions rate ( photo + chem )
!$OMP THREADPRIVATE(krpd,krate)

  INTEGER, DIMENSION(:),     ALLOCATABLE, SAVE :: nom_prod, nom_perte
  INTEGER, DIMENSION(:,:),   ALLOCATABLE, SAVE :: reactif
  INTEGER, DIMENSION(:,:),   ALLOCATABLE, SAVE :: prod
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: perte
!$OMP THREADPRIVATE(nom_prod,nom_perte,reactif,prod,perte)

  ! TEMPORARY : Dummy parameters without microphysics
  ! Here to keep the whole stuff running without modif chem. module
  ! ---------------------------------------------------------------

  INTEGER  :: utilaer(16)
  INTEGER  :: aerprod = 0
  INTEGER  :: htoh2   = 0

  ! -----------------------------------------------------------------------
  ! ***************** I. Initialisations and Firstcall ********************
  ! -----------------------------------------------------------------------

  IF (firstcall) THEN

     PRINT*, 'CHIMIE, premier appel'

     call check(nlaykim_tot,klev-15,nlrt_kim,nkim)

     ALLOCATE(r1d(131))
     ALLOCATE(ct1d(131))
     ALLOCATE(p1d(131))
     ALLOCATE(t1d(131))

     ALLOCATE(md(nlaykim_tot,nkim))
     ALLOCATE(mass(nkim))
     
     ALLOCATE(kedd(nlaykim_tot))
     
     ALLOCATE(krate(nlaykim_tot,nr_kim))
     ALLOCATE(krpd(nd_kim+1,nlrt_kim,15,nlat_actfluxes))

     ALLOCATE(nom_prod(nkim))
     ALLOCATE(nom_perte(nkim))
     ALLOCATE(reactif(5,nr_kim))
     ALLOCATE(prod(200,nkim))
     ALLOCATE(perte(2,200,nkim))
     
     ! 0. Deal with characters for C-interoperability
     ! ----------------------------------------------
     ! NB ( JVO 19 ) : Using iso_c_binding would do things in an even cleaner way !
     DO ic=1,nkim
       nomqy_c(ic) = trim(cnames(ic))//char(0) ! Add the C null terminator
     ENDDO
     nomqy_c(nkim+1)="HV"//char(0) ! For photodissociations

     ! 1. Read Vervack profile "tcp.ver", once for all
     ! -----------------------------------------------

!$OMP MASTER
     OPEN(11,FILE=TRIM(datadir)//'/tcp.ver',STATUS='old')
     READ(11,*)
     DO i=1,131
        READ(11,*) r1d(i), t1d(i), ct1d(i), p1d(i)
        ! For debug :
        ! -----------
        ! PRINT*, "tcp.ver", r1d(i), t1d(i), ct1d(i), p1d(i)
     ENDDO
     CLOSE(11)
!$OMP END MASTER
!$OMP BARRIER

     ! 2. Calculation of temp_c, densities and altitudes in planetary average
     ! ----------------------------------------------------------------------
     
     ! JVO18 : altitudes are no more calculated in firstcall, as I set kedd in pressure grid

     ! a. For GCM layers we just copy-paste ( assuming that physiq always send correct altitudes ! )

     PRINT*,'Init chemistry : pressure, density, temperature ... :'
     PRINT*,'level, press_c (mbar), nb (cm-3), temp_c (K)'
     
     IF (ngrid.NE.1) THEN
       DO l=1,klev
          temp_c(l)  = tmoy(l)                              ! K
          press_c(l) = playmoy(l)/100.                      ! mbar
          nb(l)      = 1.e-4*press_c(l) / (kbol*temp_c(l))  ! cm-3
          PRINT*, l, press_c(l), nb(l), temp_c(l)
       ENDDO
     ELSE
       DO l=1,klev
          temp_c(l)  = ctemp(1,l)                             ! K 
          press_c(l) = cplay(1,l)/100.                        ! mbar
          nb(l)      = 1.e-4*press_c(l) / (kbol*temp_c(l))  ! cm-3
          PRINT*, l, press_c(l), nb(l), temp_c(l)
       ENDDO
     ENDIF

     ! b. Extension in upper atmosphere with Vervack profile
     ! NB : Maybe the transition klev/klev+1 is harsh if T profile different from Vervack ...     

     ipres=1
     DO l=klev+1,nlaykim_tot
        press_c(l) = preskim(l-klev) / 100.0
        DO i=ipres,130
           IF ( (press_c(l).LE.p1d(i)) .AND. (press_c(l).GT.p1d(i+1)) ) THEN
              ipres=i
           ENDIF
        ENDDO
        factp = (press_c(l)-p1d(ipres)) / (p1d(ipres+1)-p1d(ipres))

        nb(l)      = exp( log(ct1d(ipres))*(1-factp) + log(ct1d(ipres+1))* factp )
        temp_c(l)  = t1d(ipres)*(1-factp) + t1d(ipres+1)*factp
        PRINT*, l , press_c(l), nb(l), temp_c(l)
     ENDDO

     ! 3. Compounds caracteristics
     ! ---------------------------
     mass(:) = 0.0
     call comp(nomqy_c,nb,temp_c,mass,md)
     PRINT*,'           Mass'
     DO ic=1,nkim
        PRINT*, nomqy_c(ic), mass(ic)
     ENDDO

     ! 4. Photodissociation rates
     ! --------------------------
     call disso(krpd,nlat_actfluxes) 

     ! 5. Init. chemical reactions with planetary average T prof.
     ! ----------------------------------------------------------

     !  NB : Chemical reactions rate are assumed to be constant within the T range of Titan's atm
     !  so we fill their krate once for all but krate for photodiss will be filled at each timestep
     
     call chimie(nomqy_c,nb,temp_c,krate,reactif, &
          nom_perte,nom_prod,perte,prod)

     ! 6. Eddy mixing coefficients (constant with time and space)
     ! ----------------------------------------------------------
     
     kedd(:) = 1.e3 ! Default value =/= zero

     ! NB : Eddy coeffs (e.g. Lavvas et al 08, Yelle et al 08) in altitude but they're rather linked to pressure
     !      Below GCM top we have dynamic mixing and for levs < nld=klev-15 the chem. solver ignores diffusion

     !! First calculate kedd for upper chemistry layers
     !DO l=klev-4,nlaykim_tot
     !   logp=-log10(press_c(l))
     !! 2E6 at 400 km ~ 10-2 mbar
     !   IF     ( logp.ge.2.0 .and. logp.le.3.0 ) THEN 
     !         kedd(l) = 2.e6 * 5.0**(logp-2.0)
     !! 1E7 at 500 km ~ 10-3 mbar
     !   ELSE IF     ( logp.ge.3.0 .and. logp.le.4.0 ) THEN 
     !         kedd(l) = 1.e7 * 3.0**(logp-3.0)
     !! 3E7 above 700 km ~ 10-4 mbar
     !   ELSEIF ( logp.gt.4.0                   ) THEN 
     !        kedd(l) = 3.e7
     !   ENDIF
     !ENDDO

     ! Kedd from (E7) in Vuitton 2019
     if (ngrid .eq. 1) then ! if 1D no dynamic mixing, we set the kedd in all column
       DO l=1,nlaykim_tot
         kedd(l) = 300.0 * ( 1.0E2 / press_c(l) )**1.5 * 3.0E7 /  &
                 ( 300.0 * ( 1.0E2 / press_c(l) )**1.5 + 3.0E7 ) 
       ENDDO
     else
       DO l=klev-4,nlaykim_tot
         ! JVO 18 : We keep the nominal profile in the GCM 5 upper layers 
         !          to have  a correct vertical mixing in the sponge layer
         kedd(l) = 300.0 * ( 1.0E2 / press_c(l) )**1.5 * 3.0E7 /  &
               ( 300.0 * ( 1.0E2 / press_c(l) )**1.5 + 3.0E7 ) 
       ENDDO
     endif
     
     if (ngrid .gt. 1) then ! not in 1D, no dynamic mixing
       ! Then adjust 10 layers profile fading to default value depending on kedd(ptop)
       DO l=klev-15,klev-5
          temp1   = ( log10(press_c(l)/press_c(klev-15)) ) / ( log10(press_c(klev-4)/press_c(klev-15)) )
          kedd(l) = 10.**( 3.0 + log10(kedd(klev-4)/1.e3) * temp1 )
       ENDDO
     endif
 
     firstcall = .FALSE.
  ENDIF  ! firstcall
  
  declin_c = declin*180./pi

  ! -----------------------------------------------------------------------
  ! *********************** II. Loop on latitudes *************************
  ! -----------------------------------------------------------------------
  
  DO ig=1,ngrid 

    IF (ig.eq.1) THEN 
        igm1=1
     ELSE
        igm1=ig-1
     ENDIF

     ! If 2D chemistry, trick to do the calculation only once per latitude band within the chunk
     ! NB1 : Will be obsolete with DYNAMICO, the chemistry will necessarly be 3D
     ! NB2 : Test of same latitude with dlat=0.1 : I think that if you run sims better than 1/10th degree then
     ! either it's with Dynamico and doesn't apply OR it is more than enough in terms of "preco / calc time" !
     ! -------------------------------------------------------------------------------------------------------

     IF ( ( moyzon_ch .AND. ( ig.EQ.1 .OR. (ABS(latitude(ig)-latitude(igm1)).GT.0.1*pi/180.0)) ) .OR. (.NOT. moyzon_ch) ) THEN

        ! 1. Compute altitude for the grid point with hydrostat. equilib.
        ! ---------------------------------------------------------------

        ! a. For GCM layers we just copy-paste
        ! JVO 19 : Now physiq always sent correct altitudes with effective g for chemistry ( even if it's not the case in physiq )

        DO l=1,klev
           rinter(l)  = (czlev(ig,l)+rad)/1000.0             ! km
           rmil(l)    = (czlay(ig,l)+rad)/1000.0             ! km
           temp_c(l)  = ctemp(ig,l)                          ! K
           phi_c(l)   = cpphi(ig,l)                          ! m2.s-2
           press_c(l) = cplay(ig,l)/100.                     ! mbar
           nb(l)      = 1.e-4*press_c(l) / (kbol*temp_c(l))  ! cm-3
        ENDDO
        rinter(klev+1)=(czlev(ig,klev+1)+rad)/1000.

        ! b. Extension in upper atmosphere with Vervack profile

        ipres=1
        DO l=klev+1,nlaykim_tot
          press_c(l) = preskim(l-klev) / 100.0
          DO i=ipres,130
              IF ( (press_c(l).LE.p1d(i)) .AND. (press_c(l).GT.p1d(i+1)) ) THEN
                ipres=i
              ENDIF
          ENDDO
          factp = (press_c(l)-p1d(ipres)) / (p1d(ipres+1)-p1d(ipres))

          nb(l)      = exp( log(ct1d(ipres))*(1-factp) + log(ct1d(ipres+1))* factp )
          temp_c(l)  = t1d(ipres)*(1-factp) + t1d(ipres+1)*factp
        ENDDO
 
        ! We build altitude with hydrostatic equilibrium on preskim grid with Vervack profile
        ! ( keeping in mind that preskim is built based on Vervack profile with dz=10km )

        DO l=klev+1,nlaykim_tot

           ! Compute geopotential on the upper grid with effective g to have correct altitudes
           ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           temp1 = 0.5*(temp_c(l-1)+temp_c(l)) ! interlayer temp
           phi_c(l) = phi_c(l-1) + r*temp1*log(press_c(l-1)/press_c(l)) ! Geopotential assuming hydrostatic equilibrium

           rmil(l)  =  ( g*rad*rad / (g*rad - ( phi_c(l) + cpphis(ig) ) ) ) / 1000.0 ! z(phi) with g varying with altitude with reference to the geoid
        ENDDO

        DO l=klev+2,nlaykim_tot
          rinter(l) = 0.5*(rmil(l-1) + rmil(l)) ! should be balanced with the intermediate pressure rather than 0.5 
        ENDDO

        ! 2. From krpd, compute krate for dissociations (declination-latitude-altitude interpolation)
        ! -------------------------------------------------------------------------------------------

        ! a. Calculate declination dependence

        if ((declin_c*10+267).lt.14.) then
           idec = 0
           dec  = 0
        else 
           if ((declin_c*10+267).gt.520.) then
              idec = 14
              dec  = 534
           else 
              idec = 1
              dec  = 27
              do while( (declin_c*10+267).ge.real(dec+20) )
                 dec  = dec+40
                 idec = idec+1
              enddo
           endif
        endif
        if ((declin_c.ge.-24.).and.(declin_c.le.24.)) then
           factdec = ( declin_c - (dec-267)/10. ) / 4.
        else
           factdec = ( declin_c - (dec-267)/10. ) / 2.7
        endif

        ! b. Calculate klat for interpolation on fixed latitudes of actinic fluxes input

        klat=1
        DO i=1,nlat_actfluxes
          IF (latitude(ig).LT.lat_actfluxes(i)) klat=i
        ENDDO
        IF (klat==nlat_actfluxes) THEN ! avoid rounding problems
          klat    = nlat_actfluxes-1
          factlat = 1.0
        ELSE
          factlat = (latitude(ig)-lat_actfluxes(klat))/(lat_actfluxes(klat+1)-lat_actfluxes(klat))
        ENDIF

        ! c. Altitude loop

        DO l=1,nlaykim_tot

           ! Calculate ialt for interpolation in altitude (krpd every 2 km)
           ialt    = int((rmil(l)-rad/1000.)/2.)+1
           factalt = (rmil(l)-rad/1000.)/2.-(ialt-1)

           ! Altitude can go above top limit of UV levels - in this case we keep the 1310km top fluxes
           IF (ialt.GT.nlrt_kim-1) THEN
             ialt     = nlrt_kim-1 ! avoid out-of-bound array
             factalt  = 1.0
           ENDIF

           DO i=1,nd_kim+1 ! nd_kim+1 is dissociation of N2 by GCR

                 krpddec   =   (   krpd(i,ialt  ,idec+1,klat)   * (1.0-factalt)                   &
                                 + krpd(i,ialt+1,idec+1,klat)   * factalt       ) * (1.0-factlat) &
                             + (   krpd(i,ialt  ,idec+1,klat+1) * (1.0-factalt)                   &
                                 + krpd(i,ialt+1,idec+1,klat+1) * factalt       ) * factlat

              if      ( factdec.lt.0. ) then 
                 krpddecm1 =   (   krpd(i,ialt  ,idec  ,klat)   * (1.0-factalt)                   &
                                 + krpd(i,ialt+1,idec  ,klat)   * factalt       ) * (1.0-factlat) &
                             + (   krpd(i,ialt  ,idec  ,klat+1) * (1.0-factalt)                   &
                                 + krpd(i,ialt+1,idec  ,klat+1) * factalt       ) * factlat
                 krate(l,i) = krpddecm1 * abs(factdec) + krpddec   * ( 1.0 + factdec)
              else if ( factdec.gt.0. ) then
                 krpddecp1 =   (   krpd(i,ialt  ,idec+2,klat)   * (1.0-factalt)                   &
                                 + krpd(i,ialt+1,idec+2,klat)   * factalt       ) * (1.0-factlat) &
                             + (   krpd(i,ialt  ,idec+2,klat+1) * (1.0-factalt)                   &
                                 + krpd(i,ialt+1,idec+2,klat+1) * factalt       ) * factlat
                 krate(l,i) = krpddecp1 * factdec      + krpddec   * ( 1.0 - factdec)
              else if ( factdec.eq.0. ) then 
                 krate(l,i) = krpddec
              endif

           ENDDO ! i=1,nd_kim+1
        ENDDO ! l=1,nlaykim_tot

        ! 3. Read composition 
        ! -------------------

        DO ic=1,nkim
           DO l=1,klev
              cqy(l,ic) = qy_c(ig,l,ic) ! advected tracers for the GCM part converted to molar frac.
           ENDDO
           
           DO l=1,nlaykim_up
              cqy(klev+l,ic) = ykim_up(ic,ig,l) ! ykim_up for the upper atm.
           ENDDO
        ENDDO

        cqy0(:,:) = cqy(:,:) ! Stores compo. before modifs

        ! 4. Call main Titan chemistry C routine
        ! --------------------------------------

        call gptitan(rinter,temp_c,nb,                  &
             nomqy_c,cqy,                               &
             dtchim,latitude(ig)*180./pi,mass,md,       &
             kedd,krate,reactif,                        &
             nom_prod,nom_perte,prod,perte,             &
             aerprod,utilaer,cmaer,cprodaer,ccsn,ccsh,  &
             htoh2,surfhaze)

        ! 5. Calculates tendencies on composition for advected tracers 
        ! ------------------------------------------------------------
        DO ic=1,nkim
           DO l=1,klev
              dqyc(ig,l,ic) = (cqy(l,ic) - cqy0(l,ic))/dtchim ! (mol/mol/s)
           ENDDO
        ENDDO

        ! 6. Update ykim_up
        ! -----------------
        DO ic=1,nkim
           DO l=1,nlaykim_up
              ykim_up(ic,ig,l) = cqy(klev+l,ic)
           ENDDO
        ENDDO
        ! NB: The full vertical composition grid will be created only for the outputs


     ELSE ! In 2D chemistry, if following grid point at same latitude, same zonal mean so don't do calculations again !
        dqyc(ig,:,:)    = dqyc(igm1,:,:) ! will be put back in 3D with longitudinal variations assuming same relative tendencies within a lat band
        ykim_up(:,ig,:) = ykim_up(:,igm1,:) ! no horizontal mixing in upper layers -> no longitudinal variations
     ENDIF

  ENDDO

END SUBROUTINE calchim
