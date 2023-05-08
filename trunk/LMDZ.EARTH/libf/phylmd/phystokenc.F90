SUBROUTINE phystokenc (nlon,nlev,pdtphys,rlon,rlat, &
     pt,pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
     pfm_therm,pentr_therm, &
     cdragh, pcoefh,yu1,yv1,ftsol,pctsrf, &
     frac_impa,frac_nucl, &
     pphis,paire,dtime,itap, &
     psh, pda, pphi, pmp, pupwd, pdnwd)
  
  USE ioipsl
  USE dimphy
  USE infotrac, ONLY : nqtot
  USE iophy
  USE control_mod
  
  IMPLICIT NONE
  
!======================================================================
! Auteur(s) FH
! Objet: Ecriture des variables pour transport offline
!
!======================================================================
  INCLUDE "dimensions.h"
  INCLUDE "tracstoke.h"
  INCLUDE "indicesol.h"
  INCLUDE "iniprint.h"
!======================================================================

! Arguments:
!
  REAL,DIMENSION(klon,klev), INTENT(IN)     :: psh   ! humidite specifique
  REAL,DIMENSION(klon,klev), INTENT(IN)     :: pda
  REAL,DIMENSION(klon,klev,klev), INTENT(IN):: pphi
  REAL,DIMENSION(klon,klev), INTENT(IN)     :: pmp
  REAL,DIMENSION(klon,klev), INTENT(IN)     :: pupwd ! saturated updraft mass flux
  REAL,DIMENSION(klon,klev), INTENT(IN)     :: pdnwd ! saturated downdraft mass flux

!   EN ENTREE:
!   ==========
!
!   divers:
!   -------
!
  INTEGER nlon ! nombre de points horizontaux
  INTEGER nlev ! nombre de couches verticales
  REAL pdtphys ! pas d'integration pour la physique (seconde)
  INTEGER itap
  INTEGER, SAVE :: physid
!$OMP THREADPRIVATE(physid)

!   convection:
!   -----------
!
  REAL pmfu(klon,klev)  ! flux de masse dans le panache montant
  REAL pmfd(klon,klev)  ! flux de masse dans le panache descendant
  REAL pen_u(klon,klev) ! flux entraine dans le panache montant
  REAL pde_u(klon,klev) ! flux detraine dans le panache montant
  REAL pen_d(klon,klev) ! flux entraine dans le panache descendant
  REAL pde_d(klon,klev) ! flux detraine dans le panache descendant
  REAL pt(klon,klev)
  REAL,ALLOCATABLE,SAVE :: t(:,:)
!$OMP THREADPRIVATE(t)
!
  REAL rlon(klon), rlat(klon), dtime
  REAL zx_tmp_3d(iim,jjm+1,klev),zx_tmp_2d(iim,jjm+1)

!   Couche limite:
!   --------------
!
  REAL cdragh(klon)          ! cdrag
  REAL pcoefh(klon,klev)     ! coeff melange CL
  REAL pcoefh_buf(klon,klev) ! coeff melange CL + cdrag
  REAL yv1(klon)
  REAL yu1(klon),pphis(klon),paire(klon)

!   Les Thermiques : (Abderr 25 11 02)
!   ---------------
  REAL, INTENT(IN) ::  pfm_therm(klon,klev+1)
  REAL pentr_therm(klon,klev)
  
  REAL,ALLOCATABLE,SAVE :: entr_therm(:,:)
  REAL,ALLOCATABLE,SAVE :: fm_therm(:,:)
!$OMP THREADPRIVATE(entr_therm)
!$OMP THREADPRIVATE(fm_therm)
!
!   Lessivage:
!   ----------
!
  REAL frac_impa(klon,klev)
  REAL frac_nucl(klon,klev)
!
! Arguments necessaires pour les sources et puits de traceur
!
  REAL ftsol(klon,nbsrf)  ! Temperature du sol (surf)(Kelvin)
  REAL pctsrf(klon,nbsrf) ! Pourcentage de sol f(nature du sol)
!======================================================================
!
  INTEGER i, k, kk
  REAL,ALLOCATABLE,SAVE :: mfu(:,:)  ! flux de masse dans le panache montant
  REAL,ALLOCATABLE,SAVE :: mfd(:,:)  ! flux de masse dans le panache descendant
  REAL,ALLOCATABLE,SAVE :: en_u(:,:) ! flux entraine dans le panache montant
  REAL,ALLOCATABLE,SAVE :: de_u(:,:) ! flux detraine dans le panache montant
  REAL,ALLOCATABLE,SAVE :: en_d(:,:) ! flux entraine dans le panache descendant
  REAL,ALLOCATABLE,SAVE :: de_d(:,:) ! flux detraine dans le panache descendant
  REAL,ALLOCATABLE,SAVE :: coefh(:,:) ! flux detraine dans le panache descendant
  
  REAL,ALLOCATABLE,SAVE :: pyu1(:)
  REAL,ALLOCATABLE,SAVE :: pyv1(:)
  REAL,ALLOCATABLE,SAVE :: pftsol(:,:)
  REAL,ALLOCATABLE,SAVE :: ppsrf(:,:)
!$OMP THREADPRIVATE(mfu,mfd,en_u,de_u,en_d,de_d,coefh)
!$OMP THREADPRIVATE(pyu1,pyv1,pftsol,ppsrf)


  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE     :: sh  
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE     :: da
  REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE   :: phi
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE     :: mp
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE     :: upwd
  REAL,DIMENSION(:,:), ALLOCATABLE,SAVE     :: dnwd
  
  REAL, SAVE :: dtcum
  INTEGER, SAVE:: iadvtr=0
!$OMP THREADPRIVATE(dtcum,iadvtr)
  REAL zmin,zmax
  LOGICAL ok_sync
  CHARACTER(len=12) :: nvar
!
!======================================================================

  iadvtr=iadvtr+1

! Dans le meme vecteur on recombine le drag et les coeff d'echange
  pcoefh_buf(:,1)      = cdragh(:)
  pcoefh_buf(:,2:klev) = pcoefh(:,2:klev)
  
  ok_sync = .TRUE.

! Initialization done only once
!======================================================================
  IF (iadvtr==1) THEN
     ALLOCATE( t(klon,klev))
     ALLOCATE( mfu(klon,klev))  
     ALLOCATE( mfd(klon,klev))  
     ALLOCATE( en_u(klon,klev)) 
     ALLOCATE( de_u(klon,klev)) 
     ALLOCATE( en_d(klon,klev)) 
     ALLOCATE( de_d(klon,klev)) 
     ALLOCATE( coefh(klon,klev)) 
     ALLOCATE( entr_therm(klon,klev))
     ALLOCATE( fm_therm(klon,klev))
     ALLOCATE( pyu1(klon))
     ALLOCATE( pyv1(klon))
     ALLOCATE( pftsol(klon,nbsrf))
     ALLOCATE( ppsrf(klon,nbsrf))
     
     ALLOCATE(sh(klon,klev)) 
     ALLOCATE(da(klon,klev)) 
     ALLOCATE(phi(klon,klev,klev)) 
     ALLOCATE(mp(klon,klev)) 
     ALLOCATE(upwd(klon,klev)) 
     ALLOCATE(dnwd(klon,klev)) 

     CALL initphysto('phystoke', dtime, dtime*istphy,dtime*istphy,physid)
     
     ! Write field phis and aire only once
     CALL histwrite_phy(physid,"phis",itap,pphis)
     CALL histwrite_phy(physid,"aire",itap,paire)
     CALL histwrite_phy(physid,"longitudes",itap,rlon)
     CALL histwrite_phy(physid,"latitudes",itap,rlat)

  END IF
  
  
! Set to zero cumulating fields
!======================================================================
  IF (MOD(iadvtr,istphy)==1.OR.istphy==1) THEN
     WRITE(lunout,*)'reinitialisation des champs cumules a iadvtr=',iadvtr
     mfu(:,:)=0.
     mfd(:,:)=0.
     en_u(:,:)=0.
     de_u(:,:)=0.
     en_d(:,:)=0.
     de_d(:,:)=0.
     coefh(:,:)=0.
     t(:,:)=0.
     fm_therm(:,:)=0.
     entr_therm(:,:)=0.
     pyv1(:)=0.
     pyu1(:)=0.
     pftsol(:,:)=0.
     ppsrf(:,:)=0.
     sh(:,:)=0.
     da(:,:)=0.
     phi(:,:,:)=0.
     mp(:,:)=0.
     upwd(:,:)=0.
     dnwd(:,:)=0.
     dtcum=0.
  ENDIF
  

! Cumulate fields at each time step
!======================================================================
  DO k=1,klev
     DO i=1,klon
        mfu(i,k)=mfu(i,k)+pmfu(i,k)*pdtphys
        mfd(i,k)=mfd(i,k)+pmfd(i,k)*pdtphys
        en_u(i,k)=en_u(i,k)+pen_u(i,k)*pdtphys
        de_u(i,k)=de_u(i,k)+pde_u(i,k)*pdtphys
        en_d(i,k)=en_d(i,k)+pen_d(i,k)*pdtphys
        de_d(i,k)=de_d(i,k)+pde_d(i,k)*pdtphys
        coefh(i,k)=coefh(i,k)+pcoefh_buf(i,k)*pdtphys
        t(i,k)=t(i,k)+pt(i,k)*pdtphys
        fm_therm(i,k)=fm_therm(i,k)+pfm_therm(i,k)*pdtphys
        entr_therm(i,k)=entr_therm(i,k)+pentr_therm(i,k)*pdtphys
        sh(i,k) = sh(i,k) + psh(i,k)*pdtphys
        da(i,k) = da(i,k) + pda(i,k)*pdtphys
        mp(i,k) = mp(i,k) + pmp(i,k)*pdtphys
        upwd(i,k) = upwd(i,k) + pupwd(i,k)*pdtphys
        dnwd(i,k) = dnwd(i,k) + pdnwd(i,k)*pdtphys
     ENDDO
  ENDDO

  DO kk=1,klev
     DO k=1,klev
        DO i=1,klon
           phi(i,k,kk) = phi(i,k,kk) + pphi(i,k,kk)*pdtphys
        END DO
     END DO
  END DO

  DO i=1,klon
     pyv1(i)=pyv1(i)+yv1(i)*pdtphys
     pyu1(i)=pyu1(i)+yu1(i)*pdtphys
  END DO
  DO k=1,nbsrf
     DO i=1,klon
        pftsol(i,k)=pftsol(i,k)+ftsol(i,k)*pdtphys
        ppsrf(i,k)=ppsrf(i,k)+pctsrf(i,k)*pdtphys
     ENDDO
  ENDDO
  
! Add time step to cumulated time
  dtcum=dtcum+pdtphys
  

! Write fields to file, if it is time to do so
!======================================================================
  IF(MOD(iadvtr,istphy)==0) THEN 

     ! normalize with time period
     DO k=1,klev
        DO i=1,klon
           mfu(i,k)=mfu(i,k)/dtcum
           mfd(i,k)=mfd(i,k)/dtcum
           en_u(i,k)=en_u(i,k)/dtcum
           de_u(i,k)=de_u(i,k)/dtcum
           en_d(i,k)=en_d(i,k)/dtcum
           de_d(i,k)=de_d(i,k)/dtcum
           coefh(i,k)=coefh(i,k)/dtcum
           t(i,k)=t(i,k)/dtcum	
           fm_therm(i,k)=fm_therm(i,k)/dtcum
           entr_therm(i,k)=entr_therm(i,k)/dtcum
           sh(i,k)=sh(i,k)/dtcum
           da(i,k)=da(i,k)/dtcum
           mp(i,k)=mp(i,k)/dtcum
           upwd(i,k)=upwd(i,k)/dtcum
           dnwd(i,k)=dnwd(i,k)/dtcum
        ENDDO
     ENDDO
     DO kk=1,klev
        DO k=1,klev
           DO i=1,klon
              phi(i,k,kk) = phi(i,k,kk)/dtcum
           END DO
        END DO
     END DO
     DO i=1,klon
        pyv1(i)=pyv1(i)/dtcum
        pyu1(i)=pyu1(i)/dtcum
     END DO
     DO k=1,nbsrf
        DO i=1,klon
           pftsol(i,k)=pftsol(i,k)/dtcum
           ppsrf(i,k)=ppsrf(i,k)/dtcum
        ENDDO
     ENDDO

     ! write fields
     CALL histwrite_phy(physid,"t",itap,t)
     CALL histwrite_phy(physid,"mfu",itap,mfu)
     CALL histwrite_phy(physid,"mfd",itap,mfd)
     CALL histwrite_phy(physid,"en_u",itap,en_u)
     CALL histwrite_phy(physid,"de_u",itap,de_u)
     CALL histwrite_phy(physid,"en_d",itap,en_d)
     CALL histwrite_phy(physid,"de_d",itap,de_d)
     CALL histwrite_phy(physid,"coefh",itap,coefh)	
     CALL histwrite_phy(physid,"fm_th",itap,fm_therm)
     CALL histwrite_phy(physid,"en_th",itap,entr_therm)
     CALL histwrite_phy(physid,"frac_impa",itap,frac_impa)
     CALL histwrite_phy(physid,"frac_nucl",itap,frac_nucl)
     CALL histwrite_phy(physid,"pyu1",itap,pyu1)
     CALL histwrite_phy(physid,"pyv1",itap,pyv1)
     CALL histwrite_phy(physid,"ftsol1",itap,pftsol(:,1))
     CALL histwrite_phy(physid,"ftsol2",itap,pftsol(:,2))
     CALL histwrite_phy(physid,"ftsol3",itap,pftsol(:,3))
     CALL histwrite_phy(physid,"ftsol4",itap,pftsol(:,4))
     CALL histwrite_phy(physid,"psrf1",itap,ppsrf(:,1))
     CALL histwrite_phy(physid,"psrf2",itap,ppsrf(:,2))
     CALL histwrite_phy(physid,"psrf3",itap,ppsrf(:,3))
     CALL histwrite_phy(physid,"psrf4",itap,ppsrf(:,4))
     CALL histwrite_phy(physid,"sh",itap,sh)
     CALL histwrite_phy(physid,"da",itap,da)
     CALL histwrite_phy(physid,"mp",itap,mp)
     CALL histwrite_phy(physid,"upwd",itap,upwd)
     CALL histwrite_phy(physid,"dnwd",itap,dnwd)


! phi
     DO k=1,klev
        IF (k<10) THEN
           WRITE(nvar,'(i1)') k
        ELSE IF (k<100) THEN
           WRITE(nvar,'(i2)') k
        ELSE
           WRITE(nvar,'(i3)') k
        END IF
        nvar='phi_lev'//trim(nvar)
        
        CALL histwrite_phy(physid,nvar,itap,phi(:,:,k))
     END DO
     
     ! Syncronize file
!$OMP MASTER
     IF (ok_sync) CALL histsync(physid)
!$OMP END MASTER
     
     
     ! Calculate min and max values for some fields (coefficients de lessivage)
     zmin=1e33
     zmax=-1e33
     DO k=1,klev
        DO i=1,klon
           zmax=MAX(zmax,frac_nucl(i,k))
           zmin=MIN(zmin,frac_nucl(i,k))
        ENDDO
     ENDDO
     WRITE(lunout,*)'------ coefs de lessivage (min et max) --------'
     WRITE(lunout,*)'facteur de nucleation ',zmin,zmax
     zmin=1e33
     zmax=-1e33
     DO k=1,klev
        DO i=1,klon
           zmax=MAX(zmax,frac_impa(i,k))
           zmin=MIN(zmin,frac_impa(i,k))
        ENDDO
     ENDDO
     WRITE(lunout,*)'facteur d impaction ',zmin,zmax
     
  ENDIF ! IF(MOD(iadvtr,istphy)==0)

END SUBROUTINE phystokenc
