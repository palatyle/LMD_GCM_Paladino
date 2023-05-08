! $Id: disvert.F90 1645 2012-07-30 16:01:50Z lguez $

SUBROUTINE disvert()

  ! Auteur : P. Le Van

#ifdef CPP_IOIPSL
  use ioipsl, only: getin
#else
  USE ioipsl_getincom, only: getin
#endif
  use new_unit_m, only: new_unit
  use assert_m, only: assert
  USE comvert_mod, ONLY: ap,bp,nivsigs,nivsig,preff,pa,presnivs,dpres,scaleheight
  USE logic_mod, ONLY: ok_strato

  IMPLICIT NONE

  include "dimensions.h"
  include "paramet.h"
  include "iniprint.h"

!-------------------------------------------------------------------------------
! Purpose: Vertical distribution functions for LMDZ.
!          Triggered by the levels number llm.
!-------------------------------------------------------------------------------
! Read    in "comvert_mod":

! pa !--- vertical coordinate is close to a PRESSURE COORDINATE FOR P
! < 0.3 * pa (relative variation of p on a model level is < 0.1 %)

! preff                      !--- REFERENCE PRESSURE                 (101325 Pa)
! Written in "comvert_mod":
! ap(llm+1), bp(llm+1)       !--- Ap, Bp HYBRID COEFFICIENTS AT INTERFACES
! aps(llm),  bps(llm)        !--- Ap, Bp HYBRID COEFFICIENTS AT MID-LAYERS
! dpres(llm)                 !--- PRESSURE DIFFERENCE FOR EACH LAYER
! presnivs(llm)              !--- PRESSURE AT EACH MID-LAYER
! scaleheight                !--- VERTICAL SCALE HEIGHT            (Earth: 8kms)
! nivsig(llm+1)              !--- SIGMA INDEX OF EACH LAYER INTERFACE
! nivsigs(llm)               !--- SIGMA INDEX OF EACH MID-LAYER
!-------------------------------------------------------------------------------
! Local variables:
  REAL sig(llm+1), dsig(llm)
  REAL sig0(llm+1), zz(llm+1)
  REAL zk, zkm1, dzk1, dzk2, z, k0, k1

  INTEGER l, unit
  REAL dsigmin
  REAL vert_scale,vert_dzmin,vert_dzlow,vert_z0low,vert_dzmid,vert_z0mid,vert_h_mid,vert_dzhig,vert_z0hig,vert_h_hig

  REAL alpha, beta, deltaz
  REAL x
  character(len=*),parameter :: modname="disvert"

  character(len=24):: vert_sampling
  ! (allowed values are "param", "tropo", "strato" and "read")

  !-----------------------------------------------------------------------

  WRITE(lunout,*) TRIM(modname)//" starts"

  ! default scaleheight is 8km for earth
  scaleheight=8.

  vert_sampling = merge("strato", "tropo ", ok_strato) ! default value
  call getin('vert_sampling', vert_sampling)
  WRITE(lunout,*) TRIM(modname)//' vert_sampling = ' // vert_sampling
  if (llm==39 .and. vert_sampling=="strato") then
     dsigmin=0.3 ! Vieille option par dÃ©faut pour CMIP5
  else
     dsigmin=1.
  endif
  call getin('dsigmin', dsigmin)
  WRITE(LUNOUT,*) trim(modname), 'Discretisation verticale DSIGMIN=',dsigmin


  select case (vert_sampling)
  case ("param")
     ! On lit les options dans sigma.def:
     OPEN(99, file='sigma.def', status='old', form='formatted')
     READ(99, *) scaleheight ! hauteur d'echelle 8.
     READ(99, *) deltaz ! epaiseur de la premiere couche 0.04
     READ(99, *) beta ! facteur d'acroissement en haut 1.3
     READ(99, *) k0 ! nombre de couches dans la transition surf
     READ(99, *) k1 ! nombre de couches dans la transition haute
     CLOSE(99)
     alpha=deltaz/(llm*scaleheight)
     write(lunout, *)trim(modname),':scaleheight, alpha, k0, k1, beta', &
                               scaleheight, alpha, k0, k1, beta

     alpha=deltaz/tanh(1./k0)*2.
     zkm1=0.
     sig(1)=1.
     do l=1, llm
        sig(l+1)=(cosh(l/k0))**(-alpha*k0/scaleheight) &
             *exp(-alpha/scaleheight*tanh((llm-k1)/k0) &
                  *beta**(l-(llm-k1))/log(beta))
        zk=-scaleheight*log(sig(l+1))

        dzk1=alpha*tanh(l/k0)
        dzk2=alpha*tanh((llm-k1)/k0)*beta**(l-(llm-k1))/log(beta)
        write(lunout, *)l, sig(l+1), zk, zk-zkm1, dzk1, dzk2
        zkm1=zk
     enddo

     sig(llm+1)=0.

     bp(: llm) = EXP(1. - 1. / sig(: llm)**2)
     bp(llmp1) = 0.

     ap = pa * (sig - bp)
  case("sigma")
     DO l = 1, llm
        x = 2*asin(1.) * (l - 0.5) / (llm + 1)
        dsig(l) = dsigmin + 7.0 * SIN(x)**2
     ENDDO
     dsig = dsig / sum(dsig)
     sig(llm+1) = 0.
     DO l = llm, 1, -1
        sig(l) = sig(l+1) + dsig(l)
     ENDDO

     bp(1)=1.
     bp(2: llm) = sig(2:llm)
     bp(llmp1) = 0.
     ap(:)=0.
  case("tropo")
     DO l = 1, llm
        x = 2*asin(1.) * (l - 0.5) / (llm + 1)
        dsig(l) = dsigmin + 7.0 * SIN(x)**2
     ENDDO
     dsig = dsig / sum(dsig)
     sig(llm+1) = 0.
     DO l = llm, 1, -1
        sig(l) = sig(l+1) + dsig(l)
     ENDDO

     bp(1)=1.
     bp(2: llm) = EXP(1. - 1. / sig(2: llm)**2)
     bp(llmp1) = 0.

     ap(1)=0.
     ap(2: llm + 1) = pa * (sig(2: llm + 1) - bp(2: llm + 1))
  case("strato")
     DO l = 1, llm
        x = 2*asin(1.) * (l - 0.5) / (llm + 1)
        dsig(l) =(dsigmin + 7. * SIN(x)**2) &
             *(0.5*(1.-tanh(1.*(x-asin(1.))/asin(1.))))**2
     ENDDO
     dsig = dsig / sum(dsig)
     sig(llm+1) = 0.
     DO l = llm, 1, -1
        sig(l) = sig(l+1) + dsig(l)
     ENDDO

     bp(1)=1.
     bp(2: llm) = EXP(1. - 1. / sig(2: llm)**2)
     bp(llmp1) = 0.

     ap(1)=0.
     ap(2: llm + 1) = pa * (sig(2: llm + 1) - bp(2: llm + 1))
  case("strato_correct")
!==================================================================
! Fredho 2014/05/18, Saint-Louis du Senegal
! Cette version de la discretisation strato est corrige au niveau
! du passage des sig aux ap, bp
! la version precedente donne un coude dans l'epaisseur des couches
! vers la tropopause
!==================================================================


     DO l = 1, llm
        x = 2*asin(1.) * (l - 0.5) / (llm + 1)
        dsig(l) =(dsigmin + 7. * SIN(x)**2) &
             *(0.5*(1.-tanh(1.*(x-asin(1.))/asin(1.))))**2
     ENDDO
     dsig = dsig / sum(dsig)
     sig0(llm+1) = 0.
     DO l = llm, 1, -1
        sig0(l) = sig0(l+1) + dsig(l)
     ENDDO
     sig=racinesig(sig0)

     bp(1)=1.
     bp(2:llm)=EXP(1.-1./sig(2: llm)**2)
     bp(llmp1)=0.

     ap(1)=0.
     ap(2:llm)=pa*(sig(2:llm)-bp(2:llm))
     ap(llm+1)=0.

  CASE("strato_custom0") 
!=======================================================
! Version Transitoire
    ! custumize strato distribution with specific alpha & beta values and function
    ! depending on llm (experimental and temporary)!
    SELECT CASE (llm)
      CASE(55)
        alpha=0.45
        beta=4.0
      CASE(63)
        alpha=0.45
        beta=5.0
      CASE(71)
        alpha=3.05
        beta=65.
      CASE(79)
        alpha=3.20
        ! alpha=2.05 ! FLOTT 79 (PLANTE)
        beta=70.
    END SELECT
    ! Or used values provided by user in def file:
    CALL getin("strato_alpha",alpha)
    CALL getin("strato_beta",beta)
    
    ! Build geometrical distribution
    scaleheight=7.
    zz(1)=0.
    IF (llm==55.OR.llm==63) THEN
      DO l=1,llm
        z=zz(l)/scaleheight
        zz(l+1)=zz(l)+0.03+z*1.5*(1.-TANH(z-0.5))+alpha*(1.+TANH(z-1.5))     &
                            +5.0*EXP((l-llm)/beta)
      ENDDO
    ELSEIF (llm==71) THEN !.OR.llm==79) THEN      ! FLOTT 79 (PLANTE)
      DO l=1,llm
        z=zz(l)
        zz(l+1)=zz(l)+0.02+0.88*TANH(z/2.5)+alpha*(1.+TANH((z-beta)/15.))
      ENDDO
    ELSEIF (llm==79) THEN
      DO l=1,llm
        z=zz(l)
        zz(l+1)=zz(l)+0.02+0.80*TANH(z/3.8)+alpha*(1+TANH((z-beta)/17.))     &
                            +0.03*TANH(z/.25)
      ENDDO
    ENDIF ! of IF (llm==55.OR.llm==63) ...
    

    ! Build sigma distribution
    sig0=EXP(-zz(:)/scaleheight)
    sig0(llm+1)=0.
!    sig=ridders(sig0)
    sig=racinesig(sig0)
    
    ! Compute ap() and bp()
    bp(1)=1.
    bp(2:llm)=EXP(1.-1./sig(2:llm)**2)
    bp(llm+1)=0.
    ap=pa*(sig-bp)

  CASE("strato_custom") 
!===================================================================
! David Cugnet, FranÃ§ois Lott, Lionel Guez, Ehouoarn Millour, Fredho
! 2014/05
! custumize strato distribution
! Al the parameter are given in km assuming a given scalehigh
    vert_scale=7.     ! scale hight
    vert_dzmin=0.02   ! width of first layer
    vert_dzlow=1.     ! dz in the low atmosphere
    vert_z0low=8.     ! height at which resolution recches dzlow
    vert_dzmid=3.     ! dz in the mid atmsophere
    vert_z0mid=70.    ! height at which resolution recches dzmid
    vert_h_mid=20.    ! width of the transition
    vert_dzhig=11.    ! dz in the high atmsophere
    vert_z0hig=80.    ! height at which resolution recches dz
    vert_h_hig=20.    ! width of the transition
!===================================================================

    call getin('vert_scale',vert_scale)
    call getin('vert_dzmin',vert_dzmin)
    call getin('vert_dzlow',vert_dzlow)
    call getin('vert_z0low',vert_z0low)
    CALL getin('vert_dzmid',vert_dzmid)
    CALL getin('vert_z0mid',vert_z0mid)
    call getin('vert_h_mid',vert_h_mid)
    call getin('vert_dzhig',vert_dzhig)
    call getin('vert_z0hig',vert_z0hig)
    call getin('vert_h_hig',vert_h_hig)

    scaleheight=vert_scale ! for consistency with further computations
    ! Build geometrical distribution
    zz(1)=0.
    DO l=1,llm
       z=zz(l)
       zz(l+1)=zz(l)+vert_dzmin+vert_dzlow*TANH(z/vert_z0low)+                &
&      (vert_dzmid-vert_dzlow)* &
&           (TANH((z-vert_z0mid)/vert_h_mid)-TANH((-vert_z0mid)/vert_h_mid)) &
&     +(vert_dzhig-vert_dzmid-vert_dzlow)*                                  &
&           (TANH((z-vert_z0hig)/vert_h_hig)-TANH((-vert_z0hig)/vert_h_hig))
    ENDDO


!===================================================================
! Comment added Fredho 2014/05/18, Saint-Louis, Senegal
! From approximate z to ap, bp, so that p=ap+bp*p0 and p/p0=exp(-z/H)
! sig0 is p/p0
! sig is an intermediate distribution introduce to estimate bp
! 1.   sig0=exp(-z/H)
! 2.   inversion of sig0=(1-pa/p0)*sig+(1-pa/p0)*exp(1-1/sig**2)
! 3.   bp=exp(1-1/sig**2)
! 4.   ap deduced from  the combination of 2 and 3 so that sig0=ap/p0+bp
!===================================================================

    sig0=EXP(-zz(:)/vert_scale)
    sig0(llm+1)=0.
    sig=racinesig(sig0)
    bp(1)=1.
    bp(2:llm)=EXP(1.-1./sig(2:llm)**2)
    bp(llm+1)=0.
    ap=pa*(sig-bp)

  case("read")
     ! Read "ap" and "bp". First line is skipped (title line). "ap"
     ! should be in Pa. First couple of values should correspond to
     ! the surface, that is : "bp" should be in descending order.
     call new_unit(unit)
     open(unit, file="hybrid.txt", status="old", action="read", &
          position="rewind")
     read(unit, fmt=*) ! skip title line
     do l = 1, llm + 1
        read(unit, fmt=*) ap(l), bp(l)
     end do
     close(unit)
     call assert(ap(1) == 0., ap(llm + 1) == 0., bp(1) == 1., &
          bp(llm + 1) == 0., "disvert: bad ap or bp values")
  case default
     call abort_gcm("disvert", 'Wrong value for "vert_sampling"', 1)
  END select

  DO l=1, llm
     nivsigs(l) = REAL(l)
  ENDDO

  DO l=1, llmp1
     nivsig(l)= REAL(l)
  ENDDO

  write(lunout, *)  trim(modname),': BP '
  write(lunout, *) bp
  write(lunout, *)  trim(modname),': AP '
  write(lunout, *) ap

  write(lunout, *) 'Niveaux de pressions approximatifs aux centres des'
  write(lunout, *)'couches calcules pour une pression de surface =', preff
  write(lunout, *) 'et altitudes equivalentes pour une hauteur d echelle de '
  write(lunout, *) scaleheight,' km'
  DO l = 1, llm
     dpres(l) = bp(l) - bp(l+1)
     presnivs(l) = 0.5 *( ap(l)+bp(l)*preff + ap(l+1)+bp(l+1)*preff )
     write(lunout, *)'PRESNIVS(', l, ')=', presnivs(l), ' Z ~ ', &
          log(preff/presnivs(l))*scaleheight &
          , ' DZ ~ ', scaleheight*log((ap(l)+bp(l)*preff)/ &
          max(ap(l+1)+bp(l+1)*preff, 1.e-10))
  ENDDO

  write(lunout, *) trim(modname),': PRESNIVS '
  write(lunout, *) presnivs

CONTAINS

!-------------------------------------------------------------------------------
!
FUNCTION ridders(sig) RESULT(sg)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Purpose: Search for s solving (Pa/Preff)*s+(1-Pa/Preff)*EXP(1-1./s**2)=sg
! Notes:   Uses Ridders' method, quite robust. Initial bracketing: 0<=sg<=1.
! Reference: Ridders, C. F. J. "A New Algorithm for Computing a Single Root of a
!       Real Continuous Function" IEEE Trans. Circuits Systems 26, 979-980, 1979
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)  :: sig(:)
  REAL              :: sg(SIZE(sig))
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: it, ns, maxit
  REAL :: c1, c2, x1, x2, x3, x4, f1, f2, f3, f4, s, xx, distrib
!-------------------------------------------------------------------------------
  ns=SIZE(sig); maxit=9999
  c1=Pa/Preff; c2=1.-c1
  DO l=1,ns
    xx=HUGE(1.)
    x1=0.0; f1=distrib(x1,c1,c2,sig(l))
    x2=1.0; f2=distrib(x2,c1,c2,sig(l))
    DO it=1,maxit
      x3=0.5*(x1+x2); f3=distrib(x3,c1,c2,sig(l))
      s=SQRT(f3**2-f1*f2);                 IF(s==0.) EXIT
      x4=x3+(x3-x1)*(SIGN(1.,f1-f2)*f3/s); IF(ABS(10.*LOG(x4-xx))<=1E-5) EXIT
      xx=x4; f4=distrib(x4,c1,c2,sig(l));  IF(f4==0.) EXIT
      IF(SIGN(f3,f4)/=f3) THEN;      x1=x3; f1=f3; x2=xx; f2=f4
      ELSE IF(SIGN(f1,f4)/=f1) THEN; x2=xx; f2=f4
      ELSE IF(SIGN(f2,f4)/=f2) THEN; x1=xx; f1=f4
      ELSE; CALL abort_gcm("ridders",'Algorithm failed (which is odd...')
      END IF
      IF(ABS(10.*LOG(ABS(x2-x1)))<=1E-5) EXIT       !--- ERROR ON SIG <= 0.01m            
    END DO
    IF(it==maxit+1) WRITE(lunout,'(a,i3)')'WARNING in ridder: failed to converg&
     &e for level ',l
    sg(l)=xx
  END DO
  sg(1)=1.; sg(ns)=0.

END FUNCTION ridders

FUNCTION racinesig(sig) RESULT(sg)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Fredho 2014/05/18
! Purpose: Search for s solving (Pa/Preff)*sg+(1-Pa/Preff)*EXP(1-1./sg**2)=s
! Notes:   Uses Newton Raphson search
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)  :: sig(:)
  REAL              :: sg(SIZE(sig))
!-------------------------------------------------------------------------------
! Local variables:
  INTEGER :: it, ns, maxit
  REAL :: c1, c2, x1, x2, x3, x4, f1, f2, f3, f4, s, xx, distrib
!-------------------------------------------------------------------------------
  ns=SIZE(sig); maxit=100
  c1=Pa/Preff; c2=1.-c1
  DO l=2,ns-1
    sg(l)=sig(l)
    DO it=1,maxit
       f1=exp(1-1./sg(l)**2)*(1.-c1)
       sg(l)=sg(l)-(c1*sg(l)+f1-sig(l))/(c1+2*f1*sg(l)**(-3))
    ENDDO
!   print*,'SSSSIG ',sig(l),sg(l),c1*sg(l)+exp(1-1./sg(l)**2)*(1.-c1)
  ENDDO
  sg(1)=1.; sg(ns)=0.

END FUNCTION racinesig




END SUBROUTINE disvert

!-------------------------------------------------------------------------------

FUNCTION distrib(x,c1,c2,x0) RESULT(res)
!
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN) :: x, c1, c2, x0
  REAL             :: res
!-------------------------------------------------------------------------------
  res=c1*x+c2*EXP(1-1/(x**2))-x0

END FUNCTION distrib



