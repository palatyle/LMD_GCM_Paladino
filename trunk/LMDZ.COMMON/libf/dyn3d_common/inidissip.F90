!
! $Id: inidissip.F90 1502 2011-03-21 16:07:54Z jghattas $
!
SUBROUTINE inidissip ( lstardis,nitergdiv,nitergrot,niterh  , &
     tetagdiv,tetagrot,tetatemp, vert_prof_dissip)
  !=======================================================================
  !   Initialization for horizontal (lateral) dissipation
  !  - in all cases, there is a multiplicative coefficient which increases
  !    the dissipation in the higher levels of the atmosphere, but there
  !    are different ways of seting the vertical profile of this coefficient
  !    (see code below).
  !  - the call to dissipation, every 'dissip_period' dynamical time step,
  !    can be imposed via 'dissip_period=...' in run.def (or via the
  !    'idissip=...' flag, but this flag should become obsolete, and is
  !    overridden by the 'dissip_period' flag). Note that setting 
  !    dissip_period=0 has the special meaning of requesting an "optimal"
  !    value for "dissip_period" (then taken as the largest possible value)
  !  - the three characteristic time scales (relative to the smallest
  !    longitudinal grid mesh size), which are privided in run.def, which 
  !    are used for the dissipations steps are:
  !     tetagdiv : time scale for the gradient of the divergence of velocities
  !     tetagrot : time scale for the curl of the curl of velocities
  !     tetatemp : time scale for the laplacian of the potential temperature
  !=======================================================================

  USE control_mod, only : dissip_period,iperiod,planet_type
  USE comvert_mod, ONLY: preff,presnivs,scaleheight,pseudoalt
  USE comconst_mod, ONLY: dtvr,dtdiss,rad,pi,dissip_zref,dissip_deltaz,		&
		dissip_factz,dissip_fac_mid,dissip_fac_up,dissip_pupstart,	&
		dissip_hdelta	
  USE logic_mod, ONLY: ok_strato

  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comdissipn.h"
  include "iniprint.h"

  LOGICAL,INTENT(in) :: lstardis
  INTEGER,INTENT(in) :: nitergdiv,nitergrot,niterh
  REAL,INTENT(in) :: tetagdiv,tetagrot,tetatemp

  integer, INTENT(in):: vert_prof_dissip ! Vertical profile of horizontal dissipation
  ! For the Earth model:
  ! Allowed values:
  ! 0: rational fraction, function of pressure
  ! 1: tanh of altitude
  ! For planets:
  ! 0: use fac_mid (read from run.def)
  ! 1: use fac_mid, fac_up, startalt, delta (hard coded in inidissip.F90)
! Local variables:
  REAL fact,zvert(llm),zz
  REAL zh(ip1jmp1),zu(ip1jmp1), gx(ip1jmp1), divgra(ip1jmp1)
  real zv(ip1jm), gy(ip1jm), deltap(ip1jmp1,llm)
  REAL ullm,vllm,umin,vmin,zhmin,zhmax
  REAL zllm

  INTEGER l,ij,idum,ii
  REAL tetamin
  REAL pseudoz
  REAL Pup
  character (len=80) :: abort_message

  REAL ran1
  logical,save :: firstcall=.true.
  real,save :: fac_mid,fac_up,startalt,delta,middle

  if (firstcall) then
    firstcall=.false.
    if ((planet_type.ne."earth").and.(vert_prof_dissip.eq.1)) then
      ! initialize values for dissipation variation along the vertical (Mars)
      fac_mid=3 ! coefficient for lower/middle atmosphere
      fac_up=30 ! coefficient for upper atmosphere
      startalt=70. ! altitude (in km) for the transition from middle to upper atm.
      delta=30.! Size (in km) of the transition region between middle/upper atm.
      
      if (ok_strato) then !overwrite defaults above with values from def file
        fac_mid=dissip_fac_mid
        fac_up=dissip_fac_up
        delta=dissip_deltaz
      endif
      
      middle=startalt+delta/2
      write(lunout,*)"inidissip: multiplicative factors in altitude:", &
        " fac_mid=",fac_mid," fac_up=",fac_up
      write(lunout,*)" transition mid/up : startalt (km) =",startalt, &
        " delta (km) =",delta
    endif
  endif !of if firstcall

  !-----------------------------------------------------------------------
  !
  !   calcul des valeurs propres des operateurs par methode iterrative:
  !   -----------------------------------------------------------------

  crot     = -1.
  cdivu    = -1.
  cdivh    = -1.

  !   calcul de la valeur propre de divgrad:
  !   --------------------------------------
  idum = 0
  DO l = 1, llm
     DO ij = 1, ip1jmp1
        deltap(ij,l) = 1.
     ENDDO
  ENDDO

  idum  = -1
  zh(1) = RAN1(idum)-.5
  idum  = 0
  DO ij = 2, ip1jmp1
     zh(ij) = RAN1(idum) -.5
  ENDDO

  CALL filtreg (zh,jjp1,1,2,1,.TRUE.,1)

  CALL minmax(iip1*jjp1,zh,zhmin,zhmax )

  IF ( zhmin .GE. zhmax  )     THEN
     write(lunout,*)'  Inidissip  zh min max  ',zhmin,zhmax
     abort_message='probleme generateur alleatoire dans inidissip'
     call abort_gcm('inidissip',abort_message,1)
  ENDIF

  zllm = ABS( zhmax )
  DO l = 1,50
     IF(lstardis) THEN
        CALL divgrad2(1,zh,deltap,niterh,divgra)
     ELSE
        CALL divgrad (1,zh,niterh,divgra)
     ENDIF

     zllm  = ABS(maxval(divgra))
     zh = divgra / zllm
  ENDDO

  IF(lstardis) THEN
     cdivh = 1./ zllm
  ELSE
     cdivh = zllm ** ( -1./niterh )
  ENDIF

  !   calcul des valeurs propres de gradiv (ii =1) et  nxgrarot(ii=2)
  !   -----------------------------------------------------------------
  write(lunout,*)'inidissip: calcul des valeurs propres'

  DO    ii = 1, 2
     !
     DO ij = 1, ip1jmp1
        zu(ij)  = RAN1(idum) -.5
     ENDDO
     CALL filtreg (zu,jjp1,1,2,1,.TRUE.,1)
     DO ij = 1, ip1jm
        zv(ij) = RAN1(idum) -.5
     ENDDO
     CALL filtreg (zv,jjm,1,2,1,.FALSE.,1)

     CALL minmax(iip1*jjp1,zu,umin,ullm )
     CALL minmax(iip1*jjm, zv,vmin,vllm )

     ullm = ABS ( ullm )
     vllm = ABS ( vllm )

     DO    l = 1, 50
        IF(ii.EQ.1) THEN
           !cccc             CALL covcont( 1,zu,zv,zu,zv )
           IF(lstardis) THEN
              CALL gradiv2( 1,zu,zv,nitergdiv,gx,gy )
           ELSE
              CALL gradiv ( 1,zu,zv,nitergdiv,gx,gy )
           ENDIF
        ELSE
           IF(lstardis) THEN
              CALL nxgraro2( 1,zu,zv,nitergrot,gx,gy )
           ELSE
              CALL nxgrarot( 1,zu,zv,nitergrot,gx,gy )
           ENDIF
        ENDIF

        zllm = max(abs(maxval(gx)), abs(maxval(gy)))
        zu = gx / zllm
        zv = gy / zllm
     end DO

     IF ( ii.EQ.1 ) THEN
        IF(lstardis) THEN
           cdivu  = 1./zllm
        ELSE
           cdivu  = zllm **( -1./nitergdiv )
        ENDIF
     ELSE
        IF(lstardis) THEN
           crot   = 1./ zllm
        ELSE
           crot   = zllm **( -1./nitergrot )
        ENDIF
     ENDIF

  end DO

  !   petit test pour les operateurs non star:
  !   ----------------------------------------

  !     IF(.NOT.lstardis) THEN
  fact    = rad*24./REAL(jjm)
  fact    = fact*fact
  write(lunout,*)'inidissip: coef u ', fact/cdivu, 1./cdivu
  write(lunout,*)'inidissip: coef r ', fact/crot , 1./crot
  write(lunout,*)'inidissip: coef h ', fact/cdivh, 1./cdivh
  !     ENDIF

  !-----------------------------------------------------------------------
  !   variation verticale du coefficient de dissipation:
  !   --------------------------------------------------
  
  if (planet_type.eq."earth") then

   if (vert_prof_dissip == 1) then
     do l=1,llm
        pseudoz=8.*log(preff/presnivs(l))
        zvert(l)=1+ &
             (tanh((pseudoz-dissip_zref)/dissip_deltaz)+1.)/2. &
             *(dissip_factz-1.)
     enddo
   else
     DO l=1,llm
        zvert(l)=1.
     ENDDO
     fact=2.
     DO l = 1, llm
        zz      = 1. - preff/presnivs(l)
        zvert(l)= fact -( fact-1.)/( 1.+zz*zz )
     ENDDO
   endif ! of if (vert_prof_dissip == 1)

  else ! other planets
  
   if (vert_prof_dissip == 0) then
! First step: going from 1 to dissip_fac_mid (in gcm.def)
!============
    DO l=1,llm
     zz      = 1. - preff/presnivs(l)
     zvert(l)= dissip_fac_mid -( dissip_fac_mid-1.)/( 1.+zz*zz )
    ENDDO

    write(lunout,*) 'Dissipation : '
    write(lunout,*) 'Multiplication de la dissipation en altitude :'
    write(lunout,*) '  dissip_fac_mid =', dissip_fac_mid

! Second step if ok_strato:  from dissip_fac_mid to dissip_fac_up (in gcm.def)
!==========================
! Utilisation de la fonction tangente hyperbolique pour augmenter
! arbitrairement la dissipation et donc la stabilite du modele a 
! partir d'une certaine altitude.

!   Le facteur multiplicatif de basse atmosphere etant deja pris 
!   en compte, il faut diviser le facteur multiplicatif de haute 
!   atmosphere par celui-ci.

    if (ok_strato) then

     Pup = dissip_pupstart*exp(-0.5*dissip_deltaz/dissip_hdelta)
     do l=1,llm
      zvert(l)= zvert(l)*(1.0+( (dissip_fac_up/dissip_fac_mid-1) &
                *(1-(0.5*(1+tanh(-6./dissip_deltaz*              &
               (-dissip_hdelta*log(presnivs(l)/Pup))  ))))  ))
     enddo 

     write(*,*) '  dissip_fac_up =', dissip_fac_up
     write(*,*) 'Transition mid /up:  Pupstart,delta =',           &
                   dissip_pupstart,'Pa', dissip_deltaz , '(km)'

    endif ! of if (ok_strato)
   elseif (vert_prof_dissip==1) then
    DO l=1,llm
     zz      = 1. - preff/presnivs(l)
!     zvert(l)= dissip_fac_mid -( dissip_fac_mid-1.)/( 1.+zz*zz )
     zvert(l)= fac_mid -( fac_mid-1.)/( 1.+zz*zz )
     
     zvert(l)= zvert(l)*(1.0+((fac_up/fac_mid-1)*    &
                (1-(0.5*(1+tanh(-6./                 &
                delta*(scaleheight*(-log(presnivs(l)/preff))-middle))))) &
                ))
    ENDDO
    write(lunout,*) "inidissip: vert_prof_disip=1, scaleheight=",scaleheight
    write(lunout,*) "           fac_mid=",fac_mid,", fac_up=",fac_up
    
   else
     write(lunout,*) 'wrong value for vert_prof_dissip:',vert_prof_dissip
     abort_message='wrong value for vert_prof_dissip'
     call abort_gcm('inidissip',abort_message,1)
   endif ! of (vert_prof_dissip == 0)
  endif ! of if (planet_type.eq."earth")


  write(lunout,*)'inidissip: Time constants for lateral dissipation'

  tetamin =  1.e+6

  DO l=1,llm
     tetaudiv(l)   = zvert(l)/tetagdiv
     tetaurot(l)   = zvert(l)/tetagrot
     tetah(l)      = zvert(l)/tetatemp

     IF( tetamin.GT. (1./tetaudiv(l)) ) tetamin = 1./ tetaudiv(l)
     IF( tetamin.GT. (1./tetaurot(l)) ) tetamin = 1./ tetaurot(l)
     IF( tetamin.GT. (1./   tetah(l)) ) tetamin = 1./    tetah(l)
  ENDDO

  ! If dissip_period=0 calculate value for dissipation period, else keep value read from gcm.def
  IF (dissip_period == 0) THEN
     dissip_period = INT( tetamin/( 2.*dtvr*iperiod) ) * iperiod
     write(lunout,*)'inidissip: tetamin dtvr iperiod dissip_period(intermed) ',tetamin,dtvr,iperiod,dissip_period
     dissip_period = MAX(iperiod,dissip_period)
  END IF

  dtdiss  = dissip_period * dtvr
  write(lunout,*)'inidissip: dissip_period=',dissip_period,' dtdiss=',dtdiss,' dtvr=',dtvr

  write(lunout,*)'pseudoZ(km)  zvert    dt(tetagdiv)   dt(tetagrot)   dt(divgrad)'
  DO l = 1,llm
     write(lunout,'(f6.1,x,4(1pe14.7,x))') &
     pseudoalt(l),zvert(l),dtdiss*tetaudiv(l),dtdiss*tetaurot(l),dtdiss*tetah(l)
     ! test if disipation is not too strong (for Explicit Euler time marching)
     if (dtdiss*tetaudiv(l).gt.1.9) then
       write(lunout,*)"STOP : lateral dissipation is too intense and will"
       write(lunout,*)"       generate instabilities in the model !"
       write(lunout,*)" You must increase tetagdiv (or increase dissip_period"
       write(lunout,*)"                             or increase day_stap)"
     endif
     if (dtdiss*tetaurot(l).gt.1.9) then
       write(lunout,*)"STOP : lateral dissipation is too intense and will"
       write(lunout,*)"       generate instabilities in the model !"
       write(lunout,*)" You must increase tetaurot (or increase dissip_period"
       write(lunout,*)"                             or increase day_stap)"
     endif
     if (dtdiss*tetah(l).gt.1.9) then
       write(lunout,*)"STOP : lateral dissipation is too intense and will"
       write(lunout,*)"       generate instabilities in the model !"
       write(lunout,*)" You must increase tetah (or increase dissip_period"
       write(lunout,*)"                          or increase day_stap)"
     endif
  ENDDO ! of DO l=1,llm

END SUBROUTINE inidissip
