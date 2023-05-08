!
! $Id: integrd_p.F 1616 2012-02-17 11:59:00Z emillour $
!
      SUBROUTINE integrd_p
     $  (  nq,vcovm1,ucovm1,tetam1,psm1,massem1,
     $     dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps0,masse,phis) !,finvmaold)
      USE parallel_lmdz, ONLY: ij_begin, ij_end, pole_nord, pole_sud,
     &                         omp_chunk
      USE control_mod, only : planet_type,force_conserv_tracer
      USE comvert_mod, ONLY: ap,bp
      USE comconst_mod, ONLY: pi
      USE logic_mod, ONLY: leapf
      USE temps_mod, ONLY: dt
      IMPLICIT NONE


c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   objet:
c   ------
c
c   Incrementation des tendances dynamiques
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "iniprint.h"

c   Arguments:
c   ----------

      integer,intent(in) :: nq ! number of tracers to handle in this routine
      real,intent(inout) :: vcov(ip1jm,llm) ! covariant meridional wind
      real,intent(inout) :: ucov(ip1jmp1,llm) ! covariant zonal wind
      real,intent(inout) :: teta(ip1jmp1,llm) ! potential temperature
      real,intent(inout) :: q(ip1jmp1,llm,nq) ! advected tracers
      real,intent(inout) :: ps0(ip1jmp1) ! surface pressure
      real,intent(inout) :: masse(ip1jmp1,llm) ! atmospheric mass
      real,intent(in) :: phis(ip1jmp1) ! ground geopotential !!! unused
      ! values at previous time step
      real,intent(inout) :: vcovm1(ip1jm,llm)
      real,intent(inout) :: ucovm1(ip1jmp1,llm)
      real,intent(inout) :: tetam1(ip1jmp1,llm)
      real,intent(inout) :: psm1(ip1jmp1)
      real,intent(inout) :: massem1(ip1jmp1,llm)
      ! the tendencies to add
      real,intent(in) :: dv(ip1jm,llm)
      real,intent(in) :: du(ip1jmp1,llm)
      real,intent(in) :: dteta(ip1jmp1,llm)
      real,intent(in) :: dp(ip1jmp1)
      real,intent(in) :: dq(ip1jmp1,llm,nq) !!! unused
!      real,intent(out) :: finvmaold(ip1jmp1,llm) !!! unused

c   Local:
c   ------

      REAL vscr( ip1jm ),uscr( ip1jmp1 ),hscr( ip1jmp1 ),pscr(ip1jmp1)
      REAL massescr( ip1jmp1,llm )
      REAL :: massratio(ip1jmp1,llm)
!      REAL finvmasse(ip1jmp1,llm)
      REAL,SAVE :: p(ip1jmp1,llmp1)
      REAL tpn,tps,tppn(iim),tpps(iim)
      REAL qpn,qps,qppn(iim),qpps(iim)
      REAL,SAVE :: deltap( ip1jmp1,llm )

      INTEGER  l,ij,iq,i,j

      REAL SSUM
      EXTERNAL SSUM
      INTEGER ijb,ije,jjb,jje
      REAL,SAVE :: ps(ip1jmp1)=0
      LOGICAL :: checksum
      INTEGER :: stop_it
c-----------------------------------------------------------------------
c$OMP BARRIER     
      if (pole_nord) THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO  l = 1,llm
          DO  ij = 1,iip1
           ucov(    ij    , l) = 0.
           uscr(     ij      ) = 0.
           ENDDO
        ENDDO
c$OMP END DO NOWAIT        
      ENDIF

      if (pole_sud) THEN
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO  l = 1,llm
          DO  ij = 1,iip1
           ucov( ij +ip1jm, l) = 0.
           uscr( ij +ip1jm   ) = 0.
          ENDDO
        ENDDO
c$OMP END DO NOWAIT      
      ENDIF

c    ............    integration  de       ps         ..............

c      CALL SCOPY(ip1jmp1*llm, masse, 1, massescr, 1)

      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO  l = 1,llm
        massescr(ijb:ije,l)=masse(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT 

c$OMP DO SCHEDULE(STATIC)
      DO 2 ij = ijb,ije
       pscr (ij)    = ps0(ij)
       ps (ij)      = psm1(ij) + dt * dp(ij)
   2  CONTINUE
c$OMP END DO  
c$OMP BARRIER
c --> ici synchro OPENMP pour ps
       
      checksum=.TRUE.
      stop_it=0

c$OMP DO SCHEDULE(STATIC)
      DO ij = ijb,ije
         IF( ps(ij).LT.0. ) THEN
           IF (checksum) stop_it=ij
           checksum=.FALSE.
         ENDIF
       ENDDO
c$OMP END DO NOWAIT 
       
        IF( .NOT. checksum ) THEN
         write(lunout,*) "integrd: negative surface pressure ",
     &                                                ps(stop_it)
         write(lunout,*) " at node ij =", stop_it
         ! since ij=j+(i-1)*jjp1 , we have
         j=modulo(stop_it,jjp1)
         i=1+(stop_it-j)/jjp1
         write(lunout,*) " lon = ",rlonv(i)*180./pi, " deg",
     &                   " lat = ",rlatu(j)*180./pi, " deg"
         write(lunout,*) " psm1(ij)=",psm1(stop_it)," dt=",dt,
     &                   " dp(ij)=",dp(stop_it)
         call abort_gcm("integrd_p", "negative surface pressure", 1)
        ENDIF

c
C$OMP MASTER
      if (pole_nord) THEN
      
        DO  ij    = 1, iim
         tppn(ij) = aire(   ij   ) * ps(  ij    )
        ENDDO
         tpn      = SSUM(iim,tppn,1)/apoln
        DO ij   = 1, iip1
         ps(   ij   )  = tpn
        ENDDO
      
      ENDIF
      
      if (pole_sud) THEN
      
        DO  ij    = 1, iim
         tpps(ij) = aire(ij+ip1jm) * ps(ij+ip1jm)
        ENDDO
         tps      = SSUM(iim,tpps,1)/apols
        DO ij   = 1, iip1
         ps(ij+ip1jm)  = tps
        ENDDO
      
      ENDIF
c$OMP END MASTER
c$OMP BARRIER
c
c  ... Calcul  de la nouvelle masse d'air au dernier temps integre t+1 ...
c

      CALL pression_p ( ip1jmp1, ap, bp, ps, p )
c$OMP BARRIER
      CALL massdair_p (     p  , masse         )

! Ehouarn : we don't use/need finvmaold and finvmasse,
!           so might as well not compute them
!c      CALL   SCOPY( ijp1llm  , masse, 1, finvmasse,  1      )
!      ijb=ij_begin
!      ije=ij_end
!      
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
!      DO  l = 1,llm
!        finvmasse(ijb:ije,l)=masse(ijb:ije,l)
!      ENDDO
!c$OMP END DO NOWAIT
!
!      jjb=jj_begin
!      jje=jj_end
!      CALL filtreg_p( finvmasse,jjb,jje, jjp1, llm, -2, 2, .TRUE., 1  )
c

c    ............   integration  de  ucov, vcov,  h     ..............

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO 10 l = 1,llm
      
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO 4 ij = ijb,ije
      uscr( ij )   =  ucov( ij,l )
      ucov( ij,l ) = ucovm1( ij,l ) + dt * du( ij,l )
   4  CONTINUE

      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1
      
      DO 5 ij = ijb,ije
      vscr( ij )   =  vcov( ij,l )
      vcov( ij,l ) = vcovm1( ij,l ) + dt * dv( ij,l )
   5  CONTINUE
      
      ijb=ij_begin
      ije=ij_end
      
      DO 6 ij = ijb,ije
      hscr( ij )    =  teta(ij,l)
      teta ( ij,l ) = tetam1(ij,l) *  massem1(ij,l) / masse(ij,l) 
     $                + dt * dteta(ij,l) / masse(ij,l)
   6  CONTINUE

c   ....  Calcul de la valeur moyenne, unique  aux poles pour  teta    ......
c
c
      IF (pole_nord) THEN
       
        DO  ij   = 1, iim
          tppn(ij) = aire(   ij   ) * teta(  ij    ,l)
        ENDDO
          tpn      = SSUM(iim,tppn,1)/apoln

        DO ij   = 1, iip1
          teta(   ij   ,l)  = tpn
        ENDDO
      
      ENDIF
      
      IF (pole_sud) THEN
       
        DO  ij   = 1, iim
          tpps(ij) = aire(ij+ip1jm) * teta(ij+ip1jm,l)
        ENDDO
          tps      = SSUM(iim,tpps,1)/apols

        DO ij   = 1, iip1
          teta(ij+ip1jm,l)  = tps
        ENDDO
      
      ENDIF
c

      IF(leapf)  THEN
c         CALL SCOPY ( ip1jmp1, uscr(1), 1, ucovm1(1, l), 1 )
c         CALL SCOPY (   ip1jm, vscr(1), 1, vcovm1(1, l), 1 )
c         CALL SCOPY ( ip1jmp1, hscr(1), 1, tetam1(1, l), 1 )
        ijb=ij_begin
        ije=ij_end
        ucovm1(ijb:ije,l)=uscr(ijb:ije)
        tetam1(ijb:ije,l)=hscr(ijb:ije)
        if (pole_sud) ije=ij_end-iip1
        vcovm1(ijb:ije,l)=vscr(ijb:ije)
      
      END IF

  10  CONTINUE
c$OMP END DO NOWAIT

c
c   .......  integration de   q   ......
c
      ijb=ij_begin
      ije=ij_end

	 if (planet_type.eq."earth") then
! Earth-specific treatment of first 2 tracers (water)
c$OMP BARRIER
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
          DO l = 1, llm
           DO ij = ijb, ije
            deltap(ij,l) =  p(ij,l) - p(ij,l+1) 
           ENDDO
          ENDDO
c$OMP END DO NOWAIT
c$OMP BARRIER

          CALL qminimum_p( q, nq, deltap )
c
c    .....  Calcul de la valeur moyenne, unique  aux poles pour  q .....
c
c$OMP BARRIER
      IF (pole_nord) THEN 
      
        DO iq = 1, nq
        
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
  
             DO ij = 1, iim
               qppn(ij) = aire(   ij   ) * q(   ij   ,l,iq)
             ENDDO
               qpn  =  SSUM(iim,qppn,1)/apoln
      
             DO ij = 1, iip1
               q(   ij   ,l,iq)  = qpn
             ENDDO    
  
          ENDDO
c$OMP END DO NOWAIT

        ENDDO
      
      ENDIF

      IF (pole_sud) THEN 
      
        DO iq = 1, nq

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
  
             DO ij = 1, iim
               qpps(ij) = aire(ij+ip1jm) * q(ij+ip1jm,l,iq)
             ENDDO
               qps  =  SSUM(iim,qpps,1)/apols 
  
             DO ij = 1, iip1
               q(ij+ip1jm,l,iq)  = qps
             ENDDO    
  
          ENDDO
c$OMP END DO NOWAIT

        ENDDO
      
      ENDIF
      
! Ehouarn: forget about finvmaold
!c         CALL  SCOPY( ijp1llm , finvmasse, 1, finvmaold, 1 )
!
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
!      DO l = 1, llm      
!        finvmaold(ijb:ije,l)=finvmasse(ijb:ije,l)        
!      ENDDO
!c$OMP END DO NOWAIT

      endif ! of if (planet_type.eq."earth")


      if (force_conserv_tracer) then
        ! Ehouarn: try to keep total amont of tracers fixed 
        ! by acounting for mass change in each cell
        massratio(ijb:ije,1:llm)=massescr(ijb:ije,1:llm)
     &                             /masse(ijb:ije,1:llm)
        do iq=1,nq
        q(ijb:ije,1:llm,iq)=q(ijb:ije,1:llm,iq)
     &                        *massratio(ijb:ije,1:llm)
        enddo
      endif ! of if (force_conserv_tracer)
c
c
c     .....   FIN  de l'integration  de   q    .......

15    continue

c$OMP DO SCHEDULE(STATIC)
      DO ij=ijb,ije  
        ps0(ij)=ps(ij)
      ENDDO
c$OMP END DO NOWAIT 

c    .................................................................


      IF( leapf )  THEN
c       CALL SCOPY (    ip1jmp1 ,  pscr   , 1,   psm1  , 1 )
c       CALL SCOPY ( ip1jmp1*llm, massescr, 1,  massem1, 1 )
c$OMP DO SCHEDULE(STATIC)
      DO ij=ijb,ije  
        psm1(ij)=pscr(ij)
      ENDDO
c$OMP END DO NOWAIT 

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
            massem1(ijb:ije,l)=massescr(ijb:ije,l)
	  ENDDO
c$OMP END DO NOWAIT	  
      END IF
c$OMP BARRIER
      RETURN
      END
