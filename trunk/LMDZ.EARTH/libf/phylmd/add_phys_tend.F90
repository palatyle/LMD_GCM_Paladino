!
! $Id: add_phys_tend.F90 1163 2009-05-20 14:11:21Z fairhead $
!
SUBROUTINE add_phys_tend (zdu,zdv,zdt,zdq,zdql,text)
!======================================================================
! Ajoute les tendances des variables physiques aux variables 
! d'etat de la dynamique t_seri, q_seri ...
! On en profite pour faire des tests sur les tendances en question.
!======================================================================


!======================================================================
! Declarations
!======================================================================

use dimphy
use phys_local_var_mod
use phys_state_var_mod
IMPLICIT none
#include "iniprint.h"

! Arguments :
!------------
REAL zdu(klon,klev),zdv(klon,klev)
REAL zdt(klon,klev),zdq(klon,klev),zdql(klon,klev)
CHARACTER*(*) text

! Local :
!--------
REAL zt,zq

INTEGER i, k,j
INTEGER jadrs(klon*klev), jbad
INTEGER jqadrs(klon*klev), jqbad
INTEGER kadrs(klon*klev)
INTEGER kqadrs(klon*klev)

integer debug_level
logical, save :: first=.true.
!$OMP THREADPRIVATE(first)
INTEGER, SAVE :: itap
!$OMP THREADPRIVATE(itap)
!======================================================================
! Initialisations

debug_level=10
     if (first) then
        itap=0
        first=.false.
     endif
! Incrementer le compteur de la physique
     itap   = itap + 1
!======================================================================
! Ajout des tendances sur le vent et l'eau liquide
!======================================================================

     u_seri(:,:)=u_seri(:,:)+zdu(:,:)
     v_seri(:,:)=v_seri(:,:)+zdv(:,:)
     ql_seri(:,:)=ql_seri(:,:)+zdql(:,:)

!======================================================================
! On ajoute les tendances de la temperature et de la vapeur d'eau
! en verifiant que ca ne part pas dans les choux
!======================================================================

      jbad=0
      jqbad=0
      DO k = 1, klev
         DO i = 1, klon
            zt=t_seri(i,k)+zdt(i,k)
            zq=q_seri(i,k)+zdq(i,k)
            IF ( zt>370. .or. zt<130. .or. abs(zdt(i,k))>50. ) then
            jbad = jbad + 1
            jadrs(jbad) = i
            kadrs(jbad) = k
            ENDIF
            IF ( zq<0. .or. zq>0.1 .or. abs(zdq(i,k))>1.e-2 ) then
            jqbad = jqbad + 1
            jqadrs(jqbad) = i
            kqadrs(jqbad) = k
            ENDIF
            t_seri(i,k)=zt
            q_seri(i,k)=zq
         ENDDO
      ENDDO

!=====================================================================================
! Impression et stop en cas de probleme important
!=====================================================================================

IF (jbad .GT. 0) THEN
      DO j = 1, jbad
         i=jadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'PLANTAGE POUR LE POINT i rlon rlat =',i,rlon(i),rlat(i),text
          print*,'l    T     dT       Q     dQ    '
          DO k = 1, klev
             write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
          ENDDO
          call print_debug_phys(i,debug_level,text)
         endif
      ENDDO
ENDIF
!
!=====================================================================================
! Impression, warning et correction en cas de probleme moins important
!=====================================================================================
IF (jqbad .GT. 0) THEN
      DO j = 1, jqbad
         i=jqadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'WARNING  : EAU POUR LE POINT i rlon rlat =',i,rlon(i),rlat(i),text
          print*,'l    T     dT       Q     dQ    '
         endif
         DO k = 1, klev
           zq=q_seri(i,k)+zdq(i,k)
           if (zq.lt.1.e-15) then
              if (q_seri(i,k).lt.1.e-15) then
               if(prt_level.ge.debug_level) THEN
                print*,' cas q_seri<1.e-15 i k q_seri zq zdq :',i,k,q_seri(i,k),zq,zdq(i,k)
               endif
               q_seri(i,k)=1.e-15
               zdq(i,k)=(1.e-15-q_seri(i,k))
              endif
           endif
!           zq=q_seri(i,k)+zdq(i,k)
!           if (zq.lt.1.e-15) then
!              zdq(i,k)=(1.e-15-q_seri(i,k))
!           endif
         ENDDO
      ENDDO
ENDIF
!

!IM ajout memes tests pour reverifier les jbad, jqbad beg
      jbad=0
      jqbad=0
      DO k = 1, klev
         DO i = 1, klon
            IF ( t_seri(i,k)>370. .or. t_seri(i,k)<130. .or. abs(zdt(i,k))>50. ) then
            jbad = jbad + 1
            jadrs(jbad) = i
!            if(prt_level.ge.debug_level) THEN
!             print*,'cas2 i k t_seri zdt',i,k,t_seri(i,k),zdt(i,k)
!            endif
            ENDIF
            IF ( q_seri(i,k)<0. .or. q_seri(i,k)>0.1 .or. abs(zdq(i,k))>1.e-2 ) then
            jqbad = jqbad + 1
            jqadrs(jqbad) = i
            kqadrs(jqbad) = k
!            if(prt_level.ge.debug_level) THEN
!             print*,'cas2 i k q_seri zdq',i,k,q_seri(i,k),zdq(i,k)
!            endif
            ENDIF
         ENDDO
      ENDDO
IF (jbad .GT. 0) THEN
      DO j = 1, jbad
         i=jadrs(j)
         k=kadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'PLANTAGE2 POUR LE POINT i itap rlon rlat txt jbad zdt t',i,itap,rlon(i),rlat(i),text,jbad, &
       &        zdt(i,k),t_seri(i,k)-zdt(i,k)
!!!       if(prt_level.ge.10.and.itap.GE.229.and.i.EQ.3027) THEN 
          print*,'l    T     dT       Q     dQ    '
          DO k = 1, klev
             write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
          ENDDO
          call print_debug_phys(i,debug_level,text)
         endif
      ENDDO
ENDIF 
!
IF (jqbad .GT. 0) THEN
      DO j = 1, jqbad
         i=jqadrs(j)
         k=kqadrs(j)
         if(prt_level.ge.debug_level) THEN
          print*,'WARNING  : EAU2 POUR LE POINT i itap rlon rlat txt jqbad zdq q zdql ql',i,itap,rlon(i),rlat(i),text,jqbad,&
       &        zdq(i,k), q_seri(i,k)-zdq(i,k), zdql(i,k), ql_seri(i,k)-zdql(i,k)
!!!       if(prt_level.ge.10.and.itap.GE.229.and.i.EQ.3027) THEN 
          print*,'l    T     dT       Q     dQ    '
          DO k = 1, klev
            write(*,'(i3,2f14.4,2e14.2)') k,t_seri(i,k),zdt(i,k),q_seri(i,k),zdq(i,k)
          ENDDO
          call print_debug_phys(i,debug_level,text)
         endif
      ENDDO
ENDIF

      CALL hgardfou(t_seri,ftsol,text)
      RETURN
      END
