!
! $Header$
!
      SUBROUTINE advect_new_p(ucov,vcov,teta,w,massebx,masseby,
     &                        du,dv,dteta)
      USE parallel_lmdz
      USE write_field_p
      USE comconst_mod, ONLY: daysec
      USE logic_mod, ONLY: conser
      IMPLICIT NONE
c=======================================================================
c
c   Auteurs:  P. Le Van , Fr. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c   *************************************************************
c   .... calcul des termes d'advection vertic.pour u,v,teta,q ...
c   *************************************************************
c        ces termes sont ajoutes a du,dv,dteta et dq .
c  Modif F.Forget 03/94 : on retire q de advect
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      include "dimensions.h"
      include "paramet.h"
      include "comgeom.h"

c   Arguments:
c   ----------

      REAL,INTENT(IN) :: vcov(ip1jm,llm)
      REAL,INTENT(IN) :: ucov(ip1jmp1,llm)
      REAL,INTENT(IN) :: teta(ip1jmp1,llm)
      REAL,INTENT(IN) :: massebx(ip1jmp1,llm)
      REAL,INTENT(IN) :: masseby(ip1jm,llm)
      REAL,INTENT(IN) :: w(ip1jmp1,llm)
      REAL,INTENT(INOUT) :: dv(ip1jm,llm)
      REAL,INTENT(INOUT) :: du(ip1jmp1,llm)
      REAL,INTENT(INOUT) :: dteta(ip1jmp1,llm)
c   Local:
c   ------

      REAL,SAVE :: dv1(ip1jm,llm),du1(ip1jmp1,llm),dteta1(ip1jmp1,llm)
      REAL,SAVE :: dv2(ip1jm,llm),du2(ip1jmp1,llm),dteta2(ip1jmp1,llm)
      REAL,SAVE :: uav(ip1jmp1,llm),vav(ip1jm,llm)
      REAL wsur2(ip1jmp1)
      REAL unsaire2(ip1jmp1), ge(ip1jmp1)
      REAL deuxjour, ww, gt, uu, vv

      INTEGER  ij,l,ijb,ije

      EXTERNAL  SSUM
      REAL      SSUM

c-----------------------------------------------------------------------
c   2. Calculs preliminaires:
c   -------------------------

      IF (conser)  THEN
         deuxjour = 2. * daysec

         DO ij   = 1, ip1jmp1
         unsaire2(ij) = unsaire(ij) * unsaire(ij)
         ENDDO
      END IF


c------------------  -yy ----------------------------------------------
c   .  Calcul de     u

c$OMP MASTER
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1

      DO ij=ijb,ije
        du2(ij,1)=0.
        du1(ij,llm)=0.
      ENDDO
      
      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1
      
      DO ij=ijb,ije
        dv2(ij,1)=0.
        dv1(ij,llm)=0.
      ENDDO
      
      ijb=ij_begin
      ije=ij_end

      DO ij=ijb,ije
        dteta2(ij,1)=0.
        dteta1(ij,llm)=0.
      ENDDO
c$OMP END MASTER

 
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)   
      DO  l=1,llm
         
         ijb=ij_begin
         ije=ij_end
         if (pole_nord) ijb=ijb+iip1
         if (pole_sud)  ije=ije-iip1
         
c         DO    ij     = iip2, ip1jmp1
c            uav(ij,l) = 0.25 * ( ucov(ij,l) + ucov(ij-iip1,l) )
c         ENDDO

c         DO    ij     = iip2, ip1jm
c            uav(ij,l) = uav(ij,l) + uav(ij+iip1,l)
c         ENDDO
         
         DO    ij     = ijb, ije
                  
           uav(ij,l)=0.25*(ucov(ij,l)+ucov(ij-iip1,l))
     .	             +0.25*(ucov(ij+iip1,l)+ucov(ij,l))
         ENDDO
         
         if (pole_nord) then
           DO      ij         = 1, iip1
              uav(ij      ,l) = 0.
           ENDDO
         endif
         
         if (pole_sud) then
           DO      ij         = 1, iip1
              uav(ip1jm+ij,l) = 0.
           ENDDO
         endif

      ENDDO
c$OMP END DO      
c      call write_field3d_p('uav',reshape(uav,(/iip1,jjp1,llm/)))
      
c------------------  -xx ----------------------------------------------
c   .  Calcul de     v
      
      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO  l=1,llm
         
         DO    ij   = ijb+1, ije
           vav(ij,l) = 0.25 * ( vcov(ij,l) + vcov(ij-1,l) )
         ENDDO
         
         DO    ij   = ijb,ije,iip1
          vav(ij,l) = vav(ij+iim,l)
         ENDDO
         
         
         DO    ij   = ijb, ije-1
          vav(ij,l) = vav(ij,l) + vav(ij+1,l)
         ENDDO
         
         DO    ij       = ijb, ije, iip1
          vav(ij+iim,l) = vav(ij,l)
         ENDDO
         
      ENDDO
c$OMP END DO
c       call write_field3d_p('vav',reshape(vav,(/iip1,jjm,llm/)))

c-----------------------------------------------------------------------
c$OMP BARRIER

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llmm1


c       ......   calcul de  - w/2.    au niveau  l+1   .......
        ijb=ij_begin
        ije=ij_end+iip1
        if (pole_sud)  ije=ij_end
      
        DO ij   = ijb, ije
          wsur2( ij ) = - 0.5 * w( ij,l+1 )
        ENDDO


c     .....................     calcul pour  du     ..................
      
        ijb=ij_begin
        ije=ij_end
        if (pole_nord) ijb=ijb+iip1
        if (pole_sud)  ije=ije-iip1
         
        DO ij = ijb ,ije-1
          ww        = wsur2 (  ij  )     + wsur2( ij+1 ) 
          uu        = 0.5 * ( ucov(ij,l) + ucov(ij,l+1) )
          du1(ij,l)  =  ww * ( uu - uav(ij, l ) )/massebx(ij, l )
          du2(ij,l+1)=  ww * ( uu - uav(ij,l+1) )/massebx(ij,l+1)
        ENDDO

c     .................    calcul pour   dv      .....................
        ijb=ij_begin
        ije=ij_end
        if (pole_sud)  ije=ij_end-iip1
      
        DO ij = ijb, ije
          ww        = wsur2( ij+iip1 )   + wsur2( ij )
          vv        = 0.5 * ( vcov(ij,l) + vcov(ij,l+1) )
          dv1(ij,l)  =  ww * (vv - vav(ij, l ) )/masseby(ij, l )
          dv2(ij,l+1)=  ww * (vv - vav(ij,l+1) )/masseby(ij,l+1)
        ENDDO

c

c     ............................................................
c     ...............    calcul pour   dh      ...................
c     ............................................................

c                       ---z
c       calcul de  - d( teta  * w )      qu'on ajoute a   dh
c                   ...............
        ijb=ij_begin
        ije=ij_end
        
        DO ij = ijb, ije
         ww            = wsur2(ij) * (teta(ij,l) + teta(ij,l+1) )
         dteta1(ij, l ) =   ww
         dteta2(ij,l+1) =   ww
        ENDDO

c ym ---> conser a voir plus tard

c      IF( conser)  THEN
c        
c        DO 17 ij = 1,ip1jmp1
c        ge(ij)   = wsur2(ij) * wsur2(ij) * unsaire2(ij)
c  17    CONTINUE
c        gt       = SSUM( ip1jmp1,ge,1 )
c        gtot(l)  = deuxjour * SQRT( gt/ip1jmp1 )
c      END IF

      ENDDO ! of DO l = 1, llmm1
c$OMP END DO

      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
        DO ij=ijb,ije-1
	  du(ij,l)=du(ij,l)+du2(ij,l)-du1(ij,l)
	ENDDO

        DO   ij   = ijb+iip1-1, ije, iip1
         du( ij, l  ) = du( ij -iim, l  )
        ENDDO 
      ENDDO
c$OMP END DO NOWAIT
      
      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
        DO ij=ijb,ije
	  dv(ij,l)=dv(ij,l)+dv2(ij,l)-dv1(ij,l)
	ENDDO
      ENDDO
c$OMP END DO NOWAIT      
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO l=1,llm
        DO ij=ijb,ije
	  dteta(ij,l)=dteta(ij,l)+dteta2(ij,l)-dteta1(ij,l)
	ENDDO
      ENDDO
c$OMP END DO NOWAIT      

      END
