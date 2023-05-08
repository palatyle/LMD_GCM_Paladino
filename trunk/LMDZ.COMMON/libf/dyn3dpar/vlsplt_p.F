      SUBROUTINE vlsplt_p(q,pente_max,masse,w,pbaru,pbarv,pdt)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c   pente_max facteur de limitation des pentes: 2 en general
c                                               0 pour un schema amont
c   pbaru,pbarv,w flux de masse en u ,v ,w
c   pdt pas de temps
c
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE mod_hallo
      USE Vampir
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"

c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
c      REAL masse(iip1,jjp1,llm),pente_max
      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
      REAL q(ip1jmp1,llm)
c      REAL q(iip1,jjp1,llm)
      REAL w(ip1jmp1,llm),pdt
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
      INTEGER ijlqmin,iqmin,jqmin,lqmin
c
      REAL zm(ip1jmp1,llm),newmasse
      REAL mu(ip1jmp1,llm)
      REAL mv(ip1jm,llm)
      REAL mw(ip1jmp1,llm+1)
      REAL zq(ip1jmp1,llm),zz
      REAL dqx(ip1jmp1,llm),dqy(ip1jmp1,llm),dqz(ip1jmp1,llm)
      REAL second,temps0,temps1,temps2,temps3
      REAL ztemps1,ztemps2,ztemps3
      REAL zzpbar, zzw
      LOGICAL testcpu
      SAVE testcpu
      SAVE temps1,temps2,temps3
      INTEGER iminn,imaxx

      REAL qmin,qmax
      DATA qmin,qmax/0.,1.e33/
      DATA testcpu/.false./
      DATA temps1,temps2,temps3/0.,0.,0./
      INTEGER ijb,ije
      type(request) :: MyRequest1
      type(request) :: MyRequest2

      call SetTag(MyRequest1,100)
      call SetTag(MyRequest2,101)
      
      zzpbar = 0.5 * pdt
      zzw    = pdt
      
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1
      
      DO l=1,llm
        DO ij = ijb,ije
            mu(ij,l)=pbaru(ij,l) * zzpbar
        ENDDO
      ENDDO
      
      ijb=ij_begin-iip1
      ije=ij_end
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1
      
      DO l=1,llm
        DO ij=ijb,ije
           mv(ij,l)=pbarv(ij,l) * zzpbar
        ENDDO
      ENDDO
      
      ijb=ij_begin
      ije=ij_end
      
      DO l=1,llm
        DO ij=ijb,ije
           mw(ij,l)=w(ij,l) * zzw
        ENDDO
      ENDDO

      DO ij=ijb,ije
         mw(ij,llm+1)=0.
      ENDDO
      
c      CALL SCOPY(ijp1llm,q,1,zq,1)
c      CALL SCOPY(ijp1llm,masse,1,zm,1)
       
       ijb=ij_begin
       ije=ij_end
       zq(ijb:ije,:)=q(ijb:ije,:)
       zm(ijb:ije,:)=masse(ijb:ije,:)
      
      
c	print*,'Entree vlx1'
c	call minmaxq(zq,qmin,qmax,'avant vlx     ')
      call vlx_p(zq,pente_max,zm,mu,ij_begin,ij_begin+2*iip1-1)
      call vlx_p(zq,pente_max,zm,mu,ij_end-2*iip1+1,ij_end)
      call VTb(VTHallo)
      call Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest1)
      call Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest1)
      call SendRequest(MyRequest1)
      call VTe(VTHallo)
      call vlx_p(zq,pente_max,zm,mu,ij_begin+2*iip1,ij_end-2*iip1)
c      call vlx_p(zq,pente_max,zm,mu,ij_begin,ij_end)
      call VTb(VTHallo)
      call WaitRecvRequest(MyRequest1)
      call VTe(VTHallo)

      
c	print*,'Sortie vlx1'
c	call minmaxq(zq,qmin,qmax,'apres vlx1    ')

c	 print*,'Entree vly1'
c      call exchange_hallo(zq,ip1jmp1,llm,2,2)
c      call exchange_hallo(zm,ip1jmp1,llm,1,1)
      
      call vly_p(zq,pente_max,zm,mv)
c	call minmaxq(zq,qmin,qmax,'apres vly1     ')
c	print*,'Sortie vly1'
      call vlz_p(zq,pente_max,zm,mw,ij_begin,ij_begin+2*iip1-1)
      call vlz_p(zq,pente_max,zm,mw,ij_end-2*iip1+1,ij_end)
      call VTb(VTHallo)
      call Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest2)
      call Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest2)
      call SendRequest(MyRequest2)
      call VTe(VTHallo)
      call vlz_p(zq,pente_max,zm,mw,ij_begin+2*iip1,ij_end-2*iip1)
      call VTb(VTHallo)
      call WaitRecvRequest(MyRequest2)
            
      call VTe(VTHallo)
      
c	call minmaxq(zq,qmin,qmax,'apres vlz     ')


      
      
      call vly_p(zq,pente_max,zm,mv)
c	call minmaxq(zq,qmin,qmax,'apres vly     ')


      call vlx_p(zq,pente_max,zm,mu,ij_begin,ij_end)
c	call minmaxq(zq,qmin,qmax,'apres vlx2    ')

	
      ijb=ij_begin
      ije=ij_end
       
      DO l=1,llm
         DO ij=ijb,ije
           q(ij,l)=zq(ij,l)
         ENDDO
      ENDDO
      
      
      DO l=1,llm
         DO ij=ijb,ije-iip1+1,iip1
            q(ij+iim,l)=q(ij,l)
         ENDDO
      ENDDO

      call WaitSendRequest(MyRequest1) 
      call WaitSendRequest(MyRequest2)
     
      RETURN
      END
      
      
      SUBROUTINE vlx_p(q,pente_max,masse,u_m,ijb_x,ije_x)

c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE Parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL u_m( ip1jmp1,llm ),pbarv( iip1,jjm,llm)
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
      INTEGER n0,iadvplus(ip1jmp1,llm),nl(llm)
c
      REAL new_m,zu_m,zdum(ip1jmp1,llm)
      REAL sigu(ip1jmp1),dxq(ip1jmp1,llm),dxqu(ip1jmp1)
      REAL zz(ip1jmp1)
      REAL adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
      REAL u_mq(ip1jmp1,llm)

      Logical extremum

      REAL      SSUM
      EXTERNAL  SSUM

      REAL z1,z2,z3

      INTEGER ijb,ije,ijb_x,ije_x
      
c   calcul de la pente a droite et a gauche de la maille

      ijb=ijb_x
      ije=ije_x
        
      if (pole_nord.and.ijb==1) ijb=ijb+iip1
      if (pole_sud.and.ije==ip1jmp1)  ije=ije-iip1
         
      IF (pente_max.gt.-1.e-5) THEN
c       IF (pente_max.gt.10) THEN

c   calcul des pentes avec limitation, Van Leer scheme I:
c   -----------------------------------------------------

c   calcul de la pente aux points u
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)         
         DO l = 1, llm
            
            DO ij=ijb,ije-1
               dxqu(ij)=q(ij+1,l)-q(ij,l)
c              IF(u_m(ij,l).lt.0.) stop'limx n admet pas les U<0'
c              sigu(ij)=u_m(ij,l)/masse(ij,l)
            ENDDO
            DO ij=ijb+iip1-1,ije,iip1
               dxqu(ij)=dxqu(ij-iim)
c              sigu(ij)=sigu(ij-iim)
            ENDDO

            DO ij=ijb,ije
               adxqu(ij)=abs(dxqu(ij))
            ENDDO

c   calcul de la pente maximum dans la maille en valeur absolue

            DO ij=ijb+1,ije
               dxqmax(ij,l)=pente_max*
     ,      min(adxqu(ij-1),adxqu(ij))
c limitation subtile
c    ,      min(adxqu(ij-1)/sigu(ij-1),adxqu(ij)/(1.-sigu(ij)))
          

            ENDDO

            DO ij=ijb+iip1-1,ije,iip1
               dxqmax(ij-iim,l)=dxqmax(ij,l)
            ENDDO

            DO ij=ijb+1,ije
#ifdef CRAY
               dxq(ij,l)=
     ,         cvmgp(dxqu(ij-1)+dxqu(ij),0.,dxqu(ij-1)*dxqu(ij))
#else
               IF(dxqu(ij-1)*dxqu(ij).gt.0) THEN
                  dxq(ij,l)=dxqu(ij-1)+dxqu(ij)
               ELSE
c   extremum local
                  dxq(ij,l)=0.
               ENDIF
#endif
               dxq(ij,l)=0.5*dxq(ij,l)
               dxq(ij,l)=
     ,         sign(min(abs(dxq(ij,l)),dxqmax(ij,l)),dxq(ij,l))
            ENDDO

         ENDDO ! l=1,llm
c$OMP END DO NOWAIT
c	print*,'Ok calcul des pentes'

      ELSE ! (pente_max.lt.-1.e-5)

c   Pentes produits:
c   ----------------
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l = 1, llm
            DO ij=ijb,ije-1
               dxqu(ij)=q(ij+1,l)-q(ij,l)
            ENDDO
            DO ij=ijb+iip1-1,ije,iip1
               dxqu(ij)=dxqu(ij-iim)
            ENDDO

            DO ij=ijb+1,ije
               zz(ij)=dxqu(ij-1)*dxqu(ij)
               zz(ij)=zz(ij)+zz(ij)
               IF(zz(ij).gt.0) THEN
                  dxq(ij,l)=zz(ij)/(dxqu(ij-1)+dxqu(ij))
               ELSE
c   extremum local
                  dxq(ij,l)=0.
               ENDIF
            ENDDO

         ENDDO
c$OMP END DO NOWAIT
      ENDIF ! (pente_max.lt.-1.e-5)

c   bouclage de la pente en iip1:
c   -----------------------------
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb+iip1-1,ije,iip1
            dxq(ij-iim,l)=dxq(ij,l)
         ENDDO
         DO ij=ijb,ije
            iadvplus(ij,l)=0
         ENDDO

      ENDDO
c$OMP END DO NOWAIT
c	 print*,'Bouclage en iip1'

c   calcul des flux a gauche et a droite

#ifdef CRAY
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb,ije-1
          zdum(ij,l)=cvmgp(1.-u_m(ij,l)/masse(ij,l),
     ,                     1.+u_m(ij,l)/masse(ij+1,l),
     ,                     u_m(ij,l))
          zdum(ij,l)=0.5*zdum(ij,l)
          u_mq(ij,l)=cvmgp(
     ,                q(ij,l)+zdum(ij,l)*dxq(ij,l),
     ,                q(ij+1,l)-zdum(ij,l)*dxq(ij+1,l),
     ,                u_m(ij,l))
          u_mq(ij,l)=u_m(ij,l)*u_mq(ij,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
#else
c   on cumule le flux correspondant a toutes les mailles dont la masse
c   au travers de la paroi pENDant le pas de temps.
c	print*,'Cumule ....'
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb,ije-1
c	print*,'masse(',ij,')=',masse(ij,l)
          IF (u_m(ij,l).gt.0.) THEN
             zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
             u_mq(ij,l)=u_m(ij,l)*(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l))
          ELSE
             zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
             u_mq(ij,l)=u_m(ij,l)*(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l))
          ENDIF
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
#endif
c	stop

c	go to 9999
c   detection des points ou on advecte plus que la masse de la
c   maille
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb,ije-1
            IF(zdum(ij,l).lt.0) THEN
               iadvplus(ij,l)=1
               u_mq(ij,l)=0.
            ENDIF
         ENDDO
      ENDDO
c$OMP END DO NOWAIT
c	print*,'Ok test 1'
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb+iip1-1,ije,iip1
          iadvplus(ij,l)=iadvplus(ij-iim,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
c	 print*,'Ok test 2'


c   traitement special pour le cas ou on advecte en longitude plus que le
c   contenu de la maille.
c   cette partie est mal vectorisee.

c  calcul du nombre de maille sur lequel on advecte plus que la maille.

      n0=0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         nl(l)=0
         DO ij=ijb,ije
            nl(l)=nl(l)+iadvplus(ij,l)
         ENDDO
         n0=n0+nl(l)
      ENDDO
c$OMP END DO NOWAIT
cym      IF(n0.gt.1) THEN
cym      IF(n0.gt.0) THEN

c      PRINT*,'Nombre de points pour lesquels on advect plus que le'
c     &       ,'contenu de la maille : ',n0
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l=1,llm
            IF(nl(l).gt.0) THEN
               iju=0
c   indicage des mailles concernees par le traitement special
               DO ij=ijb,ije
                  IF(iadvplus(ij,l).eq.1.and.mod(ij,iip1).ne.0) THEN
                     iju=iju+1
                     indu(iju)=ij
                  ENDIF
               ENDDO
               niju=iju
c              PRINT*,'niju,nl',niju,nl(l)

c  traitement des mailles
               DO iju=1,niju
                  ij=indu(iju)
                  j=(ij-1)/iip1+1
                  zu_m=u_m(ij,l)
                  u_mq(ij,l)=0.
                  IF(zu_m.gt.0.) THEN
                     ijq=ij
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)+q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m-masse(ijq,l)
                        i=mod(i-2+iim,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*
     &               (q(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq(ijq,l))
                  ELSE
                     ijq=ij+1
                     i=ijq-(j-1)*iip1
c   accumulation pour les mailles completements advectees
                     do while(-zu_m.gt.masse(ijq,l))
                        u_mq(ij,l)=u_mq(ij,l)-q(ijq,l)*masse(ijq,l)
                        zu_m=zu_m+masse(ijq,l)
                        i=mod(i,iim)+1
                        ijq=(j-1)*iip1+i
                     ENDDO
c   ajout de la maille non completement advectee
                     u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l)-
     &               0.5*(1.+zu_m/masse(ijq,l))*dxq(ijq,l))
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
c$OMP END DO NOWAIT
cym      ENDIF  ! n0.gt.0 
9999    continue


c   bouclage en latitude
c	print*,'Avant bouclage en latitude'
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        DO ij=ijb+iip1-1,ije,iip1
           u_mq(ij,l)=u_mq(ij-iim,l)
        ENDDO
      ENDDO
c$OMP END DO NOWAIT

c   calcul des tENDances
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb+1,ije
            new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+
     &      u_mq(ij-1,l)-u_mq(ij,l))
     &      /new_m
            masse(ij,l)=new_m
         ENDDO
c   ModIF Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
         DO ij=ijb+iip1-1,ije,iip1
            q(ij-iim,l)=q(ij,l)
            masse(ij-iim,l)=masse(ij,l)
         ENDDO
      ENDDO
c$OMP END DO NOWAIT
c     CALL SCOPY((jjm-1)*llm,q(iip1+iip1,1),iip1,q(iip2,1),iip1)
c     CALL SCOPY((jjm-1)*llm,masse(iip1+iip1,1),iip1,masse(iip2,1),iip1)


      RETURN
      END


      SUBROUTINE vly_p(q,pente_max,masse,masse_adv_v)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE parallel_lmdz
      USE comconst_mod, ONLY: pi

      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL masse_adv_v( ip1jm,llm)
      REAL q(ip1jmp1,llm), dq( ip1jmp1,llm)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l
c
      REAL airej2,airejjm,airescb(iim),airesch(iim)
      REAL dyq(ip1jmp1,llm),dyqv(ip1jm),zdvm(ip1jmp1,llm)
      REAL adyqv(ip1jm),dyqmax(ip1jmp1)
      REAL qbyv(ip1jm,llm)

      REAL qpns,qpsn,appn,apps,dyn1,dys1,dyn2,dys2,newmasse,fn,fs
c     REAL newq,oldmasse
      Logical extremum,first,testcpu
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
c$OMP THREADPRIVATE(temps0,temps1,temps2,temps3,temps4,temps5)
      SAVE first,testcpu
c$OMP THREADPRIVATE(first,testcpu)

      REAL convpn,convps,convmpn,convmps
      real massepn,masseps,qpn,qps
      REAL sinlon(iip1),sinlondlon(iip1)
      REAL coslon(iip1),coslondlon(iip1)
      SAVE sinlon,coslon,sinlondlon,coslondlon
c$OMP THREADPRIVATE(sinlon,coslon,sinlondlon,coslondlon)
      SAVE airej2,airejjm
c$OMP THREADPRIVATE(airej2,airejjm)
c
c
      REAL      SSUM
      EXTERNAL  SSUM

      DATA first,testcpu/.true.,.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./
      INTEGER ijb,ije

      IF(first) THEN
c         PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
         first=.false.
         do i=2,iip1
            coslon(i)=cos(rlonv(i))
            sinlon(i)=sin(rlonv(i))
            coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
            sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
         ENDDO
         coslon(1)=coslon(iip1)
         coslondlon(1)=coslondlon(iip1)
         sinlon(1)=sinlon(iip1)
         sinlondlon(1)=sinlondlon(iip1)
         airej2 = SSUM( iim, aire(iip2), 1 )
         airejjm= SSUM( iim, aire(ip1jm -iim), 1 ) 
      ENDIF

c
c	PRINT*,'CALCUL EN LATITUDE'

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1, llm
c
c   --------------------------------
c      CALCUL EN LATITUDE
c   --------------------------------

c   On commence par calculer la valeur du traceur moyenne sur le premier cercle
c   de latitude autour du pole (qpns pour le pole nord et qpsn pour
c    le pole nord) qui sera utilisee pour evaluer les pentes au pole.
      
      if (pole_nord) then
        DO i = 1, iim
          airescb(i) = aire(i+ iip1) * q(i+ iip1,l)
        ENDDO
        qpns   = SSUM( iim,  airescb ,1 ) / airej2
      endif
      
      if (pole_sud) then
        DO i = 1, iim
          airesch(i) = aire(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
        ENDDO
        qpsn   = SSUM( iim,  airesch ,1 ) / airejjm
      endif
      
      

c   calcul des pentes aux points v

      ijb=ij_begin-2*iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1
      
      DO ij=ijb,ije
         dyqv(ij)=q(ij,l)-q(ij+iip1,l)
         adyqv(ij)=abs(dyqv(ij))
      ENDDO

c   calcul des pentes aux points scalaires
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO ij=ijb,ije
         dyq(ij,l)=.5*(dyqv(ij-iip1)+dyqv(ij))
         dyqmax(ij)=min(adyqv(ij-iip1),adyqv(ij))
         dyqmax(ij)=pente_max*dyqmax(ij)
      ENDDO

c   calcul des pentes aux poles
      IF (pole_nord) THEN
        DO ij=1,iip1
           dyq(ij,l)=qpns-q(ij+iip1,l)
        ENDDO
        
        dyn1=0.
        dyn2=0.
        DO ij=1,iim
          dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
          dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
        ENDDO
        DO ij=1,iip1
          dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
         dyq(ij,l)=0.
        ENDDO
c ym tout cela ne sert pas a grand chose
      ENDIF
      
      IF (pole_sud) THEN

        DO ij=1,iip1
           dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l)-qpsn
        ENDDO

        dys1=0.
        dys2=0.

        DO ij=1,iim
          dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
          dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
        ENDDO

        DO ij=1,iip1
          dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
         dyq(ip1jm+ij,l)=0.
        ENDDO
c ym tout cela ne sert pas a grand chose
      ENDIF

c   filtrage de la derivee

c   calcul des pentes limites aux poles
c ym partie inutile
c      goto 8888
c      fn=1.
c      fs=1.
c      DO ij=1,iim
c         IF(pente_max*adyqv(ij).lt.abs(dyq(ij,l))) THEN
c            fn=min(pente_max*adyqv(ij)/abs(dyq(ij,l)),fn)
c         ENDIF
c      IF(pente_max*adyqv(ij+ip1jm-iip1).lt.abs(dyq(ij+ip1jm,l))) THEN
c         fs=min(pente_max*adyqv(ij+ip1jm-iip1)/abs(dyq(ij+ip1jm,l)),fs)
c         ENDIF
c      ENDDO
c      DO ij=1,iip1
c         dyq(ij,l)=fn*dyq(ij,l)
c         dyq(ip1jm+ij,l)=fs*dyq(ip1jm+ij,l)
c      ENDDO
c 8888    continue


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  En memoire de dIFferents tests sur la 
C  limitation des pentes aux poles.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     PRINT*,dyq(1)
C     PRINT*,dyqv(iip1+1)
C     appn=abs(dyq(1)/dyqv(iip1+1))
C     PRINT*,dyq(ip1jm+1)
C     PRINT*,dyqv(ip1jm-iip1+1)
C     apps=abs(dyq(ip1jm+1)/dyqv(ip1jm-iip1+1))
C     DO ij=2,iim
C        appn=amax1(abs(dyq(ij)/dyqv(ij)),appn)
C        apps=amax1(abs(dyq(ip1jm+ij)/dyqv(ip1jm-iip1+ij)),apps)
C     ENDDO
C     appn=min(pente_max/appn,1.)
C     apps=min(pente_max/apps,1.)
C
C
C   cas ou on a un extremum au pole
C
C     IF(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
C    &   appn=0.
C     IF(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
C    &   dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
C    &   apps=0.
C
C   limitation des pentes aux poles
C     DO ij=1,iip1
C        dyq(ij)=appn*dyq(ij)
C        dyq(ip1jm+ij)=apps*dyq(ip1jm+ij)
C     ENDDO
C
C   test
C      DO ij=1,iip1
C         dyq(iip1+ij)=0.
C         dyq(ip1jm+ij-iip1)=0.
C      ENDDO
C      DO ij=1,ip1jmp1
C         dyq(ij)=dyq(ij)*cos(rlatu((ij-1)/iip1+1))
C      ENDDO
C
C changement 10 07 96
C     IF(dyqv(ismin(iim,dyqv,1))*dyqv(ismax(iim,dyqv,1)).le.0.)
C    &   THEN
C        DO ij=1,iip1
C           dyqmax(ij)=0.
C        ENDDO
C     ELSE
C        DO ij=1,iip1
C           dyqmax(ij)=pente_max*abs(dyqv(ij))
C        ENDDO
C     ENDIF
C
C     IF(dyqv(ismax(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1)*
C    & dyqv(ismin(iim,dyqv(ip1jm-iip1+1),1)+ip1jm-iip1+1).le.0.)
C    &THEN
C        DO ij=ip1jm+1,ip1jmp1
C           dyqmax(ij)=0.
C        ENDDO
C     ELSE
C        DO ij=ip1jm+1,ip1jmp1
C           dyqmax(ij)=pente_max*abs(dyqv(ij-iip1))
C        ENDDO
C     ENDIF
C   fin changement 10 07 96
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c   calcul des pentes limitees
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1

      DO ij=ijb,ije
         IF(dyqv(ij)*dyqv(ij-iip1).gt.0.) THEN
            dyq(ij,l)=sign(min(abs(dyq(ij,l)),dyqmax(ij)),dyq(ij,l))
         ELSE
            dyq(ij,l)=0.
         ENDIF
      ENDDO

      ENDDO
c$OMP END DO NOWAIT

      ijb=ij_begin-iip1
      ije=ij_end
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       DO ij=ijb,ije
          IF(masse_adv_v(ij,l).gt.0) THEN
              qbyv(ij,l)=q(ij+iip1,l)+dyq(ij+iip1,l)*
     ,                   0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
          ELSE
              qbyv(ij,l)=q(ij,l)-dyq(ij,l)*
     ,                   0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
          ENDIF
          qbyv(ij,l)=masse_adv_v(ij,l)*qbyv(ij,l)
       ENDDO
      ENDDO
c$OMP END DO NOWAIT
      
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
         DO ij=ijb,ije
            newmasse=masse(ij,l)
     &      +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
     
            q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l))
     &         /newmasse
            masse(ij,l)=newmasse
         ENDDO
c.-. ancienne version
c        convpn=SSUM(iim,qbyv(1,l),1)/apoln
c        convmpn=ssum(iim,masse_adv_v(1,l),1)/apoln
         if (pole_nord) then
           convpn=SSUM(iim,qbyv(1,l),1)
           convmpn=ssum(iim,masse_adv_v(1,l),1)
           massepn=ssum(iim,masse(1,l),1)
           qpn=0.
           do ij=1,iim
              qpn=qpn+masse(ij,l)*q(ij,l)
           enddo
           qpn=(qpn+convpn)/(massepn+convmpn)
           do ij=1,iip1
              q(ij,l)=qpn
           enddo
         endif
         
c        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)/apols
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)/apols
         
         if (pole_sud) then
         
           convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
           convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
           masseps=ssum(iim, masse(ip1jm+1,l),1)
           qps=0.
           do ij = ip1jm+1,ip1jmp1-1
              qps=qps+masse(ij,l)*q(ij,l)
           enddo
           qps=(qps+convps)/(masseps+convmps)
           do ij=ip1jm+1,ip1jmp1
              q(ij,l)=qps
           enddo
         endif
c.-. fin ancienne version

c._. nouvelle version
c        convpn=SSUM(iim,qbyv(1,l),1)
c        convmpn=ssum(iim,masse_adv_v(1,l),1)
c        oldmasse=ssum(iim,masse(1,l),1)
c        newmasse=oldmasse+convmpn
c        newq=(q(1,l)*oldmasse+convpn)/newmasse
c        newmasse=newmasse/apoln
c        DO ij = 1,iip1
c           q(ij,l)=newq
c           masse(ij,l)=newmasse*aire(ij)
c        ENDDO
c        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
c        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
c        oldmasse=ssum(iim,masse(ip1jm-iim,l),1)
c        newmasse=oldmasse+convmps
c        newq=(q(ip1jmp1,l)*oldmasse+convps)/newmasse
c        newmasse=newmasse/apols
c        DO ij = ip1jm+1,ip1jmp1
c           q(ij,l)=newq
c           masse(ij,l)=newmasse*aire(ij)
c        ENDDO
c._. fin nouvelle version
      ENDDO
c$OMP END DO NOWAIT

      RETURN
      END
      
      
      
      SUBROUTINE vlz_p(q,pente_max,masse,w,ijb_x,ije_x)
c
c     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
c
c    ********************************************************************
c     Shema  d'advection " pseudo amont " .
c    ********************************************************************
c    q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
c     dq 	       sont des arguments de sortie pour le s-pg ....
c
c
c   --------------------------------------------------------------------
      USE Parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
c
c
c   Arguments:
c   ----------
      REAL masse(ip1jmp1,llm),pente_max
      REAL q(ip1jmp1,llm)
      REAL w(ip1jmp1,llm+1)
c
c      Local 
c   ---------
c
      INTEGER i,ij,l,j,ii
c
      REAL,SAVE :: wq(ip1jmp1,llm+1)
      REAL newmasse

      REAL,SAVE :: dzq(ip1jmp1,llm),dzqw(ip1jmp1,llm),adzqw(ip1jmp1,llm)
      REAL dzqmax
      REAL sigw

      LOGICAL testcpu
      SAVE testcpu
c$OMP THREADPRIVATE(testcpu)
      REAL temps0,temps1,temps2,temps3,temps4,temps5,second
      SAVE temps0,temps1,temps2,temps3,temps4,temps5
c$OMP THREADPRIVATE(temps0,temps1,temps2,temps3,temps4,temps5)

      REAL      SSUM
      EXTERNAL  SSUM

      DATA testcpu/.false./
      DATA temps0,temps1,temps2,temps3,temps4,temps5/0.,0.,0.,0.,0.,0./
      INTEGER ijb,ije,ijb_x,ije_x
c    On oriente tout dans le sens de la pression c'est a dire dans le
c    sens de W

#ifdef BIDON
      IF(testcpu) THEN
         temps0=second(0.)
      ENDIF
#endif

      ijb=ijb_x
      ije=ije_x

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=2,llm
         DO ij=ijb,ije
            dzqw(ij,l)=q(ij,l-1)-q(ij,l)
            adzqw(ij,l)=abs(dzqw(ij,l))
         ENDDO
      ENDDO
c$OMP END DO

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=2,llm-1
         DO ij=ijb,ije
#ifdef CRAY
            dzq(ij,l)=0.5*
     ,      cvmgp(dzqw(ij,l)+dzqw(ij,l+1),0.,dzqw(ij,l)*dzqw(ij,l+1))
#else
            IF(dzqw(ij,l)*dzqw(ij,l+1).gt.0.) THEN
                dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
            ELSE
                dzq(ij,l)=0.
            ENDIF
#endif
            dzqmax=pente_max*min(adzqw(ij,l),adzqw(ij,l+1))
            dzq(ij,l)=sign(min(abs(dzq(ij,l)),dzqmax),dzq(ij,l))
         ENDDO
      ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER
      DO ij=ijb,ije
         dzq(ij,1)=0.
         dzq(ij,llm)=0.
      ENDDO
c$OMP END MASTER
c$OMP BARRIER
#ifdef BIDON
      IF(testcpu) THEN
         temps1=temps1+second(0.)-temps0
      ENDIF
#endif
c ---------------------------------------------------------------
c   .... calcul des termes d'advection verticale  .......
c ---------------------------------------------------------------

c calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour calculer dq

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
       DO l = 1,llm-1
         do  ij = ijb,ije
          IF(w(ij,l+1).gt.0.) THEN
             sigw=w(ij,l+1)/masse(ij,l+1)
             wq(ij,l+1)=w(ij,l+1)*(q(ij,l+1)+0.5*(1.-sigw)*dzq(ij,l+1))
          ELSE
             sigw=w(ij,l+1)/masse(ij,l)
             wq(ij,l+1)=w(ij,l+1)*(q(ij,l)-0.5*(1.+sigw)*dzq(ij,l))
          ENDIF
         ENDDO
       ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER
       DO ij=ijb,ije
          wq(ij,llm+1)=0.
          wq(ij,1)=0.
       ENDDO
c$OMP END MASTER
c$OMP BARRIER

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
         DO ij=ijb,ije
            newmasse=masse(ij,l)+w(ij,l+1)-w(ij,l)
            q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l))
     &         /newmasse
            masse(ij,l)=newmasse
         ENDDO
      ENDDO
c$OMP END DO NOWAIT


      RETURN
      END
c      SUBROUTINE minmaxq(zq,qmin,qmax,comment)
c
c#include "dimensions.h"
c#include "paramet.h"

c      CHARACTER*(*) comment
c      real qmin,qmax
c      real zq(ip1jmp1,llm)

c      INTEGER jadrs(ip1jmp1), jbad, k, i


c      DO k = 1, llm
c         jbad = 0
c         DO i = 1, ip1jmp1
c         IF (zq(i,k).GT.qmax .OR. zq(i,k).LT.qmin) THEN
c            jbad = jbad + 1
c            jadrs(jbad) = i
c         ENDIF
c         ENDDO
c         IF (jbad.GT.0) THEN
c         PRINT*, comment
c         DO i = 1, jbad
cc            PRINT*, "i,k,zq=", jadrs(i),k,zq(jadrs(i),k)
c         ENDDO
c         ENDIF
c      ENDDO

c      return
c      end


      subroutine minmaxq_p(zq,qmin,qmax,comment)

#include "dimensions.h"
#include "paramet.h"

      character*20 comment
      real qmin,qmax
      real zq(ip1jmp1,llm)
      real zzq(iip1,jjp1,llm)

      integer imin,jmin,lmin,ijlmin
      integer imax,jmax,lmax,ijlmax

      integer ismin,ismax

#ifdef isminismax
      call scopy (ip1jmp1*llm,zq,1,zzq,1)

      ijlmin=ismin(ijp1llm,zq,1)
      lmin=(ijlmin-1)/ip1jmp1+1
      ijlmin=ijlmin-(lmin-1.)*ip1jmp1
      jmin=(ijlmin-1)/iip1+1
      imin=ijlmin-(jmin-1.)*iip1
      zqmin=zq(ijlmin,lmin)

      ijlmax=ismax(ijp1llm,zq,1)
      lmax=(ijlmax-1)/ip1jmp1+1
      ijlmax=ijlmax-(lmax-1.)*ip1jmp1
      jmax=(ijlmax-1)/iip1+1
      imax=ijlmax-(jmax-1.)*iip1
      zqmax=zq(ijlmax,lmax)

       if(zqmin.lt.qmin) 
c     s     write(*,9999) comment,
     s     write(*,*) comment,
     s     imin,jmin,lmin,zqmin,zzq(imin,jmin,lmin)
       if(zqmax.gt.qmax) 
c     s     write(*,9999) comment,
     s     write(*,*) comment,
     s     imax,jmax,lmax,zqmax,zzq(imax,jmax,lmax)
#endif
      return
9999  format(a20,'  q(',i3,',',i2,',',i2,')=',e12.5,e12.5)
      end



