      SUBROUTINE convadj(ngrid,nlay,nq,ptimestep,
     S                   pplay,pplev,ppopsk,
     $                   pu,pv,ph,pq,
     $                   pdufi,pdvfi,pdhfi,pdqfi,
     $                   pduadj,pdvadj,pdhadj,
     $                   pdqadj)
      IMPLICIT NONE

c=======================================================================
c
c   ajustement convectif sec
c   on ajoute les tendances pdhfi au profil pdh avant l'ajustement
c   SPECIAL VERSION : if one tracer is CO2, take into account the
c   Molecular mass variation (e.g. when CO2 condense) to trigger
c   convection  F. Forget 01/2005
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "tracer.h"


c   arguments:
c   ----------

      INTEGER ngrid,nlay
      REAL ptimestep
      REAL ph(ngrid,nlay),pdhfi(ngrid,nlay),pdhadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1),ppopsk(ngrid,nlay)
      REAL pu(ngrid,nlay),pdufi(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvfi(ngrid,nlay),pdvadj(ngrid,nlay)

c    Traceurs :
      integer nq
      real pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq)
      real pdqadj(ngrid,nlay,nq)


c   local:
c   ------

      INTEGER ig,i,l,l1,l2,jj
      INTEGER jcnt, jadrs(ngridmx)

      REAL sig(nlayermx+1),sdsig(nlayermx),dsig(nlayermx)
      REAL zu(ngridmx,nlayermx),zv(ngridmx,nlayermx)
      REAL zh(ngridmx,nlayermx)
      REAL zu2(ngridmx,nlayermx),zv2(ngridmx,nlayermx)
      REAL zh2(ngridmx,nlayermx), zhc(ngridmx,nlayermx)
      REAL zhm,zsm,zum,zvm,zalpha,zhmc
      

c    Traceurs :
      INTEGER iq,ico2
      save ico2
      REAL zq(ngridmx,nlayermx,nqmx), zq2(ngridmx,nlayermx,nqmx)
      REAL zqm(nqmx),zqco2m
      real m_co2, m_noco2, A , B
      save A, B
c    Temporaire (for diagnostic)
c      REAL diag_alpha(ngridmx)


      LOGICAL vtest(ngridmx),down,firstcall
      save firstcall
      data firstcall/.true./

      EXTERNAL SCOPY
c
c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
      IF (firstcall) THEN 
        IF(ngrid.NE.ngridmx) THEN
           PRINT*
           PRINT*,'STOP dans convadj'
           PRINT*,'ngrid    =',ngrid
           PRINT*,'ngridmx  =',ngridmx
        ENDIF
        ico2=0
        if (tracer) then
c          Prepare Special treatment if one of the tracer is CO2 gas
c           write(*,*) "convadj: nqmx=",nqmx
           do iq=1,nqmx
c           write(*,*) "convadj: iq=",iq," noms(iq)=",noms(iq)
             if (noms(iq).eq."co2") then
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/m-3)   
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/m-3)   
c               Compute A and B coefficient use to compute
c               mean molecular mass Mair defined by
c               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
c               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             end if
           enddo
           firstcall=.false.
        endif
      ENDIF

c
c-----------------------------------------------------------------------
c   detection des profils a modifier:
c   ---------------------------------
c   si le profil est a modifier
c   (i.e. ph(niv_sup) < ph(niv_inf) )
c   alors le tableau "vtest" est mis a .TRUE. ;
c   sinon, il reste a sa valeur initiale (.FALSE.)
c   cette operation est vectorisable
c   On en profite pour copier la valeur initiale de "ph"
c   dans le champ de travail "zh"


      DO l=1,nlay
         DO ig=1,ngrid
            zh(ig,l)=ph(ig,l)+pdhfi(ig,l)*ptimestep
            zu(ig,l)=pu(ig,l)+pdufi(ig,l)*ptimestep
            zv(ig,l)=pv(ig,l)+pdvfi(ig,l)*ptimestep
         ENDDO
      ENDDO

      if(tracer) then      
        DO iq =1, nq
         DO l=1,nlay
           DO ig=1,ngrid
              zq(ig,l,iq)=pq(ig,l,iq)+pdqfi(ig,l,iq)*ptimestep
           ENDDO
         ENDDO
        ENDDO
      end if


      CALL scopy(ngrid*nlay,zh,1,zh2,1)
      CALL scopy(ngrid*nlay,zu,1,zu2,1)
      CALL scopy(ngrid*nlay,zv,1,zv2,1)
      CALL scopy(ngrid*nlay*nq,zq,1,zq2,1)

      DO ig=1,ngrid
        vtest(ig)=.FALSE.
      ENDDO
c
      if (ico2.ne.0) then
c        Special case if one of the tracer is CO2 gas
c         write(*,*) "convadj: ico2=",ico2," nqmx=",nqmx
         DO l=1,nlay
           DO ig=1,ngrid
             zhc(ig,l) = zh2(ig,l)*(A*zq2(ig,l,ico2)+B)
           ENDDO
         ENDDO
       else
          CALL scopy(ngrid*nlay,zh2,1,zhc,1)
       end if

       

      DO l=2,nlay
        DO ig=1,ngrid
          IF(zhc(ig,l).LT.zhc(ig,l-1)) vtest(ig)=.TRUE.
        ENDDO
      ENDDO
c
CRAY  CALL WHENNE(ngrid, vtest, 1, 0, jadrs, jcnt)
      jcnt=0
      DO ig=1,ngrid
         IF(vtest(ig)) THEN
            jcnt=jcnt+1
            jadrs(jcnt)=ig
         ENDIF
      ENDDO


c-----------------------------------------------------------------------
c  Ajustement des "jcnt" profils instables indices par "jadrs":
c  ------------------------------------------------------------
c
      DO jj = 1, jcnt
c
          i = jadrs(jj)
c
c   Calcul des niveaux sigma sur cette colonne
          DO l=1,nlay+1
            sig(l)=pplev(i,l)/pplev(i,1)
          ENDDO
         DO l=1,nlay
            dsig(l)=sig(l)-sig(l+1)
            sdsig(l)=ppopsk(i,l)*dsig(l)
         ENDDO
          l2 = 1
c
c      -- boucle de sondage vers le haut
c
cins$     Loop
          DO
c
            l2 = l2 + 1
c
cins$       Exit
            IF (l2 .GT. nlay) EXIT
c
            IF (zhc(i, l2) .LT. zhc(i, l2-1)) THEN
c
c          -- l2 est le niveau le plus haut de la colonne instable
c
              l1 = l2 - 1
              l  = l1
              zsm = sdsig(l2)
              zhm = zh2(i, l2)
              if(ico2.ne.0) zqco2m = zq2(i,l2,ico2)
c
c          -- boucle de sondage vers le bas
c
cins$         Loop
              DO
c
                zsm = zsm + sdsig(l)
                zhm = zhm + sdsig(l) * (zh2(i, l) - zhm) / zsm
                if(ico2.ne.0) then
                  zqco2m = 
     &            zqco2m + sdsig(l) * (zq2(i,l,ico2) - zqco2m) / zsm
                  zhmc = zhm*(A*zqco2m+B)
                else 
                  zhmc = zhm
                end if
c
c            -- doit on etendre la colonne vers le bas ?
c
c_EC (M1875) 20/6/87 : AND -> AND THEN
c
                down = .FALSE.
                IF (l1 .NE. 1) THEN    !--  and then
                  IF (zhmc .LT. zhc(i, l1-1)) THEN
                    down = .TRUE.
                  END IF
                END IF
c
                IF (down) THEN
c
                  l1 = l1 - 1
                  l  = l1
c
                ELSE
c
c              -- peut on etendre la colonne vers le haut ?
c
cins$             Exit
                  IF (l2 .EQ. nlay) EXIT
c
cins$             Exit
                  IF (zhc(i, l2+1) .GE. zhmc) EXIT
c
                  l2 = l2 + 1
                  l  = l2
c
                END IF
c
cins$         End Loop
              ENDDO
c
c          -- nouveau profil : constant (valeur moyenne)
c
              zalpha=0.
              zum=0.
              zvm=0.
              do iq=1,nq
                zqm(iq) = 0.
              end do
              DO l = l1, l2
                if(ico2.ne.0) then
                  zalpha=zalpha+
     &            ABS(zhc(i,l)/(A+B*zqco2m) -zhm)*dsig(l)
                else
                  zalpha=zalpha+ABS(zh2(i,l)-zhm)*dsig(l)
                endif
                zh2(i, l) = zhm
                zum=zum+dsig(l)*zu(i,l)
                zvm=zvm+dsig(l)*zv(i,l)
                do iq=1,nq
                   zqm(iq) = zqm(iq)+dsig(l)*zq(i,l,iq)
                end do
              ENDDO
              zalpha=zalpha/(zhm*(sig(l1)-sig(l2+1)))
              zum=zum/(sig(l1)-sig(l2+1))
              zvm=zvm/(sig(l1)-sig(l2+1))
              do iq=1,nq
                 zqm(iq) = zqm(iq)/(sig(l1)-sig(l2+1))
              end do

              IF(zalpha.GT.1.) THEN
c                PRINT*,'WARNING dans convadj zalpha=',zalpha
c         STOP
                 zalpha=1.
              ELSE
c                IF(zalpha.LT.0.) STOP
                 IF(zalpha.LT.1.e-4) zalpha=1.e-4
              ENDIF

              DO l=l1,l2
                 zu2(i,l)=zu2(i,l)+zalpha*(zum-zu2(i,l))
                 zv2(i,l)=zv2(i,l)+zalpha*(zvm-zv2(i,l))
                 do iq=1,nq
                   zq2(i,l,iq)=zq2(i,l,iq)+zalpha*(zqm(iq)-zq2(i,l,iq)) 
c                  zq2(i,l,iq)=zqm(iq)
                 end do
              ENDDO
c              diag_alpha(i)=zalpha  !temporaire
 
              l2 = l2 + 1
c
            END IF
c
cins$     End Loop
          ENDDO
c
      ENDDO
c
      DO l=1,nlay
        DO ig=1,ngrid
          pdhadj(ig,l)=(zh2(ig,l)-zh(ig,l))/ptimestep
          pduadj(ig,l)=(zu2(ig,l)-zu(ig,l))/ptimestep
          pdvadj(ig,l)=(zv2(ig,l)-zv(ig,l))/ptimestep
c         pdhadj(ig,l)=0.
c         pduadj(ig,l)=0.
c         pdvadj(ig,l)=0.
        ENDDO
      ENDDO

      if(tracer) then 
        do iq=1, nq
          do  l=1,nlay
            do ig=1, ngrid
              pdqadj(ig,l,iq)=(zq2(ig,l,iq)-zq(ig,l,iq))/ptimestep 
            end do 
          end do 
        end do 
      end if
      RETURN
      END
