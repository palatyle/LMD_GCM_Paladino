










      subroutine convadj(ngrid,nlay,nq,ptimestep,
     &                   pplay,pplev,ppopsk,
     &                   pu,pv,ph,pq,
     &                   pdufi,pdvfi,pdhfi,pdqfi,
     &                   pduadj,pdvadj,pdhadj,pdqadj)

      USE tracer_h
      use comcstfi_mod, only: g
      use callkeys_mod, only: tracer,water

      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Calculates dry convective adjustment. If one tracer is CO2,
!     we take into account the molecular mass variation (e.g. when
!     CO2 condenses) to trigger convection (F. Forget 01/2005)
!
!     Authors
!     -------
!     Original author unknown.
!     Modif. 2005 by F. Forget.
!     
!==================================================================

!     ------------
!     Declarations
!     ------------


!     Arguments
!     ---------

      INTEGER ngrid,nlay
      REAL ptimestep
      REAL ph(ngrid,nlay),pdhfi(ngrid,nlay),pdhadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1),ppopsk(ngrid,nlay)
      REAL pu(ngrid,nlay),pdufi(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvfi(ngrid,nlay),pdvadj(ngrid,nlay)

!     Tracers
      integer nq
      real pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq)
      real pdqadj(ngrid,nlay,nq)


!     Local
!     -----

      INTEGER ig,i,l,l1,l2,jj
      INTEGER jcnt, jadrs(ngrid)

      REAL sig(nlay+1),sdsig(nlay),dsig(nlay)
      REAL zu(ngrid,nlay),zv(ngrid,nlay)
      REAL zh(ngrid,nlay)
      REAL zu2(ngrid,nlay),zv2(ngrid,nlay)
      REAL zh2(ngrid,nlay), zhc(ngrid,nlay)
      REAL zhm,zsm,zdsm,zum,zvm,zalpha,zhmc

!     Tracers
      INTEGER iq,ico2
      save ico2
!$OMP THREADPRIVATE(ico2)
      REAL zq(ngrid,nlay,nq), zq2(ngrid,nlay,nq)
      REAL zqm(nq),zqco2m
      real m_co2, m_noco2, A , B
      save A, B
!$OMP THREADPRIVATE(A,B)

      real mtot1, mtot2 , mm1, mm2
       integer l1ref, l2ref
      LOGICAL vtest(ngrid),down,firstcall
      save firstcall
      data firstcall/.true./
!$OMP THREADPRIVATE(firstcall)

!     for conservation test
      real masse,cadjncons

      EXTERNAL SCOPY

!     --------------
!     Initialisation
!     --------------

      IF (firstcall) THEN 
        ico2=0
        if (tracer) then
!     Prepare Special treatment if one of the tracers is CO2 gas
           do iq=1,nq
             if (noms(iq).eq."co2") then
                print*,'dont go there'
!                stop
                ico2=iq
                m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)   
                m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)   
!               Compute A and B coefficient use to compute
!               mean molecular mass Mair defined by
!               1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
!               1/Mair = A*q(ico2) + B
                A =(1/m_co2 - 1/m_noco2)
                B=1/m_noco2
             end if
           enddo
        endif
        firstcall=.false.
      ENDIF ! of IF (firstcall)

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

!     -----------------------------
!     Detection of unstable columns
!     -----------------------------
!     If ph(above) < ph(below) we set vtest=.true.

      DO ig=1,ngrid
        vtest(ig)=.false.
      ENDDO

      if (ico2.ne.0) then
!     Special case if one of the tracers is CO2 gas
         DO l=1,nlay
           DO ig=1,ngrid
             zhc(ig,l) = zh2(ig,l)*(A*zq2(ig,l,ico2)+B)
           ENDDO
         ENDDO
       else
          CALL scopy(ngrid*nlay,zh2,1,zhc,1)
       end if

!     Find out which grid points are convectively unstable
      DO l=2,nlay
        DO ig=1,ngrid
          IF (zhc(ig,l).LT.zhc(ig,l-1)) THEN
            vtest(ig)=.true.
          ENDIF
        ENDDO
      ENDDO
      
!     Make a list of them
      jcnt=0
      DO ig=1,ngrid
         IF(vtest(ig)) THEN
            jcnt=jcnt+1
            jadrs(jcnt)=ig
         ENDIF
      ENDDO


!     ---------------------------------------------------------------
!     Adjustment of the "jcnt" unstable profiles indicated by "jadrs"
!     ---------------------------------------------------------------

      DO jj = 1, jcnt   ! loop on every convective grid point

          i = jadrs(jj)
 
!     Calculate sigma in this column
          DO l=1,nlay+1
            sig(l)=pplev(i,l)/pplev(i,1)
        
          ENDDO
         DO l=1,nlay
            dsig(l)=sig(l)-sig(l+1)
            sdsig(l)=ppopsk(i,l)*dsig(l)
         ENDDO
          l2 = 1

!     Test loop upwards on l2

          DO
            l2 = l2 + 1
            IF (l2 .GT. nlay) EXIT
            IF (zhc(i, l2).LT.zhc(i, l2-1)) THEN
 
!     l2 is the highest level of the unstable column
 
              l1 = l2 - 1
              l  = l1
              zsm = sdsig(l2)
              zdsm = dsig(l2)
              zhm = zh2(i, l2)
              if(ico2.ne.0) zqco2m = zq2(i,l2,ico2)

!     Test loop downwards

              DO
                zsm = zsm + sdsig(l)
                zdsm = zdsm + dsig(l)
                zhm = zhm + sdsig(l) * (zh2(i, l) - zhm) / zsm
                if(ico2.ne.0) then
                  zqco2m = 
     &            zqco2m + dsig(l) * (zq2(i,l,ico2) - zqco2m) / zdsm
                  zhmc = zhm*(A*zqco2m+B)
                else 
                  zhmc = zhm
                end if
 
!     do we have to extend the column downwards?
 
                down = .false.
                IF (l1 .ne. 1) then    !--  and then
                  IF (zhmc.LT.zhc(i, l1-1)) then
                    down = .true.
                  END IF
                END IF
 
                ! this could be a problem...

                if (down) then
 
                  l1 = l1 - 1
                  l  = l1
 
                else
 
!     can we extend the column upwards?
 
                  if (l2 .eq. nlay) exit
 
                  if (zhc(i, l2+1) .ge. zhmc) exit

                  l2 = l2 + 1
                  l  = l2

                end if

              enddo

!     New constant profile (average value)


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
!     modifs by RDW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                zum=zum+dsig(l)*zu2(i,l)
                zvm=zvm+dsig(l)*zv2(i,l)
!                zum=zum+dsig(l)*zu(i,l)
!                zvm=zvm+dsig(l)*zv(i,l)
                do iq=1,nq
                   zqm(iq) = zqm(iq)+dsig(l)*zq2(i,l,iq)
!                   zqm(iq) = zqm(iq)+dsig(l)*zq(i,l,iq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     to conserve tracers/ KE, we must calculate zum, zvm and zqm using 
!     the up-to-date column values. If we do not do this, there are cases 
!     where convection stops at one level and starts at the next where we
!     can break conservation of stuff (particularly tracers) significantly.

                end do
              ENDDO
              zalpha=zalpha/(zhm*(sig(l1)-sig(l2+1)))
              zum=zum/(sig(l1)-sig(l2+1))
              zvm=zvm/(sig(l1)-sig(l2+1))
              do iq=1,nq
                 zqm(iq) = zqm(iq)/(sig(l1)-sig(l2+1))
              end do

              IF(zalpha.GT.1.) THEN
                 zalpha=1.
              ELSE
!                IF(zalpha.LT.0.) STOP
                 IF(zalpha.LT.1.e-4) zalpha=1.e-4
              ENDIF

              DO l=l1,l2
                 zu2(i,l)=zu2(i,l)+zalpha*(zum-zu2(i,l))
                 zv2(i,l)=zv2(i,l)+zalpha*(zvm-zv2(i,l))
                 do iq=1,nq
!                  zq2(i,l,iq)=zq2(i,l,iq)+zalpha*(zqm(iq)-zq2(i,l,iq)) 
                   zq2(i,l,iq)=zqm(iq)
                 end do
              ENDDO
              if (ico2.ne.0) then
                DO l=l1, l2
                  zhc(i,l) = zh2(i,l)*(A*zq2(i,l,ico2)+B)
                ENDDO
              end if


              l2 = l2 + 1

            END IF   ! End of l1 to l2 instability treatment
                     ! We now continue to test from l2 upwards

          ENDDO   ! End of upwards loop on l2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     check conservation
         cadjncons=0.0
         if(water)then
         do l = 1, nlay
            masse = (pplev(i,l) - pplev(i,l+1))/g
            iq    = igcm_h2o_vap
            cadjncons = cadjncons + 
     &           masse*(zq2(i,l,iq)-zq(i,l,iq))/ptimestep 
         end do
         endif

         if(cadjncons.lt.-1.e-6)then
            print*,'convadj has just crashed...'
            print*,'i  = ',i
            print*,'l1 = ',l1
            print*,'l2 = ',l2
            print*,'cadjncons        = ',cadjncons
         do l = 1, nlay
            print*,'dsig         = ',dsig(l)
         end do         
         do l = 1, nlay
            print*,'dsig         = ',dsig(l)
         end do
         do l = 1, nlay
            print*,'sig         = ',sig(l)
         end do
         do l = 1, nlay
            print*,'pplay(ig,:)         = ',pplay(i,l)
         end do
         do l = 1, nlay+1
            print*,'pplev(ig,:)         = ',pplev(i,l)
         end do
         do l = 1, nlay
            print*,'ph(ig,:)         = ',ph(i,l)
         end do
         do l = 1, nlay
            print*,'ph(ig,:)         = ',ph(i,l)
         end do
         do l = 1, nlay
            print*,'ph(ig,:)         = ',ph(i,l)
         end do
         do l = 1, nlay
            print*,'zh(ig,:)         = ',zh(i,l)
         end do
         do l = 1, nlay
            print*,'zh2(ig,:)        = ',zh2(i,l)
         end do
         do l = 1, nlay
            print*,'zq(ig,:,vap)     = ',zq(i,l,igcm_h2o_vap)
         end do
         do l = 1, nlay
            print*,'zq2(ig,:,vap)    = ',zq2(i,l,igcm_h2o_vap)
         end do
            print*,'zqm(vap)         = ',zqm(igcm_h2o_vap)
            print*,'jadrs=',jadrs

            call abort
         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ENDDO

      DO l=1,nlay
        DO ig=1,ngrid
          pdhadj(ig,l)=(zh2(ig,l)-zh(ig,l))/ptimestep
          pduadj(ig,l)=(zu2(ig,l)-zu(ig,l))/ptimestep
          pdvadj(ig,l)=(zv2(ig,l)-zv(ig,l))/ptimestep
        ENDDO
      ENDDO

      if(tracer) then 
        do iq=1, nq
          do  l=1,nlay
            DO ig=1,ngrid
              pdqadj(ig,l,iq)=(zq2(ig,l,iq)-zq(ig,l,iq))/ptimestep 
            end do 
          end do 
        end do 
      end if


!     output
!      if (ngrid.eq.1) then
!         ig=1
!         iq =1 
!         write(*,*)'**** l, pq(ig,l,iq),zq(ig,l,iq),zq2(ig,l,iq)'  
!         do l=nlay,1,-1
!           write(*,*) l, pq(ig,l,iq),zq(ig,l,iq),zq2(ig,l,iq)
!         end do
!      end if


      return
      end
