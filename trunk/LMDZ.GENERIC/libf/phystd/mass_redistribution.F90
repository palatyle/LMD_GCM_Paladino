SUBROUTINE mass_redistribution(ngrid,nlayer,nq,ptimestep,   &
                       rnat,pcapcal,pplay,pplev,pt,ptsrf,pq,pqs,     &
		       pu,pv,pdt,pdtsrf,pdq,pdu,pdv,pdmassmr,  &
		       pdtmr,pdtsrfmr,pdpsrfmr,pdumr,pdvmr,pdqmr,pdqsmr)
                                                   
       USE watercommon_h, Only: Tsat_water,RLVTT
       use surfdat_h
       use radcommon_h, only: glat
       USE tracer_h
       USE planete_mod, only: bp
       use comcstfi_mod, only: g
       USE callkeys_mod, ONLY: water
       
       IMPLICIT NONE
!=======================================================================
!   subject:
!   --------
!     Mass and momentum fluxes through sigma levels as the surface pressure is modified are also taken into account
!
!   author:   Jeremy Leconte 2012 (from F.Forget 1998)
!   ------
!
!   input:
!   ----- 
!    ngrid                 nombre de points de verticales
!                          (toutes les boucles de la physique sont au
!                          moins vectorisees sur ngrid)
!    nlayer                nombre de couches
!    pplay(ngrid,nlayer)   Pressure levels 
!    pplev(ngrid,nlayer+1) Pressure levels 
!    nq                    Number of tracers
!
!    pt(ngrid,nlayer)      temperature (en K)
!    pq(ngrid,nlayer,nq)   tracer specific concentration (kg/kg of air)
!    pu,pv (ngrid,nlayer)  wind velocity (m/s)
!
!                    
!    pdX                   physical tendency of X before mass redistribution
!
!    pdmassmr                air Mass added to the atmosphere in each layer (kg/m2/s)
!
!   output:
!   -------
!
!    pdXmr(ngrid)           physical tendency of X after mass redistribution
!
!    
!
!=======================================================================
!
!    0.  Declarations :
!    ------------------

!-----------------------------------------------------------------------
!    Arguments :
!    ---------
      INTEGER ngrid, nlayer, nq   
      REAL ptimestep
      REAL pcapcal(ngrid)
      INTEGER rnat(ngrid)      
      REAL pplay(ngrid,nlayer),pplev(ngrid,nlayer+1)
      REAL pt(ngrid,nlayer),pdt(ngrid,nlayer)
      REAL ptsrf(ngrid),pdtsrf(ngrid)
      REAL pdtmr(ngrid,nlayer)
      REAL pu(ngrid,nlayer) ,  pv(ngrid,nlayer)
      REAL pdu(ngrid,nlayer) , pdv(ngrid,nlayer)
      REAL pdmassmr(ngrid,nlayer)
      REAL pdumr(ngrid,nlayer) , pdvmr(ngrid,nlayer)
      REAL pq(ngrid,nlayer,nq),pdq(ngrid,nlayer,nq)
      REAL pqs(ngrid,nq)
      REAL pdqmr(ngrid,nlayer,nq),pdqsmr(ngrid,nq)
      REAL pdpsrfmr(ngrid),pdtsrfmr(ngrid)
!
!    Local variables :
!    -----------------

!    Boiling/sublimation
      REAL Tsat(ngrid),zmassboil(ngrid)

!    vertical reorganization of sigma levels
      REAL zzu(nlayer),zzv(nlayer)
      REAL zzq(nlayer,nq),zzt(nlayer)
!    Dummy variables      
      INTEGER n,l,ig,iq
      REAL zdtsig(ngrid,nlayer)
      REAL zmass(ngrid,nlayer),zzmass(nlayer),w(nlayer+1)
      REAL zdmass_sum(ngrid,nlayer+1)
      REAL zmflux(nlayer+1)
      REAL zq1(nlayer)
      REAL ztsrf(ngrid)
      REAL ztm(nlayer+1) 
      REAL zum(nlayer+1) , zvm(nlayer+1)
      REAL zqm(nlayer+1,nq),zqm1(nlayer+1)
      REAL sigma(nlayer+1)

!   local saved variables
      LOGICAL, SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)

!----------------------------------------------------------------------

!   Initialisation
!   --------------
!
      IF (firstcall) THEN
         firstcall=.false.	 
      ENDIF
!
!======================================================================
!    Calcul of h2o condensation  
!    ============================================================
!  
!    Used variable :
!       pdmassmr      : air Mass added to the atmosphere in each layer per unit time (kg/m2/s)
!       zdmass_sum(ngrid,l) : total air mass added to the atm above layer l per unit time (kg/m2/s)
!
!
!     Surface tracer Tendencies set to 0
!     -------------------------------------
      pdqsmr(1:ngrid,1:nq)=0.

      ztsrf(1:ngrid) = ptsrf(1:ngrid) + pdtsrf(1:ngrid)*ptimestep


      DO ig=1,ngrid
         zdmass_sum(ig,nlayer+1)=0.
         DO l = nlayer, 1, -1
           zmass(ig,l) = (pplev(ig,l)-pplev(ig,l+1))/glat(ig)
	   zdmass_sum(ig,l)= zdmass_sum(ig,l+1)+pdmassmr(ig,l)
         END DO
      END DO


      if (water) then
         do ig=1,ngrid
	    call Tsat_water(pplev(ig,1)+zdmass_sum(ig,1)*g*ptimestep,Tsat(ig))
	 enddo
#ifndef MESOSCALE
         call writediagfi(ngrid,'Tsat','saturation temperature at surface','',2,Tsat)
#endif
	 
         do ig=1,ngrid
	    if (ztsrf(ig).gt.Tsat(ig)) then
	       zmassboil(ig)=(ptsrf(ig)-Tsat(ig))*pcapcal(ig)/RLVTT/ptimestep
	       if ((zmassboil(ig)*ptimestep.gt.pqs(ig,igcm_h2o_vap)).and.(rnat(ig).eq.1)) then
	          zmassboil(ig)=pqs(ig,igcm_h2o_vap)/ptimestep
	       endif
	       zmassboil(ig)=zmassboil(ig)*0.0 !momentary, should be 1. JL12
               pdqsmr(ig,igcm_h2o_vap)=-zmassboil(ig)
	       pdtsrfmr(ig)=-zmassboil(ig)*RLVTT/pcapcal(ig)
	       ztsrf(ig)=ptsrf(ig)+pdtsrfmr(ig)*ptimestep
	    else
	       zmassboil(ig)=0.
	       pdtsrfmr(ig)=0.
	    endif
	 enddo
      endif

!     *************************
!           UPDATE SURFACE
!     *************************
!    Changing pressure at the surface:
!    """"""""""""""""""""""""""""""""""""
         
      pdpsrfmr(1:ngrid) = (zdmass_sum(1:ngrid,1)+zmassboil(1:ngrid))*g

      do ig = 1, ngrid
        IF(ABS(pdpsrfmr(ig)*ptimestep).GT.pplev(ig,1)) THEN
         PRINT*,'STOP in condens in mass_redistribution'
         PRINT*,'condensing more than total mass'
         PRINT*,'Grid point ',ig
         PRINT*,'Ps = ',pplev(ig,1)
         PRINT*,'d Ps = ',pdpsrfmr(ig)*ptimestep
         STOP
        ENDIF
      enddo ! of DO ig=1,ngrid


! ***************************************************************
!  Correction to account for redistribution between sigma or hybrid 
!  layers when changing surface pressure
!  zzx quantities have dimension (nlayer) to speed up calculation
! *************************************************************

      DO ig=1,ngrid
         zzt(1:nlayer)  = pt(ig,1:nlayer) + pdt(ig,1:nlayer) * ptimestep
         zzu(1:nlayer)  = pu(ig,1:nlayer) + pdu(ig,1:nlayer) * ptimestep
         zzv(1:nlayer)  = pv(ig,1:nlayer) + pdv(ig,1:nlayer) * ptimestep
         zzq(1:nlayer,1:nq)=pq(ig,1:nlayer,1:nq)+pdq(ig,1:nlayer,1:nq)*ptimestep ! must add the water that has fallen???

!  Mass fluxes of air through the sigma levels (kg.m-2.s-1)  (>0 when up)
!  """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
         zmflux(1) = zmassboil(ig)
         sigma(1)=1
         DO l=1,nlayer
           ! Ehouarn: shouldn't we rather compute sigma levels than use bp()?
!           sigma(l+1)=pplev(ig,l+1)/pplev(ig,1)
!           zmflux(l+1) = zmflux(l) + pdmassmr(ig,l) - &
!                        (sigma(l)-sigma(l+1))*(zdmass_sum(ig,1)+zmflux(1))
!            if (abs(zmflux(l+1)).lt.1E-13.OR.sigma(l+1).eq.0.) zmflux(l+1)=0.
           ! Ehouarn: but for now leave things as before
            zmflux(l+1) = zmflux(l) + pdmassmr(ig,l) - (bp(l)-bp(l+1))*(zdmass_sum(ig,1)+zmflux(1))
! zmflux set to 0 if very low to avoid: top layer is disappearing in v1ld  
            if (abs(zmflux(l+1)).lt.1E-13.OR.bp(l+1).eq.0.) zmflux(l+1)=0.
         END DO

! Mass of each layer
! ------------------ 
         zzmass(1:nlayer)=zmass(ig,1:nlayer)*(1.+pdpsrfmr(ig)*ptimestep/pplev(ig,1))


!  Corresponding fluxes for T,U,V,Q
!  """"""""""""""""""""""""""""""""

!        averaging operator for TRANSPORT  
!        """"""""""""""""""""""""""""""""
!        Value transfert at the surface interface when condensation
!        sublimation:
         ztm(1) = ztsrf(ig)
         zum(1) = 0.  
         zvm(1) = 0.  
         zqm(1,1:nq)=0. ! most tracer do not condense !
         if (water) zqm(1,igcm_h2o_vap)=1. ! flux is 100% h2o at surface
	 
!        Van Leer scheme:
         w(1:nlayer+1)=-zmflux(1:nlayer+1)*ptimestep
         call vl1d(nlayer,zzt,2.,zzmass,w,ztm) 
         call vl1d(nlayer,zzu,2.,zzmass,w,zum) 
         call vl1d(nlayer,zzv,2.,zzmass,w,zvm) 
         do iq=1,nq
           zq1(1:nlayer)=zzq(1:nlayer,iq)
           zqm1(1)=zqm(1,iq)
!		print*,iq
!		print*,zq1
           call vl1d(nlayer,zq1,2.,zzmass,w,zqm1)
           do l=2,nlayer
              zzq(l,iq)=zq1(l)
              zqm(l,iq)=zqm1(l)
           enddo
         enddo

!        Surface condensation affects low winds
         if (zmflux(1).lt.0) then 
            zum(1)= zzu(1) *  (w(1)/zzmass(1))
            zvm(1)= zzv(1) *  (w(1)/zzmass(1))
            if (w(1).gt.zzmass(1)) then ! ensure numerical stability
               zum(1)= (zzu(1)-zum(2))*zzmass(1)/w(1) + zum(2)
               zvm(1)= (zzv(1)-zvm(2))*zzmass(1)/w(1) + zvm(2)
            end if
         end if
                    
         ztm(nlayer+1)= zzt(nlayer) ! should not be used, but... 
         zum(nlayer+1)= zzu(nlayer)  ! should not be used, but...
         zvm(nlayer+1)= zzv(nlayer)  ! should not be used, but...
         zqm(nlayer+1,1:nq)= zzq(nlayer,1:nq)
 
!        Tendencies on T, U, V, Q 
!        """"""""""""""""""""""""
         DO l=1,nlayer
 
!           Tendencies on T
            pdtmr(ig,l) = (1/zzmass(l)) *   &
		(zmflux(l)*(ztm(l) - zzt(l))-zmflux(l+1)*(ztm(l+1)-zzt(l)))
		  !JL12 the last term in Newcondens has been set to zero because we are only dealing with redistribution here

!           Tendencies on U
            pdumr(ig,l)   = (1/zzmass(l)) *( zmflux(l)*(zum(l) - zzu(l)) - zmflux(l+1)*(zum(l+1) - zzu(l)) )

!           Tendencies on V
            pdvmr(ig,l)   = (1/zzmass(l)) *( zmflux(l)*(zvm(l) - zzv(l)) - zmflux(l+1)*(zvm(l+1) - zzv(l)) )

         END DO

!        Tendencies on Q
         do iq=1,nq
            DO l=1,nlayer
               pdqmr(ig,l,iq)= (1/zzmass(l)) *   &
		   (zmflux(l)*(zqm(l,iq)-zzq(l,iq))- zmflux(l+1)*(zqm(l+1,iq)-zzq(l,iq)) - pdmassmr(ig,l)*zzq(l,iq))
            END DO
         enddo

      END DO  ! loop on ig 

CONTAINS

! *****************************************************************
      SUBROUTINE vl1d(llm,q,pente_max,zzmass,w,qm)
!
!    
!     Operateur de moyenne inter-couche pour calcul de transport type
!     Van-Leer " pseudo amont " dans la verticale
!    q,w sont des arguments d'entree  pour le s-pg ....
!    masse : masse de la couche Dp/g
!    w : masse d'atm ``transferee'' a chaque pas de temps (kg.m-2)
!    pente_max = 2 conseillee
!
!
!   --------------------------------------------------------------------
      
      IMPLICIT NONE

!   Arguments:
!   ----------
      integer,intent(in) :: llm
      real zzmass(llm),pente_max
      REAL q(llm),qm(llm+1)
      REAL w(llm+1)
!
!      Local 
!   ---------
!
      INTEGER l
!
      real dzq(llm),dzqw(llm),adzqw(llm),dzqmax
      real sigw, Mtot, MQtot
      integer m 
!     integer ismax,ismin 


!    On oriente tout dans le sens de la pression 
!     W > 0 WHEN DOWN !!!!!!!!!!!!!

      do l=2,llm
            dzqw(l)=q(l-1)-q(l)
            adzqw(l)=abs(dzqw(l))
      enddo

      do l=2,llm-1
            if(dzqw(l)*dzqw(l+1).gt.0.) then
                dzq(l)=0.5*(dzqw(l)+dzqw(l+1))
            else
                dzq(l)=0.
            endif
            dzqmax=pente_max*min(adzqw(l),adzqw(l+1))
            dzq(l)=sign(min(abs(dzq(l)),dzqmax),dzq(l))
      enddo

         dzq(1)=0.
         dzq(llm)=0.

       do l = 1,llm-1

!         Regular scheme (transfered mass < layer mass)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(w(l+1).gt.0. .and. w(l+1).le.zzmass(l+1)) then
             sigw=w(l+1)/zzmass(l+1)
             qm(l+1)=(q(l+1)+0.5*(1.-sigw)*dzq(l+1))
          else if(w(l+1).le.0. .and. -w(l+1).le.zzmass(l)) then
             sigw=w(l+1)/zzmass(l)
             qm(l+1)=(q(l)-0.5*(1.+sigw)*dzq(l))

!         Extended scheme (transfered mass > layer mass)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else if(w(l+1).gt.0.) then
             m=l+1
             Mtot = zzmass(m)
             MQtot = zzmass(m)*q(m)
             do while ((m.lt.llm).and.(w(l+1).gt.(Mtot+zzmass(m+1))))
                m=m+1
                Mtot = Mtot + zzmass(m)
                MQtot = MQtot + zzmass(m)*q(m)
             end do
             if (m.lt.llm) then
                sigw=(w(l+1)-Mtot)/zzmass(m+1)
                qm(l+1)= (1/w(l+1))*(MQtot + (w(l+1)-Mtot)*(q(m+1)+0.5*(1.-sigw)*dzq(m+1)) )
             else
!                w(l+1) = Mtot
!                qm(l+1) = Mqtot / Mtot
                write(*,*) 'top layer is disappearing !',l,Mtot,w(l+1),qm(l+1)
		print*,zzmass
 		print*,w
		print*,q
		print*,qm
               stop
             end if
          else      ! if(w(l+1).lt.0) 
             m = l-1 
             Mtot = zzmass(m+1)
             MQtot = zzmass(m+1)*q(m+1)
             if (m.gt.0) then ! because some compilers will have problems
                              ! evaluating zzmass(0)
              do while ((m.gt.0).and.(-w(l+1).gt.(Mtot+zzmass(m))))
                m=m-1
                Mtot = Mtot + zzmass(m+1)
                MQtot = MQtot + zzmass(m+1)*q(m+1)
                if (m.eq.0) exit
              end do
             endif
             if (m.gt.0) then
                sigw=(w(l+1)+Mtot)/zzmass(m)
                qm(l+1)= (-1/w(l+1))*(MQtot + (-w(l+1)-Mtot)*(q(m)-0.5*(1.+sigw)*dzq(m))  )
             else
                qm(l+1)= (-1/w(l+1))*(MQtot + (-w(l+1)-Mtot)*qm(1))
             end if   
          end if
       enddo

!     boundary conditions (not used in newcondens !!)
!         qm(llm+1)=0.
!         if(w(1).gt.0.) then
!            qm(1)=q(1)
!         else 
!           qm(1)=0.
!         end if

       END SUBROUTINE vl1d

END SUBROUTINE mass_redistribution
