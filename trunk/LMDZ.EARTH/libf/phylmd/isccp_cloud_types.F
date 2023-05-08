!
! $Id $
!
      SUBROUTINE ISCCP_CLOUD_TYPES(
     &     debug,
     &     debugcol,
     &     npoints,
     &     sunlit,
     &     nlev,
     &     ncol,
     &     seed,
     &     pfull,
     &     phalf,
     &     qv,
     &     cc,
     &     conv,
     &     dtau_s,
     &     dtau_c,
     &     top_height,
     &     overlap,
     &     tautab,
     &     invtau,
     &     skt,
     &     emsfc_lw,
     &     at,dem_s,dem_c,
     &     fq_isccp,
     &     totalcldarea,
     &     meanptop,
     &     meantaucld,
     &     boxtau,
     &     boxptop
     &)
        

! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.
!
! This code is available without charge with the following conditions:
!
!  1. The code is available for scientific purposes and is not for 
!     commercial use.
!  2. Any improvements you make to the code should be made available 
!     to the to the authors for incorporation into a future release.
!  3. The code should not be used in any way that brings the authors 
!     or their employers into disrepute.

      implicit none

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolprint
      
!     -----
!     Input 
!     -----

      INTEGER npoints                   !  number of model points in the horizontal
c      PARAMETER(npoints=6722)
      INTEGER nlev                      !  number of model levels in column
      INTEGER ncol                      !  number of subcolumns

      INTEGER sunlit(npoints)           !  1 for day points, 0 for night time

      INTEGER seed(npoints)             !  seed value for random number generator
c                                       !  ( see Numerical Recipes Chapter 7)
c                                       !  It is recommended that the seed is set
c                                       !  to a different value for each model
c                                       !  gridbox it is called on, as it is 
c                                          !  possible that the choice of the samec 
c                                        !  seed value every time may introduce some
c                                        !  statistical bias in the results, particularly
c                                        !  for low values of NCOL.
c
      REAL pfull(npoints,nlev)                      !  pressure of full model levels (Pascals)
c                                        !  pfull(npoints,1)    is    top level of model
c                                        !  pfull(npoints,nlev) is bottom level of model

      REAL phalf(npoints,nlev+1)        !  pressure of half model levels (Pascals)
c                                        !  phalf(npoints,1)    is    top       of model
c                                        !  phalf(npoints,nlev+1) is the surface pressure

      REAL qv(npoints,nlev)             !  water vapor specific humidity (kg vapor/ kg air)
c                                        !         on full model levels

      REAL cc(npoints,nlev)             !  input cloud cover in each model level (fraction) 
c                                        !  NOTE:  This is the HORIZONTAL area of each
c                                        !         grid box covered by clouds

      REAL conv(npoints,nlev)           !  input convective cloud cover in each model level (fraction) 
c                                        !  NOTE:  This is the HORIZONTAL area of each
c                                        !         grid box covered by convective clouds

      REAL dtau_s(npoints,nlev)         !  mean 0.67 micron optical depth of stratiform
c                                        !  clouds in each model level
c                                        !  NOTE:  this the cloud optical depth of only the
c                                        !         cloudy part of the grid box, it is not weighted
c                                        !         with the 0 cloud optical depth of the clear
c                                        !         part of the grid box

      REAL dtau_c(npoints,nlev)         !  mean 0.67 micron optical depth of convective
c                                        !  clouds in each
c                                        !  model level.  Same note applies as in dtau_s.

      INTEGER overlap                   !  overlap type
                                        
!  1=max
                                        
!  2=rand
!  3=max/rand

      INTEGER top_height                !  1 = adjust top height using both a computed
c                                        !  infrared brightness temperature and the visible
c                                        !  optical depth to adjust cloud top pressure. Note
c                                        !  that this calculation is most appropriate to compare
c                                        !  to ISCCP data during sunlit hours.
c                                        !  2 = do not adjust top height, that is cloud top
c                                        !  pressure is the actual cloud top pressure
c                                        !  in the model
c                                        !  3 = adjust top height using only the computed
c                                        !  infrared brightness temperature. Note that this
c                                        !  calculation is most appropriate to compare to ISCCP
c                                        !  IR only algortihm (i.e. you can compare to nighttime
c                                        !  ISCCP data with this option)

      REAL tautab(0:255)                !  ISCCP table for converting count value to 
c                                        !  optical thickness

      INTEGER invtau(-20:45000)         !  ISCCP table for converting optical thickness 
c                                        !  to count value
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL skt(npoints)                 !  skin Temperature (K)
      REAL emsfc_lw                     !  10.5 micron emissivity of surface (fraction)                                            
      REAL at(npoints,nlev)                   !  temperature in each model level (K)
      REAL dem_s(npoints,nlev)                !  10.5 micron longwave emissivity of stratiform
c                                        !  clouds in each
c                                        !  model level.  Same note applies as in dtau_s.
      REAL dem_c(npoints,nlev)                  !  10.5 micron longwave emissivity of convective
c                                        !  clouds in each
c                                        !  model level.  Same note applies as in dtau_s.
cIM reg.dyn BEG
      REAL t1, t2
!     REAL w(npoints)                   !vertical wind at 500 hPa
!     LOGICAL pct_ocean(npoints)        !TRUE if oceanic point, FALSE otherway
!     INTEGER iw(npoints) , nw
!     REAL wmin, pas_w
!     INTEGER k, l, iwmx
!     PARAMETER(wmin=-100.,pas_w=10.,iwmx=30)
!     REAL fq_dynreg(7,7,iwmx)
!     REAL nfq_dynreg(7,7,iwmx) 
!     LOGICAL pctj(7,7,iwmx)
cIM reg.dyn END
!     ------
!     Output
!     ------

      REAL fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
c                                        !  each of the 49 ISCCP D level cloud types

      REAL totalcldarea(npoints)        !  the fraction of model grid box columns
c                                        !  with cloud somewhere in them.  This should
c                                        !  equal the sum over all entries of fq_isccp
        
        
c      ! The following three means are averages over the cloudy areas only.  If no
c      ! clouds are in grid box all three quantities should equal zero.        
                                        
      REAL meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
c                                        !  in cloud top pressure.
                                        
      REAL meantaucld(npoints)          !  mean optical thickness 
c                                        !  linear averaging in albedo performed.
      
      REAL boxtau(npoints,ncol)         !  optical thickness in each column
      
      REAL boxptop(npoints,ncol)        !  cloud top pressure (mb) in each column
                                        
                                                                                                                        
!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
c                                        ! Equivalent of BOX in original version, but
c                                        ! indexed by column then row, rather than
c                                        ! by row then column

      REAL tca(npoints,0:nlev) ! total cloud cover in each model level (fraction)
c                                        ! with extra layer of zeroes on top
c                                        ! in this version this just contains the values input
c                                        ! from cc but with an extra level
      REAL cca(npoints,nlev)         ! convective cloud cover in each model level (fraction)
c                                        ! from conv 

      REAL threshold(npoints,ncol)   ! pointer to position in gridbox
      REAL maxocc(npoints,ncol)      ! Flag for max overlapped conv cld
      REAL maxosc(npoints,ncol)      ! Flag for max overlapped strat cld

      REAL boxpos(npoints,ncol)      ! ordered pointer to position in gridbox

      REAL threshold_min(npoints,ncol) ! minimum value to define range in with new threshold
c                                        ! is chosen

      REAL dem(npoints,ncol),bb(npoints)     !  working variables for 10.5 micron longwave 
c                                        !  emissivity in part of
c                                        !  gridbox under consideration

      REAL ran(npoints)                   ! vector of random numbers
      REAL ptrop(npoints)
      REAL attrop(npoints)
      REAL attropmin (npoints)
      REAL atmax(npoints)
      REAL atmin(npoints)
      REAL btcmin(npoints)
      REAL transmax(npoints)

      INTEGER i,j,ilev,ibox,itrop(npoints)
      INTEGER ipres(npoints)
      INTEGER itau(npoints),ilev2
      INTEGER acc(nlev,ncol)
      INTEGER match(npoints,nlev-1)
      INTEGER nmatch(npoints)
      INTEGER levmatch(npoints,ncol)
      
c      !variables needed for water vapor continuum absorption
      real fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      real taumin(npoints)
      real dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
      real press(npoints), dpress(npoints), atmden(npoints)
      real rvh20(npoints), wk(npoints), rhoave(npoints)
      real rh20s(npoints), rfrgn(npoints)
      real tmpexp(npoints),tauwv(npoints)
      
      character*1 cchar(6),cchar_realtops(6)
      integer icycle
      REAL tau(npoints,ncol)
      LOGICAL box_cloudy(npoints,ncol)
      REAL tb(npoints,ncol)
      REAL ptop(npoints,ncol)
      REAL emcld(npoints,ncol)
      REAL fluxtop(npoints,ncol)
      REAL trans_layers_above(npoints,ncol)
      real isccp_taumin,fluxtopinit(npoints),tauir(npoints)
      real meanalbedocld(npoints) 
      REAL albedocld(npoints,ncol)
      real boxarea
      integer debug       ! set to non-zero value to print out inputs
c                          ! with step debug
      integer debugcol    ! set to non-zero value to print out column
c                          ! decomposition with step debugcol

      integer index1(npoints),num1,jj
      real rec2p13,tauchk

      character*10 ftn09
      
      DATA isccp_taumin / 0.3 /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

      tauchk = -1.*log(0.9999999)
      rec2p13=1./2.13

      ncolprint=0

cIM
c     PRINT*,' isccp_cloud_types npoints=',npoints
c
c      if ( debug.ne.0 ) then
c          j=1
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write(6,'(a10)') 'debug='
c          write(6,'(8I10)') debug
c          write(6,'(a10)') 'debugcol='
c          write(6,'(8I10)') debugcol
c          write(6,'(a10)') 'npoints='
c          write(6,'(8I10)') npoints
c          write(6,'(a10)') 'nlev='
c          write(6,'(8I10)') nlev
c          write(6,'(a10)') 'ncol='
c          write(6,'(8I10)') ncol
c          write(6,'(a10)') 'top_height='
c          write(6,'(8I10)') top_height
c          write(6,'(a10)') 'overlap='
c          write(6,'(8I10)') overlap
c          write(6,'(a10)') 'emsfc_lw='
c          write(6,'(8f10.2)') emsfc_lw
c          write(6,'(a10)') 'tautab='
c          write(6,'(8f10.2)') tautab
c          write(6,'(a10)') 'invtau(1:100)='
c          write(6,'(8i10)') (invtau(i),i=1,100)
c          do j=1,npoints,debug
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write(6,'(a10)') 'sunlit='
c          write(6,'(8I10)') sunlit(j)
c          write(6,'(a10)') 'seed='
c          write(6,'(8I10)') seed(j)
c          write(6,'(a10)') 'pfull='
c          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
c          write(6,'(a10)') 'phalf='
c          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
c          write(6,'(a10)') 'qv='
c          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
c          write(6,'(a10)') 'cc='
c          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
c          write(6,'(a10)') 'conv='
c          write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
c          write(6,'(a10)') 'dtau_s='
c          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
c          write(6,'(a10)') 'dtau_c='
c          write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
c          write(6,'(a10)') 'skt='
c          write(6,'(8f10.2)') skt(j)
c          write(6,'(a10)') 'at='
c          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
c          write(6,'(a10)') 'dem_s='
c          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
c          write(6,'(a10)') 'dem_c='
c          write(6,'(8f10.2)') (dem_c(j,i),i=1,nlev)
c          enddo
c      endif

!     ---------------------------------------------------!

!     assign 2d tca array using 1d input array cc

      do j=1,npoints
        tca(j,0)=0
      enddo
  
      do ilev=1,nlev
        do j=1,npoints
          tca(j,ilev)=cc(j,ilev)
        enddo
      enddo

!     assign 2d cca array using 1d input array conv

      do ilev=1,nlev
cIM pas besoin        do ibox=1,ncol
          do j=1,npoints
            cca(j,ilev)=conv(j,ilev)
          enddo
cIM        enddo
      enddo

cIM
!     do j=1, iwmx
!     do l=1, 7
!     do k=1, 7
!       fq_dynreg(k,l,j) =0. 
!       nfq_dynreg(k,l,j) =0. 
!      enddo !k
!     enddo !l
!     enddo !j
cIM
cIM
c      if (ncolprint.ne.0) then
c        do j=1,npoints,1000
c        write(6,'(a10)') 'j='
c        write(6,'(8I10)') j
c        write (6,'(a)') 'seed:'
c        write (6,'(I3.2)') seed(j)

c        write (6,'(a)') 'tca_pp_rev:'
c        write (6,'(8f5.2)') 
c     &   ((tca(j,ilev)),
c     &      ilev=1,nlev)

c        write (6,'(a)') 'cca_pp_rev:'
c        write (6,'(8f5.2)') 
c     &   ((cca(j,ilev),ibox=1,ncolprint),ilev=1,nlev)
c        enddo
c      endif

      if (top_height .eq. 1 .or. top_height .eq. 3) then 

      do j=1,npoints 
          ptrop(j)=5000.
          atmin(j) = 400.
          attropmin(j) = 400.
          atmax(j) = 0.
          attrop(j) = 120.
          itrop(j) = 1
      enddo 

      do 12 ilev=1,nlev
        do j=1,npoints 
         if (pfull(j,ilev) .lt. 40000. .and.
     &          pfull(j,ilev) .gt.  5000. .and.
     &          at(j,ilev) .lt. attropmin(j)) then
                ptrop(j) = pfull(j,ilev)
                attropmin(j) = at(j,ilev)
                attrop(j) = attropmin(j)
                itrop(j)=ilev
           end if
           if (at(j,ilev) .gt. atmax(j)) atmax(j)=at(j,ilev)
           if (at(j,ilev) .lt. atmin(j)) atmin(j)=at(j,ilev)
        enddo
12    continue

      end if

!     -----------------------------------------------------!

!     ---------------------------------------------------!

cIM
c     do 13 ilev=1,nlev
cnum1=0
c       do j=1,npoints
c           if (cc(j,ilev) .lt. 0. .or. cc(j,ilev) .gt. 1.) then
c        num1=num1+1
c        index1(num1)=j
c           end if
c       enddo
c       do jj=1,num1   
c        j=index1(jj)
c               write(6,*)  ' error = cloud fraction less than zero'
c        write(6,*) ' or '
c               write(6,*)  ' error = cloud fraction greater than 1'
c        write(6,*) 'value at point ',j,' is ',cc(j,ilev)
c        write(6,*) 'level ',ilev
c               STOP
c       enddo
cnum1=0
c       do j=1,npoints
c           if (conv(j,ilev) .lt. 0. .or. conv(j,ilev) .gt. 1.) then
c        num1=num1+1
c        index1(num1)=j
c           end if
c       enddo
c       do jj=1,num1   
c        j=index1(jj)
c               write(6,*)  
c    &           ' error = convective cloud fraction less than zero'
c        write(6,*) ' or '
c               write(6,*)  
c    &           ' error = convective cloud fraction greater than 1'
c        write(6,*) 'value at point ',j,' is ',conv(j,ilev)
c        write(6,*) 'level ',ilev
c               STOP
c       enddo

cnum1=0
c       do j=1,npoints
c           if (dtau_s(j,ilev) .lt. 0.) then
c        num1=num1+1
c        index1(num1)=j
c           end if
c       enddo
c       do jj=1,num1   
c        j=index1(jj)
c               write(6,*)  
c    &           ' error = stratiform cloud opt. depth less than zero'
c        write(6,*) 'value at point ',j,' is ',dtau_s(j,ilev)
c        write(6,*) 'level ',ilev
c               STOP
c       enddo
cnum1=0
c       do j=1,npoints
c           if (dtau_c(j,ilev) .lt. 0.) then
c        num1=num1+1
c        index1(num1)=j
c           end if
c       enddo
c       do jj=1,num1   
c        j=index1(jj)
c               write(6,*)  
c    &           ' error = convective cloud opt. depth less than zero'
c        write(6,*) 'value at point ',j,' is ',dtau_c(j,ilev)
c        write(6,*) 'level ',ilev
c               STOP
c       enddo

cnum1=0
c       do j=1,npoints
c           if (dem_s(j,ilev) .lt. 0. .or. dem_s(j,ilev) .gt. 1.) then
c        num1=num1+1
c        index1(num1)=j
c           end if
c       enddo
c       do jj=1,num1   
c        j=index1(jj)
c               write(6,*)  
c    &           ' error = stratiform cloud emissivity less than zero'
c        write(6,*)'or'
c               write(6,*)  
c    &           ' error = stratiform cloud emissivity greater than 1'
c        write(6,*) 'value at point ',j,' is ',dem_s(j,ilev)
c        write(6,*) 'level ',ilev
c               STOP
c       enddo

cnum1=0
c       do j=1,npoints
c           if (dem_c(j,ilev) .lt. 0. .or. dem_c(j,ilev) .gt. 1.) then
c        num1=num1+1
c        index1(num1)=j
c           end if
c       enddo
c       do jj=1,num1   
c        j=index1(jj)
c               write(6,*)  
c    &           ' error = convective cloud emissivity less than zero'
c        write(6,*)'or'
c               write(6,*)  
c    &           ' error = convective cloud emissivity greater than 1'
c               write (6,*) 
c    &          'j=',j,'ilev=',ilev,'dem_c(j,ilev) =',dem_c(j,ilev) 
c               STOP
c       enddo
c13    continue


      do ibox=1,ncol
        do j=1,npoints 
          boxpos(j,ibox)=(ibox-.5)/ncol
        enddo
      enddo

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero

      do ilev=1,nlev
        do ibox=1,ncol
          do j=1,npoints
            frac_out(j,ibox,ilev)=0.0
          enddo
        enddo
      enddo

cIM
c      if (ncolprint.ne.0) then
c        write (6,'(a)') 'frac_out_pp_rev:'
c          do j=1,npoints,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(8f5.2)') 
c     &     ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

c          enddo
c        write (6,'(a)') 'ncol:'
c        write (6,'(I3)') ncol
c      endif
c      if (ncolprint.ne.0) then
c        write (6,'(a)') 'last_frac_pp:'
c          do j=1,npoints,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(8f5.2)') (tca(j,0))
c          enddo
c      endif

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information 
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud
      
      !loop over vertical levels
      DO 200 ilev = 1,nlev
                                  
!     Initialise threshold

        IF (ilev.eq.1) then
            ! If max overlap 
            IF (overlap.eq.1) then
              ! select pixels spread evenly
              ! across the gridbox
              DO ibox=1,ncol
                do j=1,npoints
                  threshold(j,ibox)=boxpos(j,ibox)
                enddo
              enddo
            ELSE
              DO ibox=1,ncol
                call ran0_vec(npoints,seed,ran)
                ! select random pixels from the non-convective
                ! part the gridbox ( some will be converted into
                ! convective pixels below )
                do j=1,npoints
                  threshold(j,ibox)=
     &            cca(j,ilev)+(1-cca(j,ilev))*ran(j)
                enddo
              enddo
            ENDIF
cIM
c            IF (ncolprint.ne.0) then
c              write (6,'(a)') 'threshold_nsf2:'
c                do j=1,npoints,1000
c                write(6,'(a10)') 'j='
c                write(6,'(8I10)') j
c                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
c                enddo
c            ENDIF
        ENDIF

c        IF (ncolprint.ne.0) then
c            write (6,'(a)') 'ilev:'
c            write (6,'(I2)') ilev
c        ENDIF

        DO ibox=1,ncol

          ! All versions
          do j=1,npoints
            if (boxpos(j,ibox).le.cca(j,ilev)) then
cIM REAL           maxocc(j,ibox) = 1
              maxocc(j,ibox) = 1.0
            else
cIM REAL           maxocc(j,ibox) = 0
              maxocc(j,ibox) = 0.0
            end if
          enddo

          ! Max overlap
          if (overlap.eq.1) then 
            do j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
cIM REAL           maxosc(j,ibox)=1
              maxosc(j,ibox)=1.0
            enddo
          endif

          ! Random overlap
          if (overlap.eq.2) then 
            do j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
cIM REAL           maxosc(j,ibox)=0
              maxosc(j,ibox)=0.0
            enddo
          endif

          ! Max/Random overlap
          if (overlap.eq.3) then 
            do j=1,npoints
              threshold_min(j,ibox)=max(cca(j,ilev),
     &          min(tca(j,ilev-1),tca(j,ilev)))
              if (threshold(j,ibox)
     &          .lt.min(tca(j,ilev-1),tca(j,ilev))
     &          .and.(threshold(j,ibox).gt.cca(j,ilev))) then
cIM REAL                maxosc(j,ibox)= 1
                   maxosc(j,ibox)= 1.0
              else
cIM REAL                 maxosc(j,ibox)= 0
                   maxosc(j,ibox)= 0.0
              end if
            enddo
          endif
    
          ! Reset threshold 
          call ran0_vec(npoints,seed,ran)
           
          do j=1,npoints
            threshold(j,ibox)=
              !if max overlapped conv cloud
     &        maxocc(j,ibox) * (                                       
     &            boxpos(j,ibox)                                               
     &        ) +                                                      
              !else
     &        (1-maxocc(j,ibox)) * (                                   
                  !if max overlapped strat cloud
     &            (maxosc(j,ibox)) * (                                 
                      !threshold=boxpos
     &                threshold(j,ibox)                                        
     &            ) +                                                  
                  !else
     &            (1-maxosc(j,ibox)) * (                               
                      !threshold_min=random[thrmin,1]
     &                threshold_min(j,ibox)+
     &                  (1-threshold_min(j,ibox))*ran(j)  
     &           ) 
     &        )
          enddo

        ENDDO ! ibox

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
             do j=1,npoints 
               if (tca(j,ilev).gt.threshold(j,ibox)) then
cIM REAL             frac_out(j,ibox,ilev)=1
               frac_out(j,ibox,ilev)=1.0
               else
cIM REAL             frac_out(j,ibox,ilev)=0
               frac_out(j,ibox,ilev)=0.0
               end if               
             enddo
           ENDDO

!           Code to partition boxes into startiform and convective parts
!           goes here

           DO ibox=1,ncol
             do j=1,npoints 
                if (threshold(j,ibox).le.cca(j,ilev)) then
                    ! = 2 IF threshold le cca(j)
cIM REAL                  frac_out(j,ibox,ilev) = 2 
                    frac_out(j,ibox,ilev) = 2.0 
                else
                    ! = the same IF NOT threshold le cca(j) 
                    frac_out(j,ibox,ilev) = frac_out(j,ibox,ilev)
                end if
             enddo
           ENDDO

!         Set last_frac to tca at this level, so as to be tca 
!         from last level next time round

cIM
c          if (ncolprint.ne.0) then

c            do j=1,npoints ,1000
c            write(6,'(a10)') 'j='
c            write(6,'(8I10)') j
c            write (6,'(a)') 'last_frac:'
c            write (6,'(8f5.2)') (tca(j,ilev-1))
    
c            write (6,'(a)') 'cca:'
c            write (6,'(8f5.2)') (cca(j,ilev),ibox=1,ncolprint)
    
c            write (6,'(a)') 'max_overlap_cc:'
c            write (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)
    
c            write (6,'(a)') 'max_overlap_sc:'
c            write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
    
c            write (6,'(a)') 'threshold_min_nsf2:'
c            write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
    
c            write (6,'(a)') 'threshold_nsf2:'
c            write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
    
c            write (6,'(a)') 'frac_out_pp_rev:'
c            write (6,'(8f5.2)') 
c     &       ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
c            enddo
c          endif

200   CONTINUE    !loop over nlev

!
!     ---------------------------------------------------!

      
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      do 15 ibox=1,ncol
        do j=1,npoints 
            tau(j,ibox)=0.
            albedocld(j,ibox)=0.
            boxtau(j,ibox)=0.
            boxptop(j,ibox)=0.
            box_cloudy(j,ibox)=.false.
        enddo
15    continue

      !compute total cloud optical depth for each column     
      do ilev=1,nlev
            !increment tau for each of the boxes
            do ibox=1,ncol
              do j=1,npoints 
cIM REAL              if (frac_out(j,ibox,ilev).eq.1) then
                 if (frac_out(j,ibox,ilev).eq.1.0) then
                        tau(j,ibox)=tau(j,ibox)
     &                     + dtau_s(j,ilev)
                 endif
cIM REAL              if (frac_out(j,ibox,ilev).eq.2) then
                 if (frac_out(j,ibox,ilev).eq.2.0) then
                        tau(j,ibox)=tau(j,ibox)
     &                     + dtau_c(j,ilev)
                 end if
              enddo
            enddo ! ibox
      enddo ! ilev
cIM
c          if (ncolprint.ne.0) then

c              do j=1,npoints ,1000
c                write(6,'(a10)') 'j='
c                write(6,'(8I10)') j
c                write(6,'(i2,1X,8(f7.2,1X))') 
c     &          ilev,
c     &          (tau(j,ibox),ibox=1,ncolprint)
c              enddo
c          endif 
!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      if (top_height .eq. 1 .or. top_height .eq. 3) then


        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644
        wtmh20 = 18.01534
        Navo = 6.023E+23
        grav = 9.806650E+02
        pstd = 1.013250E+06
        t0 = 296.
cIM
c        if (ncolprint .ne. 0) 
c     &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
        do 125 ilev=1,nlev
          do j=1,npoints 
               !press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10.
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
               !atmden = g/cm2 = kg/m2 / 10 
               atmden(j) = dpress(j)/grav
               rvh20(j) = qv(j,ilev)*wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
               tmpexp(j) = exp(-0.02*(at(j,ilev)-t0))
               tauwv(j) = wk(j)*1.e-20*( 
     &           (0.0224697*rh20s(j)*tmpexp(j)) + 
     &                (3.41817e-7*rfrgn(j)) )*0.98
               dem_wv(j,ilev) = 1. - exp( -1. * tauwv(j))
          enddo
cIM
c               if (ncolprint .ne. 0) then
c               do j=1,npoints ,1000
c               write(6,'(a10)') 'j='
c               write(6,'(8I10)') j
c               write(6,'(i2,1X,3(f8.3,3X))') ilev,
c     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100.),
c     &           tauwv(j),dem_wv(j,ilev)
c               enddo
c               endif
125     continue

        !initialize variables
        do j=1,npoints 
          fluxtop_clrsky(j) = 0.
          trans_layers_above_clrsky(j)=1.
        enddo

        do ilev=1,nlev
          do j=1,npoints 
 
            ! Black body emission at temperature of the layer

                bb(j)=1 / ( exp(1307.27/at(j,ilev)) - 1. )
                !bb(j)= 5.67e-8*at(j,ilev)**4

                ! increase TOA flux by flux emitted from layer
                ! times total transmittance in layers above

                fluxtop_clrsky(j) = fluxtop_clrsky(j) 
     &            + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j) 
            
                ! update trans_layers_above with transmissivity
                ! from this layer for next time around loop

                trans_layers_above_clrsky(j)=
     &            trans_layers_above_clrsky(j)*(1.-dem_wv(j,ilev))
                   

          enddo 
cIM  
c            if (ncolprint.ne.0) then
c             do j=1,npoints ,1000
c              write(6,'(a10)') 'j='
c              write(6,'(8I10)') j
c              write (6,'(a)') 'ilev:'
c              write (6,'(I2)') ilev
    
c              write (6,'(a)') 
c     &        'emiss_layer,100.*bb(j),100.*f,total_trans:'
c              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100.*bb(j),
c     &             100.*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
c             enddo   
c            endif

        enddo   !loop over level
        
        do j=1,npoints 
          !add in surface emission
          bb(j)=1/( exp(1307.27/skt(j)) - 1. )
          !bb(j)=5.67e-8*skt(j)**4

          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j) 
     &     * trans_layers_above_clrsky(j)
        enddo

cIM
c        if (ncolprint.ne.0) then
c        do j=1,npoints ,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(a)') 'id:'
c          write (6,'(a)') 'surface'

c          write (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
c          write (6,'(4(f7.2,1X))') emsfc_lw,100.*bb(j),
c     &      100.*fluxtop_clrsky(j),
c     &       trans_layers_above_clrsky(j)
c        enddo
c        endif
    

        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------


cIM
c        if (ncolprint.ne.0) then

c        do j=1,npoints ,1000
c            write(6,'(a10)') 'j='
c            write(6,'(8I10)') j
c            write (6,'(a)') 'ts:'
c            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
    
c            write (6,'(a)') 'ta_rev:'
c            write (6,'(8f7.2)') 
c     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)

c        enddo
c        endif 
        !loop over columns 
        do ibox=1,ncol
          do j=1,npoints
            fluxtop(j,ibox)=0.
            trans_layers_above(j,ibox)=1.
          enddo
        enddo

        do ilev=1,nlev
              do j=1,npoints 
                ! Black body emission at temperature of the layer

                bb(j)=1 / ( exp(1307.27/at(j,ilev)) - 1. )
                !bb(j)= 5.67e-8*at(j,ilev)**4
              enddo

            do ibox=1,ncol
              do j=1,npoints 

                ! emissivity for point in this layer
cIM REAL             if (frac_out(j,ibox,ilev).eq.1) then
                if (frac_out(j,ibox,ilev).eq.1.0) then
                dem(j,ibox)= 1. - 
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_s(j,ilev)) )
cIM REAL             else if (frac_out(j,ibox,ilev).eq.2) then
                else if (frac_out(j,ibox,ilev).eq.2.0) then
                dem(j,ibox)= 1. - 
     &             ( (1. - dem_wv(j,ilev)) * (1. -  dem_c(j,ilev)) )
                else
                dem(j,ibox)=  dem_wv(j,ilev)
                end if
                

                ! increase TOA flux by flux emitted from layer
                ! times total transmittance in layers above

                fluxtop(j,ibox) = fluxtop(j,ibox) 
     &            + dem(j,ibox) * bb(j)
     &            * trans_layers_above(j,ibox) 
            
                ! update trans_layers_above with transmissivity
                ! from this layer for next time around loop

                trans_layers_above(j,ibox)=
     &            trans_layers_above(j,ibox)*(1.-dem(j,ibox))

              enddo ! j
            enddo ! ibox

cIM
c            if (ncolprint.ne.0) then
c              do j=1,npoints,1000
c              write (6,'(a)') 'ilev:'
c              write (6,'(I2)') ilev
    
c              write(6,'(a10)') 'j='
c              write(6,'(8I10)') j
c              write (6,'(a)') 'emiss_layer:'
c              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
        
c              write (6,'(a)') '100.*bb(j):'
c              write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
        
c              write (6,'(a)') '100.*f:'
c              write (6,'(8f7.2)') 
c     &         (100.*fluxtop(j,ibox),ibox=1,ncolprint)
        
c              write (6,'(a)') 'total_trans:'
c              write (6,'(8f7.2)') 
c     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
c              enddo
c          endif

        enddo ! ilev


          do j=1,npoints 
            !add in surface emission
            bb(j)=1/( exp(1307.27/skt(j)) - 1. )
            !bb(j)=5.67e-8*skt(j)**4
          end do

        do ibox=1,ncol
          do j=1,npoints 

            !add in surface emission

            fluxtop(j,ibox) = fluxtop(j,ibox) 
     &         + emsfc_lw * bb(j) 
     &         * trans_layers_above(j,ibox) 
            
          end do
        end do

cIM
c        if (ncolprint.ne.0) then

c          do j=1,npoints ,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j
c          write (6,'(a)') 'id:'
c          write (6,'(a)') 'surface'

c          write (6,'(a)') 'emiss_layer:'
c          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') '100.*bb(j):'
c          write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
    
c          write (6,'(a)') '100.*f:'
c          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
c          end do
c        endif
    
        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature

        !account for partially transmitting cloud recompute flux 
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the 
        !adjustment performed in this section
        !
        !If it turns out that the cloud brightness temperature
        !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
        !
        !Note that this is discussed on pages 85-87 of 
        !the ISCCP D level documentation (Rossow et al. 1996)
           
          do j=1,npoints  
            !compute minimum brightness temperature and optical depth
            btcmin(j) = 1. /  ( exp(1307.27/(attrop(j)-5.)) - 1. ) 
          enddo 
        do ibox=1,ncol
          do j=1,npoints  
            transmax(j) = (fluxtop(j,ibox)-btcmin(j))
     &                /(fluxtop_clrsky(j)-btcmin(j))
            !note that the initial setting of tauir(j) is needed so that
            !tauir(j) has a realistic value should the next if block be
            !bypassed
            tauir(j) = tau(j,ibox) * rec2p13
            taumin(j) = -1. * log(max(min(transmax(j),0.9999999),0.001))

          enddo 

          if (top_height .eq. 1) then
            do j=1,npoints  
              if (transmax(j) .gt. 0.001 .and. 
     &          transmax(j) .le. 0.9999999) then
                fluxtopinit(j) = fluxtop(j,ibox)
                tauir(j) = tau(j,ibox) *rec2p13
              endif
            enddo
            do icycle=1,2
              do j=1,npoints  
                if (tau(j,ibox) .gt. (tauchk            )) then 
                if (transmax(j) .gt. 0.001 .and. 
     &            transmax(j) .le. 0.9999999) then
                  emcld(j,ibox) = 1. - exp(-1. * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -   
     &              ((1.-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.E-06,
     &              (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= 1307.27
     &              / (log(1. + (1./fluxtop(j,ibox))))
                  if (tb(j,ibox) .gt. 260.) then
                    tauir(j) = tau(j,ibox) / 2.56
                  end if                         
                end if
                end if
              enddo
            enddo
                
          endif
        
          do j=1,npoints
            if (tau(j,ibox) .gt. (tauchk            )) then 
                !cloudy box
                tb(j,ibox)= 1307.27/ (log(1. + (1./fluxtop(j,ibox))))
                if (top_height.eq.1.and.tauir(j).lt.taumin(j)) then
                         tb(j,ibox) = attrop(j) - 5. 
                         tau(j,ibox) = 2.13*taumin(j)
                end if
            else
                !clear sky brightness temperature
                tb(j,ibox) = 1307.27/(log(1.+(1./fluxtop_clrsky(j))))
            end if
          enddo ! j
        enddo ! ibox

cIM
c        if (ncolprint.ne.0) then

c          do j=1,npoints,1000
c          write(6,'(a10)') 'j='
c          write(6,'(8I10)') j

c          write (6,'(a)') 'attrop:'
c          write (6,'(8f7.2)') (attrop(j))
    
c          write (6,'(a)') 'btcmin:'
c          write (6,'(8f7.2)') (btcmin(j))
    
c          write (6,'(a)') 'fluxtop_clrsky*100:'
c          write (6,'(8f7.2)') 
c     &      (100.*fluxtop_clrsky(j))

c          write (6,'(a)') '100.*f_adj:'
c          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'transmax:'
c          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'tau:'
c          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'emcld:'
c          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'total_trans:'
c          write (6,'(8f7.2)') 
c     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'total_emiss:'
c          write (6,'(8f7.2)') 
c     &          (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'total_trans:'
c          write (6,'(8f7.2)') 
c     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
    
c          write (6,'(a)') 'ppout:'
c          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
c          enddo ! j
c        endif

      end if

!     ---------------------------------------------------!

!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      do 30 ibox=1,ncol
        !segregate according to optical thickness
        if (top_height .eq. 1 .or. top_height .eq. 3) then  
          !find level whose temperature
          !most closely matches brightness temperature
          do j=1,npoints 
            nmatch(j)=0
          enddo
          do 29 ilev=1,nlev-1
!cdir nodep
            do j=1,npoints 
              if ((at(j,ilev)   .ge. tb(j,ibox) .and. 
     &          at(j,ilev+1) .lt. tb(j,ibox)) .or.
     &          (at(j,ilev) .le. tb(j,ibox) .and. 
     &          at(j,ilev+1) .gt. tb(j,ibox))) then 
   
                nmatch(j)=nmatch(j)+1
                if(abs(at(j,ilev)-tb(j,ibox)) .lt.
     &            abs(at(j,ilev+1)-tb(j,ibox))) then
                  match(j,nmatch(j))=ilev
                else
                  match(j,nmatch(j))=ilev+1
                end if
              end if                        
            enddo
29        continue

          do j=1,npoints 
            if (nmatch(j) .ge. 1) then
              ptop(j,ibox)=pfull(j,match(j,nmatch(j)))
              levmatch(j,ibox)=match(j,nmatch(j))   
            else
              if (tb(j,ibox) .lt. atmin(j)) then
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              end if
              if (tb(j,ibox) .gt. atmax(j)) then
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              end if                                
            end if
          enddo ! j

        else ! if (top_height .eq. 1 .or. top_height .eq. 3) 
 
          do j=1,npoints     
            ptop(j,ibox)=0.
          enddo
          do ilev=1,nlev
            do j=1,npoints     
              if ((ptop(j,ibox) .eq. 0. )
cIM  &           .and.(frac_out(j,ibox,ilev) .ne. 0)) then
     &           .and.(frac_out(j,ibox,ilev) .ne. 0.0)) then
                ptop(j,ibox)=pfull(j,ilev)
                levmatch(j,ibox)=ilev
              end if
            end do
          end do
        end if                            
          
        do j=1,npoints
          if (tau(j,ibox) .le. (tauchk            )) then
            ptop(j,ibox)=0.
            levmatch(j,ibox)=0      
          endif 
        enddo

30    continue
              
!
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy 
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.  
!

      !compute isccp frequencies

      !reset frequencies
      do 38 ilev=1,7
      do 38 ilev2=1,7
        do j=1,npoints ! 
             fq_isccp(j,ilev,ilev2)=0.
        enddo
38    continue

      !reset variables need for averaging cloud properties
      do j=1,npoints 
        totalcldarea(j) = 0.
        meanalbedocld(j) = 0.
        meanptop(j) = 0.
        meantaucld(j) = 0.
      enddo ! j

      boxarea = 1./real(ncol)
     
              !determine optical depth category
cIM       do 39 j=1,npoints
cIM       do ibox=1,ncol
        do 39 ibox=1,ncol
          do j=1,npoints

cIM
c         CALL CPU_time(t1)
cIM

          if (tau(j,ibox) .gt. (tauchk            )
     &      .and. ptop(j,ibox) .gt. 0.) then
              box_cloudy(j,ibox)=.true.
          endif

cIM
c         CALL CPU_time(t2)
c         print*,'IF tau t2 - t1',t2 - t1

c         CALL CPU_time(t1)
cIM

          if (box_cloudy(j,ibox)) then

              ! totalcldarea always diagnosed day or night
              totalcldarea(j) = totalcldarea(j) + boxarea

              if (sunlit(j).eq.1) then

                ! tau diagnostics only with sunlight

                boxtau(j,ibox) = tau(j,ibox)

                !convert optical thickness to albedo
                  albedocld(j,ibox)
     &            =real(invtau(min(nint(100.*tau(j,ibox)),45000)))
            
                !contribute to averaging
                meanalbedocld(j) = meanalbedocld(j) 
     &            +albedocld(j,ibox)*boxarea

            endif

          endif
          
cIM
c         CALL CPU_time(t2)
c         print*,'IF box_cloudy t2 - t1',t2 - t1
          
c         CALL CPU_time(t1)
cIM BEG 
cIM           !convert ptop to millibars
              ptop(j,ibox)=ptop(j,ibox) / 100.
            
cIM           !save for output cloud top pressure and optical thickness
              boxptop(j,ibox) = ptop(j,ibox)
cIM END

cIM BEG
              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0

              if (tau(j,ibox) .lt. isccp_taumin) then
                  itau(j)=1
              else if (tau(j,ibox) .ge. isccp_taumin
     &                                    
     &          .and. tau(j,ibox) .lt. 1.3) then
                itau(j)=2
              else if (tau(j,ibox) .ge. 1.3 
     &          .and. tau(j,ibox) .lt. 3.6) then
                itau(j)=3
              else if (tau(j,ibox) .ge. 3.6 
     &          .and. tau(j,ibox) .lt. 9.4) then
                  itau(j)=4
              else if (tau(j,ibox) .ge. 9.4 
     &          .and. tau(j,ibox) .lt. 23.) then
                  itau(j)=5
              else if (tau(j,ibox) .ge. 23. 
     &          .and. tau(j,ibox) .lt. 60.) then
                  itau(j)=6
              else if (tau(j,ibox) .ge. 60.) then
                  itau(j)=7
              end if

              !determine cloud top pressure category
              if (    ptop(j,ibox) .gt. 0.  
     &          .and.ptop(j,ibox) .lt. 180.) then
                  ipres(j)=1
              else if(ptop(j,ibox) .ge. 180.
     &          .and.ptop(j,ibox) .lt. 310.) then
                  ipres(j)=2
              else if(ptop(j,ibox) .ge. 310.
     &          .and.ptop(j,ibox) .lt. 440.) then
                  ipres(j)=3
              else if(ptop(j,ibox) .ge. 440.
     &          .and.ptop(j,ibox) .lt. 560.) then
                  ipres(j)=4
              else if(ptop(j,ibox) .ge. 560.
     &          .and.ptop(j,ibox) .lt. 680.) then
                  ipres(j)=5
              else if(ptop(j,ibox) .ge. 680.
     &          .and.ptop(j,ibox) .lt. 800.) then
                  ipres(j)=6
              else if(ptop(j,ibox) .ge. 800.) then
                  ipres(j)=7
              end if 
cIM END

          if (sunlit(j).eq.1 .or. top_height .eq. 3) then 

cIM         !convert ptop to millibars
cIM           ptop(j,ibox)=ptop(j,ibox) / 100.
            
cIM         !save for output cloud top pressure and optical thickness
cIM             boxptop(j,ibox) = ptop(j,ibox)
    
            if (box_cloudy(j,ibox)) then
            
              meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea

cIM             !reset itau(j), ipres(j)
cIM           itau(j) = 0
cIM           ipres(j) = 0

c             if (tau(j,ibox) .lt. isccp_taumin) then
c                 itau(j)=1
c             else if (tau(j,ibox) .ge. isccp_taumin
c    &                                    
c    &          .and. tau(j,ibox) .lt. 1.3) then
c               itau(j)=2
c             else if (tau(j,ibox) .ge. 1.3 
c    &          .and. tau(j,ibox) .lt. 3.6) then
c               itau(j)=3
c             else if (tau(j,ibox) .ge. 3.6 
c    &          .and. tau(j,ibox) .lt. 9.4) then
c                 itau(j)=4
c             else if (tau(j,ibox) .ge. 9.4 
c    &          .and. tau(j,ibox) .lt. 23.) then
c                 itau(j)=5
c             else if (tau(j,ibox) .ge. 23. 
c    &          .and. tau(j,ibox) .lt. 60.) then
c                 itau(j)=6
c             else if (tau(j,ibox) .ge. 60.) then
c                 itau(j)=7
c             end if

c             !determine cloud top pressure category
c             if (    ptop(j,ibox) .gt. 0.  
c    &          .and.ptop(j,ibox) .lt. 180.) then
c                 ipres(j)=1
c             else if(ptop(j,ibox) .ge. 180.
c    &          .and.ptop(j,ibox) .lt. 310.) then
c                 ipres(j)=2
c             else if(ptop(j,ibox) .ge. 310.
c    &          .and.ptop(j,ibox) .lt. 440.) then
c                 ipres(j)=3
c            else if(ptop(j,ibox) .ge. 440.
c    &          .and.ptop(j,ibox) .lt. 560.) then
c                 ipres(j)=4
c             else if(ptop(j,ibox) .ge. 560.
c    &          .and.ptop(j,ibox) .lt. 680.) then
c                 ipres(j)=5
c             else if(ptop(j,ibox) .ge. 680.
c    &          .and.ptop(j,ibox) .lt. 800.) then
c                 ipres(j)=6
c             else if(ptop(j,ibox) .ge. 800.) then
c                 ipres(j)=7
c             end if 

              !update frequencies
              if(ipres(j) .gt. 0.and.itau(j) .gt. 0) then
              fq_isccp(j,itau(j),ipres(j))=
     &          fq_isccp(j,itau(j),ipres(j))+ boxarea
              end if

cIM calcul stats regime dynamique BEG
!             iw(j) = int((w(j)-wmin)/pas_w) +1
!             pctj(itau(j),ipres(j),iw(j))=.FALSE.
!             !update frequencies W500
!             if (pct_ocean(j)) then
!             if (ipres(j) .gt. 0.and.itau(j) .gt. 0) then
!             if (iw(j) .gt. int(wmin).and.iw(j) .le. iwmx) then
c             print*,' ISCCP iw=',iw(j),j
!             fq_dynreg(itau(j),ipres(j),iw(j))=
!    &          fq_dynreg(itau(j),ipres(j),iw(j))+
!    &          boxarea
c    &          fq_isccp(j,itau(j),ipres(j))
!             pctj(itau(j),ipres(j),iw(j))=.TRUE.
c             nfq_dynreg(itau(j),ipres(j),iw(j))=
c    &          nfq_dynreg(itau(j),ipres(j),iw(j))+1.
!              end if
!             end if
!             end if
cIM calcul stats regime dynamique END
            end if !IM boxcloudy

          end if !IM sunlit
                       
cIM
c         CALL CPU_time(t2)
c         print*,'IF sunlit boxcloudy t2 - t1',t2 - t1
cIM
        enddo !IM ibox/j


cIM ajout stats s/ W500 BEG
cIM ajout stats s/ W500 END

c             if(j.EQ.6722) then
c             print*,' ISCCP',w(j),iw(j),ipres(j),itau(j)
c             endif

!     if (pct_ocean(j)) then
!     if (ipres(j) .gt. 0.and.itau(j) .gt. 0) then
!     if (iw(j) .gt. int(wmin).and.iw(j) .le. iwmx) then
!     if(pctj(itau(j),ipres(j),iw(j))) THEN 
!         nfq_dynreg(itau(j),ipres(j),iw(j))=
!    &    nfq_dynreg(itau(j),ipres(j),iw(j))+1.
c         if(itau(j).EQ.4.AND.ipres(j).EQ.2.AND.
c    &    iw(j).EQ.10) then
c         PRINT*,' isccp AVANT',
c    &    nfq_dynreg(itau(j),ipres(j),iw(j)),
c    &    fq_dynreg(itau(j),ipres(j),iw(j))
c         endif
!     endif
!     endif
!     endif
!     endif
39    continue !IM j/ibox
      
      !compute mean cloud properties
      do j=1,npoints 
       if (totalcldarea(j) .gt. 0.) then
         meanptop(j) = meanptop(j) / totalcldarea(j)
         if (sunlit(j).eq.1) then
           meanalbedocld(j)=meanalbedocld(j) / totalcldarea(j)
           meantaucld(j)=tautab(min(255,max(1,nint(meanalbedocld(j)))))
         end if
       end if
      enddo ! j
!
cIM ajout stats s/ W500 BEG
!     do nw = 1, iwmx
!     do l = 1, 7
!     do k = 1, 7
!       if (nfq_dynreg(k,l,nw).GT.0.) then
!       fq_dynreg(k,l,nw) = fq_dynreg(k,l,nw)/nfq_dynreg(k,l,nw)
c        if(k.EQ.4.AND.l.EQ.2.AND.nw.EQ.10) then
c        print*,' isccp APRES',nfq_dynreg(k,l,nw),
c    &   fq_dynreg(k,l,nw)
c        endif
!       else
!        if(fq_dynreg(k,l,nw).NE.0.) then
!        print*,'nfq_dynreg = 0 tau,pc,nw',k,l,nw,fq_dynreg(k,l,nw)
!        endif 
c        fq_dynreg(k,l,nw) = -1.E+20
c        nfq_dynreg(k,l,nw) = 1.E+20
!       end if 
!     enddo !k
!     enddo !l
!     enddo !nw
cIM ajout stats s/ W500 END
!     ---------------------------------------------------!

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
!cIM
c      if (debugcol.ne.0) then
!     
c         do j=1,npoints,debugcol

c            !produce character output
c            do ilev=1,nlev
c              do ibox=1,ncol
c                   acc(ilev,ibox)=0
c              enddo
c            enddo

c            do ilev=1,nlev
c              do ibox=1,ncol
c                   acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
c                   if (levmatch(j,ibox) .eq. ilev) 
c     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
c              enddo
c            enddo

             !print test

c          write(ftn09,11) j
c11        format('ftn09.',i4.4)
c         open(9, FILE=ftn09, FORM='FORMATTED')

c             write(9,'(a1)') ' '
c                    write(9,'(10i5)')
c     &                  (ilev,ilev=5,nlev,5)
c             write(9,'(a1)') ' '
             
c             do ibox=1,ncol
c               write(9,'(40(a1),1x,40(a1))')
c     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) 
c     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
c             end do
c            close(9)
c
cIM
c             if (ncolprint.ne.0) then
c               write(6,'(a1)') ' '
c                    write(6,'(a2,1X,5(a7,1X),a50)') 
c     &                  'ilev',
c     &                  'pfull','at',
c     &                  'cc*100','dem_s','dtau_s',
c     &                  'cchar'

!               do 4012 ilev=1,nlev
!                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
!                   write(6,'(i2,1X,5(f7.2,1X),50(a1))') 
!     &                  ilev,
!     &                  pfull(j,ilev)/100.,at(j,ilev),
!     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
!     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
!4012           continue
c               write (6,'(a)') 'skt(j):'
c               write (6,'(8f7.2)') skt(j)
                                      
c               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
              
c               write (6,'(a)') 'tau:'
c               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    
c               write (6,'(a)') 'tb:'
c               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
    
c               write (6,'(a)') 'ptop:'
c               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
c             endif 
    
c        enddo
       
c      end if 

      return
      end 
