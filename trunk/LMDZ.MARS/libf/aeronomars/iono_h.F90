MODULE iono_h

      IMPLICIT NONE

      character*1,save,allocatable :: o1d_eq(:)
      character*1,save,allocatable :: ho2_eq(:)
      character*1,save,allocatable :: oh_eq(:)
      character*1,save,allocatable :: h_eq(:)
      character*1,save,allocatable :: n2d_eq(:)
      character*1,save,allocatable :: no2_eq(:)
      character*1,save,allocatable :: o3_eq(:)
      character*1,save,allocatable :: no_eq(:)
      character*1,save,allocatable :: cplus_eq(:)
      character*1,save,allocatable :: coplus_eq(:)
      character*1,save,allocatable :: oplus_eq(:)
      character*1,save,allocatable :: n2plus_eq(:)
      character*1,save,allocatable :: hplus_eq(:)
      character*1,save,allocatable :: co2plus_eq(:)
      character*1,save,allocatable :: o2plus_eq(:)
      character*1,save,allocatable :: noplus_eq(:)
      character*1,save,allocatable :: nplus_eq(:)
      character*1,save,allocatable :: hco2plus_eq(:)
      real*8,save,allocatable ::   tauco2(:,:)
      real*8,save,allocatable ::   tauo2(:,:)
      real*8,save,allocatable ::   tauo3p(:,:)
      real*8,save,allocatable ::   tauco(:,:)
      real*8,save,allocatable ::   tauh(:,:)
      real*8,save,allocatable ::   tauoh(:,:)
      real*8,save,allocatable ::   tauho2(:,:)
      real*8,save,allocatable ::   tauh2(:,:)
      real*8,save,allocatable ::   tauh2o(:,:)
      real*8,save,allocatable ::   tauo1d(:,:)
      real*8,save,allocatable ::   tauh2o2(:,:)
      real*8,save,allocatable ::   tauo3(:,:)
      real*8,save,allocatable ::   taun(:,:)
      real*8,save,allocatable ::   tauno(:,:)
      real*8,save,allocatable ::   taun2(:,:)
      real*8,save,allocatable ::   taun2d(:,:)
      real*8,save,allocatable ::   tauno2(:,:)
      real*8,save,allocatable ::   tauco2plus(:,:)
      real*8,save,allocatable ::   tauoplus(:,:)
      real*8,save,allocatable ::   tauo2plus(:,:)
      real*8,save,allocatable ::   taucoplus(:,:)
      real*8,save,allocatable ::   taucplus(:,:)
      real*8,save,allocatable ::   taunplus(:,:) 
      real*8,save,allocatable ::   taunoplus(:,:)
      real*8,save,allocatable ::   taun2plus(:,:)
      real*8,save,allocatable ::   tauhplus(:,:)
      real*8,save,allocatable ::   tauhco2plus(:,:)

      CONTAINS

         SUBROUTINE allocate_param_iono(nlayer,nreact)

           IMPLICIT NONE

           INTEGER :: nreact
           INTEGER :: nlayer

           allocate(o1d_eq(nlayer))
           allocate(ho2_eq(nlayer))
           allocate(oh_eq(nlayer))
           allocate(h_eq(nlayer))
           allocate(n2d_eq(nlayer))
           allocate(no2_eq(nlayer))
           allocate(o3_eq(nlayer))
           allocate(no_eq(nlayer))
           allocate(cplus_eq(nlayer))
           allocate(coplus_eq(nlayer))
           allocate(oplus_eq(nlayer))
           allocate(n2plus_eq(nlayer))
           allocate(hplus_eq(nlayer))
           allocate(co2plus_eq(nlayer))
           allocate(o2plus_eq(nlayer))
           allocate(noplus_eq(nlayer))
           allocate(nplus_eq(nlayer))
           allocate(hco2plus_eq(nlayer))
           allocate(tauco2(nreact,nlayer))
           allocate(tauo2(nreact,nlayer))
           allocate(tauo3p(nreact,nlayer))
           allocate(tauco(nreact,nlayer))
           allocate(tauh(nreact,nlayer))
           allocate(tauoh(nreact,nlayer))
           allocate(tauho2(nreact,nlayer))
           allocate(tauh2(nreact,nlayer))
           allocate(tauh2o(nreact,nlayer))
           allocate(tauo1d(nreact,nlayer))
           allocate(tauh2o2(nreact,nlayer))
           allocate(tauo3(nreact,nlayer))
           allocate(taun(nreact,nlayer))
           allocate(tauno(nreact,nlayer))
           allocate(taun2(nreact,nlayer))
           allocate(taun2d(nreact,nlayer))
           allocate(tauno2(nreact,nlayer))
           allocate(tauco2plus(nreact,nlayer))
           allocate(tauoplus(nreact,nlayer))
           allocate(tauo2plus(nreact,nlayer))
           allocate(taucoplus(nreact,nlayer))
           allocate(taucplus(nreact,nlayer))
           allocate(taunplus(nreact,nlayer)) 
           allocate(taunoplus(nreact,nlayer))
           allocate(taun2plus(nreact,nlayer))
           allocate(tauhplus(nreact,nlayer))
           allocate(tauhco2plus(nreact,nlayer))
         
         END SUBROUTINE allocate_param_iono


!***********************************************************************
      function temp_elect(zkm,tt,origin)

!     Computes the electronic temperature, either from Viking (origin=1) 
!     or MAVEN (origin=2) 

!***********************************************************************
      
!     Arguments 	

      real 	      tt        ! Temperature 
      real            zkm       !  Altitude in km
      integer         origin    ! Viking (origin=1) or MAVEN (origin=2)

! local variables: 
      real          temp_elect     ! electronic temperatures
      real          zhanson(9),tehanson(9)
      real          incremento
      integer       ii, i1, i2

      zhanson(1:9) = (/ 120.,130.,150.,175.,200.,225.,250.,275.,300. /)
      tehanson(2:9) = (/ 200.,300.,500.,1250.,2000.,2200.,2400.,2500. /)
      tehanson(1) = tt

      if(origin.eq.1) then
         if ( zkm .le. 120. ) then
            temp_elect = tt 
         else if(zkm .ge.300.) then
            temp_elect=tehanson(9)
         else
            do ii=9,2,-1 
               if ( zkm .lt. zhanson(ii) ) then 
                  i1 = ii - 1
                  i2 = ii 
               endif
            enddo
            incremento=(tehanson(i2)-tehanson(i1))/(zhanson(i2)-zhanson(i1))
            temp_elect = tehanson(i1) + (zkm-zhanson(i1)) * incremento
	endif
      else if(origin.eq.2) then
         !MAVEN measured electronic temperature (Ergun et al., GRL 2015)
         !Note that the Langmuir probe is not sensitive below ~500K, so 
         !electronic temperatures in the lower thermosphere (<150 km) may
         !be overestimated by this formula
         if(zkm.le.120) then
            temp_elect = tt
         else 
            temp_elect=((3140.+120.)/2.)+((3140.-120.)/2.)*tanh((zkm-241.)/60.)
         endif
      else
         write(*,*)'Error in function temp_elect:'
         write(*,*)'Origin must be either 1 or 2'
         write(*,*)'Using neutral temperature instead'
         temp_elect = tt
      endif

      return

      end function temp_elect

END MODULE iono_h
