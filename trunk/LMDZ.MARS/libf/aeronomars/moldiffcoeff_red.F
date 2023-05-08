      subroutine moldiffcoeff_red(dij,indic,gcmind,ncompdiff2)

       use tracer_mod, only: nqmx, noms, mmol
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing molecular diffusion coefficients
c   following Nair 94 (pg 131)
c   author:  MAC 2002
c   ------
c
c=======================================================================
#include "callkeys.h"
#include "diffusion.h"

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
c       integer,parameter :: ncompmoldiff = 12
	integer ncompdiff2
	integer gcmind(ncompdiff2)
      real dij(ncompdiff2,ncompdiff2)
      integer indic(nqmx)

c    Local variables:
c    ---------------
      INTEGER nq, n, nn, i,iq
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the molecular diffusion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real :: dijh2co,dijh2n2,dijh2co2,dijh2o2,dijho,dijref
!	integer :: i_h2,i_h,i_o
	integer :: g_h2,g_h,g_o
!      integer,parameter :: i_o   = 1
!      integer,parameter :: i_n2   = 2
!      integer,parameter :: i_co   = 3
!      integer,parameter :: i_ar  = 4
!      integer,parameter :: i_h2   = 5
!      integer,parameter :: i_h    = 6
!      integer,parameter :: i_o2   = 7
!      integer,parameter :: i_oh  = 8
!      integer,parameter :: i_ho2  = 9
!      integer,parameter :: i_h2o = 10
!      integer,parameter :: i_h2o2  = 11
!      integer,parameter :: i_o1d   = 12
!      integer,parameter :: i_o3   = 13
!      integer,parameter :: i_n    = 13
!      integer,parameter :: i_no   = 14
!      integer,parameter :: i_no2  = 15
!      integer,parameter :: i_n2d  = 17
!      integer,parameter :: i_oplus = 18
!      integer,parameter :: i_co2    = 16
!      integer,parameter :: i_oplus = 17
!      integer,parameter :: i_hplus = 18

! Tracer indexes in the GCM:
!      integer,save :: g_co2=0
!      integer,save :: g_co=0
!      integer,save :: g_o=0
!      integer,save :: g_o1d=0
!      integer,save :: g_o2=0
!      integer,save :: g_o3=0
!      integer,save :: g_h=0
!      integer,save :: g_h2=0
!      integer,save :: g_oh=0
!      integer,save :: g_ho2=0
!      integer,save :: g_h2o2=0
!      integer,save :: g_n2=0
!      integer,save :: g_ar=0
!      integer,save :: g_h2o=0
!      integer,save :: g_n=0
!      integer,save :: g_no=0
!      integer,save :: g_no2=0
!      integer,save :: g_n2d=0
!      integer,save :: g_oplus=0
!      integer,save :: g_hplus=0

!      integer,save :: gcmind(ncompdiff)

      real dnh
      logical,save :: firstcall=.true.
      logical,parameter :: outputcoeffs=.false. ! to output 'coeffs.dat' file,
                                                ! set outputcoeffs=.true.

! Initializations at first call (and some sanity checks)
      if (firstcall) then
        ! identify the indexes of the tracers we'll need
!        g_co2=igcm_co2
!        if (g_co2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no CO2 tracer !!!"
!          stop
!        endif
!        g_n2=igcm_n2
!        if (g_n2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no N2 tracer !!!"
!          stop
!        endif
!        g_ar=igcm_ar
!        if (g_ar.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no Ar tracer !!!"
!          stop
!        endif        
!        g_h2=igcm_h2
!        if (g_h2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no H2 tracer !!!"
!          stop
!        endif
!        g_h=igcm_h
!        if (g_h.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no H tracer !!!"
!          stop
!        endif
!        g_co=igcm_co
!        if (g_co.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no CO tracer !!!"
!          stop
!        endif
!        g_o2=igcm_o2
!        if (g_o2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no O2 tracer !!!"
!          stop
!        endif
!        g_oh=igcm_oh
!        if (g_oh.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no OH tracer !!!"
!          stop
!        endif
!        g_ho2=igcm_ho2
!        if (g_ho2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no HO2 tracer !!!"
!          stop
!        endif
!        g_h2o=igcm_h2o_vap
!        if (g_h2o.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no H2O tracer !!!"
!          stop
!        endif
!        g_h2o2=igcm_h2o2
!        if (g_h2o2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no H2O2 tracer !!!"
!          stop
!        endif
!        g_o1d=igcm_o1d
!        if (g_h.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no O1d tracer !!!"
!          stop
!        endif
!        g_o3=igcm_o3
!        if (g_o3.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no O3 tracer !!!"
!          stop
!        endif
!        g_n=igcm_n
!        if (g_n.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no N tracer !!!"
!          stop
!        endif
!        g_no=igcm_no
!        if (g_no.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no NO tracer !!!"
!          stop
!        endif
!        g_no2=igcm_no2
!        if (g_no2.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no NO2 tracer !!!"
!          stop
!        endif
!        g_n2d=igcm_n2d
!        if (g_n2d.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no N2(D) tracer !!!"
!          stop
!        endif
!        g_oplus=igcm_oplus
!        if (g_oplus .eq. 0) then
!        write(*,*) "moldiffcoeff: Error; no Oplus tracer !!!"
!        stop 
!        endif
!	g_hplus=igcm_hplus
!        if (g_hplus .eq. 0) then
!        write(*,*) "moldiffcoeff: Error; no Hplus tracer !!!"
!        stop
!        endif
!        g_o=igcm_o
!        if (g_o.eq.0) then
!          write(*,*) "moldiffcoeff: Error; no O tracer !!!"
!          stop
!        endif
       
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    fill array to relate local indexes to gcm indexes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!        gcmind(i_co2)  =   g_co2
!        gcmind(i_n2)  =   g_n2
!        gcmind(i_ar)  =   g_ar
!        gcmind(i_h2) =   g_h2
!        gcmind(i_h)  =   g_h
!        gcmind(i_co)   =   g_co
!        gcmind(i_o2) =   g_o2
!        gcmind(i_oh)=   g_oh
!        gcmind(i_ho2)  =   g_ho2
!        gcmind(i_h2o) = g_h2o
!        gcmind(i_h2o2)= g_h2o2
!        gcmind(i_o1d) = g_o1d
!        gcmind(i_o3) = g_o3
!        gcmind(i_n)= g_n
!        gcmind(i_no) = g_no
!        gcmind(i_no2) = g_no2
!        gcmind(i_n2d) = g_n2d
!        gcmind(i_oplus) =  g_oplus
!        gcmind(i_hplus) = g_hplus
!        gcmind(i_o)   =   g_o

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        firstcall= .false.
      endif ! of if (firstcall)

	dijh2co = 0.0000651
	dijh2n2 = 0.0000674
	dijh2o2 = 0.0000697
	dijh2co2 = 0.0000550
	dijho = 0.000114

!      dij(i_h2,i_co)   = 0.0000651
!      dij(i_h2,i_n2)   = 0.0000674
!      dij(i_h2,i_o2)   = 0.0000697
!      dij(i_h2,i_co2)  = 0.0000550
!      dij(i_h2,i_h2)   = 0.0
!      dij(i_h2,i_h)    = 0.0
!      dij(i_h2,i_h2o)  = 0.0	!0003
!      dij(i_h2,i_h2o2) = 0.0	!0003
!      dij(i_h2,i_o3)   = 0.0	!0003
!      dij(i_h2,i_o)    = 0.0
!      dij(i_h2,i_ar)   = 0.0
!      dij(i_h2,i_n)    = 0.0

!c      dij(i_h,i_o)     = 0.0000144
!      dij(i_h,i_o)     = 0.000114

! find h2, h and o index in gcm
! these species are used to define the diffusion coefficients

	do n=1,nqmx
	if (noms(n) .eq. 'h2') g_h2=n
	if (noms(n) .eq. 'h') g_h=n
	if (noms(n) .eq. 'o') g_o=n
	enddo
	print*,'moldiffcoeff_red: gh2',g_h2,g_h,g_o

       print*,'moldiffcoeff_red: COEFF CALC'

      do n=1,ncompdiff2
	dijref=0.
	if (noms(gcmind(n)) .eq. 'co') dijref=dijh2co
	if (noms(gcmind(n)) .eq. 'n2') dijref=dijh2n2
	if (noms(gcmind(n)) .eq. 'o2') dijref=dijh2o2
	if (noms(gcmind(n)) .eq. 'co2') dijref=dijh2co2
!	print*,'test',n,dijref
        if (dijref .gt. 0.0) then
          do nn=n,ncompdiff2
            dij(nn,n)=dijref
     &                  *sqrt(mmol(g_h2)/mmol(gcmind(nn)))
            if(n.eq.nn) dij(nn,n)=1.0
            dij(n,nn)=dij(nn,n)
          enddo 
        endif
        if (dijref .eq. 0.0) then
	dijref=dijho
          dnh=dijref*sqrt(mmol(g_o)/mmol(gcmind(n)))
          do nn=n,ncompdiff2
            dij(nn,n)=dnh*sqrt(mmol(g_h)/mmol(gcmind(nn)))
            if(n.eq.nn) dij(nn,n)=1.0
            dij(n,nn)=dij(nn,n)
          enddo 
        endif
      enddo 

      if (outputcoeffs) then
       ! output coefficients in 'coeffs.dat' file
       open(56,file='coeffs.dat',status='unknown')
       do n=1,ncompdiff2
        do nn=n,ncompdiff2
          write(56,*) n,nn,dij(n,nn)	!*1.e5/1.381e-23/(273**1.75)
        enddo
       enddo
       close(56)
      endif


      return   
      end 
      
