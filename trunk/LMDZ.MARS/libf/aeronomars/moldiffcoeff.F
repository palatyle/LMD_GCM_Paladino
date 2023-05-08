      subroutine moldiffcoeff(dij)

      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_o1d,
     &                      igcm_o2, igcm_o3, igcm_h, igcm_h2, igcm_oh,
     &                      igcm_ho2, igcm_h2o2, igcm_n2, igcm_ar,
     &                      igcm_h2o_vap, mmol

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

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
       integer,parameter :: ncompmoldiff = 14
      real dij(ncompmoldiff,ncompmoldiff)

c    Local variables:
c    ---------------
      INTEGER nq, n, nn, i,iq
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the molecular diffusion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Atomic oxygen must always be the LAST species of the list as
c it is the dominant species at high altitudes. 
      integer,parameter :: i_co   = 1
      integer,parameter :: i_n2   = 2
      integer,parameter :: i_o2   = 3
      integer,parameter :: i_co2  = 4
      integer,parameter :: i_h2   = 5
      integer,parameter :: i_h    = 6
      integer,parameter :: i_oh   = 7
      integer,parameter :: i_ho2  = 8
      integer,parameter :: i_h2o  = 9
      integer,parameter :: i_h2o2 = 10
      integer,parameter :: i_o1d  = 11
      integer,parameter :: i_o3   = 12
      integer,parameter :: i_ar   = 13
      integer,parameter :: i_o    = 14

! Tracer indexes in the GCM:
      integer,save :: g_co2=0
      integer,save :: g_co=0
      integer,save :: g_o=0
      integer,save :: g_o1d=0
      integer,save :: g_o2=0
      integer,save :: g_o3=0
      integer,save :: g_h=0
      integer,save :: g_h2=0
      integer,save :: g_oh=0
      integer,save :: g_ho2=0
      integer,save :: g_h2o2=0
      integer,save :: g_n2=0
      integer,save :: g_ar=0
      integer,save :: g_h2o=0

      integer,save :: gcmind(ncompmoldiff)

      real dnh
      logical,save :: firstcall=.true.

! Initializations at first call (and some sanity checks)
      if (firstcall) then
        ! identify the indexes of the tracers we'll need
        g_co2=igcm_co2
        if (g_co2.eq.0) then
          write(*,*) "moldiffcoeff: Error; no CO2 tracer !!!"
          stop
        endif
        g_co=igcm_co
        if (g_co.eq.0) then
          write(*,*) "moldiffcoeff: Error; no CO tracer !!!"
          stop
        endif
        g_o=igcm_o
        if (g_o.eq.0) then
          write(*,*) "moldiffcoeff: Error; no O tracer !!!"
          stop
        endif
        g_o1d=igcm_o1d
        if (g_o1d.eq.0) then
          write(*,*) "moldiffcoeff: Error; no O1D tracer !!!"
          stop
        endif
        g_o2=igcm_o2
        if (g_o2.eq.0) then
          write(*,*) "moldiffcoeff: Error; no O2 tracer !!!"
          stop
        endif
        g_o3=igcm_o3
        if (g_o3.eq.0) then
          write(*,*) "moldiffcoeff: Error; no O3 tracer !!!"
          stop
        endif
        g_h=igcm_h
        if (g_h.eq.0) then
          write(*,*) "moldiffcoeff: Error; no H tracer !!!"
          stop
        endif
        g_h2=igcm_h2
        if (g_h2.eq.0) then
          write(*,*) "moldiffcoeff: Error; no H2 tracer !!!"
          stop
        endif
        g_oh=igcm_oh
        if (g_oh.eq.0) then
          write(*,*) "moldiffcoeff: Error; no OH tracer !!!"
          stop
        endif
        g_ho2=igcm_ho2
        if (g_ho2.eq.0) then
          write(*,*) "moldiffcoeff: Error; no HO2 tracer !!!"
          stop
        endif
        g_h2o2=igcm_h2o2
        if (g_h2o2.eq.0) then
          write(*,*) "moldiffcoeff: Error; no H2O2 tracer !!!"
          stop
        endif
        g_n2=igcm_n2
        if (g_n2.eq.0) then
          write(*,*) "moldiffcoeff: Error; no N2 tracer !!!"
          stop
        endif
        g_ar=igcm_ar
        if (g_ar.eq.0) then
          write(*,*) "moldiffcoeff: Error; no AR tracer !!!"
          stop
        endif
        g_h2o=igcm_h2o_vap
        if (g_h2o.eq.0) then
          write(*,*) "moldiffcoeff: Error; no water vapor tracer !!!"
          stop
        endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    fill array to relate local indexes to gcm indexes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        gcmind(i_co)  =   g_co
        gcmind(i_n2)  =   g_n2
        gcmind(i_o2)  =   g_o2
        gcmind(i_co2) =   g_co2
        gcmind(i_h2)  =   g_h2
        gcmind(i_h)   =   g_h
        gcmind(i_oh)  =   g_oh
        gcmind(i_ho2) =   g_ho2
        gcmind(i_h2o) =   g_h2o
        gcmind(i_h2o2)=   g_h2o2
        gcmind(i_o1d) =   g_o1d
        gcmind(i_o3)  =   g_o3
        gcmind(i_o)   =   g_o
        gcmind(i_ar)   =  g_ar
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        firstcall= .false.
      endif ! of if (firstcall)


      dij(i_h2,i_co)   = 0.0000651
      dij(i_h2,i_n2)   = 0.0000674
      dij(i_h2,i_o2)   = 0.0000697
      dij(i_h2,i_co2)  = 0.0000550
      dij(i_h2,i_h2)   = 0.0
      dij(i_h2,i_h)    = 0.0
      dij(i_h2,i_oh)   = 0.0	!0003
      dij(i_h2,i_ho2)  = 0.0	!0003
      dij(i_h2,i_h2o)  = 0.0	!0003
      dij(i_h2,i_h2o2) = 0.0	!0003
      dij(i_h2,i_o1d)  = 0.0
      dij(i_h2,i_o3)   = 0.0	!0003
      dij(i_h2,i_o)    = 0.0
      dij(i_h2,i_ar)   = 0.0

c      dij(i_h,i_o)     = 0.0000144
      dij(i_h,i_o)     = 0.000114

       print*,'moldiffcoef: COEFF CALC'
       open(56,file='coeffs.dat',status='unknown')
      do n=1,ncompmoldiff
        if (dij(i_h2,n).gt.0.0) then
          do nn=n,ncompmoldiff
            dij(nn,n)=dij(i_h2,n)
     &                  *sqrt(mmol(g_h2)/mmol(gcmind(nn)))
            if(n.eq.nn) dij(nn,n)=1.0
            dij(n,nn)=dij(nn,n)
          enddo 
        endif
        if (dij(i_h2,n).eq.0.0) then
          dnh=dij(i_h,i_o)*sqrt(mmol(g_o)/mmol(gcmind(n)))
          do nn=n,ncompmoldiff
            dij(nn,n)=dnh*sqrt(mmol(g_h)/mmol(gcmind(nn)))
            if(n.eq.nn) dij(nn,n)=1.0
            dij(n,nn)=dij(nn,n)
          enddo 
        endif
      enddo 

      do n=1,ncompmoldiff
        do nn=n,ncompmoldiff
          write(56,*) n,nn,dij(n,nn)	!*1.e5/1.381e-23/(273**1.75)
        enddo
      enddo
      close(56)


      return   
      end 
