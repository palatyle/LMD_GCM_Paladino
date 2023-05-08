      subroutine moldiffcoeff(dij)

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
#include "dimensions.h"
#include "dimphys.h"
#include "callkeys.h"
#include "comdiurn.h"
#include "chimiedata.h"
#include "tracer.h"
#include "conc.h"

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER nq, n, nn, i
      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2, i_n2, i_o3,
     $        i_o1d, i_o, i_h, i_oh, i_ho2, i_ar
      integer g_co2, g_co, g_o2, g_h2, g_h2o, g_h2o2,
     $        g_o1d, g_o, g_h, g_oh, g_ho2, g_o3, g_n2, g_ar
      integer gcmind(ncomptot)

      real dnh
      real dij(ncomptot,ncomptot)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      g_co2      =  nqchem_min
      g_co       =  nqchem_min + 1
      g_o        =  nqchem_min + 2
      g_o1d      =  nqchem_min + 3
      g_o2       =  nqchem_min + 4
      g_o3       =  nqchem_min + 5
      g_h        =  nqchem_min + 6
      g_h2       =  nqchem_min + 7
      g_oh       =  nqchem_min + 8
      g_ho2      =  nqchem_min + 9
      g_h2o2     =  nqchem_min + 10
      g_n2       =  nqchem_min + 11
      g_ar       =  nqchem_min + 12
      g_h2o      =  nqmx

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the molecular diffusion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Atomic oxygen must always be the LAST species of the list as 
c it is the dominant species at high altitudes.  
      i_co   = 1
      i_n2   = 2
      i_o2   = 3
      i_co2  = 4
      i_h2   = 5
      i_h    = 6
      i_oh   = 7
      i_ho2  = 8
      i_h2o  = 9
      i_h2o2 = 10
      i_o1d  = 11
      i_o3   = 12
      i_ar   = 13
      i_o    = 14
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     array to relate local indexes to gcm indexes
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

       print*,' COEFF CALC'
       open(56,file='coeffs.dat',status='unknown')
      do n=1,ncomptot
        if (dij(i_h2,n).gt.0.0) then
          do nn=n,ncomptot
            dij(nn,n)=dij(i_h2,n)
     &                  *sqrt(mmol(g_h2)/mmol(gcmind(nn)))
            if(n.eq.nn) dij(nn,n)=1.0
            dij(n,nn)=dij(nn,n)
          enddo 
        endif
        if (dij(i_h2,n).eq.0.0) then
          dnh=dij(i_h,i_o)*sqrt(mmol(g_o)/mmol(gcmind(n)))
          do nn=n,ncomptot
            dij(nn,n)=dnh*sqrt(mmol(g_h)/mmol(gcmind(nn)))
            if(n.eq.nn) dij(nn,n)=1.0
            dij(n,nn)=dij(nn,n)
          enddo 
        endif
      enddo 

      do n=1,ncomptot
        do nn=n,ncomptot
          write(56,*) n,nn,dij(n,nn)	!*1.e5/1.381e-23/(273**1.75)
        enddo
      enddo
      close(56)


      return   
      end 
