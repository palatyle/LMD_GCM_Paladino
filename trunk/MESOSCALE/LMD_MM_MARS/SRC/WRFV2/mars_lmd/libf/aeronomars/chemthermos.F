      SUBROUTINE chemthermos(ig,lswitch,zycol,ztemp,zdens,zpress,
     $                       zlocal,zenit,ptimestep,zday)

       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing chemical variations in the thermosphere
c
c   author:  MAC July 2003
c   ------
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comdiurn.h"
#include "param.h"
#include "param_v3.h"
#include "chimiedata.h"
#include "conc.h"
c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER lswitch,ig

      REAL zday,zycol(nlayermx,nqmx)
      REAL ptimestep
      real zenit
      real ztemp(nlayermx)
      real zdens(nlayermx)
      real zpress(nlayermx)			! in mbar
c
c    Local variables :
c    -----------------
      INTEGER nlayer,l,nesptherm
      parameter (nesptherm = 11)
      real tmean
      real aux1(nlayermx),aux2(nlayermx)
      real zlocal(nlayermx)
      real rm(nlayermx,nesptherm) 		!number density (cm-3)

      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2, 
     $        i_o1d, i_o, i_h, i_oh, i_ho2
      integer g_co2, g_co, g_o2, g_h2, g_h2o, g_h2o2, 
     $        g_o1d, g_o, g_h, g_oh, g_ho2, g_o3, g_n2

      logical firstcall
      save firstcall
      data firstcall /.true./

c      if (firstcall) then
c        call param_read
c        firstcall= .false.
c      endif

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
      g_h2o      =  nqmx

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the thermospheric chemistry
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2  = 1
      i_o2   = 2
      i_o    = 3
      i_co   = 4
      i_h    = 5
      i_oh   = 6
      i_ho2  = 7
      i_h2   = 8
      i_h2o  = 9
      i_h2o2 = 10
      i_o1d  = 11
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      nlayer=nlayermx
c      zlocal(1)=0.00625
c      do l=2,nlayer
c        tmean=ztemp(l)
c        if(ztemp(l).ne.ztemp(l-1))
c     &    tmean=(ztemp(l)-ztemp(l-1))/log(ztemp(l)/ztemp(l-1))
c        zlocal(l)= zlocal(l-1)-log(zpress(l)/zpress(l-1))
c     &                         *Rnew(ig,l-1)*tmean/g/1000.
c      enddo

      do l=1,nlayer
        aux1(l)=0.
        aux2(l)=0.
        rm(l,i_co2)  = zycol(l,g_co2)  *zdens(l)  
        rm(l,i_co)   = zycol(l,g_co)   *zdens(l) 
        rm(l,i_o)    = zycol(l,g_o)    *zdens(l)
        rm(l,i_o1d)  = zycol(l,g_o1d)  *zdens(l) 
        rm(l,i_o2)   = zycol(l,g_o2)   *zdens(l)
        rm(l,i_h)    = zycol(l,g_h)    *zdens(l) 
        rm(l,i_h2)   = zycol(l,g_h2)   *zdens(l) 
        rm(l,i_oh)   = zycol(l,g_oh)   *zdens(l) 
        rm(l,i_ho2)  = zycol(l,g_ho2)  *zdens(l)
        rm(l,i_h2o2) = zycol(l,g_h2o2) *zdens(l) 
        rm(l,i_h2o)  = zycol(l,g_h2o)  *zdens(l) 
      enddo
         
      call flujo(solarcondate+zday/365.)

      call jthermcalc
     $ (rm(1,i_co2),rm(1,i_o2),rm(1,i_o),rm(1,i_h2),rm(1,i_h2o),
     &  rm(1,i_h2o2),aux1,aux2,ztemp,nlayermx,zlocal,
     &  solarcondate+zday/365.,zenit)

      call paramfoto(lswitch,zdens,ztemp,ptimestep/3600.,zenit,
     &             nlayer,rm(1,i_co2),rm(1,i_o2),rm(1,i_o),
     &             rm(1,i_co),rm(1,i_h),rm(1,i_oh),rm(1,i_ho2),
     &             rm(1,i_h2),rm(1,i_h2o),rm(1,i_h2o2),rm(1,i_o1d))

      do l=lswitch,nlayer
        zycol(l,g_co2)  = max(rm(l,i_co2)  / zdens(l) , 1.e-30)
        zycol(l,g_co)   = max(rm(l,i_co)   / zdens(l) , 1.E-30)
        zycol(l,g_o2)   = max(rm(l,i_o2)   / zdens(l) , 1.e-30)
        zycol(l,g_h2)   = max(rm(l,i_h2)   / zdens(l) , 1.e-30)
        zycol(l,g_h)    = max(rm(l,i_h)    / zdens(l) , 1.e-30)
        zycol(l,g_oh)   = max(rm(l,i_oh)   / zdens(l) , 1.e-30)
        zycol(l,g_ho2)  = max(rm(l,i_ho2)  / zdens(l) , 1.e-30)
        zycol(l,g_h2o)  = max(rm(l,i_h2o)  / zdens(l) , 1.e-30)
        zycol(l,g_h2o2) = max(rm(l,i_h2o2) / zdens(l) , 1.e-30)
        zycol(l,g_o1d)  = max(rm(l,i_o1d)  / zdens(l) , 1.e-30)
        zycol(l,g_o)    = max(rm(l,i_o)    / zdens(l) , 1.e-30)
       enddo	!nlayer


      return
      end 
