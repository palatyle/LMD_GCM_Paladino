      SUBROUTINE euvheat(pt,pdt,pplev,pplay,zzlay,dist_sol,
     $     mu0,ptimestep,ptime,zday,pq,pdq,pdteuv)

       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Computing heating rate due to EUV absorption
c
c   author:  MAC 2002
c   ------
c
c   input:
c   ----- 
c   dist_sol		  sun-Mars distance (AU)
c   mu0(ngridmx)           
c   pplay(ngrid,nlayer)   pressure at middle of layers (Pa)
c
c   output:
c   -------
c
c   pdteuv(ngrid,nlayer)      Heating rate (K/s)
c
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
#include "tracer.h"
#include "conc.h"
c-----------------------------------------------------------------------
c    Input/Output
c    ------------

      REAL pplay(ngridmx,nlayermx)
      REAL zzlay(ngridmx,nlayermx)
      real pplev(ngridmx,nlayermx+1)
      REAL pt(ngridmx,nlayermx)
      REAL pdt(ngridmx,nlayermx)
      real zday
      REAL dist_sol
      real mu0(ngridmx)
      real pq(ngridmx,nlayermx,nqmx)
      real pdq(ngridmx,nlayermx,nqmx)
      real ptimestep,ptime


      REAL pdteuv(ngridmx,nlayermx)
c
c    Local variables :
c    -----------------
      integer nespeuv
      parameter (nespeuv=6)

      INTEGER l,iq,ig,n,inif, ngrid, nlayer
      real rm(nlayermx,nespeuv)				! number density (cm-3)
      real zq(ngridmx,nlayermx,nqmx)
      real zt(ngridmx,nlayermx)
      real zlocal(nlayermx)
      real zenit
      real aux1(nlayermx)
      real aux2(nlayermx)
      real jtot(nlayermx)
      real dens								! amu/cm-3
      real tx(nlayermx)
      real tmean
      
      integer i_co2, i_o2, i_h2, i_h2o, i_h2o2, i_o
      integer g_co2, g_co, g_o2, g_h2, g_h2o, g_h2o2,
     $        g_o1d, g_o, g_h, g_oh, g_ho2, g_o3, g_n2


      logical firstcall
      save firstcall
      data firstcall /.true./

      if (firstcall) then
        if ((nqchem_min+nespeuv).gt.nqmx) then
          print*,'******* Dimension problem in EUVHEAT ********'
          STOP
        endif

        firstcall= .false.
        print*,'EUV',nlayer,nlayermx,ngrid,ngridmx   
      endif

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
c     tracer numbering in the EUV heating
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2  = 1
      i_o2   = 2
      i_o    = 3
      i_h2   = 4
      i_h2o  = 5
      i_h2o2 = 6
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do l=1,nlayermx
        do ig=1,ngridmx
          do iq=nqchem_min,nqmx
            zq(ig,l,iq)=pq(ig,l,iq)+pdq(ig,l,iq)*ptimestep
          enddo
          zt(ig,l)=pt(ig,l)+pdt(ig,l)*ptimestep
        enddo
      enddo

      call flujo(solarcondate)
c     call flujo(solarcondate+zday/365.) ! version with EUV time evolution

      do ig=1,ngridmx

        zenit=acos(mu0(ig))*180./acos(-1.)
 
        do l=1,nlayermx
          dens=pplay(ig,l)/(rnew(ig,l)*zt(ig,l)) / 1.66e-21
          rm(l,i_co2)  = zq(ig,l,g_co2)  * dens / mmol(g_co2)
          rm(l,i_o2)   = zq(ig,l,g_o2)   * dens / mmol(g_o2)
          rm(l,i_o)    = zq(ig,l,g_o)    * dens / mmol(g_o)
          rm(l,i_h2)   = zq(ig,l,g_h2)   * dens / mmol(g_h2)
          rm(l,i_h2o)  = zq(ig,l,g_h2o)  * dens / mmol(g_h2o)
          rm(l,i_h2o2) = zq(ig,l,g_h2o2) * dens / mmol(g_h2o2)
        enddo

c        zlocal(1)=-log(pplay(ig,1)/pplev(ig,1))
c     &            *Rnew(ig,1)*zt(ig,1)/g
        zlocal(1)=zzlay(ig,1)
        zlocal(1)=zlocal(1)/1000.
        tx(1)=zt(ig,1)

        do l=2,nlayermx
          tx(l)=zt(ig,l)
          tmean=tx(l)
          if(tx(l).ne.tx(l-1))
     &    tmean=(tx(l)-tx(l-1))/log(tx(l)/tx(l-1))
c          zlocal(l)= zlocal(l-1)
c     &   -log(pplay(ig,l)/pplay(ig,l-1))*Rnew(ig,l-1)*tmean/g/1000.
        zlocal(l)=zzlay(ig,l)/1000.
        enddo

        call hrtherm (rm(1,i_co2),rm(1,i_o2),rm(1,i_o),
     &                rm(1,i_h2),rm(1,i_h2o),rm(1,i_h2o2),
     &                aux1,aux2,tx,nlayermx,zlocal,
     &                solarcondate+zday/365.,zenit,jtot)

        do l=1,nlayermx
          pdteuv(ig,l)=0.16*jtot(l)/10.
     &            /(cpnew(ig,l)*pplay(ig,l)/(rnew(ig,l)*zt(ig,l)))
     &                  *(1.52/dist_sol)**2 

        enddo	

      enddo  !ngrid

      return
      end 
