      SUBROUTINE concentrations(pplay,pt,pdt,pq,pdq,ptimestep)
                                             
       IMPLICIT NONE
c=======================================================================

c CALCULATION OF MEAN MOLECULAR MASS, Cp, Akk and R
c
c mmean(ngridmx,nlayermx)			amu
c cpnew(ngridmx,nlayermx)			J/kg/K
c rnew(ngridmx,nlayermx)			J/kg/K
c akknew(ngridmx,nlayermx)			coefficient of thermal concduction
c
c=======================================================================
c    0.  Declarations :
c    ------------------
c
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comdiurn.h"
#include "chimiedata.h"
#include "tracer.h"
#include "conc.h"

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER ngrid,nlayer

      REAL pplay(ngridmx,nlayermx)
      REAL pt(ngridmx,nlayermx)
      REAL pdt(ngridmx,nlayermx)
      real pq(ngridmx,nlayermx,nqmx)
      REAL pdq(ngridmx,nlayermx,nqmx)
      REAL ptimestep

c    Local variables :
c    -----------------
      INTEGER iq,l,ig,ll,n,k,nq
      integer gind(ncomptot)
      real ni(ncomptot)
      real nt, ntot
      real q2(ngridmx,nlayermx,ncomptot)
      real zt(ngridmx,nlayermx)
      real q2tot(ngridmx,nlayermx)
      real aki(ncomptot)
      real cpi(ncomptot)

      integer i_co2, i_co, i_o2, i_h2, i_h2o, i_h2o2,
     $        i_o1d, i_o, i_h, i_oh, i_ho2, i_n2, i_o3, i_ar

c      save    aki, cpi

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering for the thermal conduction and
c     specific heat coefficients
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      i_co2  = 1
      i_co   = 2
      i_o    = 3
      i_o1d  = 4
      i_o2   = 5
      i_o3   = 6
      i_h    = 7
      i_h2   = 8
      i_oh   = 9
      i_ho2  = 10
      i_h2o2 = 11
      i_n2   = 12
      i_ar   = 13
      i_h2o  = 14

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracer numbering in the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      gind(i_co2)   =  nqchem_min			! co2
      gind(i_co)    =  nqchem_min + 1			! co
      gind(i_o)     =  nqchem_min + 2			! o
      gind(i_o1d)   =  nqchem_min + 3			! o1d
      gind(i_o2)    =  nqchem_min + 4			! o2
      gind(i_o3)    =  nqchem_min + 5			! o3
      gind(i_h)     =  nqchem_min + 6			! h
      gind(i_h2)    =  nqchem_min + 7			! h2
      gind(i_oh)    =  nqchem_min + 8			! oh
      gind(i_ho2)   =  nqchem_min + 9			! ho2
      gind(i_h2o2)  =  nqchem_min + 10			! h2o2
      gind(i_n2)    =  nqchem_min + 11			! n2
      gind(i_ar)    =  nqchem_min + 12			! ar
      gind(i_h2o)   =  nqmx				! h2o

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Thermal conductivity and specific heat coefficients
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      aki(i_co2)  = 3.072e-4
      aki(i_co)   = 4.87e-4
      aki(i_o)    = 7.59e-4
      aki(i_o1d)  = 7.59e-4		!?
      aki(i_o2)   = 5.68e-4
      aki(i_o3)   = 3.00e-4		!?
      aki(i_h)    = 0.0
      aki(i_h2)   = 36.314e-4
      aki(i_oh)   = 7.00e-4		!?
      aki(i_ho2)  = 0.0
      aki(i_h2o2) = 0.0
      aki(i_n2)   = 5.6e-4
      aki(i_ar)   = 0.0                 !?
      aki(i_h2o)  = 0.0

      cpi(i_co2)   = 0.834e3
      cpi(i_co)    = 1.034e3
      cpi(i_o)     = 1.3e3
      cpi(i_o1d)   = 1.3e3		!?
      cpi(i_o2)    = 0.9194e3
      cpi(i_o3)    = 0.800e3	!?
      cpi(i_h)     = 20.780e3
      cpi(i_h2)    = 14.266e3
      cpi(i_oh)    = 1.045e3
      cpi(i_ho2)   = 1.065e3	!?
      cpi(i_h2o2)  = 1.000e3	!?
      cpi(i_n2)    = 1.034e3
      cpi(i_ar)    = 1.000e3    !?
      cpi(i_h2o)   = 1.870e3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      nlayer=nlayermx
      ngrid=ngridmx
      nq=nqmx

      DO l=1,nlayer
        DO ig=1,ngrid
          DO n=1,ncomptot
            q2(ig,l,n)=max(1.e-30,
     .         pq(ig,l,gind(n))+pdq(ig,l,gind(n))*ptimestep)
          ENDDO
          zt(ig,l)=pt(ig,l) +pdt(ig,l)*ptimestep
        ENDDO
      ENDDO
  
      do l=1,nlayermx
        do ig=1,ngridmx
          ntot=pplay(ig,l)/(1.381e-23*zt(ig,l))*1.e-6  ! in #/cm3
          cpnew(ig,l)=0.
          akknew(ig,l)=0.
          mmean(ig,l)=0.
          q2tot(ig,l)=0.
          nt=0.
          do n=1,ncomptot
            ni(n)=0.0 
            do k=1,ncomptot
              if(k.ne.n) ni(n)=ni(n)+q2(ig,l,k)/mmol(gind(k))
            enddo
            ni(n)=ntot/(1.+mmol(gind(n))/q2(ig,l,n)*ni(n))
            cpnew(ig,l)=cpnew(ig,l)+ni(n)*cpi(n)
            akknew(ig,l)=akknew(ig,l)+ni(n)*aki(n)
            mmean(ig,l)=mmean(ig,l)+q2(ig,l,n)/mmol(gind(n))
            q2tot(ig,l)=q2tot(ig,l)+q2(ig,l,n)
c        if(ig.eq.1.and.l.eq.1) write(*,*)'q2tot(ig,l)',n,q2tot(ig,l) 
            if(cpi(n) .ne. 0.0) nt=nt+ni(n)
          enddo 

c     print*,"concentrations rep 3",l,ig,nt,mmean(ig,l),zt(ig+1,l)

          cpnew(ig,l)=cpnew(ig,l)/nt
          akknew(ig,l)=akknew(ig,l)/nt
          mmean(ig,l)=1/mmean(ig,l)    ! in amu
          rnew(ig,l)=8.314/mmean(ig,l)*1.e3     ! J/kg/K		
        enddo
c        print*,l,mmean(1,l),cpnew(1,l),rnew(1,l)
c         write(228,*),l,pplay(1,l),ntot
      enddo
      return
      end 
