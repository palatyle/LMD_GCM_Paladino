      subroutine inichim_newstart(pq,ps,sig,nq,latfi,lonfi,airefi,
     $ thermo)

      implicit none

c=======================================================================
c
c   subject:
c   --------
c
c  Initialization of the composition for use in the building of a new start file
c  using newstart.F
c
c   Author:   Sebastien Lebonnois (08/11/2002)
c   -------
c Revised 07/2003 by Monica Angelats-i-Coll to use input files
c
c   Arguments:
c   ----------
c
c    pq(iip1,jjp1,llm,nqmx)  Advected fields, ie chemical species here
c    sig = sigma                vertical coordinate (interface of layers)
c    nq: permet de ne pas modifier l'eau (et la glace d'eau si iceparty)
c
c=======================================================================

c    Declarations :
c    --------------

#include "dimensions.h"
#include "dimphys.h"
#include "paramet.h"
#include "chimiedata.h"
#include "tracer.h"
#include "comcstfi.h"
#include "comdiurn.h"
#include "callkeys.h"
#include "temps.h"
#include "datafile.h"

c Arguments :
c -----------

      real    ps(iip1,jjp1)
      real    pq(iip1,jjp1,llm,nqmx)     ! Advected fields, ie chemical species
      real    sig(llm+1)                 ! vertical coordinate (interface of layers)
      integer nq                         !    =nqmx 
                                         ! or =nqmx-1 if H2O kept
                                         ! or =nqmx-2 if H2O and ice kept
      REAL latfi(ngridmx),lonfi(ngridmx),airefi(ngridmx)
      integer thermo ! flag for thermosphere initialisation only

c Local variables :
c -----------------

      integer  iq,i,j,l,n, nspecini,inil,iqmax
      INTEGER aa(nqmx)
      REAL qtot(iip1,jjp1,llm)
      real densconc,zgcm,zfile(252),vmrint,nt, mmean
      real nxf(252),zfilemin(252),sigfile(252)
      real vmrf(252,14),nf(252,14)
      real tfile(252),zzfile(252)
      real ni(14),niold(14)
      real mmolf(14)    ! molecular mass in amu (species in files)
      data mmolf/44.,40.,28.,32.,28.,16.,2.,   ! majors
     &           1.,17.,33.,18.,34.,16.,48./   ! minors
      character*3 tmp ! temporary variable
      integer ierr,lnblnk
      external lnblnk


c-----------------------------------------------------------------------
c    Input/Output
c    ------------

       call inichim_readcallphys()

c Dimension test:
c ---------------

       if (iceparty) then
         if ((nqchem_min+ncomp+1).ne.nqmx) then
            print*,'********* Dimension problem! ********'
            print*,"nqchem_min+ncomp+1).ne.nqmx"
            print*,"ncomp: ",ncomp
            print*,"nqchem_min: ",nqchem_min
            print*,"nqmx: ",nqmx
            print*,'Change ncomp in chimiedata.h'
            STOP
         endif
       else
         if ((nqchem_min+ncomp).ne.nqmx) then
            print*,'********* Dimension problem! ********'
            print*,"nqchem_min+ncomp).ne.nqmx"
            print*,"ncomp: ",ncomp
            print*,"nqchem_min: ",nqchem_min
            print*,"nqmx: ",nqmx
            print*,'Change ncomp in chimiedata.h'
            STOP
         endif
       endif

c noms and mmol vectors:
c ----------------------

       if (iceparty) then
          do iq=nqchem_min,nqmx-2
            noms(iq) = nomchem(iq-nqchem_min+1)
            mmol(iq) = mmolchem(iq-nqchem_min+1)
          enddo
          noms(nqmx-1) = 'ice'
          mmol(nqmx-1) = 18.
          noms(nqmx)   = 'h2o'
          mmol(nqmx)   = 18.
       else
          do iq=nqchem_min,nqmx-1
            noms(iq) = nomchem(iq-nqchem_min+1)
            mmol(iq) = mmolchem(iq-nqchem_min+1)
          enddo
          noms(nqmx)   = 'h2o'
          mmol(nqmx)   = 18.
       endif
       if (nqchem_min.gt.1) then
          do iq=1,nqchem_min-1
            noms(iq) = 'dust'
            mmol(iq) = 100.
          enddo
       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tracers numbering in the gcm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     co2      =  nqchem_min
c     co       =  nqchem_min + 1         
c     o        =  nqchem_min + 2
c     o(1d)    =  nqchem_min + 3
c     o2       =  nqchem_min + 4
c     o3       =  nqchem_min + 5
c     h        =  nqchem_min + 6
c     h2       =  nqchem_min + 7
c     oh       =  nqchem_min + 8
c     ho2      =  nqchem_min + 9
c     h2o2     =  nqchem_min + 10
c     n2       =  nqchem_min + 11
c     ar       =  nqchem_min + 12
c     h2o      =  nqmx
c
c----------------------------------------------------------------------
c Major species in files:
c
c    1=CO2 
c    2=AR
c    3=N2  
c    4=O2  
c    5=CO  
c    6=O   
c    7=H2
c
c Minor species in files:
c
c    1=H  
c    2=OH 
c    3=HO2 
c    4=H2O 
c    5=H2O2 
c    6=O1D
c    7=O3
c
c----------------------------------------------------------------------
c Major species:
        aa(nqchem_min)      = 1
        aa(nqchem_min + 12) = 2
        aa(nqchem_min + 11) = 3
        aa(nqchem_min + 4)  = 4
        aa(nqchem_min + 1)  = 5
        aa(nqchem_min + 2)  = 6
        aa(nqchem_min + 7)  = 7
c Minor species:
        aa(nqchem_min + 6)  = 8
        aa(nqchem_min + 8)  = 9
        aa(nqchem_min + 9)  = 10
        aa(nqmx)            = 11
        aa(nqchem_min + 10) = 12
        aa(nqchem_min + 3)  = 13 
        aa(nqchem_min + 5)  = 14
c----------------------------------------------------------------------
      open(210, iostat=ierr,file=
     & datafile(1:lnblnk(datafile))//'/atmosfera_LMD_may.dat')
      if (ierr.ne.0) then
        write(*,*)'Error : cannot open file atmosfera_LMD_may.dat '
        write(*,*)'(in aeronomars/inichim.F)'
        write(*,*)'It should be in :', datafile(1:lnblnk(datafile)),'/'
        write(*,*)'1) You can change this directory address in '
        write(*,*)'   file phymars/datafile.h'
        write(*,*)'2) If necessary atmosfera_LMD_may.dat (and others)'
        write(*,*)'   can be obtained online on:'
        write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
         STOP
      endif
      open(220, iostat=ierr,file=
     & datafile(1:lnblnk(datafile))//'/atmosfera_LMD_min.dat')
      if (ierr.ne.0) then
        write(*,*)'Error : cannot open file atmosfera_LMD_min.dat '
        write(*,*)'(in aeronomars/inichim.F)'
        write(*,*)'It should be in :', datafile(1:lnblnk(datafile)),'/'
        write(*,*)'1) You can change this directory address in '
        write(*,*)'   file phymars/datafile.h'
        write(*,*)'2) If necessary atmosfera_LMD_min.dat (and others)'
        write(*,*)'   can be obtained online on:'
        write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
         STOP
      endif

      read(210,*) tmp
      read(220,*) tmp
        do l=1,252
          read(210,112) tfile(l),zfile(l),nxf(l),(vmrf(l,n),n=1,7)
          zfile(l)=zfile(l)*100.              !pressure (Pa)
          sigfile(l)=zfile(l)/zfile(1)       
        enddo
        do l=1,252
          read(220,113)zfilemin(l),(vmrf(l,n),n=8,14)
          zfilemin(l)=zfilemin(l)*1000.       !height (m)
          do n=1,14
            nf(l,n)=vmrf(l,n)*nxf(l)
          enddo
        enddo
        close(210)
        close(220)
c flag thermo set to 1 or 0 in newstart 
c inil=33 for initialisation of thermosphere only        
c inil=1 for initialisation of all layers        
        if (thermo.eq.1) then
        inil=33
        else
        inil=1
        endif

      do i=1,iip1
       do j=1,jjp1
         do l=inil,llm

c          zgcm=sig(l)
          zgcm=sig(l)*ps(i,j)
          densconc=0.
          nt=0.

          do n=1,14
            call intrplf(zgcm,vmrint,zfile,nf(1,n),252)
c            call intrplf(zgcm,vmrint,sigfile,nf(1,n),252)
            ni(n)=vmrint

            densconc=densconc+ni(n)*mmolf(n)
            nt=nt+ni(n)
          enddo

          mmean=densconc/nt                     ! in amu

          if (nq.ne.nqmx) then
            do n=nqchem_min,nq-1
            pq(i,j,l,n)=ni(aa(n))/nt*mmol(n)/mmean
            if(i.eq.iip1) pq(i,j,l,n)=pq(1,j,l,n) 
            enddo
            if (iceparty .and. nq.gt.nqmx-2) then
              pq(i,j,l,nq) = 0.
              if(i.eq.iip1) pq(i,j,l,nq)=pq(1,j,l,nq) 
            else
              pq(i,j,l,nq)=ni(aa(nq))/nt*mmol(nq)/mmean
              if(i.eq.iip1) pq(i,j,l,nq)=pq(1,j,l,nq) 
            endif
          endif
          if (nq.eq.nqmx) then
            do n=nqchem_min,nq-2
              pq(i,j,l,n)=ni(aa(n))/nt*mmol(n)/mmean
              if(i.eq.iip1) pq(i,j,l,n)=pq(1,j,l,n) 
            enddo
            if (iceparty) then
              pq(i,j,l,nq-1) = 0.
              if(i.eq.iip1) pq(i,j,l,nq-1)=pq(1,j,l,nq-1) 
              pq(i,j,l,nq)=ni(aa(nq))/nt*mmol(nq)/mmean
              if(i.eq.iip1) pq(i,j,l,nq)=pq(1,j,l,nq) 
            else
              do n=nq-1,nq
                pq(i,j,l,n)=ni(aa(n))/nt*mmol(n)/mmean
                if(i.eq.iip1) pq(i,j,l,n)=pq(1,j,l,n) 
              enddo
            endif
          endif
          enddo     !nlayer
         enddo      !ngrid_j
        enddo       !ngrid_i

cccccccccccccccccccccccccccccc      
c  make sure that sum of q = 1      
c dominent species is = 1 - sum(all other species)      
cccccccccccccccccccccccccccccc      
      iqmax=nqchem_min
      
      if ((nqmx-nqchem_min).gt.10) then
       do l=1,llm
        do j=1,jjp1
          do i=1,iip1
           do iq=nqchem_min,nqmx
            if ( (pq(i,j,l,iq).gt.pq(i,j,l,iqmax)).and.
     .           (noms(iq).ne."ice") )  then 
              iqmax=iq
            endif
           enddo
           pq(i,j,l,iqmax)=1.
           qtot(i,j,l)=0
           do iq=nqchem_min,nqmx
            if ( (iq.ne.iqmax).and.
     .           (noms(iq).ne."ice") ) then        
              pq(i,j,l,iqmax)=pq(i,j,l,iqmax)-pq(i,j,l,iq)        
            endif
           enddo !iq
           do iq=nqchem_min,nqmx
            if (noms(iq).ne."ice") then
             qtot(i,j,l)=qtot(i,j,l)+pq(i,j,l,iq)
            endif
c            if (i.eq.1.and.j.eq.1.and.l.eq.1) write(*,*) 'qtot(i,j,l)',
c     $    qtot(i,j,l)
           enddo !iq
            if (i.eq.1.and.j.eq.1) write(*,*) 'qtot(i,j,l)',
     $    qtot(i,j,l)
          enddo !i   
         enddo !j   
       enddo !l  
      endif
ccccccccccccccccccccccccccccccc

66      format(i2,6(1x,e11.5))
112     format(7x,f7.3,3x,e12.6,3x,e12.6,7(3x,e12.6))
113     format(1x,f6.2,7(3x,e12.6))

      end
