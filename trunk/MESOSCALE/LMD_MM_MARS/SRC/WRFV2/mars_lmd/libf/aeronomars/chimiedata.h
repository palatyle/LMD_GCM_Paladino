c ------------------------------------
c ---  Data for chemical routines  ---
c ------------------------------------
c    Author: Sebastien Lebonnois
c            Franck Lefevre
c    
c  Update after water clouds included
c ------------------------------------
c
      integer    nd, nozo, nr, nsza, ntemp
c
      parameter (nd    = 11)
      parameter (nozo  = 7)
      parameter (nr    = nd + 28)
      parameter (nsza  = 27)
      parameter (ntemp = 4)
c
      real       kb
      parameter (kb = 1.38e-23)
c
c     factor: molecular gas constant/(boltzmann*gravity)
c
      real       factor
      parameter (factor = 192./(kb*1.e4*3.72))
c
      common/chimiedata/reactifs,nprod,nperte,prod,perte,jphot,
     $                  colairtab,table_ozo
c
c     used in obsolete scheme A :
      integer    reactifs(nr,5),nprod(nqmx),nperte(nqmx)
      integer    prod(nqmx,100),perte(nqmx,100,2)
c
c     VERSION: 28/05/2004  (update: Ar)
c
c     photodissociation rates matrix: version dependent
c     for this version, file jmars.20030707
c
      real       jphot(ntemp,nsza,0:200,nozo,nd)
      real       colairtab(0:200)
      real       szatab(nsza)
      real       table_ozo(nozo)
c
      data       szatab/0.,  5., 10., 15., 20., 25.,
     $                 30., 35., 40., 45., 50., 55.,
     $                 60., 65., 70., 75., 80., 82.,
     $                 84., 86., 88., 90., 91., 92.,
     $                 93., 94., 95./
c
c     compounds taken into account: version dependent
c     DOES NOT INCLUDE h2o (always iq=nqmx)
c            NOR water ice (if iceparty, iq=nqmx-1)
c
      integer       ncomp
      parameter    (ncomp=13)
      character*10  nomchem(ncomp)
      data          nomchem/"co2","co","o","o1d","o2","o3","h","h2",
     $                   "oh","ho2","h2o2", "n2", "ar"/
c
      real  mmolchem(ncomp)
      data  mmolchem/44.,28.,16.,16.,32.,48.,1.,2.,17.,33.,34.,28.,40./

