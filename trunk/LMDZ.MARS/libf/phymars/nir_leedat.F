c***********************************************************************
      subroutine NIR_leedat                              
                                                
c 	reads parameters for NIR NLTE calculation    
                                                
c 	nov 2011    fgg+malv    first version                
c***********************************************************************

      use datafile_mod, only: datadir

      implicit none                                  
                                                
      include 'nirdata.h'
                                                
                                                
c local variables                               

      integer 	ind                      

                              
c***********************************************************************

      open(43,file=trim(datadir)//'/NIRcorrection_feb2011.dat',
     $       status='old')         
      do ind=1,9
         read(43,*)
      enddo
      
      do ind=1,npres
         read(43,*)pres1d(ind),corgcm(ind),oco21d(ind),p1999(ind),
     $        alfa(ind)
                                !Tabulated pression to Pa
         pres1d(ind)=pres1d(ind)*100.
      enddo
      close(43)

      return

      end
