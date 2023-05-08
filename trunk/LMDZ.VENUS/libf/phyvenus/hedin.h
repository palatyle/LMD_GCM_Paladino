
      integer musize,zsize
      parameter (zsize=23)        ! pour z de 100 a 210 km
      parameter (musize=10)	  ! pour sza de 1 a -1
      common /hedin/ o_hedin,co2_hedin,n2_hedin,co_hedin,pres_hedin,    &
     &			mu_hedin,n_hedin
      		

      real  :: co2_hedin(musize,zsize),co_hedin(musize,zsize)
      real  :: n2_hedin(musize,zsize),o_hedin(musize,zsize)
      real  :: n_hedin(musize,zsize)
      real  :: pres_hedin(musize,zsize), mu_hedin(musize)
    
