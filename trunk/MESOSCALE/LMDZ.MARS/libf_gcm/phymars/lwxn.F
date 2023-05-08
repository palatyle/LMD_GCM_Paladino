      subroutine lwxn ( ig0,kdlon,kflev
     .                , dp
     .                , aer_t,co2_u,co2_up)

c----------------------------------------------------------------------
c     LWXN   computes transmission function and exchange coefficiants
c                        for neighbours layers
c                          (co2 / aerosols)
c                       (bands 1 and 2 of co2) 
c----------------------------------------------------------------------
c  
c              |---|---|---|---|---|---|---|---|
c   kflev+1    |   |   |   |   |   |   |   | 0 |  (space)
c              |---|---|---|---|---|---|---|---|
c    kflev     |   |   |   |   |   |***| 0 |   |
c              |---|---|---|---|---|---|---|---|
c     ...      |   |   |   |   |***| 0 |***|   |
c              |---|---|---|---|---|---|---|---|
c      4       |   |   |   |***| 0 |***|   |   | 
c              |---|---|---|---|---|---|---|---|
c      3       |   |   |***| 0 |***|   |   |   |
c              |---|---|---|---|---|---|---|---|
c      2       |   |***| 0 |***|   |   |   |   |
c              |---|---|---|---|---|---|---|---|
c      1       |   | 0 |***|   |   |   |   |   |
c              |---|---|---|---|---|---|---|---|
c      0       | 0 |   |   |   |   |   |   |   |  (ground)
c              |---|---|---|---|---|---|---|---|
c                0   1   2   3   4  ...  k |k+1 
c             (ground)                    (space)
c  
c  (*)  xi computed in this subroutine
c----------------------------------------------------------------------
c                                                                       
c  *********************************************************** nj=1
c
c
c                                       sublayer    1                  
c
c
c                                ----------------------------- nj=2
c
c                                       sublayer    2                  
c
c    - - - -  LAYER j - - - - -  ----------------------------- nj=3
c
c                                       sublayer    3                  
c                                ----------------------------- nj=4
c                                       sublayer  ncouche                  
c  *********************************************************** ni=nj=ncouche+1
c                                       sublayer  ncouche                  
c                                ----------------------------- ni=4
c                                       sublayer    3                  
c
c    - - - -  LAYER i - - - - -  ----------------------------- ni=3
c
c                                       sublayer    2                  
c                     
c                                ----------------------------- ni=2
c
c
c                                       sublayer    1                  
c
c
c  *********************************************************** ni=1
c                                                                       
c-----------------------------------------------------------------------
c ATTENTION AUX UNITES:
c le facteur 10*g fait passer des kg m-2 aux g cm-2
c-----------------------------------------------------------------------

      implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
 
#include "yomlw.h"
#include "callkeys.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c                                                            inputs:
c                                                            -------
      integer ig0
      integer kdlon     ! part of ngrid
      integer kflev     ! part of nalyer

      real dp (ndlo2,kflev)              ! layer pressure thickness (Pa)

      real aer_t (ndlo2,nuco2,kflev+1)   ! transmission (aer)
      real co2_u (ndlo2,nuco2,kflev+1)   ! absorber amounts (co2)
      real co2_up (ndlo2,nuco2,kflev+1)  ! idem scaled by the pressure (co2)

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer ja,jl,jk,nlmd,ni,nj

      integer nmax
      parameter (nmax=50)               ! max: 50 sublayers

      real cn (nmax), cb (nmax)

      real zu_layer_i (ndlon,nuco2)
      real zup_layer_i (ndlon,nuco2)
      real zt_aer_layer_i (ndlon,nuco2)

      real zu_layer_j (ndlon,nuco2)
      real zup_layer_j (ndlon,nuco2)
      real zt_aer_layer_j (ndlon,nuco2)

      real zu (ndlon,nuco2)
      real zup (ndlon,nuco2)
      real zt_co2 (ndlon,nuco2)
      real zt_aer (ndlon,nuco2)

      real zu_i (ndlon,nuco2,nmax+1)
      real zup_i (ndlon,nuco2,nmax+1)
      real zu_j (ndlon,nuco2,nmax+1)
      real zup_j (ndlon,nuco2,nmax+1)
      real zt_aer_i (ndlon,nuco2,nmax+1)
      real zt_aer_j (ndlon,nuco2,nmax+1)

      real trans (ndlon,nuco2,nmax+1,nmax+1)
      real ksi (ndlon,nuco2,nflev-1)

c----------------------------------------------------------------------
c         0.3   Initialisation
c               --------------

      jk=ncouche+1
      do ja = 1 ,nuco2
        do jl = 1 , kdlon
          zu_i (jl,ja,jk) = 0.
          zup_i (jl,ja,jk) = 0.
          zu_j (jl,ja,jk) = 0.
          zup_j (jl,ja,jk) = 0.
          zt_aer_i (jl,ja,jk) = 1.
          zt_aer_j (jl,ja,jk) = 1.
        enddo
      enddo

      if (linear) then

      do nlmd = 1 ,ncouche
        cn(nlmd)=(1.0/ncouche)
        cb(nlmd)=(ncouche-nlmd+0.5)/ncouche
c       print*,nlmd,cb(nlmd),cn(nlmd)
      enddo

      else

      do nlmd = 1 ,ncouche-1
        cn(nlmd)=(1-alphan)*alphan**(nlmd-1)
        cb(nlmd)=0.5*(1+alphan)*alphan**(nlmd-1)
      enddo
      cn(ncouche)=alphan**(ncouche-1)
      cb(ncouche)=0.5*alphan**(ncouche-1)

      endif

c test
      if (nmax .LT. ncouche) then
        print*,'!!!!! ATTENTION !!!!! '
        print*,'probleme dans lwxn.F'
        print*,' nmax=',nmax,'  < ncouche=',ncouche
        call exit(1)
      endif

c----------------------------------------------------------------------
      do jk = 1 , nlaylte-1
c----------------------------------------------------------------------
c         1.0   (co2) amount and (aer) transmission for all sublayers  
c               ----------------------------------------------------

        do ja = 1 , nuco2
          do jl = 1 , kdlon

c                                                        layer i (down)
c                                                        -------------
      zu_layer_i(jl,ja) =  co2_u(jl,ja,jk)  - co2_u(jl,ja,jk+1)
      zup_layer_i(jl,ja) = co2_up(jl,ja,jk) - co2_up(jl,ja,jk+1)
      zt_aer_layer_i(jl,ja) = aer_t(jl,ja,jk)
     .                       / aer_t(jl,ja,jk+1)
                                  
      do nlmd=1,ncouche
        zu_i(jl,ja,nlmd)=cn(nlmd)*zu_layer_i(jl,ja)
        zup_i(jl,ja,nlmd)=cn(nlmd)*zup_layer_i(jl,ja)
        zt_aer_i(jl,ja,nlmd)=zt_aer_layer_i(jl,ja)**cn(nlmd)
      enddo

c                                                          layer j (up)
c                                                          ------------
      zu_layer_j(jl,ja) =  co2_u(jl,ja,jk+1)  - co2_u(jl,ja,jk+2)
      zup_layer_j(jl,ja) = co2_up(jl,ja,jk+1) - co2_up(jl,ja,jk+2)
      zt_aer_layer_j(jl,ja) = aer_t(jl,ja,jk+1)
     .                       / aer_t(jl,ja,jk+2)
        
      do nlmd=1,ncouche
        zu_j(jl,ja,nlmd)=cn(nlmd)*zu_layer_j(jl,ja)
        zup_j(jl,ja,nlmd)=cn(nlmd)*zup_layer_j(jl,ja)
        zt_aer_j(jl,ja,nlmd)=zt_aer_layer_j(jl,ja)**cn(nlmd)
      enddo

          enddo
        enddo

c----------------------------------------------------------------------
c         2.0   transmissions between all sublayers
c               ------------------------------------

        do ni = 1 ,ncouche+1

            do ja = 1 ,nuco2
              do jl = 1 , kdlon
      zu(jl,ja)=0.    
      zup(jl,ja)=0.    
      zt_aer(jl,ja)=1.    
                        
      do nlmd=ni,ncouche+1
        zu(jl,ja)=zu(jl,ja)+zu_i(jl,ja,nlmd)    
        zup(jl,ja)=zup(jl,ja)+zup_i(jl,ja,nlmd)    
        zt_aer(jl,ja)=zt_aer(jl,ja)*zt_aer_i(jl,ja,nlmd)    
      enddo
              enddo
            enddo

      call lwtt(kdlon,zu,zup,nuco2,zt_co2)
                         
        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            trans(jl,ja,ni,ncouche+1)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo
                   
c on ajoute la couche J
            do ja = 1 ,nuco2
              do jl = 1 , kdlon
      zu(jl,ja)=zu(jl,ja)+zu_layer_j(jl,ja) 
      zup(jl,ja)=zup(jl,ja)+zup_layer_j(jl,ja) 
      zt_aer(jl,ja)=zt_aer(jl,ja)*zt_aer_layer_j(jl,ja) 
              enddo
            enddo
                     
      call lwtt(kdlon,zu,zup,nuco2,zt_co2)
                         
        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            trans(jl,ja,ni,1)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo
                   
        enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        do nj = 1 ,ncouche+1

            do ja = 1 ,nuco2
              do jl = 1 , kdlon
      zu(jl,ja)=0.    
      zup(jl,ja)=0.    
      zt_aer(jl,ja)=1.    
                        
      do nlmd=nj,ncouche+1
        zu(jl,ja)=zu(jl,ja)+zu_j(jl,ja,nlmd)    
        zup(jl,ja)=zup(jl,ja)+zup_j(jl,ja,nlmd)    
        zt_aer(jl,ja)=zt_aer(jl,ja)*zt_aer_j(jl,ja,nlmd)    
      enddo
              enddo
            enddo
                     
      call lwtt(kdlon,zu,zup,nuco2,zt_co2)
                         
        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            trans(jl,ja,ncouche+1,nj)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo

c on ajoute la couche I
            do ja = 1 ,nuco2
              do jl = 1 , kdlon
      zu(jl,ja)=zu(jl,ja)+zu_layer_i(jl,ja) 
      zup(jl,ja)=zup(jl,ja)+zup_layer_i(jl,ja) 
      zt_aer(jl,ja)=zt_aer(jl,ja)*zt_aer_layer_i(jl,ja) 
              enddo
            enddo
                     
      call lwtt(kdlon,zu,zup,nuco2,zt_co2)
                         
        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            trans(jl,ja,1,nj)=zt_co2(jl,ja)*zt_aer(jl,ja)
          enddo
        enddo
                   
        enddo
        
c----------------------------------------------------------------------
c         3.0   global exchange coefficiant between neigthbours
c               -----------------------------------------------
        
        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            ksi(jl,ja,jk) = 0.
          enddo
        enddo

        do ni = 1 ,ncouche
          do ja = 1 ,nuco2
            do jl = 1 , kdlon

      ksi(jl,ja,jk)=ksi(jl,ja,jk) +
     .             ( trans(jl,ja,ni+1,ncouche+1)
     .             - trans(jl,ja,ni,ncouche+1)
     .             - trans(jl,ja,ni+1,1)
     .             + trans(jl,ja,ni,1) )
     .             * (cb(ni)*dp(jl,jk)) * 2    
     .             /  (dp(jl,jk) + dp(jl,jk+1))       !!!!!!!!!!!!!!!!!!!

            enddo
          enddo
        enddo

        do nj = 1 ,ncouche
          do ja = 1 ,nuco2
            do jl = 1 , kdlon

      ksi(jl,ja,jk)=ksi(jl,ja,jk) +
     .             ( trans(jl,ja,ncouche+1,nj+1)
     .             - trans(jl,ja,1,nj+1)
     .             - trans(jl,ja,ncouche+1,nj)
     .             + trans(jl,ja,1,nj) )
     .             * (cb(nj)*dp(jl,jk+1)) * 2 
     .             /  (dp(jl,jk) + dp(jl,jk+1))       !!!!!!!!!!!!!!!!!!!

            enddo
          enddo
        enddo

        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            xi(ig0+jl,ja,jk,jk+1) = ksi(jl,ja,jk)
     .                            + xi_emis(ig0+jl,ja,jk)

c                                                        ksi reciprocity
c                                                        ---------------
            xi(ig0+jl,ja,jk+1,jk) = xi(ig0+jl,ja,jk,jk+1)
          enddo
        enddo

c----------------------------------------------------------------------
c         4.0   Special treatment for ground
c               ----------------------------


      if (jk .EQ. 1) then

        do ja = 1 ,nuco2
          do jl = 1 , kdlon
            xi_ground(ig0+jl,ja)=0.
          enddo
        enddo

        do ni = 1 ,ncouche
          do ja = 1 ,nuco2
              do jl = 1 , kdlon

      xi_ground(ig0+jl,ja) = xi_ground(ig0+jl,ja)
     .                     + ( trans(jl,ja,ni+1,ncouche+1)
     .                        -trans(jl,ja,ni,ncouche+1))
     .                     * 2 * cb(ni)
            enddo
          enddo
        enddo

      endif

c----------------------------------------------------------------------
      enddo    !  boucle sur jk
c----------------------------------------------------------------------
      return
      end



