      SUBROUTINE inicondens(ny,press,temp,nomy,yc)

c=======================================================================
c   initialisation des profils de saturation des traceurs 
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

      use dimphy
      use mod_grid_phy_lmdz, only: nbp_lev
      IMPLICIT NONE
#include "YOMCST.h"

c    Arguments :
c    -----------
      INTEGER ny
      real press(nbp_lev),temp(nbp_lev)  ! pressure in mbar !
      REAL yc(nbp_lev,ny)
      character*10 nomy(ny+1)
      
c    Local variables :
c    -----------------
      INTEGER  l,ic
      real     sy,x
       
      do ic=1,ny
       print*, 'traceur CH(', ic, ')=', nomy(ic),'------------'
           do l=1,nbp_lev

c  Par defaut, yc est a 1 c'est a dire qu'on ne condense pas
              yc(l,ic)=1.

              if(nomy(ic).eq."CH4") then
                 if (temp(l).lt.90.65) then
                   yc(l,ic)= 
     s             10.0**(4.42507e0 - ( ( ( 1165560.7e0 / TEMP(l) -
     s             115352.19e0 ) / TEMP(l) + 4055.6016e0 ) / TEMP(l)
     s            + 453.92414e0 ) / TEMP(l) ) / PRESS(l) * 1013.25e0
                 else
                   yc(l,ic)= 
     s             10.0**(3.901408e0 - ( ( 154567.02e0 / TEMP(l) -
     s             1598.8512e0 ) / TEMP(l) + 437.54809e0 ) / TEMP(l))
     s             / PRESS(l) * 1013.25e0;
                 endif
c maintient a 1.4% minimum               
                 if (yc(l,ic).lt.0.014) yc(l,ic)=0.014
              endif

              if(nomy(ic).eq."C2H2") then
                 yc(l,ic)= 
     s            10.0**(6.09748e0-1644.1e0/TEMP(l)+7.42346e0
     s           * alog10(1.0e3/TEMP(l)) ) / PRESS(l)*1013.25e0/760.0
              endif

              if(nomy(ic).eq."C2H4") then
                 if (temp(l).lt.89.0) then
                   yc(l,ic)= 
     s              10.0**(1.5477e0 + (1.0e0/TEMP(l) - 0.011e0)
     s              *(16537.0e0*(1.0e0/TEMP(l) - 0.011e0) - 1038.1e0)) 
     s              / PRESS(l) * 1.01325e0 / 760.0
                 elseif (temp(l).lt.104.0) then
                   yc(l,ic)= 
     s              10.0**(8.724e0 - 901.6e0/(TEMP(l) - 2.555e0) )
     s              / PRESS(l) * 1013.25e0 / 760.0
                 elseif (temp(l).lt.120.0) then
                   yc(l,ic)= 
     s              10.0**(50.79e0 - 1703.0e0/TEMP(l) - 17.141e0 *
     s              alog10(TEMP(l)) ) / PRESS(l) * 1013.25e0 / 760.0
                 elseif (temp(l).lt.155.0) then
                   yc(l,ic)= 
     s              10.0**(6.74756e0 - 585.0e0/(TEMP(l) - 18.16e0) ) 
     s              / PRESS(l) * 1013.25e0 / 760.0
                 endif
              endif

              if(nomy(ic).eq."C2H6") then
                 if (temp(l).lt.90.) then
                   yc(l,ic)= 
     s              10.0**(10.01e0-1085.0e0/(TEMP(l)-0.561e0) ) 
     s                  / PRESS(l) * 1013.25e0 / 760.0e0
                 else
                   yc(l,ic)= 
     s              10.0**(5.9366e0 - 1086.17e0/TEMP(l) + 3.83464e0 *
     s              alog10(1.0e3/TEMP(l)) ) / PRESS(l)*1013.25e0/760.0
                 endif
              endif

              if((nomy(ic).eq."CH3CCH")
     s       .or.(nomy(ic).eq."CH2CCH2")) then
                 yc(l,ic)= 
     s            10.0**(2.8808e0 - 4.5e0*(249.9e0 - TEMP(l))
     s                            /(1.15e0*TEMP(l) - 37.485e0) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C3H6")  then
                 yc(l,ic)= 
     s            10.0**(7.4463e0 - 1028.5654e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C3H8")  then
                 yc(l,ic)= 
     s            10.0**(7.217e0 - 994.30251e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if((nomy(ic).eq."C4H2")
     s       .or.(nomy(ic).eq."C4H2s")) then
                 yc(l,ic)= 
     s            10.0**(96.26781e0 - 4651.872e0/TEMP(l) - 31.68595e0
     s            *alog10(TEMP(l)) ) / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C4H4")  then
                 yc(l,ic)= 
     s            1.0e3 * exp(9.3898e0 - 2203.57/(TEMP(l)-43.15e0) ) 
     s            / PRESS(l)
              endif

              if(nomy(ic).eq."C4H6")  then
                 yc(l,ic)= 
     s            10.0**(2.8808e0 - 4.6e0*(262.3e0 - TEMP(l))
     s                            /(1.15e0*TEMP(l) - 39.345e0) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C4H10")  then
                 yc(l,ic)= 
     s            10.0**(8.446e0 - 1461.2e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C6H2")  then
                 yc(l,ic)= 
     s            10.0**(4.666e0 - 4956e0/TEMP(l) + 25.845e0 *
     s            alog10(1.0e3/TEMP(l)) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C8H2")  then
                 yc(l,ic)= 
     s            10.0**(3.95e0 - 6613e0/TEMP(l) + 35.055e0 *
     s            alog10(1.0e3/TEMP(l)) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."AC6H6")  then
                 x = 1.0e0 - TEMP(l) / 562.2e0
                 yc(l,ic)= 
     s            48.9e3 * exp( ( 1.33213 * x**1.5 - 6.98273 * x
     s                           - x**3 * (2.62863 + 3.33399 * x**3) )
     s                         * 562.2e0/TEMP(l) ) / PRESS(l)
              endif

              if(nomy(ic).eq."HCN")  then
                 yc(l,ic)= 
     s            10.0**(8.6165e0 - 1516.5e0/(TEMP(l) - 26.2e0) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."CH3CN")  then
                 yc(l,ic)= 
     s            10.0**(8.458e0 - 1911.7e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C2H3CN")  then
                 yc(l,ic)= 
     s            10.0**(9.3051e0 - 2782.21/(TEMP(l) - 51.15e0) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."NCCN")  then
                 yc(l,ic)= 
     s            10.0**(7.454e0 - 1832e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."HC3N")  then
                 yc(l,ic)= 
     s            10.0**(7.7446e0 - 1453.5609e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

              if(nomy(ic).eq."C4N2")  then
                 yc(l,ic)= 
     s            10.0**(8.269e0 - 2155.0e0/TEMP(l) )
     s                 / PRESS(l) * 1013.25e0 / 760.0e0
              endif

           enddo
       enddo

       print*, 'inicondens end'
       
      RETURN
      END
