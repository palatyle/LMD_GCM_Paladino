










      subroutine forceWCfn(ngrid,nlayer,nq,pplev,pt,dq,dqs)

      USE tracer_h
      use comcstfi_mod, only: g

      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Force tracer conservation in a column for a given pair of 
!     delta q, delta q_s
!
!     Authors
!     ------- 
!     R. Wordsworth
!     
!==================================================================

      INTEGER ngrid,nlayer,nq

      real masse, Wtot, Wdiff

      real pplev(ngrid,nlayer+1)
      real pt(ngrid)

      real dqs(ngrid,nq) 
      real dq(ngrid,nlayer,nq)

      integer iq, ig, ilay

      do iq=1,nq 
        do ig=1,ngrid
           Wtot = 0.0
           do ilay=1,nlayer
              masse = (pplev(ig,ilay) - pplev(ig,ilay+1))/g
              Wtot  = Wtot + masse*dq(ig,ilay,iq)
           enddo
           Wdiff = Wtot + dqs(ig,iq)
         
           dqs(ig,iq) = dqs(ig,iq) - Wdiff
        enddo
      enddo

      end 

