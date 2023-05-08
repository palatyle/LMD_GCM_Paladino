      MODULE calldrag_noro_mod
      
      IMPLICIT NONE
      
      CONTAINS
      
      SUBROUTINE calldrag_noro(ngrid,nlayer,ptimestep,
     &                 pplay,pplev,pt,pu,pv,pdtgw,pdugw,pdvgw) 



       use surfdat_h, only: zstd, zsig, zgam, zthe
       use dimradmars_mod, only: ndomainsz
       use drag_noro_mod, only: drag_noro
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Subroutine designed to call SUBROUTINE drag_noro 
c   Interface for sub-grid scale orographic scheme
c   The purpose of this subroutine is
c      1) Make some initial calculation at first call
c      2) Split the calculation in several sub-grid
c        ("sub-domain") to save memory and
c        be able run on a workstation at high resolution
c        The sub-grid size is defined in dimradmars_mod.
c
c   author:   
c   ------
c           Christophe Hourdin/ Francois Forget
c
c   changes:
c   -------
c   > J.-B. Madeleine 10W12
c   This version uses the variable's splitting, which can be usefull
c     when performing very high resolution simulation like LES.
c
c   input:
c   ----- 
c   ngrid                 number of gridpoint of horizontal grid
c   nlayer                Number of layer
c   ptimestep             Physical timestep (s)
c   pplay(ngrid,nlayer)    pressure (Pa) in the middle of each layer
c   pplev(ngrid,nlayer+1)  pressure (Pa) at boundaries of each layer
c   pt(ngrid,nlayer)       atmospheric temperature  (K)
c   pu(ngrid,nlayer)       zonal wind (m s-1)
c   pv(ngrid,nlayer)       meridional wind (m s-1)
c
c   output:
c   -------
c   pdtgw(ngrid,nlayer)    Temperature trend (K.s-1)
c   pdugw(ngrid,nlayer)    zonal wind trend  (m.s-2)
c   pdvgw(ngrid,nlayer)    meridional wind trend  (m.s-2)
c
c
c
c
c
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c

c-----------------------------------------------------------------------
c    Input/Output
c    ------------
      INTEGER ngrid,nlayer  

      real ptimestep

      REAL pplev(ngrid,nlayer+1),pplay(ngrid,nlayer)
      REAL pt(ngrid,nlayer), pu(ngrid,nlayer),pv(ngrid,nlayer)
      REAL pdtgw(ngrid,nlayer), pdugw(ngrid,nlayer),pdvgw(ngrid,nlayer)


c
c    Local variables :
c    -----------------

      REAL sigtest(nlayer+1)
      INTEGER igwd,igwdim,itest(ngrid)

      INTEGER :: ndomain
!      parameter (ndomain = (ngrid-1) / ndomainsz + 1)

      INTEGER l,ig
      INTEGER jd,ig0,nd

      REAL zulow(ngrid),zvlow(ngrid)
      REAL zustr(ngrid),zvstr(ngrid)

      REAL zplev(ndomainsz,nlayer+1)
      REAL zplay(ndomainsz,nlayer)
      REAL zt(ndomainsz,nlayer)
      REAL zu(ndomainsz,nlayer)
      REAL zv(ndomainsz,nlayer)
      INTEGER zidx(ndomainsz)
      REAL zzdtgw(ndomainsz,nlayer)
      REAL zzdugw(ndomainsz,nlayer)
      REAL zzdvgw(ndomainsz,nlayer)

      logical ll


c   local saved variables
c   ---------------------

      LOGICAL firstcall
      DATA firstcall/.true./
      SAVE firstcall


c----------------------------------------------------------------------

c     Initialisation
c     --------------

      IF (firstcall) THEN

         do l=1,nlayer+1
           sigtest(l)=pplev(1,l)/pplev(1,1)
         enddo
         call sugwd(nlayer,sigtest)

         if (ngrid .EQ. 1) then
           if (ndomainsz .NE. 1) then
             print*
             print*,'ATTENTION !!!'
             print*,'pour tourner en 1D, meme pour drag_noro '
             print*,'fixer ndomainsz=1 dans phymars/dimradmars_mod'
             print*
             call exit(1)
           endif
         endif

         firstcall=.false.
      END IF

      !! AS: moved out of firstcall to allow nesting+evoluting horiz domain
      ndomain = (ngrid-1) / ndomainsz + 1

c     Starting loop on sub-domain
c     ----------------------------

      DO jd=1,ndomain
        ig0=(jd-1)*ndomainsz
        if (jd.eq.ndomain) then
          nd=ngrid-ig0
        else
          nd=ndomainsz
        endif

c       Detecting points concerned by the scheme
c       ----------------------------------------

        igwd=0
        DO ig=ig0+1,ig0+nd
          itest(ig)=0
          ll=zstd(ig).gt.50.0
          IF(ll) then
            itest(ig)=1
            igwd=igwd+1
            zidx(igwd)=ig - ig0
          ENDIF
        ENDDO
        IGWDIM=MAX(1,IGWD)

c       Spliting input variable in sub-domain input variables
c       ---------------------------------------------------

        do l=1,nlayer+1
          do ig = 1,nd
           zplev(ig,l) = pplev(ig0+ig,l)
          enddo
        enddo

        do l=1,nlayer
          do ig = 1,nd
            zplay(ig,l) = pplay(ig0+ig,l)
            zt(ig,l) = pt(ig0+ig,l)
            zu(ig,l) = pu(ig0+ig,l)
            zv(ig,l) = pv(ig0+ig,l)
          enddo
        enddo

c       Calling gravity wave and subgrid scale topo parameterization
c       -------------------------------------------------------------

        call drag_noro (nd,nlayer,ptimestep,zplay,zplev,
     e        zstd(ig0+1),zsig(ig0+1),zgam(ig0+1),zthe(ig0+1),
     e        igwd,igwdim,zidx,itest(ig0+1),
     e        zt, zu, zv,
     s        zulow(ig0+1),zvlow(ig0+1),zustr(ig0+1),zvstr(ig0+1),
     s        zzdtgw,zzdugw,zzdvgw)

c       Un-spliting output variable from sub-domain input variables
c       ------------------------------------------------------------
c       (and devide by ptimestep -> true tendancies)

        do l=1,nlayer
         do ig = 1,nd
          pdtgw(ig0+ig,l) = zzdtgw(ig,l)/ptimestep
          pdugw(ig0+ig,l) = zzdugw(ig,l)/ptimestep
          pdvgw(ig0+ig,l) = zzdvgw(ig,l)/ptimestep
         enddo
        enddo

      ENDDO         !   (boucle jd=1, ndomain)

      END SUBROUTINE calldrag_noro
      
      END MODULE calldrag_noro_mod

