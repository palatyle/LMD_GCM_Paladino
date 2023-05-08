      SUBROUTINE surface_nature(ngrid,nq,obliquit,qsurf,qsurfliquid
     &   ,qsurfsnow,rnat,oceanarea)

      USE surfdat_h
      USE comsoil_h
      USE geometry_mod, ONLY: cell_area
      USE tracer_h

      IMPLICIT none

!==================================================================
!     
!     Purpose
!     -------
!     Defines a few things
!     
!     Authors
!     ------- 
!     B. Charnay (2010)
!     
!     Called by
!     ---------
!     physiq.F
!     
!     Calls
!     -----
!     none
!
!     Notes
!     -----
!     rnat is terrain type: 0-ocean; 1-continent; 2-continental ice
!     
!==================================================================

        integer ngrid,nq

	REAL qsurf(ngrid,nq),ps(ngrid)
	REAL qsurfliquid(ngrid)
	REAL qsurfsnow(ngrid)
	INTEGER iq, ig
	INTEGER rnat(ngrid)
 	REAL oceanarea
	REAL obliquit

	do ig=1,ngrid
           rnat(ig)=1
           dryness(ig)=1        !(coefficient for evaporation)
           if (inertiedat(ig,1).gt.1E4) then
              rnat(ig)=0
           end if
	end do

! surface of all the oceans
        
        oceanarea=0.
	do ig=1,ngrid
           if (rnat(ig).eq.0)then
              oceanarea=oceanarea+cell_area(ig)
           end if
	enddo

        return
        end

