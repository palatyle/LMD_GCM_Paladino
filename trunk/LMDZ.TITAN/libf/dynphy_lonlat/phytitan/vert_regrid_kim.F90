SUBROUTINE vert_regrid_kim(nq,q)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Purpose : * Calculates the zonally averaged upper chemistry fields according
  !           to the new pressure grid - based on interp_vert.
  !           * In case the GCM top is lowered we interpolate upper fields
  !           between the former ones and the GCM top layer zonally averaged.
  !           * In case the GCM top is highered we interpolate tracers on the
  !           GCM grid uppermost layers between the upper chemsitry fields and
  !           the chem. tracers below.
  !
  ! Author : Jan Vatant d'Ollone (2018)
  ! ------
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE comchem_h, ONLY: nkim, nlaykim_up, preskim, ykim_up
  USE comchem_newstart_h
  USE tracer_h
  USE comvert_mod, ONLY: aps
  
  IMPLICIT NONE
  
  INCLUDE "netcdf.inc"

  INCLUDE "dimensions.h"
 
  ! -----------------
  ! Declarations
  ! -----------------
  INTEGER, INTENT(IN) :: nq ! Total number of advected fields (tracers)
 
  REAL,DIMENSION(iim+1,jjm+1,llm,nq), INTENT(INOUT) :: q   ! Advected fields (kg/kg) on 3D dyn. grid

  REAL, DIMENSION(:,:), ALLOCATABLE   :: avg_qtop  ! Zonally averaged q (mol/mol) on top layer

  REAL :: coef, ykimlon

  INTEGER :: ngridmx
 
  INTEGER :: ng0, isup

  INTEGER :: ln, lo, ilay, ichem, ilon, ilat

  LOGICAL :: lowered

  ! -----------------------------
  ! 0. Get useful size of arrays
  ! -----------------------------

  ngridmx = size(ykim_up,DIM=2)

  ! ------------------------------------------------------------
  ! 1. Compute zonal mean of last layer for every chem and lat
  ! and convert it to molar mixing fraction and to physics grid 
  ! Preliminary, only in case ceiling has been lowered
  ! ------------------------------------------------------------

  lowered = .FALSE.

  IF ( preskim(1) .GT. preskimold(1) ) THEN

    lowered = .TRUE.

    ALLOCATE(avg_qtop(nkim,ngridmx))

    DO ichem=1,nkim
      avg_qtop(ichem,:)=0.0
      DO ilon=1,iim
        avg_qtop(ichem,1)=avg_qtop(ichem,1)+q(ilon,1,llm,chimi_indx(ichem))
        avg_qtop(ichem,ngridmx)=avg_qtop(ichem,ngridmx)+q(ilon,jjm+1,llm,chimi_indx(ichem))
        DO ilat=2,jjm
          ng0 = iim*(ilat-2)+1
          avg_qtop(ichem,ng0+1:ng0+iim)=avg_qtop(ichem,ng0+1:ng0+iim)+q(ilon,ilat,llm,chimi_indx(ichem))
        ENDDO
      ENDDO
      ! mass -> molar mixing ratio to be comparable to ykim_up later
      avg_qtop(ichem,:)=avg_qtop(ichem,:)/rat_mmol(chimi_indx(ichem)) 
    ENDDO

  avg_qtop(:,:) = avg_qtop(:,:) / real(iim)

  ENDIF

  ! ------------------------------
  ! 2. Process upper chem. fields
  ! ------------------------------

  DO ln=1,nlaykim_up

    ! Standard case between 2 old upper pressure grid points

    DO lo=1,nlaykimold-1

      IF ( ( preskim(ln) .LE. preskimold(lo)   ) .AND. &
           ( preskim(ln) .GT. preskimold(lo+1) ) ) THEN

        coef = (preskim(ln)-preskimold(lo)) / (preskimold(lo+1)-preskimold(lo))
        ykim_up(:,:,ln) = (1.0-coef)*ykim_up_oldv(:,:,lo) + coef*ykim_up_oldv(:,:,lo+1)

      ENDIF

    ENDDO

    ! Special cases

    IF ( lowered )  THEN

    ! If the ceiling of GCM has been lowered we interpolate a zonal mean between new GCM top and former upper fields
    ! NB : We could have kept in memory the former tracers at this altitude and zonally averaged them
    ! but it's certainly useless as we're in newstart and the fields will re-equilibrate

      IF ( preskim(ln) .GT. preskimold(1) ) THEN
          DO ichem=1,nkim
            coef = ( preskim(ln)-aps(llm) ) / ( preskimold(1)-aps(llm) )
            ykim_up(ichem,:,ln) = (1.0-coef)*avg_qtop(ichem,:) + coef*ykim_up_oldv(ichem,:,1)
          ENDDO
      ENDIF

    ENDIF

    IF ( preskim(ln) .LE. preskimold(nlaykimold) ) THEN ! upper ceiling at 1300km can have slight variations
      ykim_up(:,:,ln) = ykim_up_oldv(:,:,lo) 
    ENDIF  
  
  ENDDO ! do ln=1,nlaykim_up

  ! ----------------------------------------------------------------------------
  ! 3. Correct the 3D advected chem. tracer fields if model ceiling is highered
  ! In this case we interpolate between tracers below and upper_chemistry values
  ! Doing this we convert via rat_mmol ykim_up from molar to mass mixing ratio
  ! ----------------------------------------------------------------------------
 
  IF ( preskim(1) .LT. preskimold(1) ) THEN
 
    ! We just want to process the concerned upper layers of the GCM
    DO ilay=1,llm
      IF ( aps(ilay) .LT. preskimold(1) ) THEN
        isup = ilay-1
        EXIT
      ENDIF
    ENDDO

    DO ilay=isup+1,llm
  
      coef = ( aps(ilay) - preskim(1) ) / ( aps(isup) - preskim(1) )

      DO ichem=1,nkim

        ! We need to convert ykim_up on phys grid to q on dyn grid
        ! so we deal with mono-gridpoints for North and South Poles   

        q(:,1,ilay,chimi_indx(ichem)) = (1.0-coef)*ykim_up(ichem,1,1)*rat_mmol(chimi_indx(ichem)) &
                                      + coef*q(:,1,isup,chimi_indx(ichem))
        q(:,jjm+1,ilay,chimi_indx(ichem)) = (1.0-coef)*ykim_up(ichem,ngridmx,1)*rat_mmol(chimi_indx(ichem)) &
                                          + coef*q(:,jjm+1,isup,chimi_indx(ichem))

          DO ilat=2,jjm

            ng0 = iim*(ilat-2)+1
           
            DO ilon=2,iim
              ! ykim_up and q are shifted one to the other on longitudinal grid
              ykimlon = 0.5*(ykim_up(ichem,ng0+ilon-1,1)+ykim_up(ichem,ng0+ilon,1)) * rat_mmol(chimi_indx(ichem))

              q(ilon,ilat,ilay,chimi_indx(ichem)) = (1.0-coef)*ykimlon + coef*q(ilon,ilat,isup,chimi_indx(ichem))
            ENDDO

            ! Periodicity on longitude at 180 and -180

            ykimlon = 0.5*(ykim_up(ichem,ng0+1,1)+ykim_up(ichem,ng0+iim,1)) * rat_mmol(chimi_indx(ichem))

            q(1,ilat,ilay,chimi_indx(ichem)) = (1.0-coef)*ykimlon + coef*q(1,ilat,isup,chimi_indx(ichem))
            q(iim+1,ilat,ilay,chimi_indx(ichem)) = q(1,ilat,ilay,chimi_indx(ichem))

          ENDDO ! do ilat=2,jjm
      ENDDO ! do ichem=1,nkim
    ENDDO ! do ilay=1,llm

  ENDIF
 
END SUBROUTINE vert_regrid_kim
