      SUBROUTINE soil_tifeedback(ngrid,nsoil,icecover,newtherm_i)

      use tracer_mod, only: nqmx, igcm_h2o_ice, rho_ice
      use comsoil_h, only: layer, inertiedat
      use surfdat_h, only: watercaptag, inert_h2o_ice
      IMPLICIT NONE

c=======================================================================
c   Description :
c       Surface water ice / Thermal inertia feedback.
c
c   When surface water-ice is thick enough, this routine creates a new
c   soil thermal inertia with three different layers :
c   - One layer of surface water ice (the thickness is given
c     by the variable icecover (in kg of ice per m2) and the thermal
c     inertia is prescribed by inert_h2o_ice (see surfdat_h)); 
c   - A transitional layer of mixed thermal inertia;
c   - A last layer of regolith below the ice cover whose thermal inertia
c     is equal to inertiedat.
c
c   To use the model :
c       SET THE tifeedback LOGICAL TO ".true." in callphys.def.
c
c   Author: J.-B. Madeleine Mars 2008 - Updated November 2012
c=======================================================================

c Local variables
c ---------------

      INTEGER :: ig                     ! Grid point (ngrid)
      INTEGER :: ik                     ! Grid point (nsoil)
      INTEGER :: iref                   ! Ice/Regolith boundary index
      INTEGER :: ngrid                  ! Number of horizontal grid points 
      INTEGER :: nsoil                  ! Number of soil layers 
      REAL :: icedepth                  ! Ice cover thickness (m)

c Inputs
c ------

      REAL icecover(ngrid,nqmx)         ! tracer on the surface (kg.m-2)
                                        ! (iq=igcm_h2o_ice) is surface
                                        ! water ice
c Outputs
c -------

      REAL newtherm_i(ngrid,nsoil)    ! New soil thermal inertia

c Initialization
c --------------

      newtherm_i(1:ngrid,1:nsoil) = 0

c Creating the new soil thermal inertia table
c -------------------------------------------
      DO ig=1,ngrid
c       Calculating the ice cover thickness
        icedepth=icecover(ig,igcm_h2o_ice)/rho_ice
c       If the ice cover is too thick or watercaptag=true,
c         the entire column is changed :
        IF ((icedepth.ge.layer(nsoil)).or.(watercaptag(ig))) THEN
          DO ik=1,nsoil
               newtherm_i(ig,ik)=inert_h2o_ice
          ENDDO 
c       We neglect the effect of a very thin ice cover :
        ELSE IF (icedepth.lt.layer(1)) THEN
          DO ik=1,nsoil
               newtherm_i(ig,ik)=inertiedat(ig,ik)
          ENDDO 
        ELSE
c         Ice/regolith boundary index :
          iref=1
c         Otherwise, we find the ice/regolith boundary:
          DO ik=1,nsoil-1
              IF ((icedepth.ge.layer(ik)).and.
     &        (icedepth.lt.layer(ik+1))) THEN
                  iref=ik+1
                  EXIT
              ENDIF
          ENDDO
c         And we change the thermal inertia:
          DO ik=1,iref-1
            newtherm_i(ig,ik)=inert_h2o_ice
          ENDDO
c         Transition (based on the equations of thermal conduction):
          newtherm_i(ig,iref)=sqrt( (layer(iref)-layer(iref-1)) /
     &        ( ((icedepth-layer(iref-1))/inert_h2o_ice**2) +
     &        ((layer(iref)-icedepth)/inertiedat(ig,ik)**2) ) )
c         Underlying regolith:
          DO ik=iref+1,nsoil
              newtherm_i(ig,ik)=inertiedat(ig,ik)
          ENDDO
        ENDIF ! icedepth
      ENDDO ! ig

c=======================================================================
      RETURN
      END
