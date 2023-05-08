      subroutine aaam_bud (iam,nlon,nlev,rjour,rsec,
     i                   rea,rg,ome,      
     i                   plat,plon,phis,
     i                   dragu,liftu,clu,
     i                   dragv,liftv,clv,
     i                   p, u, v)
c
      use dimphy
      use mod_grid_phy_lmdz, only: nbp_lon, nbp_lat, klon_glo
      implicit none
c======================================================================
c Auteur(s): F.Lott (LMD/CNRS) date: 20031020
c Object: Compute different terms of the axial AAAM Budget.
C No outputs, every AAM quantities are written on the IAM
C File. 
C WARNING: Only valid for regular rectangular grids.
C REMARK: CALL DANS PHYSIQ AFTER lift_noro:
C        CALL aaam_bud (27,klon,klev,rjourvrai,gmtime,
C    C               ra,rg,romega,
C    C               rlat,rlon,pphis,
C    C               zustrdr,zustrli,zustrcl,
C    C               zvstrdr,zvstrli,zvstrcl,
C    C               paprs,u,v)
C
C======================================================================
c Explicit Arguments:
c ==================
c iam-----input-I-File number where AAMs and torques are written
c                 It is a formatted file that has been opened
c                 in physiq.F
c nlon----input-I-Total number of horizontal points that get into physics
c nlev----input-I-Number of vertical levels
c rjour---input-R-Jour compte depuis le debut de la simu (run.def)
c rsec----input-R-Seconde de la journee
c rea-----input-R-Earth radius
c rg------input-R-gravity constant
c ome-----input-R-Earth rotation rate
c plat ---input-R-Latitude en degres
c plon ---input-R-Longitude en degres
c phis ---input-R-Geopotential at the ground
c dragu---input-R-orodrag stress (zonal)
c liftu---input-R-orolift stress (zonal)
c clu-----input-R-Boundary layer stress (zonal)
c dragv---input-R-orodrag stress (Meridional)
c liftv---input-R-orolift stress (Meridional)
c clv-----input-R-Boundary layer stress (Meridional)
c p-------input-R-Pressure (Pa) at model half levels
c u-------input-R-Horizontal wind (m/s)
c v-------input-R-Meridional wind (m/s)
c
c
c Implicit Arguments:
c ===================
c
c nbp_lon--common-I: Number of longitude intervals
c nbp_lat-1--common-I: Number of latitude intervals
c klon-common-I: Number of points seen by the physics
c                nbp_lon*(nbp_lat-1-1)+2 for instance
c klev-common-I: Number of vertical layers
c======================================================================
c Local Variables:
c ================
c dlat-----R: Latitude increment (Radians)
c dlon-----R: Longitude increment (Radians)
c raam  ---R: Wind AAM (3 Components, 1 & 2 Equatoriales; 3 Axiale)
c oaam  ---R: Mass AAM (3 Components, 1 & 2 Equatoriales; 3 Axiale)
c tmou-----R: Resolved Mountain torque (3 components)
c tsso-----R: Parameterised Moutain drag torque (3 components)
c tbls-----R: Parameterised Boundary layer torque (3 components)
c
c LOCAL ARRAY:
c ===========
c zs    ---R: Topographic height
c ps    ---R: Surface Pressure  
c ub    ---R: Barotropic wind zonal
c vb    ---R: Barotropic wind meridional
c zlat  ---R: Latitude in radians
c zlon  ---R: Longitude in radians
c======================================================================

c
c ARGUMENTS
c
      INTEGER iam,nlon,nlev
      REAL rjour,rsec,rea,rg,ome
      REAL plat(nlon),plon(nlon),phis(nlon)
      REAL dragu(nlon),liftu(nlon),clu(nlon)             
      REAL dragv(nlon),liftv(nlon),clv(nlon)             
      REAL p(nlon,nlev+1), u(nlon,nlev), v(nlon,nlev)
c
c Variables locales:
c
      INTEGER i,j,k,l
      REAL xpi,hadley,hadday
      REAL dlat,dlon
      REAL raam(3),oaam(3),tmou(3),tsso(3),tbls(3)
      integer iax


      REAL ZS(801,401),PS(801,401)
      REAL UB(801,401),VB(801,401)
      REAL SSOU(801,401),SSOV(801,401)
      REAL BLSU(801,401),BLSV(801,401)
      REAL ZLON(801),ZLAT(401)
C
C  PUT AAM QUANTITIES AT ZERO:
C
      if(nbp_lon+1.gt.801.or.nbp_lat.gt.401)then
      print *,' Pb de dimension dans aaam_bud'
      stop
      endif

      xpi=acos(-1.)
      hadley=1.e18
      hadday=1.e18*1.e7
      IF (klon_glo.EQ.1) THEN
        dlat=xpi
      ELSE
        dlat=xpi/float(nbp_lat-1)
      ENDIF
      dlon=2.*xpi/float(nbp_lon) 
      
      do iax=1,3
      oaam(iax)=0.
      raam(iax)=0.
      tmou(iax)=0.
      tsso(iax)=0.
      tbls(iax)=0.
      enddo

C MOUNTAIN HEIGHT, PRESSURE AND BAROTROPIC WIND:

C North pole values (j=1):
 
      l=1

        ub(1,1)=0.
        vb(1,1)=0.
        do k=1,nlev
          ub(1,1)=ub(1,1)+u(l,k)*(p(l,k)-p(l,k+1))/rg
          vb(1,1)=vb(1,1)+v(l,k)*(p(l,k)-p(l,k+1))/rg
        enddo

          zlat(1)=plat(l)*xpi/180.

        do i=1,nbp_lon+1

          zs(i,1)=phis(l)/rg
          ps(i,1)=p(l,1)
          ub(i,1)=ub(1,1)                             
          vb(i,1)=vb(1,1)                             
          ssou(i,1)=dragu(l)+liftu(l)
          ssov(i,1)=dragv(l)+liftv(l)
          blsu(i,1)=clu(l)
          blsv(i,1)=clv(l)

        enddo


      do j = 2,nbp_lat-1

C Values at Greenwich (Periodicity)

      zs(nbp_lon+1,j)=phis(l+1)/rg
      ps(nbp_lon+1,j)=p(l+1,1)
          ssou(nbp_lon+1,j)=dragu(l+1)+liftu(l+1)
          ssov(nbp_lon+1,j)=dragv(l+1)+liftv(l+1)
          blsu(nbp_lon+1,j)=clu(l+1)
          blsv(nbp_lon+1,j)=clv(l+1)
      zlon(nbp_lon+1)=-plon(l+1)*xpi/180.
      zlat(j)=plat(l+1)*xpi/180.

      ub(nbp_lon+1,j)=0.
      vb(nbp_lon+1,j)=0.
         do k=1,nlev
         ub(nbp_lon+1,j)=ub(nbp_lon+1,j)+u(l+1,k)*
     &                                   (p(l+1,k)-p(l+1,k+1))/rg
         vb(nbp_lon+1,j)=vb(nbp_lon+1,j)+v(l+1,k)*
     &                                   (p(l+1,k)-p(l+1,k+1))/rg
         enddo
      

      do i=1,nbp_lon

      l=l+1
      zs(i,j)=phis(l)/rg
      ps(i,j)=p(l,1)
          ssou(i,j)=dragu(l)+liftu(l)
          ssov(i,j)=dragv(l)+liftv(l)
          blsu(i,j)=clu(l)
          blsv(i,j)=clv(l)
      zlon(i)=plon(l)*xpi/180.

      ub(i,j)=0.
      vb(i,j)=0.
         do k=1,nlev
         ub(i,j)=ub(i,j)+u(l,k)*(p(l,k)-p(l,k+1))/rg
         vb(i,j)=vb(i,j)+v(l,k)*(p(l,k)-p(l,k+1))/rg
         enddo

      enddo

      enddo


C South Pole

      l=l+1
      ub(1,nbp_lat)=0.
      vb(1,nbp_lat)=0.
      do k=1,nlev
         ub(1,nbp_lat)=ub(1,nbp_lat)+u(l,k)*(p(l,k)-p(l,k+1))/rg
         vb(1,nbp_lat)=vb(1,nbp_lat)+v(l,k)*(p(l,k)-p(l,k+1))/rg
      enddo
      zlat(nbp_lat)=plat(l)*xpi/180.

      do i=1,nbp_lon+1
      zs(i,nbp_lat)=phis(l)/rg
      ps(i,nbp_lat)=p(l,1)
          ssou(i,nbp_lat)=dragu(l)+liftu(l)
          ssov(i,nbp_lat)=dragv(l)+liftv(l)
          blsu(i,nbp_lat)=clu(l)
          blsv(i,nbp_lat)=clv(l)
      ub(i,nbp_lat)=ub(1,nbp_lat)                               
      vb(i,nbp_lat)=vb(1,nbp_lat)                                
      enddo

C
C  MOMENT ANGULAIRE 
C
        DO j=1,nbp_lat-1    
        DO i=1,nbp_lon

           raam(1)=raam(1)-rea**3*dlon*dlat*0.5*
     c    (cos(zlon(i  ))*sin(zlat(j  ))*cos(zlat(j  ))*ub(i  ,j  )
     c    +cos(zlon(i  ))*sin(zlat(j+1))*cos(zlat(j+1))*ub(i  ,j+1))
     c                    +rea**3*dlon*dlat*0.5*
     c    (sin(zlon(i  ))*cos(zlat(j  ))*vb(i  ,j  )
     c    +sin(zlon(i  ))*cos(zlat(j+1))*vb(i  ,j+1))

           oaam(1)=oaam(1)-ome*rea**4*dlon*dlat/rg*0.5*
     c   (cos(zlon(i  ))*cos(zlat(j  ))**2*sin(zlat(j  ))*ps(i  ,j  )
     c   +cos(zlon(i  ))*cos(zlat(j+1))**2*sin(zlat(j+1))*ps(i  ,j+1))

           raam(2)=raam(2)-rea**3*dlon*dlat*0.5*
     c    (sin(zlon(i  ))*sin(zlat(j  ))*cos(zlat(j  ))*ub(i  ,j  )
     c    +sin(zlon(i  ))*sin(zlat(j+1))*cos(zlat(j+1))*ub(i  ,j+1))
     c                    -rea**3*dlon*dlat*0.5*
     c    (cos(zlon(i  ))*cos(zlat(j  ))*vb(i  ,j  )
     c    +cos(zlon(i  ))*cos(zlat(j+1))*vb(i  ,j+1))

           oaam(2)=oaam(2)-ome*rea**4*dlon*dlat/rg*0.5*
     c   (sin(zlon(i  ))*cos(zlat(j  ))**2*sin(zlat(j  ))*ps(i  ,j  )
     c   +sin(zlon(i  ))*cos(zlat(j+1))**2*sin(zlat(j+1))*ps(i  ,j+1))

           raam(3)=raam(3)+rea**3*dlon*dlat*0.5*
     c           (cos(zlat(j))**2*ub(i,j)+cos(zlat(j+1))**2*ub(i,j+1))

           oaam(3)=oaam(3)+ome*rea**4*dlon*dlat/rg*0.5*
     c        (cos(zlat(j))**3*ps(i,j)+cos(zlat(j+1))**3*ps(i,j+1))

        ENDDO
        ENDDO

C
C COUPLE DES MONTAGNES:
C

        DO j=1,nbp_lat-1
        DO i=1,nbp_lon
           tmou(1)=tmou(1)-rea**2*dlon*0.5*sin(zlon(i))
     c  *(zs(i,j)-zs(i,j+1))
     c  *(cos(zlat(j+1))*ps(i,j+1)+cos(zlat(j))*ps(i,j)) 
           tmou(2)=tmou(2)+rea**2*dlon*0.5*cos(zlon(i))
     c  *(zs(i,j)-zs(i,j+1))
     c  *(cos(zlat(j+1))*ps(i,j+1)+cos(zlat(j))*ps(i,j)) 
        ENDDO
        ENDDO
           
        DO j=2,nbp_lat-1 
        DO i=1,nbp_lon
           tmou(1)=tmou(1)+rea**2*dlat*0.5*sin(zlat(j))
     c  *(zs(i+1,j)-zs(i,j))
     c  *(cos(zlon(i+1))*ps(i+1,j)+cos(zlon(i))*ps(i,j))
           tmou(2)=tmou(2)+rea**2*dlat*0.5*sin(zlat(j))
     c  *(zs(i+1,j)-zs(i,j))
     c  *(sin(zlon(i+1))*ps(i+1,j)+sin(zlon(i))*ps(i,j))
           tmou(3)=tmou(3)-rea**2*dlat*0.5*
     c  cos(zlat(j))*(zs(i+1,j)-zs(i,j))*(ps(i+1,j)+ps(i,j))
        ENDDO
        ENDDO
C
C COUPLES DES DIFFERENTES FRICTION AU SOL:
C
        l=1
        DO j=2,nbp_lat-1
        DO i=1,nbp_lon
        l=l+1
           tsso(1)=tsso(1)-rea**3*cos(zlat(j))*dlon*dlat*
     c     ssou(i,j)          *sin(zlat(j))*cos(zlon(i))
     c                    +rea**3*cos(zlat(j))*dlon*dlat*
     c     ssov(i,j)          *sin(zlon(i))

           tsso(2)=tsso(2)-rea**3*cos(zlat(j))*dlon*dlat*
     c     ssou(i,j)          *sin(zlat(j))*sin(zlon(i))
     c                    -rea**3*cos(zlat(j))*dlon*dlat*
     c     ssov(i,j)          *cos(zlon(i))

           tsso(3)=tsso(3)+rea**3*cos(zlat(j))*dlon*dlat*
     c     ssou(i,j)          *cos(zlat(j))

           tbls(1)=tbls(1)-rea**3*cos(zlat(j))*dlon*dlat*
     c     blsu(i,j)          *sin(zlat(j))*cos(zlon(i))
     c                    +rea**3*cos(zlat(j))*dlon*dlat*
     c     blsv(i,j)          *sin(zlon(i))

           tbls(2)=tbls(2)-rea**3*cos(zlat(j))*dlon*dlat*
     c     blsu(i,j)          *sin(zlat(j))*sin(zlon(i))
     c                    -rea**3*cos(zlat(j))*dlon*dlat*
     c     blsv(i,j)          *cos(zlon(i))

           tbls(3)=tbls(3)+rea**3*cos(zlat(j))*dlon*dlat*
     c     blsu(i,j)          *cos(zlat(j))

        ENDDO
        ENDDO
            

c     write(*,*) 'AAM',rsec,
c    c      raam(3)/hadday,oaam(3)/hadday,
c    c      tmou(3)/hadley,tsso(3)/hadley,tbls(3)/hadley

      write(iam,100)rjour+rsec/1.e7,
c    c      raam(1)/hadday,oaam(1)/hadday,
c    c      tmou(1)/hadley,tsso(1)/hadley,tbls(1)/hadley,
c    c      raam(2)/hadday,oaam(2)/hadday,
c    c      tmou(2)/hadley,tsso(2)/hadley,tbls(2)/hadley,
     c      raam(3)/hadday,oaam(3)/hadday,
     c      tmou(3)/hadley,tsso(3)/hadley,tbls(3)/hadley 
c100   format(F12.5,15(1x,F12.5))
100   format(F12.5,5(1x,F12.5))

      write(iam+1,*)((zs(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((ps(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((ub(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((vb(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((ssou(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((ssov(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((blsu(i,j),i=1,nbp_lon),j=1,nbp_lat)
      write(iam+1,*)((blsv(i,j),i=1,nbp_lon),j=1,nbp_lat)

      RETURN
      END
