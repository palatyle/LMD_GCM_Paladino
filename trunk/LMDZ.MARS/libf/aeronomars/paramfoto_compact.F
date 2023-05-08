c**********************************************************************

      subroutine paramfoto_compact
     $(ig,nlayer,chemthermod,lswitch,tx,timestep,zenit,zx,rm,nesptherm)

c     Main thermospheric photochemistry routine.
  
c     may 2008      FGG+MALV,GG             
c**********************************************************************

      use iono_h
      use param_v4_h
      implicit none

c     arguments

      integer lswitch,ig,nesptherm,chemthermod,nlayer
      real zdens(nlayer)
      real tx(nlayer)
      real zenit
      real zx(nlayer)
      real rm(nlayer,nesptherm)
      real timestep


c     local variables

      real*8  deltat,timefrac_sec
      real*8  co2xnew,o2xnew,o3pxnew,coxnew,hxnew,ohxnew
      real*8  ho2xnew,h2xnew
      real*8  h2o2xnew,o1dxnew,o3xnew,h2oxnew
      real*8  noxnew, nxnew, n2xnew, n2dxnew,no2xnew
      real*8  oplusxnew, o2plusxnew,co2plusxnew
      real*8  nplusxnew, n2plusxnew,noplusxnew, hplusxnew
      real*8  electxnew,coplusxnew, cplusxnew, hco2plusxnew

      real*8  co2xoutput,o2xoutput,o3pxoutput,coxoutput
      real*8  ho2xoutput,h2xoutput,hxoutput,ohxoutput
      real*8  h2o2xoutput,o1dxoutput,o3xoutput,h2oxoutput
      real*8  nxoutput,noxoutput,n2xoutput,n2dxoutput,no2xoutput
      real*8  co2plusxoutput,coplusxoutput,oplusxoutput,o2plusxoutput
      real*8  cplusxoutput,noplusxoutput,n2plusxoutput,hplusxoutput
      real*8  electxoutput,nplusxoutput,hco2plusxoutput
      real*8  electxoutput_timemarching, electxoutput_neutrality 

      real*8  co2xinput,o2xinput,o3pxinput,coxinput
      real*8  ho2xinput,h2xinput,hxinput,ohxinput
      real*8  h2o2xinput,o1dxinput,o3xinput,h2oxinput
      real*8  nxinput,noxinput,n2xinput,n2dxinput,no2xinput
      real*8  co2plusxinput,coplusxinput,oplusxinput,o2plusxinput
      real*8  cplusxinput,noplusxinput,n2plusxinput,hplusxinput
      real*8  electxinput,nplusxinput,hco2plusxinput

      real*8  co2xini,o2xini,o3pxini,coxini
      real*8  ho2xini,h2xini,hxini,ohxini
      real*8  h2o2xini,o1dxini,o3xini,h2oxini
      real*8  nxini,noxini,n2xini,n2dxini,no2xini
      real*8  co2plusxini,coplusxini,oplusxini,o2plusxini
      real*8  cplusxini,noplusxini,n2plusxini,hplusxini
      real*8  electxini,nplusxini,hco2plusxini

      real*8  dco2x,do2x,do3px,dcox,dhx,dohx,dho2x,dh2x
      real*8  dh2ox,dh2o2x,do1dx,do3x,dnx,dnox,dn2x,dn2dx,dno2x
      real*8  dco2plusx,dcoplusx,doplusx,do2plusx 
      real*8  dcplusx,dnoplusx,dn2plusx,dhplusx,dhco2plusx
      real*8  delectx,dnplusx

      real*8  jdistot8(nabs,nlayer)
      real*8  jdistot8_b(nabs,nlayer)
      real*8  jion8(nabs,nlayer,4)
      real*8  tx8
	  
      

      real*8  alfa_laststep, IonMostAbundant

      real*8  tmin(nlayer)
      real*8  fmargin1,critere

      integer compmin(nlayer)
      integer i,j,k
      integer numpasos
      integer n_comp_en_EQ(nlayer), paso

! Tracer indexes in the thermospheric chemistry:
!!! ATTENTION. These values have to be identical to those in chemthermos.F90
!!! If the values are changed there, the same has to be done here  !!!
!      integer,parameter :: i_co2=1
!      integer,parameter :: i_o2=2
!      integer,parameter :: i_o=3
!      integer,parameter :: i_co=4
!      integer,parameter :: i_h=5
!      integer,parameter :: i_oh=6
!      integer,parameter :: i_ho2=7
!      integer,parameter :: i_h2=8
!      integer,parameter :: i_h2o=9
!      integer,parameter :: i_h2o2=10
!      integer,parameter :: i_o1d=11
!      integer,parameter :: i_o3=12
!      integer,parameter :: i_n2=13
!      integer,parameter :: i_n=14
!      integer,parameter :: i_no=15
!      integer,parameter :: i_n2d=16
!      integer,parameter :: i_no2=17
      integer,parameter :: i_co2  =  1
      integer,parameter :: i_co   =  2
      integer,parameter :: i_o    =  3
      integer,parameter :: i_o1d  =  4
      integer,parameter :: i_o2   =  5
      integer,parameter :: i_o3   =  6
      integer,parameter :: i_h    =  7
      integer,parameter :: i_h2   =  8
      integer,parameter :: i_oh   =  9
      integer,parameter :: i_ho2  = 10
      integer,parameter :: i_h2o2 = 11
      integer,parameter :: i_h2o  = 12
      integer,parameter :: i_n    = 13
      integer,parameter :: i_n2d  = 14
      integer,parameter :: i_no   = 15
      integer,parameter :: i_no2  = 16
      integer,parameter :: i_n2   = 17
      integer,parameter :: i_co2plus=18
      integer,parameter :: i_oplus=19
      integer,parameter :: i_o2plus=20
      integer,parameter :: i_coplus=21
      integer,parameter :: i_cplus=22
      integer,parameter :: i_nplus=23
      integer,parameter :: i_noplus=24
      integer,parameter :: i_n2plus=25
      integer,parameter :: i_hplus=26
      integer,parameter :: i_hco2plus=27
      integer,parameter :: i_elec=28

c     formats

c**********************************************************************


c     external timestep
      timefrac_sec=dble(timestep)

C     Start: altitude loop
      do i=nlayer,lswitch,-1
c     Temperature and concentrations to real*8
         tx8=dble(tx(i))
         co2xini=dble(rm(i,i_co2))
         o2xini=dble(rm(i,i_o2))
         o3pxini=dble(rm(i,i_o))
         coxini=dble(rm(i,i_co))
         hxini=dble(rm(i,i_h))
         ohxini=dble(rm(i,i_oh))
         ho2xini=dble(rm(i,i_ho2))
         h2xini=dble(rm(i,i_h2))
         h2oxini=dble(rm(i,i_h2o))
         h2o2xini=dble(rm(i,i_h2o2))
         o1dxini=dble(rm(i,i_o1d))
         !Only if O3, N or ion chemistry requested
         if(chemthermod.ge.1) o3xini=dble(rm(i,i_o3))
         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
            n2xini=dble(rm(i,i_n2))
            nxini=dble(rm(i,i_n))
            noxini=dble(rm(i,i_no))
            n2dxini=dble(rm(i,i_n2d))
            no2xini=dble(rm(i,i_no2))
         endif
         !Only if ion chemistry requested
         if(chemthermod.eq.3) then
            co2plusxini=dble(rm(i,i_co2plus))
            oplusxini=dble(rm(i,i_oplus))
            o2plusxini=dble(rm(i,i_o2plus))
            coplusxini=dble(rm(i,i_coplus))
            cplusxini=dble(rm(i,i_cplus))
            nplusxini=dble(rm(i,i_nplus))
            n2plusxini=dble(rm(i,i_n2plus))
            noplusxini=dble(rm(i,i_noplus))
            hplusxini=dble(rm(i,i_hplus))
            hco2plusxini=dble(rm(i,i_hco2plus))
            electxini=dble(rm(i,i_elec))
         endif

         !Calculation of photodissociation and photoionization rates
         !from photoabsorption rates and ionization-to-dissociation 
         !branching ratios
         call phdisrate(ig,nlayer,chemthermod,zenit,i)   
         ! Conversion to double precision
         do j=1,nabs
	    jdistot8(j,i) = dble(jdistot(j,i))
            jdistot8_b(j,i) = dble(jdistot_b(j,i))
            do k=1,4
               jion8(j,i,k)=dble(jion(j,i,k))
            enddo
         end do

         !Reaction rates
         call getch( ig, chemthermod,tx8, zx(i))
		 
         !Lifetimes and temporal integration 
         call lifetimes(ig,i,nlayer,chemthermod,zenit,zx,
     $        jdistot8,jdistot8_b,jion8,
     $        tmin(i),compmin(i), 
     $        n_comp_en_EQ(i),co2xini,o2xini,o3pxini,coxini,hxini,
     $        ohxini,ho2xini,h2xini,h2oxini,h2o2xini,o1dxini,o3xini,
     $        n2xini,nxini,noxini,no2xini,n2dxini,co2plusxini,oplusxini,
     $        o2plusxini,coplusxini,cplusxini,nplusxini,noplusxini,
     $        n2plusxini,hplusxini,hco2plusxini,electxini )

         !Calculation of the internal timestep and revision of the 
         !validity of the photochemical equilibrium approximation 
         !for each species

! JYC criteria added to avoid instabilities in (H) + (O+) <-> (H+) + (O) reactions when H+ is important
	fmargin1=5
        !Only if ion chemistry requested
        if(chemthermod.eq.3) then
           critere=hplusxini/(o3pxini+hxini+h2xini)
           if (critere .gt. 5d-4) then
              fmargin1=2000.*critere
              if (fmargin1 .gt. 50.) fmargin1=50
           endif
        endif   !Of chemthermod.eq.3

         call timemarching ( ig,i,nlayer,chemthermod,n_comp_en_EQ,
     .       compmin,tmin,timefrac_sec, deltat,fmargin1)

         !Number of timesteps
         numpasos = int( timefrac_sec / deltat ) 
         alfa_laststep = 1.d0 + timefrac_sec/deltat - dble(numpasos)
         do paso=1,numpasos  

            !Concentrations at the first step
            if(paso.eq.1) then
               co2xinput=co2xini
               o2xinput=o2xini
               o3pxinput=o3pxini
               coxinput=coxini
               hxinput=hxini
               ohxinput=ohxini
               ho2xinput=ho2xini
               h2xinput=h2xini
               h2oxinput=h2oxini
               h2o2xinput=h2o2xini
               o1dxinput=o1dxini
               o3xinput=o3xini
               nxinput=nxini
               noxinput=noxini
               n2xinput=n2xini
               n2dxinput=n2dxini
               no2xinput=no2xini
               !
	       co2plusxinput = co2plusxini
	       oplusxinput   = oplusxini
               o2plusxinput  = o2plusxini
               coplusxinput  = coplusxini
               cplusxinput   = cplusxini
               nplusxinput   = nplusxini
               n2plusxinput  = n2plusxini
               noplusxinput  = noplusxini
               hplusxinput   = hplusxini
               hco2plusxinput= hco2plusxini
               electxinput   = electxini
            else
               !Concentrations for the new step
               co2xinput=co2xinput+dco2x
               o2xinput=o2xinput+do2x
               o3pxinput=o3pxinput+do3px
               coxinput=coxinput+dcox
               hxinput=hxinput+dhx
               ohxinput=ohxinput+dohx
               ho2xinput=ho2xinput+dho2x
               h2xinput=h2xinput+dh2x
               h2oxinput=h2oxinput+dh2ox
               h2o2xinput=h2o2xinput+dh2o2x
               o1dxinput=o1dxinput+do1dx
               !Only if O3, N or ion chemistry requested
               if(chemthermod.ge.1) o3xinput=o3xinput+do3x
               !Only if N or ion chemistry requested
               if(chemthermod.ge.2) then
                  nxinput=nxinput+dnx
                  noxinput=noxinput+dnox
                  n2xinput=n2xinput+dn2x
                  n2dxinput=n2dxinput+dn2dx
                  no2xinput=no2xinput+dno2x
               endif
	       !Only if ion chemistry requested
               if(chemthermod.eq.3) then
                  co2plusxinput = co2plusxinput + dco2plusx
                  oplusxinput   = oplusxinput   + doplusx
                  o2plusxinput  = o2plusxinput  + do2plusx
                  coplusxinput  = coplusxinput  + dcoplusx
                  cplusxinput   = cplusxinput   + dcplusx
                  nplusxinput   = nplusxinput   + dnplusx
                  n2plusxinput  = n2plusxinput  + dn2plusx
                  noplusxinput  = noplusxinput  + dnoplusx
                  hplusxinput   = hplusxinput   + dhplusx
                  hco2plusxinput= hco2plusxinput+ dhco2plusx
                  electxinput   = electxinput   + delectx
               endif

            end if
            !Calculation of productions and losses
            call prodsandlosses (ig,i,nlayer,chemthermod,zenit, zx,
     &                 jdistot8, jdistot8_b, jion8,
     &                              co2xinput, o2xinput, o3pxinput,
     &                              coxinput,  h2xinput, o3xinput,
     &                              h2oxinput, nxinput,  noxinput, 
     &                              h2o2xinput, n2xinput, 
     &                           o1dxinput, ohxinput,  ho2xinput,
     &                           hxinput,   n2dxinput, no2xinput,
     &                 co2plusxinput,  o2plusxinput, coplusxinput,
     &                 oplusxinput,    cplusxinput,  noplusxinput,
     &                 n2plusxinput,   hplusxinput,  nplusxinput,
     &                 hco2plusxinput,electxinput )


            !New abundances, implicit scheme for the timemarching

            !First, for the 11 species that can not be in PE
            !( CO2, O2, O3P, CO, H2, H2O, H2O2, O3, N, NO, N2 )

            call implicito ( ig, co2xoutput,  ' CO2', 
     &             co2xinput, Pco2tot(i), Lco2tot(i), deltat )
            call implicito ( ig, o2xoutput,   '  O2', 
     &             o2xinput, Po2tot(i), Lo2tot(i), deltat )
            call implicito ( ig, o3pxoutput,  ' O3P', 
     &             o3pxinput, Po3ptot(i), Lo3ptot(i), deltat )
            call implicito ( ig, coxoutput,   '  CO', 
     &             coxinput, Pcotot(i), Lcotot(i), deltat )
            call implicito ( ig, h2xoutput,   '  H2', 
     &             h2xinput, Ph2tot(i), Lh2tot(i), deltat )
            call implicito ( ig, h2oxoutput,  ' H2O', 
     &             h2oxinput, Ph2otot(i), Lh2otot(i), deltat )
            call implicito ( ig, h2o2xoutput, 'H2O2', 
     &             h2o2xinput, Ph2o2tot(i), Lh2o2tot(i), deltat )
            !only if O3, N or ion chemistry requested
            if(chemthermod.ge.1) 
     $           call implicito ( ig, o3xoutput,   '  O3', 
     &             o3xinput, Po3tot(i), Lo3tot(i), deltat )
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               call implicito ( ig, nxoutput,    '   N', 
     &              nxinput, Pntot(i), Lntot(i), deltat )
               call implicito ( ig, noxoutput,   '  NO', 
     &              noxinput, Pnotot(i), Lnotot(i), deltat )
               call implicito ( ig, n2xoutput,   '  N2', 
     &              n2xinput, Pn2tot(i), Ln2tot(i), deltat )
            endif


            !Second, 6+10 species that can be in PE, but are not
            ! 6 neutral , O1D, OH, HO2, H, N2D, NO2
            if(o1d_eq(i).eq.'N') then
               call implicito ( ig, o1dxoutput,   ' O1D', 
     &             o1dxinput, Po1dtot(i), Lo1dtot(i), deltat )
            end if
            if(oh_eq(i).eq.'N') then
               call implicito ( ig, ohxoutput,   '  OH', 
     &             ohxinput, Pohtot(i), Lohtot(i), deltat )
            end if
            if(ho2_eq(i).eq.'N') then
               call implicito ( ig, ho2xoutput,   ' HO2', 
     &             ho2xinput, Pho2tot(i), Lho2tot(i), deltat )
            end if
            if(h_eq(i).eq.'N') then
               call implicito ( ig, hxoutput,   '   H', 
     &             hxinput, Phtot(i), Lhtot(i), deltat )
            end if
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if(n2d_eq(i).eq.'N') then
                  call implicito ( ig, n2dxoutput,   ' N2D', 
     &                 n2dxinput, Pn2dtot(i), Ln2dtot(i), deltat )
               end if
               if(no2_eq(i).eq.'N') then
                  call implicito ( ig, no2xoutput,   ' NO2', 
     &                 no2xinput, Pno2tot(i), Lno2tot(i), deltat )
               end if
            endif

            ! 9 ions (all of them) and electrons
            !Only if ion chemistry requested
            if(chemthermod.ge.3) then
               if(n2plus_eq(i).eq.'N') then
                  call implicito ( ig, n2plusxoutput,   ' N2+', 
     &                 n2plusxinput,Pn2plustot(i),Ln2plustot(i),deltat)
               end if
               if(cplus_eq(i).eq.'N') then
                  call implicito ( ig, cplusxoutput,   '  C+', 
     &                 cplusxinput,Pcplustot(i),Lcplustot(i),deltat)
               end if
               if(coplus_eq(i).eq.'N') then
                  call implicito ( ig, coplusxoutput,   ' CO+', 
     &                 coplusxinput,Pcoplustot(i),Lcoplustot(i),deltat)
               end if
               if(co2plus_eq(i).eq.'N') then
                  call implicito ( ig, co2plusxoutput,   'CO2+', 
     &               co2plusxinput,Pco2plustot(i),Lco2plustot(i),deltat)
               end if
               if(oplus_eq(i).eq.'N') then
                  call implicito ( ig, oplusxoutput,   '  O+', 
     &                 oplusxinput,Poplustot(i),Loplustot(i),deltat)
               end if
               if(hplus_eq(i).eq.'N') then 
                  call implicito ( ig, hplusxoutput,   '  H+', 
     &                 hplusxinput,Phplustot(i),Lhplustot(i),deltat)
               end if            
               if(o2plus_eq(i).eq.'N') then
                  call implicito ( ig, o2plusxoutput,   ' O2+', 
     &                 o2plusxinput,Po2plustot(i),Lo2plustot(i),deltat)
               end if
               if(noplus_eq(i).eq.'N') then
                  call implicito ( ig, noplusxoutput,   ' NO+', 
     &                 noplusxinput,Pnoplustot(i),Lnoplustot(i),deltat)
               end if
               if(nplus_eq(i).eq.'N') then
                  call implicito ( ig, nplusxoutput,   '  N+', 
     &                 nplusxinput,Pnplustot(i),Lnplustot(i),deltat)
               end if
               if(hco2plus_eq(i).eq.'N') then
                  call implicito ( ig, hco2plusxoutput,   'CO2+', 
     &               hco2plusxinput,Phco2plustot(i),Lhco2plustot(i),
     $                 deltat)
               end if
           ! elect
            call implicito ( ig, electxoutput_timemarching,   'elec', 
     &              electxinput,Pelecttot(i),Lelecttot(i),deltat)
            endif !Of chemthermod.eq.3


            !Third, those species (among the 16 that can be in PE) that are in PE
            call EF_oscilacion
     &           ( ig,i,nlayer, paso,chemthermod,zenit, zx,
     &           jdistot8, jdistot8_b,jion8,
     &           deltat,
     $           co2xoutput,     co2xinput,
     $           o2xoutput,      o2xinput,
     $           o3pxoutput,     o3pxinput,
     $           coxoutput,      coxinput,
     $           h2xoutput,      h2xinput,
     $           h2oxoutput,     h2oxinput,
     $           h2o2xoutput,    h2o2xinput,
     $           o3xoutput,      o3xinput,
     $           nxoutput,       nxinput,
     $           noxoutput,      noxinput,
     $           n2xoutput,      n2xinput,
     &           o1dxoutput, o1dxinput, 
     &           ohxoutput,  ohxinput, 
     &           ho2xoutput, ho2xinput, 
     &           hxoutput,   hxinput, 
     &           n2dxoutput,  n2dxinput, 
     &           no2xoutput, no2xinput, 
     &           co2plusxoutput, co2plusxinput, 
     &           o2plusxoutput,  o2plusxinput, 
     &           coplusxoutput,  coplusxinput, 
     &           oplusxoutput,   oplusxinput, 
     &           cplusxoutput,   cplusxinput, 
     &           noplusxoutput,  noplusxinput, 
     &           n2plusxoutput,  n2plusxinput, 
     &           hplusxoutput,   hplusxinput, 
     &           nplusxoutput,   nplusxinput,
     $           hco2plusxoutput,hco2plusxinput,
     &           electxoutput,   electxinput, 
     &           electxoutput_timemarching )

            !Electrons given by the condition of global neutrality
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               electxoutput = o2plusxoutput +
     @              co2plusxoutput +
     @              coplusxoutput +
     @              oplusxoutput +
     @              cplusxoutput +
     @              n2plusxoutput +
     @              nplusxoutput +
     @              noplusxoutput +
     @              hplusxoutput +
     $              hco2plusxoutput
               electxoutput_neutrality = electxoutput
                                !
               IonMostAbundant = o2plusxoutput
               IonMostAbundant = max( co2plusxoutput, IonMostAbundant)
               IonMostAbundant = max( coplusxoutput, IonMostAbundant)
               IonMostAbundant = max( oplusxoutput, IonMostAbundant)
               IonMostAbundant = max( cplusxoutput, IonMostAbundant)
               IonMostAbundant = max( n2plusxoutput, IonMostAbundant)
               IonMostAbundant = max( noplusxoutput, IonMostAbundant)
               IonMostAbundant = max( nplusxoutput, IonMostAbundant)
               IonMostAbundant = max( hplusxoutput, IonMostAbundant)
               IonMostAbundant = max( hco2plusxoutput, IonMostAbundant)
               IonMostAbundant = IonMostAbundant / electxoutput
            endif !Of chemthermod.eq.3

            !Concentration changes for this time step
            dco2x=co2xoutput-co2xinput
            do2x=o2xoutput-o2xinput
            do3px=o3pxoutput-o3pxinput
            dcox=coxoutput-coxinput
            dhx=hxoutput-hxinput
            dohx=ohxoutput-ohxinput
            dho2x=ho2xoutput-ho2xinput
            dh2x=h2xoutput-h2xinput
            dh2ox=h2oxoutput-h2oxinput
            dh2o2x=h2o2xoutput-h2o2xinput
            do1dx=o1dxoutput-o1dxinput
            !Only if O3, N or ion chemistry requested
            if(chemthermod.ge.1) do3x=o3xoutput-o3xinput
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               dnx=nxoutput-nxinput
               dnox=noxoutput-noxinput
               dn2x=n2xoutput-n2xinput
               dn2dx=n2dxoutput-n2dxinput
               dno2x=no2xoutput-no2xinput
            endif
	    !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               dco2plusx=co2plusxoutput-co2plusxinput
               do2plusx=o2plusxoutput-o2plusxinput
               doplusx=oplusxoutput-oplusxinput
               dcoplusx=coplusxoutput-coplusxinput
               dcplusx=cplusxoutput-cplusxinput
               dnplusx=nplusxoutput-nplusxinput
               dn2plusx=n2plusxoutput-n2plusxinput
               dnoplusx=noplusxoutput-noplusxinput
               dhplusx=hplusxoutput-hplusxinput
               dhco2plusx=hco2plusxoutput-hco2plusxinput
               delectx=electxoutput- electxinput
            endif
            if(paso.eq.numpasos) then		
               !Final concentrations after last time step
               co2xnew = co2xinput +  dco2x * alfa_laststep
               if(co2xnew.lt.0)co2xnew=1.e-30
               o2xnew = o2xinput +  do2x * alfa_laststep
               if(o2xnew.lt.0)o2xnew=1.e-30
               o3pxnew = o3pxinput +  do3px * alfa_laststep
               if(o3pxnew.lt.0)o3pxnew=1.e-30
               coxnew = coxinput +  dcox * alfa_laststep
               if(coxnew.lt.0)coxnew=1.e-30
               hxnew = hxinput +  dhx * alfa_laststep
               if(hxnew.lt.0)hxnew=1.e-30
               ohxnew = ohxinput +  dohx * alfa_laststep
               if(ohxnew.lt.0)ohxnew=1.e-30
               ho2xnew = ho2xinput +  dho2x * alfa_laststep
               if(ho2xnew.lt.0)ho2xnew=1.e-30
               h2xnew = h2xinput +  dh2x * alfa_laststep
               if(h2xnew.lt.0)h2xnew=1.e-30
               h2oxnew = h2oxinput +  dh2ox * alfa_laststep
               if(h2oxnew.lt.0)h2oxnew=1.e-30
               h2o2xnew = h2o2xinput +  dh2o2x * alfa_laststep
               if(h2o2xnew.lt.0)h2o2xnew=1.e-30
               o1dxnew = o1dxinput +  do1dx * alfa_laststep
               if(o1dxnew.lt.0)o1dxnew=1.e-30
               !Only if O3, N or ion chemistry requested
               if(chemthermod.ge.1) then
                  o3xnew = o3xinput +  do3x * alfa_laststep
                  if(o3xnew.lt.0)o3xnew=1.e-30
               endif
               !Only if N or ion chemistry requested
               if(chemthermod.ge.2) then
                  nxnew = nxinput +  dnx * alfa_laststep
                  if(nxnew.lt.0)nxnew=1.e-30
                  noxnew = noxinput +  dnox * alfa_laststep
                  if(noxnew.lt.0)noxnew=1.e-30
                  n2xnew = n2xinput +  dn2x * alfa_laststep
                  if(n2xnew.lt.0)n2xnew=1.e-30
                  n2dxnew = n2dxinput +  dn2dx * alfa_laststep
                  if(n2dxnew.lt.0)n2dxnew=1.e-30
                  no2xnew = no2xinput +  dno2x * alfa_laststep
                  if(no2xnew.lt.0)no2xnew=1.e-30
               endif
               !Only if ion chemistry requested
               if(chemthermod.ge.3) then
                  co2plusxnew = co2plusxinput+dco2plusx*alfa_laststep
                  if(co2plusxnew.lt.0)co2plusxnew=1.e-30
                  o2plusxnew = o2plusxinput+do2plusx*alfa_laststep
                  if(o2plusxnew.lt.0)o2plusxnew=1.e-30
                  oplusxnew = oplusxinput+doplusx*alfa_laststep
                  if(oplusxnew.lt.0)oplusxnew=1.e-30
                  coplusxnew = coplusxinput+dcoplusx*alfa_laststep
                  if(coplusxnew.lt.0)coplusxnew=1.e-30
                  nplusxnew = nplusxinput +dnplusx*alfa_laststep
                  if(nplusxnew.lt.0)nplusxnew=1.e-30
                  n2plusxnew = n2plusxinput+dn2plusx*alfa_laststep
                  if(n2plusxnew.lt.0)n2plusxnew=1.e-30
                  noplusxnew = noplusxinput+dnoplusx*alfa_laststep
                  if(noplusxnew.lt.0)noplusxnew=1.e-30
                  hplusxnew = hplusxinput+dhplusx*alfa_laststep
                  if(hplusxnew.lt.0)hplusxnew=1.e-30
                  cplusxnew = cplusxinput+dcplusx*alfa_laststep
                  if(cplusxnew.lt.0)cplusxnew=1.e-30
                  hco2plusxnew = hco2plusxinput+dhco2plusx*alfa_laststep
                  if(hco2plusxnew.lt.0)hco2plusxnew=1.e-30
                  electxnew = electxinput+delectx*alfa_laststep
                  if(electxnew.lt.0)electxnew=1.e-30
               endif    !Of chemthermod.ge.3
            endif       !Of paso.eq.numpasos


         end do   
         
         !New concentrations to be returned
         rm(i,i_co2)     = real(co2xnew)
         rm(i,i_o2)      = real(o2xnew)
         rm(i,i_o)       = real(o3pxnew)
         rm(i,i_co)      = real(coxnew)
         rm(i,i_h)       = real(hxnew)
         rm(i,i_oh)      = real(ohxnew)
         rm(i,i_ho2)     = real(ho2xnew)
         rm(i,i_h2)      = real(h2xnew)
         rm(i,i_h2o)     = real(h2oxnew)
         rm(i,i_h2o2)    = real(h2o2xnew)
         rm(i,i_o1d)     = real(o1dxnew)
         !Only if O3, N or ion chemistry requested
         if(chemthermod.ge.1) 
     $        rm(i,i_o3) = real(o3xnew)
         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
            rm(i,i_n)    = real(nxnew)
            rm(i,i_n2)   = real(n2xnew)
            rm(i,i_no)   = real(noxnew)
            rm(i,i_n2d)  = real(n2dxnew)
            rm(i,i_no2)  = real(no2xnew)
         endif
         !Only if ion chemistry requested
         if(chemthermod.eq.3) then
            rm(i,i_co2plus) = real(co2plusxnew)
            rm(i,i_oplus)   = real(oplusxnew)
            rm(i,i_o2plus)  = real(o2plusxnew)
            rm(i,i_coplus)  = real(coplusxnew)
            rm(i,i_cplus)   = real(cplusxnew)
            rm(i,i_nplus)   = real(nplusxnew)
            rm(i,i_n2plus)  = real(n2plusxnew)
            rm(i,i_noplus)  = real(noplusxnew)
            rm(i,i_hplus)   = real(hplusxnew)
            rm(i,i_hco2plus)= real(hco2plusxnew)
            rm(i,i_elec)    = real(electxnew)
         endif
      end do     
cccccc End altitude loop 

      return


      end


c**********************************************************************
c**********************************************************************

      subroutine implicito ( ig,c_output, text, 
     &                      c_input, Prod, Loss, tstep )

  
c Given the productions and losses, calculates new concentrations using
c an implicit integration scheme. Checks if there are negative values and
c avoids underflows.

c     jul 2008      MA        Version en subrutina
c**********************************************************************

      implicit none

c     arguments
c
      integer   ig
      real*8    c_output                      ! O.
      real*8    c_input                            ! I.
      real*8    tstep                              ! I. 
      real*8    Prod                               ! I.
      real*8    Loss                               ! I.
      character*4 text                             ! I.

ccccccccccccccc CODE STARTS 

         c_output = (c_input + Prod * tstep) / (1.d0 + Loss * tstep)

         ! Stop is negative prods, losses, concentrations or times
         !
         if ( c_output.lt.0.d0 ) then 
              write(*,*)  text//' < 0 !!!'
              write (*,*) '  Terms of the implicit equation: '
              write (*,*) '   c_input =', c_input
              write (*,*) '   Prod = ', Prod
              write (*,*) '   Loss = ', Loss
              write (*,*) '   tstep = ', tstep
              write (*,*) '   c_output =', c_output
              write (*,*) '   ig = ', ig
              stop ' Stop at IMPLICIT , PHCHEM ' 
         endif

         ! Avoid underflow
         !
         if ( c_output.lt.1.d-30) c_output=1.d-30


         return
c END
         end



c***********************************************************************
      function ionsec_nplus (zenit, alt)                      

c       Calculates the N+ production by photoelectrons, following
c       Nicholson et al. 2009

c       FGG    sep 2010   first version
c***********************************************************************
                                                
      implicit none             

! Arguments 
      real*8  ionsec_nplus
      real    zenit
      real    alt

! Local variables 
      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
      real*8 b0,b1,b2,b3,b4
      real*8 altaux
      real*8 zenit_rad

!!!!!!! Program starts

      zenit_rad=dble(zenit*3.141592/180.)

      if(zenit.le.90.) then
         altaux=dble(alt)+
     $        15.*cos(zenit_rad)-40.*sqrt(cos(zenit_rad))+25.
      else
         altaux=dble(alt)
      endif

      if(altaux.gt.108.4) then
         a0  = 1.139925703d3
         a1  = -4.742398258d1
         a2  = 8.404232989d-1
         a3  = -8.108229906d-3
         a4  = 4.420892285d-5
         a5  = -1.154901432d-7
         a6  = -3.514073816d-11
         a7  = 8.790819159d-13
         a8  = -1.320788149d-16
         a9  = -8.202233732d-18
         a10 = -1.256480521d-22
         a11 = 1.329336168e-22
         a12 = -4.403185142d-25
         a13 = 6.098474897d-28
         a14 = -3.256951018d-31
         ionsec_nplus = a0 + a1*altaux + a2*altaux**2 + a3*altaux**3 +
     $        a4*altaux**4 + a5*altaux**5 + a6*altaux**6 + a7*altaux**7
     $        + a8*altaux**8 + a9*altaux**9 + a10*altaux**10 +
     $        a11*altaux**11 + a12*altaux**12 + a13*altaux**13 +
     $        a14*altaux**14
         ionsec_nplus = 10**(ionsec_nplus-2.)
      elseif(altaux.gt.80..and.altaux.le.108.4) then 
         b0 = 6.346190854d4
         b1 = -2.623253212d3
         b2 = 4.050319629d1
         b3 = -2.767987276d-1
         b4 = 7.064439029d-4
         ionsec_nplus = b0 + b1*altaux + b2*altaux**2 + b3*altaux**3 +
     $        b4*altaux**4
      else
         ionsec_nplus=0.d0
      endif
      if(ionsec_nplus.gt.100.d0.or.ionsec_nplus.lt.0.d0)
     $     ionsec_nplus=0.d0
      
      
      return                                         

      end


c***********************************************************************
      function ionsec_n2plus (zenit, alt)                      

c     N2+ production by photoelectrons, following Nicholson et al. 2009

c     FGG    sep 2010   first version
c***********************************************************************
                                                
      implicit none             

! Arguments 
      real*8  ionsec_n2plus
      real    zenit
      real    alt

! Local variables 

      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
      real*8 b0,b1,b2,b3,b4,b5
      real*8 altaux
      real*8 zenit_rad

!!!!!!! Program starts

      zenit_rad=dble(zenit*3.141592/180.)
      if(zenit.le.90.) then
         altaux=dble(alt)+
     $        15.*cos(zenit_rad)-40.*sqrt(cos(zenit_rad))+25.
      else
         altaux=dble(alt)
      endif


      if(altaux.gt.108.4) then
         a0 = 9.843804026d2
         a1 = -3.978696855d1
         a2 = 7.028448262d-1
         a3 = -7.11195117d-3
         a4 = 4.545683986d-5
         a5 = -1.905046447d-7
         a6 = 5.240068127d-10
         a7 = -9.130399894d-13
         a8 = 9.151792207d-16
         a9 = -4.023230948d-19
         ionsec_n2plus = a0 + a1*altaux + a2*altaux**2 + a3*altaux**3 +
     $        a4*altaux**4 + a5*altaux**5 + a6*altaux**6 +
     $        a7*altaux**7 + a8*altaux**8 + a9*altaux**9
         ionsec_n2plus = 10**(ionsec_n2plus-2.)
      elseif(altaux.gt.80..and.altaux.le.108.4) then
         b0 = 5.146111566d4
         b1 = -1.771736158d3
         b2 = 1.811156914d1
         b3 = 3.583722498d-3
         b4 = -9.891151731d-4
         b5 = 3.994474097d-6
         ionsec_n2plus = b0 + b1*altaux + b2*altaux**2 + b3*altaux**3 +
     $        b4*altaux**4 + b5*altaux**5
      else
         ionsec_n2plus = 0.d0
      endif
      if(ionsec_n2plus.gt.100.d0.or.ionsec_n2plus.lt.0.d0)
     $     ionsec_n2plus=0.d0
      
      return                                         

      end  



c***********************************************************************
      function ionsec_oplus (zenit, alt)                      

c     O+ production by photoelectrons, from Nicholson et al. 2009

c     FGG    sep 2010   first version
c***********************************************************************
                                                
      implicit none             

! Arguments 
      real*8  ionsec_oplus
      real    zenit
      real    alt

! Local variables 
      
      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
      real*8 b0,b1,b2,b3,b4,b5,b6
      real*8 altaux
      real*8 zenit_rad

!!!!!!! Program starts

      zenit_rad=dble(zenit*3.141592/180.)

      if(zenit.le.90.) then
         altaux=dble(alt)+
     $        15.*cos(zenit_rad)-40.*sqrt(cos(zenit_rad))+25.
      else
         altaux=dble(alt)
      endif

      if(altaux.gt.112.9) then
         a0  = 6.453740731d2 
         a1  = -2.547798991d1
         a2  = 4.384613636d-1
         a3  = -4.288387072d-3
         a4  = 2.596437447d-5
         a5  = -9.750300694d-8
         a6  = 1.986722344d-10
         a7  = -2.293667699d-14
         a8  = -1.080547019d-15
         a9  = 3.171787989d-18
         a10 = -4.566493384d-21
         a11 = 3.472393897d-24
         a12 = -1.115699979d-27
         ionsec_oplus = a0 + a1*altaux + a2*altaux**2 + a3*altaux**3 +
     $        a4*altaux**4 + a5*altaux**5 + a6*altaux**6 + a7*altaux**7
     $        + a8*altaux**8 + a9*altaux**9 + a10*altaux**10 +
     $        a11*altaux**11 +a12*altaux**12
         ionsec_oplus = 10**(ionsec_oplus-2.)
      elseif(altaux.gt.80..and.altaux.le.112.9) then 
         b0 = -5.934881676d5
         b1 = 3.546095705d4
         b2 = -8.806801303d2
         b3 = 1.163735173d1
         b4 = -8.62992438d-2
         b5 = 3.40547333d-4
         b6 = -5.587037506d-7
         ionsec_oplus = b0 + b1*altaux + b2*altaux**2 + b3*altaux**3 +
     $        b4*altaux**4 + b5*altaux**5 + b6*altaux**6
      else
         ionsec_oplus=0.d0   
      endif

      if(ionsec_oplus.gt.100.d0.or.ionsec_oplus.lt.0.d0)
     $     ionsec_oplus=0.d0
	
      return                                         

      end



c***********************************************************************
      function ionsec_coplus (zenit, alt)                      

c     CO+ production by photoelectrons from Nicholson et al. 2009

c     FGG    sep 2010   first version
c***********************************************************************
                                                
      implicit none             
      
! Arguments 
      real*8  ionsec_coplus
      real    zenit
      real    alt

! Local variables 

      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
      real*8 b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12
      real*8 altaux
      real*8 zenit_rad

!!!!!!! Program starts

      zenit_rad=dble(zenit*3.141592/180.)

      if(zenit.le.90.) then
         altaux=dble(alt)+
     $        15.*cos(zenit_rad)-40.*sqrt(cos(zenit_rad))+25.
      else
         altaux=dble(alt)
      endif

      if(altaux.gt.110.6) then
         a0  = 7.33258229d2
         a1  = -2.919984139d1
         a2  = 5.079651482d-1
         a3  = -5.057170037d-3
         a4  = 3.178156709d-5
         a5  = -1.309076957d-7
         a6  = 3.53858799d-10
         a7  = -6.060315762d-13
         a8  = 5.973573923d-16
         a9  = -2.584454326d-19
         ionsec_coplus = a0 + a1*altaux + a2*altaux**2 + a3*altaux**3 +
     $        a4*altaux**4 + a5*altaux**5 + a6*altaux**6 + a7*altaux**7
     $        + a8*altaux**8 + a9*altaux**9 
         ionsec_coplus = 10**(ionsec_coplus-2.)
      elseif(altaux.gt.80..and.altaux.le.110.6) then
         b0  = -1.165107657d6
         b1  = 4.315606169d4
         b2  = -3.480483017d2
         b3  = -3.831253024d0
         b4  = 4.33316742d-2
         b5  = 2.43075291d-4
         b6  = -7.835732322d-8
         b7  = -3.80859833d-8
         b8  = -1.947628467d-10
         b9  = 3.440753726d-12
         b10 = 2.336227916d-14
         b11 = -3.575877198d-16
         b12 = 1.030801684d-18
         ionsec_coplus = b0 + b1*altaux + b2*altaux**2 + b3*altaux**3 +
     $        b4*altaux**4 + b5*altaux**5 + b6*altaux**6 +
     $        b7*altaux**7 + b8*altaux**8 + b9*altaux**9 +
     $        b10*altaux**10 + b11*altaux**11 + b12*altaux**12
      else
         ionsec_coplus=0.d0
      endif
      if(ionsec_coplus.gt.100..or.ionsec_coplus.lt.0.d0)
     $     ionsec_coplus=0.d0
      
      return                                         

      end  
      


c***********************************************************************
      function ionsec_co2plus (zenit, alt)                      
      
c     CO2+ production by photoelectrons, from Nicholson et al. 2009
      
c     FGG    sep 2010   first version
c***********************************************************************
                                                
      implicit none             
      
!     Arguments 
      real*8  ionsec_co2plus
      real    zenit
      real    alt
      
!     Local variables 
      
      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
      real*8 b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
      real*8 altaux
      real*8 zenit_rad

!!!!!!!Program starts
      
      zenit_rad=dble(zenit*3.141592/180.)
      
      if(zenit.le.90.) then
         altaux=dble(alt)+
     $        15.*cos(zenit_rad)-40.*sqrt(cos(zenit_rad))+25.
      else
         altaux=dble(alt)
       endif

      if(altaux.gt.112.9) then
         a0  = 8.64096192d2
         a1  = -3.471713960d1
         a2  = 6.072614479d-1
         a3  = -6.050002721d-3
         a4  = 3.779639483d-5
         a5  = -1.533626303d-7
         a6  = 4.032987841d-10
         a7  = -6.602964344d-13
         a8  = 6.067681784d-16
         a9  = -2.356829271d-19
         ionsec_co2plus = a0 + a1*altaux + a2*altaux**2 + a3*altaux**3 +
     $        a4*altaux**4 + a5*altaux**5 + a6*altaux**6 +
     $        a7*altaux**7 + a8*altaux**8 + a9*altaux**9 
         ionsec_co2plus = 10**(ionsec_co2plus-2.)
      elseif(altaux.ge.80..and.altaux.le.112.9) then
         b0  = 1.159404818d6
         b1  = -5.617523193d4
         b2  = 8.535640078d2
         b3  = -5.956128524d-1
         b4  = -8.532255532d-2
         b5  = 1.217829692d-4
         b6  = 9.641535217d-6
         b7  = -4.502450788d-8
         b8  = -4.9684920146d-10
         b9  = 4.830572346d-12
         b10 = -1.168635127d-14
         ionsec_co2plus = b0 + b1*altaux + b2*altaux**2 + b3*altaux**3 +
     $        b4*altaux**4 + b5*altaux**5 + b6*altaux**6 +
     $        b7*altaux**7 + b8*altaux**8 + b9*altaux**9 +
     $        b10*altaux**10
      else
         ionsec_co2plus = 0.d0
      endif
      if(ionsec_co2plus.gt.100.d0.or.ionsec_co2plus.lt.0.d0)
     $     ionsec_co2plus=0.d0
      
      return                                         

      end  
      
      
c***********************************************************************
      function ionsec_o2plus (zenit, alt)                      

c     O2+ production by photoelectrons, from Nicholson et al. 2009

c     FGG    sep 2010   first version
c***********************************************************************
                                                
      implicit none             

! Arguments 
      real*8  ionsec_o2plus
      real    zenit
      real    alt
      
!     Local variables 

      real*8 a0,a1,a2,a3,a4,a5,a6,a7
      real*8 b0,b1,b2,b3,b4,b5,b6,b7,b8
      real*8 altaux
      real*8 zenit_rad
      
!!!!!!! Program starts

      zenit_rad=dble(zenit*3.141592/180.)
      
      if(zenit.le.90.) then
         altaux=dble(alt)+
     $        15.*cos(zenit_rad)-40.*sqrt(cos(zenit_rad))+25.
      else
         altaux=dble(alt)
      endif

      if(altaux.gt.112.9) then
         a0  = 7.265142765d2
         a1  = -2.714716832d1
         a2  = 4.315022800d-1
         a3  = -3.774025573d-3
         a4  = 1.962771814d-5
         a5  = -6.076128732d-8
         a6  = 1.037835637d-10
         a7  = -7.552930040d-14
         ionsec_o2plus = a0 + a1*altaux + a2*altaux**2 + a3*altaux**3 +
     $        a4*altaux**4 + a5*altaux**5 + a6*altaux**6 + a7*altaux**7
         ionsec_o2plus = 10**(ionsec_o2plus-2.)
      elseif(altaux.gt.80..and.altaux.le.112.9) then
         b0  = 3.622091694d6
         b1  = -1.45914419d5
         b2  = 1.764604914d3
         b3  = 1.186141031d0
         b4  = -1.331821089d-1
         b5  = -3.661686584d-4
         b6  = 1.934372959d-5
         b7  = -1.307294421d-7
         b8  = 2.846288872d-10
         ionsec_o2plus = b0 + b1*altaux + b2*altaux**2 + b3*altaux**3 +
     $        b4*altaux**4 + b5*altaux**5 + b6*altaux**6 + b7*altaux**7
     $        + b8*altaux**8 
      else
         ionsec_o2plus = 0.d0
      endif
      if(ionsec_o2plus.gt.100.d0.or.ionsec_o2plus.lt.0.d0)
     $     ionsec_o2plus=0.d0
      

      return                                         

      end  




c**********************************************************************
c**********************************************************************

      subroutine phdisrate(ig,nlayer,chemthermod,zenit,i)

c     Calculates photoionization and photodissociation rates from the
c     photoabsorption rates calculated in jthermcalc_e107 and the 
c     ionization/dissociation branching ratios in param_read_e107

c     apr 2002       fgg           first version
c**********************************************************************

      use param_v4_h, only: ninter,nabs,
     .    jfotsout,fluxtop,
     .    jion,jdistot,jdistot_b,
     .    efdisco2,efdiso2,efdish2o,
     .    efdish2o2,efdish2,efdiso3,
     .    efdiso,efdisn,efdish,
     .    efdisno,efdisn2,efdisno2,
     .    efdisco,efionco2,efiono2,efionn2,
     .    efionco,efiono3p,efionn,
     .    efionno,efionh

      implicit none

c     arguments

      integer         i                !altitude
      integer         ig,chemthermod,nlayer
      real            zenit

c     local variables

      integer         inter,iz,j
      real            lambda
      real            jdis(nabs,ninter,nlayer)
      character*1     dn


c**********************************************************************

c     photodissociation and photoionization rates initialization

      jion(:,i,:)    = 0.d0
      jdistot(:,i)   = 0.d0
      jdistot_b(:,i) = 0.d0

!      jion(1,i,1) = 0.d0          ! CO2 channel 1  ( --> CO2+ + e- )
!      jion(1,i,2) = 0.d0          ! CO2 channel 2  (  --> O+ + CO + e- )
!      jion(1,i,3) = 0.d0          ! CO2 channel 3  (  --> CO+ + O + e- )
!      jion(1,i,4) = 0.d0          ! CO2 channel 4  (  --> C+ + O2 + e- )
!      jion(2,i,1) = 0.d0          ! O2 (only one ionization channel)
!      jion(3,i,1) = 0.d0          ! O3P (only one ionization channel)
!      jion(4,i,1) = 0.d0          ! H2O (no ionization)
!      jion(5,i,1) = 0.d0          ! H2 (no ionization)
!      jion(6,i,1) = 0.d0          ! H2O2 (no ionization)
!      jion(7,i,1) = 0.d0          ! O3 (no ionization)
!      jion(8,i,1) = 0.d0          ! N2 channel 1 ( --> n2+ + e- )
!      jion(8,i,2) = 0.d0          ! N2 channel 2 ( --> n+ + n + e- )
!      jion(9,i,1) = 0.d0          ! N ( --> N+ + e- )
!      jion(10,i,1)= 0.d0          ! We do not mind its ionization
!      jion(11,i,1) = 0.d0          ! CO channel 1 ( --> co+ + e- )
!      jion(11,i,2) = 0.d0          ! CO channel 2 ( --> c+ + o + e- )
!      jion(11,i,3) = 0.d0          ! CO channel 3 ( --> o+ + c + e- )
!      jion(12,i,1) = 0.d0         ! H ( --> H+ + e- )
!      jion(13,i,1) = 0.d0
!      do j=1,nabs
!         jdistot(j,i) = 0.
!         jdistot_b(j,i) = 0.
!      end do

      if(zenit.gt.140.) then
         dn='n'
         else
         dn='d'
      end if

      if(dn.eq.'n') return

c     photodissociation and photoionization rates for each species

      do inter=1,ninter
c     CO2
         jdis(1,inter,i) = jfotsout(inter,1,i) * fluxtop(inter) 
     $        * efdisco2(inter)
         if(inter.gt.29.and.inter.le.32) then
            jdistot(1,i) = jdistot(1,i) + jdis(1,inter,i)
         else if(inter.le.29) then
            jdistot_b(1,i) = jdistot_b(1,i) + jdis(1,inter,i)
         end if
         jion(1,i,1)=jion(1,i,1) + 
     $        jfotsout(inter,1,i)*fluxtop(inter)*efionco2(inter,1)
         jion(1,i,2)=jion(1,i,2) + 
     $        jfotsout(inter,1,i)*fluxtop(inter)*efionco2(inter,2)
         jion(1,i,3)=jion(1,i,3) + 
     $        jfotsout(inter,1,i)*fluxtop(inter)*efionco2(inter,3)
         jion(1,i,4)=jion(1,i,4) + 
     $        jfotsout(inter,1,i)*fluxtop(inter)*efionco2(inter,4)


c     O2
         jdis(2,inter,i) = jfotsout(inter,2,i) * fluxtop(inter)
     $        * efdiso2(inter)
         if(inter.ge.31) then
            jdistot(2,i) = jdistot(2,i) + jdis(2,inter,i)
         else if(inter.eq.30) then
            jdistot(2,i)=jdistot(2,i)+0.02*jdis(2,inter,i)
            jdistot_b(2,i)=jdistot_b(2,i)+0.98*jdis(2,inter,i)
         else if(inter.lt.31) then
            jdistot_b(2,i) = jdistot_b(2,i) + jdis(2,inter,i)
         end if
         jion(2,i,1)=jion(2,i,1) +
     $        jfotsout(inter,2,i) * fluxtop(inter) * efiono2(inter,1)
         jion(2,1,2)=jion(2,1,2) + 
     $        jfotsout(inter,2,i) * fluxtop(inter) * efiono2(inter,2)
!(1.-efdiso2(inter))
         
c     O3P
         jion(3,i,1)=jion(3,i,1) +
     $        jfotsout(inter,3,i) * fluxtop(inter) * efiono3p(inter)
         
c     H2O
         jdis(4,inter,i) = jfotsout(inter,4,i) * fluxtop(inter)
     $        * efdish2o(inter)
         jdistot(4,i) = jdistot(4,i) + jdis(4,inter,i)
         
c     H2
         jdis(5,inter,i) = jfotsout(inter,5,i) * fluxtop(inter)
     $        * efdish2(inter)
         jdistot(5,i) = jdistot(5,i) + jdis(5,inter,i)
         
c     H2O2
         jdis(6,inter,i) = jfotsout(inter,6,i) * fluxtop(inter)
     $        * efdish2o2(inter) 
         jdistot(6,i) = jdistot(6,i) + jdis(6,inter,i)
         
         !Only if O3, N or ion chemistry requested
         if(chemthermod.ge.1) then
c     O3
            jdis(7,inter,i) = jfotsout(inter,7,i) * fluxtop(inter)
     $           * efdiso3(inter)
            if(inter.eq.34) then
               jdistot(7,i) = jdistot(7,i) + jdis(7,inter,i)
            else if (inter.eq.35) then
               jdistot(7,i) = jdistot(7,i) + 0.997 * jdis(7,inter,i)
               jdistot_b(7,i) = jdistot_b(7,i) + 0.003 * jdis(7,inter,i)
            else if (inter.eq.36) then
               jdistot_b(7,i) = jdistot_b(7,i) + jdis(7,inter,i)
            endif
         endif   !Of chemthermod.ge.1

         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
c     N2
            jdis(8,inter,i) = jfotsout(inter,8,i) * fluxtop(inter)
     $           * efdisn2(inter)
            jdistot(8,i) = jdistot(8,i) + jdis(8,inter,i)
            jion(8,i,1) = jion(8,i,1) + jfotsout(inter,8,i) * 
     $           fluxtop(inter) * efionn2(inter,1)
            jion(8,i,2) = jion(8,i,2) + jfotsout(inter,8,i) *
     $           fluxtop(inter) * efionn2(inter,2)
            
c     N
            jion(9,i,1) = jion(9,i,1) + jfotsout(inter,9,i) *
     $           fluxtop(inter) * efionn(inter)
            
c     NO
            jdis(10,inter,i) = jfotsout(inter,10,i) * fluxtop(inter)
     $           * efdisno(inter)
            jdistot(10,i) = jdistot(10,i) + jdis(10,inter,i)
            jion(10,i,1) = jion(10,i,1) + jfotsout(inter,10,i) *
     $           fluxtop(inter) * efionno(inter)
            
c     NO2
            jdis(13,inter,i) = jfotsout(inter,13,i) * fluxtop(inter)
     $           * efdisno2(inter)
            jdistot(13,i) = jdistot(13,i) + jdis(13,inter,i)

         endif   !Of chemthermod.ge.2
            
         !Only if ion chemistry requested
         if(chemthermod.eq.3) then
c     CO
            jdis(11,inter,i) = jfotsout(inter,11,i) * fluxtop(inter)
     $           * efdisco(inter)
            jdistot(11,i) = jdistot(11,i) + jdis(11,inter,i)
            jion(11,i,1) = jion(11,i,1) + jfotsout(inter,11,i) *
     $           fluxtop(inter) * efionco(inter,1)
            jion(11,i,2) = jion(11,i,2) + jfotsout(inter,11,i) *
     $           fluxtop(inter) * efionco(inter,2)
            jion(11,i,3) = jion(11,i,3) + jfotsout(inter,11,i) *
     $           fluxtop(inter) * efionco(inter,3) 

c     H
            jion(12,i,1) = jion(12,i,1) + jfotsout(inter,12,i) * 
     $           fluxtop(inter) * efionh(inter)
         endif    !Of chemthermod.eq.3


      end do 


      return
      

      end




c**********************************************************************
c***************************************************************************
	subroutine getch (ig,chemthermod,tt,zkm)


c     Reaction rates. The parameters rcoef are read from the 
c     chemthermos_reactionrates.def file
 

c***************************************************************************

        use param_v4_h, only: rcoef,
     .  ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,ch13,ch14,ch15,ch18,
     .  ch19,ch20,ch21,ch22,ch23,ch24,ch30,ch31,ch32,ch33,ch34,
     .  ch35,ch36,ch37,ch38,ch39,ch40,ch41,ch42,ch43,ch45,
     .  ch46,ch47,ch48,ch49,ch50,ch55,ch56,ch57,ch58,ch59,ch62,
     .  ch63,ch64,ch65,ch66,ch67,ch68,ch69,ch70,ch71,
     .  ch72,ch73,ch74,ch75,ch76,ch85,ch86,ch87



	implicit none

c Arguments 	
	integer       ig,chemthermod
	real*8 	      tt         ! Temperature 
	real          zkm	!  Altitude in km

c local variables: 
	real*8        tcte   
	real*8        t_elect ! electronic temperatures
	real*8        val ! valores de alturas corresp a t_elect
	real*8        zhanson(9),tehanson(9)
	real*8        incremento
	integer       ii, i1, i2
	
c**************************************************************************

	tcte = tt
	
!        goto 151
	!Electronic temperatures
        ! (Hanson et al. 1977) approx. from Mars p. 107 	   
	zhanson(1) = 120. 
	zhanson(2) = 130. 
	zhanson(3) = 150. 
	zhanson(4) = 175. 
	zhanson(5) = 200. 
	zhanson(6) = 225. 
	zhanson(7) = 250. 
	zhanson(8) = 275. 
	zhanson(9) = 300.
	tehanson(1) = tt
	tehanson(2) = 200.
	tehanson(3) = 300.
	tehanson(4) = 500.
	tehanson(5) = 1250.
	tehanson(6) = 2000.
	tehanson(7) = 2200.
	tehanson(8) = 2400.
	tehanson(9) = 2500.
	if ( zkm .le. 120. ) then
	   t_elect = tt 
	else if(zkm .ge.300.) then
	   t_elect=tehanson(9)
	else
	   do ii=9,2,-1 
	      if ( zkm .lt. zhanson(ii) ) then 
		 i1 = ii - 1
		 i2 = ii 
	      endif
	   enddo
	   incremento = ( tehanson(i2)-tehanson(i1) ) / 
     $          ( zhanson(i2)-zhanson(i1) )
	   t_elect = tehanson(i1) + (zkm-zhanson(i1)) * incremento

!           t_elect = t_elect * 2.
	endif
! 151    continue
        !MAVEN measured electronic temperature (Ergun et al., GRL 2015)
!        t_elect=((3140.+120.)/2.)+((3140.-120.)/2.)*tanh((zkm-241.)/60.)

	!Initializations
	ch2  = 0.d0
	ch3  = 0.0
	ch4  = 0.d0
	ch5  = 0.d0
	ch7  = 0.d0
	ch9  = 0.d0
	ch10 = 0.d0
	ch11 = 0.d0
	ch13 = 0.d0
	ch14 = 0.d0
	ch15 = 0.d0
	ch18 = 0.d0
	ch19 = 0.d0
	ch20 = 0.d0
	ch21 = 0.d0
	ch22 = 0.d0
	ch23 = 0.d0
	ch24 = 0.d0
	ch30 = 0.d0
	ch31 = 0.d0
	ch32 = 0.d0
	ch33 = 0.d0
	ch34 = 0.d0
	ch35 = 0.d0
	ch36 = 0.d0
	ch37 = 0.d0
	ch38 = 0.d0
	ch39 = 0.d0
	ch40 = 0.d0
	ch41 = 0.d0
	ch42 = 0.d0
	ch43 = 0.d0
	ch45 = 0.d0
	ch46 = 0.d0
	ch47 = 0.d0
	ch48 = 0.d0
	ch49 = 0.d0
	ch50 = 0.d0
	ch55 = 0.d0
	ch56 = 0.d0
	ch57 = 0.d0
	ch58 = 0.d0
	ch59 = 0.d0
	ch62 = 0.d0
	ch63 = 0.d0
	ch64 = 0.d0
	ch65 = 0.d0
	ch66 = 0.d0
	ch67 = 0.d0
	ch68 = 0.d0
	ch69 = 0.d0
	ch70 = 0.d0
	ch71 = 0.d0
	ch72 = 0.d0
	ch73 = 0.d0
	ch74 = 0.d0
	ch75 = 0.d0
	ch76 = 0.d0
	ch85 = 0.d0
	ch86 = 0.d0


	!Reaction rates
!ch2: h + o2 + co2 --> ho2 + co2
        ! JPL 2003 (low pressure limit)*2.5
!	ch2 = 1.425d-31 * (tcte / 300.)**(-1.6d0)
        ! JPL 2011 (low pressure limit)*2.5
!        ch2 = 1.1e-31 * (tcte / 300.)**(-1.3)
        ch2=rcoef(1,1)*((tcte/300.)**rcoef(1,2))*exp(rcoef(1,3)/tcte)

!ch3: o + ho2 --> oh + o2
        ! JPL 2011: 
!	ch3 = 3.0d-11 * exp(200.d0 / tcte)
        ch3=rcoef(2,1)*((tcte/300.)**rcoef(2,2))*exp(rcoef(2,3)/tcte)

ch4: co + oh --> co2 + h
        !Nair et al, 1994:
	!ch4 = 3.2d-13 * exp(-300.d0 / tcte)
        !mccabe et al., grl, 28, 3135, 2001
        !ch4 = 1.57d-13 + 3.54d-33*concco2 
        !JPL 2011 (low pressure limit):
!        ch4 = 1.5d-13 * (tcte/300.)**0.6
        ch4=rcoef(3,1)*((tcte/300.)**rcoef(3,2))*exp(rcoef(3,3)/tcte)

ch5: ho2 + ho2 --> h2o2 + o2
        !JPL 2003:
	!ch5 = 2.3d-13 * exp(600.d0 / tcte)
        !JPL 2011:
!        ch5 = 3.0d-13 * exp(460.d0 / tcte)
        ch5=rcoef(4,1)*((tcte/300.)**rcoef(4,2))*exp(rcoef(4,3)/tcte)

ch7: oh + ho2 --> h2o + o2
        !JPL 2011:
!	ch7 = 4.8d-11 * exp(250.d0 / tcte)
        ch7=rcoef(5,1)*((tcte/300.)**rcoef(5,2))*exp(rcoef(5,3)/tcte)

ch9: o(1d) + h2o --> 2oh
        !JPL 2003:
	!ch9 = 2.2d-10
        !JPL 2011:
!        ch9 = 1.63d-10 * exp(60.d0 / tcte)
        ch9=rcoef(6,1)*((tcte/300.)**rcoef(6,2))*exp(rcoef(6,3)/tcte)

ch10: o + o + co2 --> o2 + co2
        !JPL 1990:
!	ch10 = 1.1d-27 * (tcte **(-2.0d0)) !Estandard en el 1-D
        !Tsang and Hampson, 1986:
!	ch10 = 1.3d-34 * exp(900.d0/tcte)
        ch10=rcoef(7,1)*((tcte/300.)**rcoef(7,2))*exp(rcoef(7,3)/tcte)

ch11: o + oh --> o2 + h
        !JPL 2003:
	!ch11 = 2.2d-11 * exp(120.d0 / tcte)
        !JPL 2011:
!        ch11 = 1.8d-11 * exp(180.d0 /tcte)
        ch11=rcoef(8,1)*((tcte/300.)**rcoef(8,2))*exp(rcoef(8,3)/tcte)

ch13: h + ho2 --> h2 + o2
        !JPL 2003:
	!ch13 = 6.5d-12
        !JPL 2011:
!        ch13 = 6.9d-12
        ch13=rcoef(9,1)*((tcte/300.)**rcoef(9,2))*exp(rcoef(9,3)/tcte)

ch14: o(1d) + h2 --> h + oh
        !JPL 2003:
	!ch14 = 1.1d-10
        !JPL 2011:
!        ch14 = 1.2d-10
        ch14=rcoef(10,1)*((tcte/300.)**rcoef(10,2))*
     $       exp(rcoef(10,3)/tcte)

ch15: oh + h2 --> h + h2o
        !JPL 2003:
	!ch15 = 5.5d-12 * exp (-2000.d0 / tcte)
        !JPL 2011:
!        ch15 = 2.8d-12 * exp (-1800.d0 / tcte)
        ch15=rcoef(11,1)*((tcte/300.)**rcoef(11,2))*
     $       exp(rcoef(11,3)/tcte)

ch18: oh + h2o2 --> h2o + ho2
        !JPL 2003:
	!ch18 = 2.9d-12 * exp (-160.d0 / tcte)
        !JPL 2011:
!        ch18 = 1.8d-12
        ch18=rcoef(12,1)*((tcte/300.)**rcoef(12,2))*
     $       exp(rcoef(12,3)/tcte)

ch19: o(1d) + co2 --> o + co2
        !JPL 2003:
	!ch19 = 7.4d-11 * exp(120.d0 / tcte)
        !JPL 2011:
!        ch19 = 7.5d-11 * exp(115.d0 / tcte)        
        ch19=rcoef(13,1)*((tcte/300.)**rcoef(13,2))*
     $       exp(rcoef(13,3)/tcte)

ch20: o(1d) + o2 --> o + o2
        !JPL 2003:
	!ch20 = 3.2d-11 * exp (70.d0 / tcte)
        !JPL 2011:
!        ch20 = 3.3d-11 * exp(55.d0 / tcte)
        ch20=rcoef(14,1)*((tcte/300.)**rcoef(14,2))*
     $       exp(rcoef(14,3)/tcte)

ch21: o + o2 + co2 --> o3 + co2
        !JPL 2011 * 2.5:
!	ch21 = 1.5d-33 * ((tcte / 300.d0) ** (-2.4d0))
        ch21=rcoef(15,1)*((tcte/300.)**rcoef(15,2))*
     $       exp(rcoef(15,3)/tcte)
	

	!Only if O3, N or ion chemistry requested
	if(chemthermod.ge.1) then

ch22: o3 + h --> o2 + oh
           !JPL 2011:
!	   ch22 = 1.4d-10 * exp (-470.d0 / tcte)
           ch22=rcoef(16,1)*((tcte/300.)**rcoef(16,2))*
     $       exp(rcoef(16,3)/tcte)

ch23: o3 + oh --> ho2 + o2
           !JPL 2011:
!	   ch23 = 1.7d-12 * exp (-940.d0 / tcte)
           ch23=rcoef(17,1)*((tcte/300.)**rcoef(17,2))*
     $       exp(rcoef(17,3)/tcte)

ch24: o3 + ho2 --> oh + 2o2
           !JPL 2011:
!	   ch24 = 1.0d-14 * exp (-490.d0 / tcte)
           ch24=rcoef(18,1)*((tcte/300.)**rcoef(18,2))*
     $       exp(rcoef(18,3)/tcte)

	endif

	!Only if N or ion chemistry requested
	if(chemthermod.ge.2) then
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c N chemistry
ch30: n + no --> n2 + o
           !JPL 2011:
!	   ch30 = 2.1d-11 * exp (100.d0 / tcte)
           ch30=rcoef(19,1)*((tcte/300.)**rcoef(19,2))*
     $       exp(rcoef(19,3)/tcte)

ch31: n2 + o(1d) --> n2 + o
           !JPL 2011:
!	   ch31 = 2.15d-11 * exp (110.d0 / tcte)
           ch31=rcoef(20,1)*((tcte/300.)**rcoef(20,2))*
     $       exp(rcoef(20,3)/tcte)

ch32: n + o2 --> no + o
           !JPL 2011:
!	   ch32 = 1.5d-11 * exp (-3600.d0 / tcte)
           ch32=rcoef(21,1)*((tcte/300.)**rcoef(21,2))*
     $       exp(rcoef(21,3)/tcte)

ch33: n + oh --> no + h
           !Atkinson et al., 1989 (usado en Nair et al., 1994)
!	   ch33 = 3.8d-11 * exp (85.d0 / tcte)
           ch33=rcoef(22,1)*((tcte/300.)**rcoef(22,2))*
     $       exp(rcoef(22,3)/tcte)

ch34: n + o3 --> no + o2
           !JPL 2011 (it is an upper limit):
!	   ch34 = 1.0d-16
           ch34=rcoef(23,1)*((tcte/300.)**rcoef(23,2))*
     $       exp(rcoef(23,3)/tcte)

ch35: n + ho2 --> no + oh
           !Brune et al., 1983 (from Nair et al., 1994)
!	   ch35 = 2.2d-11
           ch35=rcoef(24,1)*((tcte/300.)**rcoef(24,2))*
     $       exp(rcoef(24,3)/tcte)

ch36: n(2d) + o --> n + o
           !Fell et al., 1990 (from Nair et al., 1994)
	   !ch36 = 6.9d-13
           !Herron, 1999:
!           ch36 = 3.3d-12 * exp(-260.d0 / tcte)
           ch36=rcoef(25,1)*((tcte/300.)**rcoef(25,2))*
     $       exp(rcoef(25,3)/tcte)

ch37: n(2d) + n2 --> n + n2
           !Herron, 1999:
           !Coincides with Nair et al., 1994:
!	   ch37 = 1.7d-14
           ch37=rcoef(26,1)*((tcte/300.)**rcoef(26,2))*
     $       exp(rcoef(26,3)/tcte)

ch38: n(2d) + co2 --> no + co
           !Pipper et al., 1989 (from Nair et al., 1994):
	   !ch38 = 3.5d-13
           !Herron, 1999:
!           ch38 = 3.6d-13
           ch38=rcoef(27,1)*((tcte/300.)**rcoef(27,2))*
     $       exp(rcoef(27,3)/tcte)

ch39: no + ho2 --> no2+oh
           !JPL 2006:
	   !ch39 = 3.5d-12 * exp (250.d0 / tcte)
           !JPL 2011:
!           ch39 = 3.3d-12 * exp(270.d0 / tcte)
           ch39=rcoef(28,1)*((tcte/300.)**rcoef(28,2))*
     $       exp(rcoef(28,3)/tcte)

ch40: o + no + co2 --> no2 + co2
           !JPL 2011 * 2.5 (low pressure limit)
!	   ch40 = 2.5d0 * 9.0d-32 * ((tcte / 300.d0) ** (-1.5d0))
           ch40=rcoef(29,1)*((tcte/300.)**rcoef(29,2))*
     $       exp(rcoef(29,3)/tcte)

ch41: o + no2 --> no + o2
           !JPL 2011:
!	   ch41 = 5.1d-12 * exp (210.d0 / tcte)
           ch41=rcoef(30,1)*((tcte/300.)**rcoef(30,2))*
     $       exp(rcoef(30,3)/tcte)

ch42: no + o3 --> no2 + o2
           !JPL 2011:
!	   ch42 = 3.0d-12 * exp (-1500.d0 / tcte)
           ch42=rcoef(31,1)*((tcte/300.)**rcoef(31,2))*
     $       exp(rcoef(31,3)/tcte)

ch43: h + no2 --> no + oh
           !JPL 2011:
!	   ch43 = 4.0d-10 * exp (-340.d0 / tcte)
           ch43=rcoef(32,1)*((tcte/300.)**rcoef(32,2))*
     $       exp(rcoef(32,3)/tcte)

ch45: n + o --> no
           !Du and Dalgarno, 1990
!	   ch45 = 2.8d-17 * ((300.d0 / tcte) ** 0.5)
           ch45=rcoef(33,1)*((tcte/300.)**rcoef(33,2))*
     $       exp(rcoef(33,3)/tcte)
	   
	endif    !of if(chemthermod.ge.2)


	!Only if ion chemistry requested
	if(chemthermod.eq.3) then
c
c Ionosphere
c
ch46:  co2+ +  O2 --> O2+ + CO2
           !Moffat et al., 2005 (given by GG):
	   !ch46 = 5.0d-11
           !Copp et al., 1982:
           !ch46 = 5.3d-11
           !Aninich 1993 (from Fox and Sung 2001):
!           ch46 = 5.5d-11 * (300.d0/t_elect)**0.82
           ch46=rcoef(34,1)*((t_elect/300.)**rcoef(34,2))*
     $       exp(rcoef(34,3)/t_elect)
           

ch47: CO2+ + O -->  O+  +  CO2
           !Original (incorrect) value (corresponds to ch48):
	   !ch47 = 1.64d-10
           !Fehsenfeld et al., 1970 (from UMIST, 
           !Fox and Sung 2001, Krasnopolsky 2002):
!           ch47 = 9.6d-11
           ch47=rcoef(35,1)*((t_elect/300.)**rcoef(35,2))*
     $       exp(rcoef(35,3)/t_elect)

ch48: CO2+ + O --> O2+  + CO
           !Original (incorrect) value (corresponds to ch47):
	   !ch48 = 9.6d-11
           !Fehsenfeld et al., 1970 (from UMIST,
           !Fox and Sung 2001, Krasnopolsky 2002):
!           ch48 = 1.64d-10
           ch48=rcoef(36,1)*((t_elect/300.)**rcoef(36,2))*
     $       exp(rcoef(36,3)/t_elect)

ch49: O2+ + elect --> O + O
           !Alge et al., 1983:
           !Here we do not divide into reaction producing different
           !O atomic states. O + O(1d) seems to be the dominant products
           !(see Fox and Sung 2002). We should consider dividing 
           !into two different reactions
!	   ch49 = 2.0d-7*(300.d0/t_elect)**(0.7d0)
           ch49=rcoef(37,1)*((t_elect/300.)**rcoef(37,2))*
     $       exp(rcoef(37,3)/t_elect)

ch50: O+  + CO2 --> O2+  + CO
           !Adams et al., 1980 (from UMIST):
!	   ch50 = 9.4d-10
           !Anicich 1993 (from Fox and Sung 2001):
           !ch50 = 1.1d-9
           ch50=rcoef(38,1)*((t_elect/300.)**rcoef(38,2))*
     $       exp(rcoef(38,3)/t_elect)

ch55: CO2+ +  e ---->   CO  +	O
           !Mitchell, 1990 (from UMIST):
!	   ch55 = 3.8d-7*(300.d0/t_elect)**(0.5d0)
           !Gougousi et al., 1997 (from Fox and Sung 2001):
           !ch55 = 3.5d-7*(300.d0/t_elect)**0.5d0
           ch55=rcoef(39,1)*((t_elect/300.)**rcoef(39,2))*
     $       exp(rcoef(39,3)/t_elect)
        
ch56: O+ + CO2  --->  O2   +   CO+
           !Original, Kim et al., 1989:
	   !ch56 = 9.4d-10
           !It does not appear in any other paper. Its presence in 
           !Kim et al., 1989 is probably a confusion with ch50.
!           ch56 = 0.d0
           ch56=rcoef(40,1)*((t_elect/300.)**rcoef(40,2))*
     $       exp(rcoef(40,3)/t_elect)

ch57: CO+ + CO2  ---> CO2+ + CO
           !Adams et al., 1978 (from UMIST):
!	   ch57 = 1.0d-9
           !Anicich 1993 (from Fox and Sung 2001):
           !ch57 = 1.1d-9
           ch57=rcoef(41,1)*((t_elect/300.)**rcoef(41,2))*
     $       exp(rcoef(41,3)/t_elect)


ch58: CO+ + O   --->   O+  + CO
           !Fenhsenfeld et al. 1970 (from UMIST, F&S2001, K2002):
!	   ch58 = 1.4d-10
           ch58=rcoef(42,1)*((t_elect/300.)**rcoef(42,2))*
     $       exp(rcoef(42,3)/t_elect)
  
ch59: C+ +  CO2 --->  CO+  + CO    !!!!  NEW  !!!
           !Fahey et al., 1981 (from UMIST, F&S2001, K2002):
!	   ch59 = 1.1d-9
           ch59=rcoef(43,1)*((t_elect/300.)**rcoef(43,2))*
     $       exp(rcoef(43,3)/t_elect)

ch62:  CO2+  +  NO   -->  NO+  +  CO2
           !Copp et al., 1982 (from UMIST):
!	   ch62 = 1.2d-10
           !Anicich 1993 (from Fox and Sung 2001):
           !ch62 = 1.23d-10
           ch62=rcoef(44,1)*((t_elect/300.)**rcoef(44,2))*
     $       exp(rcoef(44,3)/t_elect)

ch63:  CO2+  +  N    -->  NO  + CO+
           !Kim et al., 1989:
	   !ch63 = 1.0d-11
           !Scott et al., 1998 (from Fox and Sung 2001):
!           ch63 = 3.4d-10
           ch63=rcoef(45,1)*((t_elect/300.)**rcoef(45,2))*
     $       exp(rcoef(45,3)/t_elect)

ch64:   O2+  +  NO    -->  NO+  + O2
           !Middey and Vigiano 1999 (from Fox and Sung 2001):
!	   ch64 = 4.5d-10
           !Aninich 1993 (from UMIST):
           !ch64 = 4.6d-10
           ch64=rcoef(46,1)*((t_elect/300.)**rcoef(46,2))*
     $       exp(rcoef(46,3)/t_elect)

ch65:  O2+  +  N2    -->  NO+  + NO
           !Original from GG, Moffat 2005:
	   !ch65 = 1.0d-16
           !Ferguson 1973 (from Fox and Sung 2001):
!           ch65 = 1.0d-15
           ch65=rcoef(47,1)*((t_elect/300.)**rcoef(47,2))*
     $       exp(rcoef(47,3)/t_elect)

ch66:  O2+  +  N    -->  NO+  + O
           !Kim et al., 1989:
	   !ch66 = 1.2d-10
           !Scott et al., 1998 (from Fox and Sung 2001):
!           ch66 = 1.0d-10
           !Goldan et al., 1966 (from UMIST):
           !ch66 = 1.8d-10
           ch66=rcoef(48,1)*((t_elect/300.)**rcoef(48,2))*
     $       exp(rcoef(48,3)/t_elect)


ch67:   O+  +  N2    -->  NO+  + N
           !Moffat 2005:
	   !ch67 = 1.2d-12 * (300.d0/t_elect)**(0.41d0)
           !Hierl et al. 1997 (from Fox and Sung 2001):
!           ch67 = 1.2d-12 * (300.d0/t_elect)**0.45d0
           !Adams et al., 1980 (from UMIST):
           !ch67=2.42d-12 * (300.d0/t_elec)**(-0.21)*exp(44./t_elec)
           ch67=rcoef(49,1)*((t_elect/300.)**rcoef(49,2))*
     $       exp(rcoef(49,3)/t_elect)

ch68:  N2+  +  CO2    -->  CO2+  + N2
           !Adams et al. 1980 (from UMIST):
	   !ch68 = 7.7d-10
           !Dotan et al. 2000 (from F&S2001):
!           ch68 = 9.0d-10 * (300./t_elect)**0.23
           ch68=rcoef(50,1)*((t_elect/300.)**rcoef(50,2))*
     $       exp(rcoef(50,3)/t_elect)

ch69:   N2+  +  O3p    -->  NO+  + N
           !McFarland et al., 1974 (from UMIST):
	   !ch69 = 1.3d-10
           !Scott et al. 1999 (from F&S2001):
!           ch69 = 1.33d-10 * (300./t_elect)**0.44
           ch69=rcoef(51,1)*((t_elect/300.)**rcoef(51,2))*
     $       exp(rcoef(51,3)/t_elect)

ch70:  N2+  +  CO    -->  N2  + CO+
           !Adams et al., 1980 (from UMIST):
!	   ch70 = 7.4d-11
           !Frost et al., 1998 (from F&S2001):
           !ch70 = 7.6d-11
           ch70=rcoef(52,1)*((t_elect/300.)**rcoef(52,2))*
     $       exp(rcoef(52,3)/t_elect)

ch71:  N2+  +  e-    -->  N  + N
           !Moffat 2005
	   !ch71 = 3.5d-7 * (300.d0/t_elect)**(0.5d0)
           !Peterson et al. 1998 (from UMIST):
!           ch71 = 1.7d-7 * (300.d0/t_elect)**0.3
           !Zipf 1980+Kella et al 1996 (from F&S2001):
           !ch71 = 2.2d-7 * (300.d0/t_elect)**0.39
           ch71=rcoef(53,1)*((t_elect/300.)**rcoef(53,2))*
     $       exp(rcoef(53,3)/t_elect)


ch72:  N2+  +  O3p    -->  O+  + N2
           !Moffat 2005:
	   !ch72 = 4.1d-10
           !McFarland et al. 1974 (From UMIST):
           !ch72 = 1.1d-11
           !Scott et al., 1999 (from F&S2001):
!           ch72 = 7.0d-12 * (300.d0/t_elect)**0.23
           ch72=rcoef(54,1)*((t_elect/300.)**rcoef(54,2))*
     $       exp(rcoef(54,3)/t_elect)
           

ch73  CO+  +  H    -->  H+  + CO
           !Scott et al., 1997 (from F&S2001):
!	   ch73 = 4.0d-10
           !Federer et al. 1984 (from UMIST):
           !ch73 = 7.5d-10
           ch73=rcoef(55,1)*((t_elect/300.)**rcoef(55,2))*
     $       exp(rcoef(55,3)/t_elect)


ch74:   O+  +  H    -->  H+  + O
           !Krasnopolsky 2002:
	   !ch74 = 5.7d-10 * (tt/300.d0)**(0.36d0)
           !Stancil et al. 1999 (from UMIST):
!           ch74 = 5.66d-10*(tt/300.)**0.36*exp(8.6/tt)
           !Aninich 1993 (from F&S2001):
           !ch74 = 6.4e-10
           ch74=rcoef(56,1)*((tcte/300.)**rcoef(56,2))*
     $       exp(rcoef(56,3)/tcte)

ch75:  NO+  +  e-  -->  N  +  O
           !Mitchel 1990 (from UMIST):
!	   ch75 = 4.3d-7 * (300.d0/t_elect)**(0.37d0)
           !Vejby-Christensen et al. 1996 (from F&S2001):
           !ch75=4.0d-7 * (300.d0/t_elect)**0.5d0
           ch75=rcoef(57,1)*((t_elect/300.)**rcoef(57,2))*
     $       exp(rcoef(57,3)/t_elect)

ch76:   H+  +  O3p  -->  O+  +  H
           !Krasnopolsky et al. 2002:
	   !ch76 = 7.3d-10 * (tt/300.d0)**(0.23d0) * exp(-226./tt)
           !Stancil et al. 1999 (from UMIST):
!           ch76 = 6.86e-10* (tt/300.)**0.26*exp(-224.3/tt)
           ch76=rcoef(58,1)*((tcte/300.)**rcoef(58,2))*
     $       exp(rcoef(58,3)/tcte)

ch85:  N+  +  CO2    -->  CO2+  + N
           !Krasnopolsky 2002:
	   !ch85 = 1.d-9
           !Adams et al. 1980 (from UMIST):
!           ch85 = 7.5d-10
           !Aninich et al. 1993 (from F&S2001):
           !ch85 = 9.2d-10
           ch85=rcoef(59,1)*((t_elect/300.)**rcoef(59,2))*
     $       exp(rcoef(59,3)/t_elect)

ch86:  H2 + CO2+  --> H + HCO2+
           !Scott et al. 1998 (from F&S2001 and K2002):
	   !ch86 = 8.7d-10
           !Copp et al. 1982 (from UMIST):
!           ch86 = 9.5d-10
           ch86=rcoef(60,1)*((t_elect/300.)**rcoef(60,2))*
     $       exp(rcoef(60,3)/t_elect)


c     h87:  HCO2+ + e -> CO2 + H
           !Krasnopolsky 2002:
!           ch87 = 3.4d-7 * (300.d0/t_elect)**(0.5d0)
           !UMIST 2012: the reactions has 3 different sets of products: CO2+H,
           !CO+O+H (dominante) y CO+OH. Habria que tener esto en cuenta
           ch87=rcoef(61,1)*((t_elect/300.)**rcoef(61,2))*
     $       exp(rcoef(61,3)/t_elect)

	   
	endif    !Of if(chemthermod.eq.3)

	return 
	end



c**********************************************************************
c**********************************************************************

      subroutine lifetimes
     &   ( ig,i,nlayer,chemthermod,zenit,zx,jdistot8, jdistot8_b, jion8,
     $     xtmin, xcompmin, xn_comp_en_EQ,
     $     co2xini,o2xini,o3pxini,coxini,hxini,ohxini,ho2xini,h2xini,
     $     h2oxini,h2o2xini,o1dxini,o3xini,n2xini,nxini,noxini,no2xini,
     $     n2dxini,co2plusxini,oplusxini,o2plusxini,coplusxini,
     $     cplusxini,nplusxini,noplusxini,n2plusxini,hplusxini,
     $     hco2plusxini,electxini)

  
c Calculates the lifetime of each species at each time step (itime) 
c and each altitude (i), and the minimum of them (tmin)
c It also computes the number of species in PE
c 
c
c     jul 2008      MA        Version en subrutina
c         2009      FGG       Adaptation to GCM
c**********************************************************************

      use iono_h
      use param_v4_h

      implicit none

      include 'callkeys.h'

c     arguments
c
      integer   i,ig,nlayer                               ! I. Layer  
      integer   chemthermod
      real      zenit
      real      zx(nlayer)
      real*8    jdistot8(nabs,nlayer)                  ! I.
      real*8    jdistot8_b(nabs,nlayer)                ! I.
      real*8    jion8(nabs,nlayer,4)                ! I.

      real*8    xtmin                         ! O. 
      integer   xcompmin                      ! O.
      integer   xn_comp_en_EQ                 ! O.

      real*8    co2xini,o2xini,o3pxini,coxini
      real*8    ho2xini,h2xini,hxini,ohxini
      real*8    h2o2xini,o1dxini,o3xini,h2oxini
      real*8    nxini,noxini,n2xini,n2dxini,no2xini
      real*8    co2plusxini,coplusxini,oplusxini,o2plusxini
      real*8    cplusxini,noplusxini,n2plusxini,hplusxini
      real*8    electxini,nplusxini,hco2plusxini

c     local variables
c
      integer   j 


      external  ionsec_nplus
      real*8    ionsec_nplus

      external  ionsec_n2plus
      real*8    ionsec_n2plus

      external  ionsec_oplus
      real*8    ionsec_oplus
     
      external  ionsec_coplus
      real*8    ionsec_coplus

      external  ionsec_co2plus
      real*8    ionsec_co2plus

      external  ionsec_o2plus
      real*8    ionsec_o2plus



ccccccccccccccc CODE STARTS 

         !Initialization of lifetimes
         do j = 1, nreact
            tauco2(j,i)  = 1.d30
            tauo2(j,i)   = 1.d30
            tauo3p(j,i)  = 1.d30
            tauco(j,i)   = 1.d30
            tauh(j,i)    = 1.d30
            tauoh(j,i)   = 1.d30
            tauho2(j,i)  = 1.d30
            tauh2(j,i)   = 1.d30
            tauh2o(j,i)  = 1.d30
            tauo1d(j,i)  = 1.d30
            tauh2o2(j,i) = 1.d30
            tauo3(j,i)   = 1.d30
            taun2(j,i)   = 1.d30
            taun(j,i)    = 1.d30
            tauno(j,i)   = 1.d30
            taun2d(j,i)  = 1.d30
            tauno2(j,i)  = 1.d30
            tauco2plus(j,i) = 1.d30
            tauoplus(j,i)   = 1.d30
	    tauo2plus(j,i)  = 1.d30
            taucoplus(j,i)  = 1.d30
	    taucplus(j,i)   = 1.d30
	    taunplus(j,i)   = 1.d30
	    taun2plus(j,i)  = 1.d30
            taunoplus(j,i)  = 1.d30
	    tauhplus(j,i)   = 1.d30
            tauhco2plus(j,i)= 1.d30
         end do

         !Lifetime of each species in each reaction
         if(jdistot8(1,i).gt.1.d-30) tauco2(1,i) = 1.d0 / jdistot8(1,i)
		 
         if(ch2*o2xini*co2xini.gt.1.d-30) 
     $      tauh(2,i) = 1.d0 / (ch2 * o2xini * co2xini)  
         if(ch2*hxini*co2xini.gt.1.d-30)
     $      tauo2(2,i) = 1.d0 / (ch2 * hxini * co2xini)
 
         if(ch3*o3pxini.gt.1.d-30) tauho2(3,i) = 1.d0 / 
     $        (ch3 * o3pxini)
         if(ch3*ho2xini.gt.1.d-30) tauo3p(3,i) = 1.d0 / 
     $        (ch3 * ho2xini)
 
         if(ch4*coxini.gt.1.d-30) tauoh(4,i) = 1.d0 / 
     $        (ch4 * coxini)
         if(ch4*ohxini.gt.1.d-30) tauco(4,i) = 1.d0 / 
     $        (ch4 * ohxini)
 
         if(ch5*ho2xini.gt.1.d-30)tauho2(5,i)=1.d0 / 
     $        (2.d0*ch5*ho2xini)

		 
	 if(jdistot8(6,i).gt.1.d-30) tauh2o2(6,i) = 1.d0 / jdistot8(6,i)

         if(ch7*ohxini.gt.1.d-30) tauho2(7,i) = 1.d0 / 
     $        (ch7 * ohxini)
         if(ch7*ho2xini.gt.1.d-30) tauoh(7,i) = 1.d0 / 
     $        (ch7 * ho2xini)
 
         if(jdistot8(4,i).gt.1.d-30) tauh2o(8,i) = 1.d0 / jdistot8(4,i)

	 if(ch9*o1dxini.gt.1.d-30) tauh2o(9,i) = 1.d0 / 
     $        (ch9 * o1dxini)
         if(ch9*h2oxini.gt.1.d-30) tauo1d(9,i) = 1.d0 / 
     $        (ch9 * h2oxini)

         if(ch10*o3pxini*co2xini.gt.1.d-30) 
     $     tauo3p(10,i) = 1.d0 / 
     $        (2.d0 * ch10 * o3pxini * co2xini)

         if(ch11*o3pxini.gt.1.d-30) tauoh(11,i)=1.d0 / 
     $        (ch11 * o3pxini)
         if(ch11*ohxini.gt.1.d-30) tauo3p(11,i) = 1.d0 / 
     $        (ch11 * ohxini)

         if(jdistot8(2,i).gt.1.d-30) tauo2(12,i) = 1.d0 / jdistot8(2,i)

         if(ch13*hxini.gt.1.d-30) tauho2(13,i) = 1.d0 / 
     $        (ch13 * hxini)
         if(ch13*ho2xini.gt.1.d-30) tauh(13,i) = 1.d0 / 
     $        (ch13 * ho2xini)

         if(ch14*o1dxini.gt.1.d-30) tauh2(14,i) = 1.d0 / 
     $        (ch14 * o1dxini)
         if(ch14*h2xini.gt.1.d-30) tauo1d(14,i) = 1.d0 / 
     $        (ch14 * h2xini)

         if(ch15*ohxini.gt.1.d-30) tauh2(15,i) = 1.d0 / 
     $        (ch15 * ohxini)
         if(ch15*h2xini.gt.1.d-30) tauoh(15,i) = 1.d0 / 
     $        (ch15 * h2xini)

         if(jdistot8_b(1,i).gt.1.d-30) tauco2(16,i)=1.d0/jdistot8_b(1,i)

         if(jdistot8_b(2,i).gt.1.d-30) tauo2(17,i)=1.d0/jdistot8_b(2,i)

         if(ch18*ohxini.gt.1.d-30) tauh2o2(18,i)=1.d0 / 
     $        (ch18 * ohxini)
         if(ch18*h2o2xini.gt.1.d-30) tauoh(18,i)=1.d0 / 
     $        (ch18 * h2o2xini)

         if(ch19*co2xini.gt.1.d-30)tauo1d(19,i)=1.d0 / 
     $        (ch19 * co2xini)

         if(ch20*o2xini.gt.1.d-30)tauo1d(20,i)= 1.d0 / 
     $        (ch20 * o2xini)

         if(ch21*o2xini*co2xini.gt.1.d-30)tauo3p(21,i)= 1.d0 / 
     $        (ch21 * o2xini * co2xini)
         if(ch21*o3pxini*co2xini.gt.1.d-30) tauo2(21,i) = 1.d0 / 
     $        (ch21 * o3pxini * co2xini)

         !Only if O3, N or ion chemistry requested
         if(chemthermod.ge.1) then
            if(ch22*hxini.gt.1.d-30) tauo3(22,i) = 1.d0 / 
     $           (ch22 * hxini)
            if(ch22*o3xini.gt.1.d-30) tauh(22,i) = 1.d0 / 
     $           (ch22 * o3xini)
            
            if(ch23*ohxini.gt.1.d-30) tauo3(23,i) = 1.d0 / 
     $           (ch23 * ohxini)
            if(ch23*o3xini.gt.1.d-30) tauoh(23,i) = 1.d0 / 
     $           (ch23 * o3xini)
            
            if(ch24 * ho2xini.gt.1.d-30)tauo3(24,i)= 1.d0 / 
     $           (ch24 * ho2xini)
            if(ch24 * o3xini.gt.1.d-30)tauho2(24,i)= 1.d0 / 
     $           (ch24 * o3xini)
            
            if(jdistot8(7,i).gt.1.d-30) tauo3(25,i)=1.d0 /
     $           jdistot8(7,i)

            if(jdistot8_b(7,i).gt.1.d-30) tauo3(26,i)= 1.d0 /
     $           jdistot8_b(7,i)

         endif    !Of chemthermod.ge.1

         if(jdistot8(5,i).gt.1.d-30) tauh2(27,i)= 1.d0/jdistot8(5,i)
         
         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
            if(jdistot8(8,i).gt.1.d-30) taun2(28,i) = 1.d0 / 
     $           jdistot8(8,i)

            if(jdistot8(10,i).gt.1.d-30) tauno(29,i) = 1.d0 /
     $           jdistot8(10,i)

            if(ch30 * noxini.gt.1.d-30) taun(30,i) = 1.d0 / 
     $           (ch30 * noxini)
            if(ch30 * nxini.gt.1.d-30) tauno(30,i) = 1.d0 / 
     $           (ch30 * nxini)

            if(ch31 * o1dxini.gt.1.d-30) taun2(31,i) = 1.d0 / 
     $           (ch31 * o1dxini)
            if(ch31 * n2xini.gt.1.d-30) tauo1d(31,i) = 1.d0 /
     $           (ch31 * n2xini)

            if(ch32 * o2xini.gt.1.d-30) taun(32,i) = 1.d0 /
     $           (ch32 * o2xini)
            if(ch32 * nxini.gt.1.d-30) tauo2(32,i) = 1.d0 /
     $           (ch32 * nxini)
            
            if(ch33 * ohxini.gt.1.d-30) taun(33,i) = 1.d0 /
     $           (ch33 * ohxini)
            if(ch33 * nxini.gt.1.d-30) tauoh(33,i) = 1.d0 /
     $           (ch33 * nxini)

            if(ch34 * o3xini.gt.1.d-30) taun(34,i) = 1.d0 /
     $           (ch34 * o3xini)
            if(ch34 * nxini.gt.1.d-30) tauo3(34,i) = 1.d0 /
     $           (ch34 * nxini)
            
            if(ch35 * ho2xini.gt.1.d-30) taun(35,i) = 1.d0 /
     $           (ch35 * ho2xini)
            if(ch35 * nxini.gt.1.d-30) tauho2(35,i) = 1.d0 /
     $           (ch35 * nxini)
            
            if(ch36 * o3pxini.gt.1.d-30) taun2d(36,i) = 1.d0 /
     $           (ch36 * o3pxini)
            if(ch36 * n2dxini.gt.1.d-30) tauo3p(36,i) = 1.d0 /
     $           (ch36 * n2dxini)
            
            if(ch37 * n2xini.gt.1.d-30) taun2d(37,i) = 1.d0 /
     $           (ch37 * n2xini)
            if(ch37 * n2dxini.gt.1.d-30) taun2(37,i) = 1.d0 /
     $           (ch37 * n2dxini)
            
            if(ch38 * co2xini.gt.1.d-30) taun2d(38,i) = 1.d0 /
     $           (ch38 * co2xini)
            if(ch38 * n2dxini.gt.1.d-30) tauco2(38,i) = 1.d0 /
     $           (ch38 * n2dxini)
            
            if(ch39 * ho2xini.gt.1.d-30) tauno(39,i) = 1.d0 /
     $           (ch39 * ho2xini)
            if(ch39 * noxini.gt.1.d-30) tauho2(39,i) = 1.d0 /
     $           (ch39 * noxini)
            
            if(ch40 * noxini * co2xini.gt.1.d-30) tauo3p(40,i) = 1.d0 /
     $           (ch40 * noxini * co2xini)
            if(ch40 * o3pxini * co2xini.gt.1.d-30) tauno(40,i) = 1.d0 /
     $           (ch40 * o3pxini * co2xini)
            
            if(ch41 * no2xini.gt.1.d-30) tauo3p(41,i) = 1.d0 /
     $           (ch41 * no2xini)
            if(ch41 * o3pxini.gt.1.d-30) tauno2(41,i) = 1.d0 /
     $           (ch41 * o3pxini)
            
            if(ch42 * noxini.gt.1.d-30) tauo3(42,i) = 1.d0 /
     $           (ch42 * noxini)
            if(ch42 * o3xini.gt.1.d-30) tauno(42,i) = 1.d0 /
     $           (ch42 * o3xini)
            
            if(ch43 * no2xini.gt.1.d-30) tauh(43,i) = 1.d0 /
     $           (ch43 * no2xini)
            if(ch43 * hxini.gt.1.d-30) tauno2(43,i) = 1.d0 /
     $           (ch43 * hxini)

            if(jdistot8(13,i).gt.1.d-30) tauno2(44,i) = 1.d0 / 
     $           jdistot8(13,i)

            if(ch45 * nxini.gt.1.d-30) tauo3p(45,i) = 1.d0 / 
     $           (ch45 * nxini)
            if(ch45 * o3pxini.gt.1.d-30) taun(45,i) = 1.d0 /
     $           (ch45 * o3pxini)

         endif    !Of chemthermod.ge.2

c>>>>>>>>>>>>>>>>>>>>>>>>  IONOSPHERE >>>>>>>>>>>>>>>>

         !Only if ion chemistry requested
         if(chemthermod.eq.3) then
            if(ch46 * co2plusxini .gt.1.d-30) tauo2(46,i) = 
     @           1.d0/(ch46*co2plusxini)
            if(ch46 * o2xini .gt.1.d-30) tauco2plus(46,i) = 
     @           1.d0/(ch46*o2xini)

            if ( ch47*o3pxini .gt. 1.d-30 ) tauco2plus(47,i) = 
     @           1.d0/( ch47*o3pxini )
            if ( ch47*co2plusxini .gt. 1.d-30 ) tauo3p(47,i) = 
     @           1.d0/( ch47*co2plusxini )
            
            if ( ch48*o3pxini .gt. 1.d-30 ) tauco2plus(48,i) =
     @           1.d0/(ch48*o3pxini)	   
            if ( ch48*co2plusxini .gt. 1.d-30 ) tauo3p(48,i) =
     @           1.d0/(ch48*co2plusxini)	 

            if ( ch49*electxini .gt. 1.d-30 ) tauo2plus(49,i) =
     @           1.d0/(ch49*electxini)	 

            if ( ch50*co2xini  .gt. 1.d-30 ) tauoplus(50,i) =
     @           1.d0/(ch50*co2xini)
            if ( ch50*oplusxini .gt. 1.d-30 ) tauco2(50,i) =
     @           1.d0/(ch50*oplusxini)

            if ( jion8(1,i,1).gt.1.d-30 ) tauco2(51,i) = 
     $           1.d0 / jion8(1,i,1)

            if ( jion8(1,i,2).gt.1.d-30 ) tauco2(52,i) = 
     $           1.d0 / jion8(1,i,2)

            if ( jion8(1,i,3).gt.1.d-30 ) tauco2(53,i) = 
     $           1.d0 / jion8(1,i,3)

            if ( jion8(1,i,4).gt.1.d-30 ) tauco2(54,i) = 
     $           1.d0 / jion8(1,i,4) 
		
            if ( ch55*electxini .gt. 1.d-30 ) tauco2plus(55,i) =
     @           1.d0/(ch55*electxini)	 

            if ( ch56*oplusxini .gt. 1.d-30 ) tauco2(56,i) =
     @           1.d0/(ch56*oplusxini)
            if ( ch56*co2xini .gt. 1.d-30 ) tauoplus(56,i) =
     @           1.d0/(ch56*co2xini)	
    
            if ( ch57*coplusxini .gt. 1.d-30 ) tauco2(57,i) =
     @           1.d0/(ch57*coplusxini)	 
            if ( ch57*co2xini .gt. 1.d-30 ) taucoplus(57,i) =
     @           1.d0/(ch57*co2xini)

            if ( ch58*coplusxini .gt. 1.d-30 ) tauo3p(58,i) =
     @           1.d0/(ch58*coplusxini)	
            if ( ch58*o3pxini .gt. 1.d-30 ) taucoplus(58,i) =
     @           1.d0/(ch58*o3pxini)

            if ( ch59*cplusxini .gt. 1.d-30 ) tauco2(59,i) =
     @           1.d0/(ch59*cplusxini)	
            if ( ch59*co2xini .gt. 1.d-30 ) taucplus(59,i) =
     @           1.d0/(ch59*co2xini)	   	 	 
         
            if ( jion8(2,i,1).gt.1.d-30 ) tauo2(60,i) = 
     $           1.d0 / jion8(2,i,1)

            if ( jion8(3,i,1).gt.1.d-30 ) tauo3p(61,i) = 
     $           1.d0 / jion8(3,i,1)  

            if ( ch62*co2plusxini .gt. 1.d-30 ) tauno(62,i) =
     @           1.d0/(ch62*co2plusxini)	
            if ( ch62*noxini .gt. 1.d-30 ) tauco2plus(62,i) =
     @           1.d0/(ch62*noxini)	   	 	 
 
            if ( ch63*co2plusxini .gt. 1.d-30 ) taun(63,i) =
     @           1.d0/(ch63*cplusxini)	
            if ( ch63*nxini .gt. 1.d-30 ) tauco2plus(63,i) =
     @           1.d0/(ch63*nxini)	   	 	 
     
            if ( ch64*o2plusxini .gt. 1.d-30 ) tauno(64,i) =
     @           1.d0/(ch64*o2plusxini)	
            if ( ch64*noxini .gt. 1.d-30 ) tauo2plus(64,i) =
     @           1.d0/(ch64*noxini)	   	 	 
   
            if ( ch65*o2plusxini .gt. 1.d-30 ) taun2(65,i) =
     @           1.d0/(ch65*o2plusxini)	
            if ( ch65*n2xini .gt. 1.d-30 ) tauo2plus(65,i) =
     @           1.d0/(ch65*n2xini)	   	 	 
   
            if ( ch66*o2plusxini .gt. 1.d-30 ) taun(66,i) =
     @           1.d0/(ch66*o2plusxini)	
            if ( ch66*nxini .gt. 1.d-30 ) tauo2plus(66,i) =
     @           1.d0/(ch66*nxini)	   	 	 
   
            if ( ch67*oplusxini .gt. 1.d-30 ) taun2(67,i) =
     @           1.d0/(ch67*oplusxini)	
            if ( ch67*n2xini .gt. 1.d-30 ) tauoplus(67,i) =
     @           1.d0/(ch67*n2xini)	   	 	 

            if ( ch68*n2plusxini .gt. 1.d-30 ) tauco2(68,i) =
     @           1.d0/(ch68*n2plusxini)
            if ( ch68*co2xini .gt. 1.d-30 ) taun2plus(68,i) =
     @           1.d0/(ch68*co2xini)	   	 	 
    
            if ( ch69*cplusxini .gt. 1.d-30 ) tauco2(69,i) =
     @           1.d0/(ch69*cplusxini)	
            if ( ch69*co2xini .gt. 1.d-30 ) taucplus(69,i) =
     @           1.d0/(ch69*co2xini)	   	 	 
            
            if ( ch70*n2plusxini .gt. 1.d-30 ) tauco(70,i) =
     @           1.d0/(ch70*n2plusxini)	
            if ( ch70*coxini .gt. 1.d-30 ) taun2plus(70,i) =
     @           1.d0/(ch70*coxini)	   	 	 
            
            if ( ch71*electxini .gt. 1.d-30 ) taun2plus(71,i) =
     @           1.d0/(ch71*electxini)	
            
            if ( ch72*n2plusxini .gt. 1.d-30 ) tauo3p(72,i) =
     @           1.d0/(ch72*n2plusxini)	
            if ( ch72*o3pxini .gt. 1.d-30 ) taun2plus(72,i) =
     @           1.d0/(ch72*o3pxini)	   	 	 
     	 	 
            if ( ch73*coplusxini .gt. 1.d-30 ) tauh(73,i) =
     @           1.d0/(ch73*coplusxini)	
            if ( ch73*hxini .gt. 1.d-30 ) taucoplus(73,i) =
     @           1.d0/(ch73*hxini)	   	 	 
         
            if ( ch74*oplusxini .gt. 1.d-30 ) tauh(74,i) =
     @           1.d0/(ch74*oplusxini)	
            if ( ch74*hxini .gt. 1.d-30 ) tauoplus(74,i) =
     @           1.d0/(ch74*hxini)	
 
            if ( ch75*electxini .gt. 1.d-30 ) taunoplus(75,i) =
     @           1.d0/(ch75*electxini)	   	 	 
            
            if ( ch76*hplusxini .gt. 1.d-30 ) tauo3p(76,i) =
     @           1.d0/(ch76*hplusxini)	
            if ( ch76*o3pxini .gt. 1.d-30 ) tauhplus(76,i) =
     @           1.d0/(ch76*o3pxini)	
	 
            if( jion8(11,i,1).gt.1.d-30 )  tauco(77,i) = 
     $           1.d0 / jion8(11,i,1)

            if( jion8(11,i,2).gt.1.d-30 )  tauco(78,i) = 
     $           1.d0 / jion8(11,i,2)  

         !if( jion8(11,i,3).gt.1.d-30 ) tauco(79,i) = 
!     $        1.d0 / jion8(11,i,3) 

            if( jion8(10,i,1).gt.1.d-30 )  tauno(80,i) = 
     $           1.d0 / jion8(10,i,1) 
         
            if( jion8(8,i,1).gt.1.d-30 )   taun2(81,i)  = 
     $           1.d0 / jion8(8,i,1) 

            if( jion8(8,i,2).gt.1.d-30 )   taun2(82,i)  = 
     $           1.d0 / jion8(8,i,2)
		 
            if( jion8(12,i,1).gt.1.d-30 )  tauh(83,i)  = 
     $           1.d0 / jion8(12,i,1)

            if( jion8(9,i,1).gt.1.d-30 )   taun(84,i)  = 
     $           1.d0 / jion8(9,i,1)

            if ( ch85*nplusxini .gt. 1.d-30 ) tauco2(85,i) =
     @           1.d0/(ch85*nplusxini)	
            if ( ch85*co2xini .gt. 1.d-30 ) taunplus(85,i) =
     @           1.d0/(ch85*co2xini)

            if ( ch86*co2plusxini .gt. 1.d-30) tauh2(86,i) =
     $           1.d0/(ch86*co2plusxini)
            if ( ch86*h2xini .gt. 1.d-30) tauco2plus(86,i) =
     $           1.d0/(ch86*h2xini)
         
            if ( ch87*electxini .gt. 1.d-30) tauhco2plus(87,i) =
     $                  1.d0/(ch87*electxini)

            if ( jion8(9,i,1)*ionsec_nplus(zenit,zx(i)).gt.1.d-30) 
     $           taun(88,i) = 1.d0 / 
     $           (jion8(9,i,1)*ionsec_nplus(zenit,zx(i)))

            if ( jion8(8,i,1)*ionsec_n2plus(zenit,zx(i)).gt.1.d-30) 
     $           taun2(89,i) = 1.d0 / 
     $           (jion8(8,i,1)*ionsec_n2plus(zenit,zx(i))) 
            
            if ( jion8(3,i,1)*ionsec_oplus(zenit,zx(i)).gt.1.d-30) 
     $           tauo3p(90,i) = 1.d0 /
     $           (jion8(3,i,1)*ionsec_oplus(zenit,zx(i)))
            
            if (jion8(11,i,1)*ionsec_coplus(zenit,zx(i)).gt.1.d-30) 
     $           tauco(91,i) = 1.d0 /
     $           (jion8(11,i,1)*ionsec_coplus(zenit,zx(i)))

            if (jion8(1,i,1)*ionsec_co2plus(zenit,zx(i)).gt.1.d-30) 
     $           tauco2(92,i) = 1.d0 /
     $           (jion8(1,i,1)*ionsec_co2plus(zenit,zx(i)))

            if ( jion8(2,i,1)*ionsec_o2plus(zenit,zx(i)).gt.1.d-30) 
     $           tauo2(93,i) = 1.d0 /
     $           (jion8(2,i,1)*ionsec_o2plus(zenit,zx(i)))

         endif                  !Of chemthermod.eq.3
c>>>>>>>>>>>>>>>>>>>>>>>>
	 
         !Minimum lifetime for each species
         tminco2(i)  = 1.d30
         tmino2(i)   = 1.d30
         tmino3p(i)  = 1.d30
         tminco(i)   = 1.d30
         tminh(i)    = 1.d30
         tminoh(i)   = 1.d30
         tminho2(i)  = 1.d30
         tminh2(i)   = 1.d30
         tminh2o(i)  = 1.d30
         tmino1d(i)  = 1.d30
         tminh2o2(i) = 1.d30
         tmino3(i)   = 1.d30
	 tminn(i)    = 1.d30
         tminno(i)   = 1.d30
	 tminn2(i) = 1.d30
         tminn2d(i) = 1.d30
         tminno2(i) = 1.d30
         tminco2plus(i) = 1.d30
	 tminoplus(i)   = 1.d30
	 tmino2plus(i)  = 1.d30
         tmincoplus(i)  = 1.d30
	 tmincplus(i)   = 1.d30
	 tminnplus(i)   = 1.d30
	 tminn2plus(i)  = 1.d30
         tminnoplus(i)  = 1.d30
	 tminhplus(i)   = 1.d30
         tminhco2plus(i)=1.d30

         do j=1,nreact
            tminco2(i)  = min(tminco2(i),tauco2(j,i))
            tmino2(i)   = min(tmino2(i),tauo2(j,i))
            tmino3p(i)  = min(tmino3p(i),tauo3p(j,i))
            tminco(i)   = min(tminco(i),tauco(j,i))
            tminh(i)    = min(tminh(i),tauh(j,i))
            tminoh(i)   = min(tminoh(i),tauoh(j,i))
            tminho2(i)  = min(tminho2(i),tauho2(j,i))
            tminh2(i)   = min(tminh2(i),tauh2(j,i))
            tminh2o(i)  = min(tminh2o(i),tauh2o(j,i))
            tmino1d(i)  = min(tmino1d(i),tauo1d(j,i))
            tminh2o2(i) = min(tminh2o2(i),tauh2o2(j,i))
            tmino3(i)   = min(tmino3(i),tauo3(j,i))
            tminn(i)    = min(tminn(i),taun(j,i))
            tminno(i)   = min(tminno(i),tauno(j,i))
	    tminn2(i)   = min(tminn2(i),taun2(j,i))
	    tminn2d(i)  = min(tminn2d(i),taun2d(j,i))
            tminno2(i)  = min(tminno2(i),tauno2(j,i))
            !
            tminco2plus(i) = min(tminco2plus(i),tauco2plus(j,i))
	    tminoplus(i)   = min(tminoplus(i),tauoplus(j,i))
            tmino2plus(i)  = min(tmino2plus(i),tauo2plus(j,i))
            tmincoplus(i)  = min(tmincoplus(i),taucoplus(j,i))
            tmincplus(i)   = min(tmincplus(i),taucplus(j,i))
	    tminnplus(i)   = min(tminnplus(i),taunplus(j,i))
            tminn2plus(i)  = min(tminn2plus(i),taun2plus(j,i))
            tminnoplus(i)  = min(tminnoplus(i),taunoplus(j,i))
            tminhplus(i)   = min(tminhplus(i),tauhplus(j,i))
            tminhco2plus(i)= min(tminhco2plus(i),tauhco2plus(j,i))
         end do

       !!! Minimum lifetime for all species

         xtmin = 1.d30

         ! Neutrals that can not be in PE

         xtmin = min(xtmin,tminco2(i))
         if(xtmin.eq.tminco2(i)) xcompmin=1
         xtmin = min(xtmin,tmino2(i))
         if(xtmin.eq.tmino2(i)) xcompmin=2
         xtmin = min(xtmin,tmino3p(i))
         if(xtmin.eq.tmino3p(i)) xcompmin=3
         xtmin = min(xtmin,tminco(i))
         if(xtmin.eq.tminco(i)) xcompmin=4

         xtmin = min(xtmin,tminh2(i))
         if(xtmin.eq.tminh2(i)) xcompmin=8
         xtmin = min(xtmin,tminh2o(i))
         if(xtmin.eq.tminh2o(i)) xcompmin=9
  
         xtmin = min(xtmin,tminh2o2(i))
         if(xtmin.eq.tminh2o2(i)) xcompmin=11
         xtmin = min(xtmin,tmino3(i))	
         if(xtmin.eq.tmino3(i)) xcompmin=12
         xtmin = min(xtmin,tminn(i))
         if(xtmin.eq.tminn(i)) xcompmin=13
         xtmin = min(xtmin,tminno(i))
         if(xtmin.eq.tminno(i)) xcompmin=14
         xtmin = min(xtmin,tminn2(i))
         if(xtmin.eq.tminn2(i)) xcompmin=15

         ! Neutrals that can be in PE
         !       [ O1D , OH , HO2 , H , N2D , NO2 ]

         xn_comp_en_EQ = 0
         
         ! H
         h_eq(i)='Y'
         xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! OH
         oh_eq(i)='Y'
         xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! HO2
         ho2_eq(i)='Y'
         xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! O1D
         o1d_eq(i)='Y'
         xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! O3
         o3_eq(i)='N'
!         o3_eq(i)='Y'
!         xn_comp_en_EQ = xn_comp_en_EQ + 1


         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
         ! N2D
            n2d_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! NO2
            no2_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1


         ! NO
            no_eq(i)='N'

!            no_eq(i)='Y'
!            xn_comp_en_EQ = xn_comp_en_EQ + 1

         endif   !Of chemthermod.ge.2

         ! Ions
         !Only if ion chemistry requested

         if(chemthermod.eq.3) then
         ! C+
            cplus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! CO+
            coplus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! O+
!            oplus_eq(i)='N'
         oplus_eq(i)='Y'
         xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! N2+
            n2plus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! H+
!            hplus_eq(i)='N'

         hplus_eq(i)='Y'
         xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! CO2+
            co2plus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! O2+
            o2plus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! NO+
            noplus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! N+
            nplus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1

         ! HCO2+
            hco2plus_eq(i)='Y'
            xn_comp_en_EQ = xn_comp_en_EQ + 1
         endif  !Of chemthermod.eq.3

         return
c END
         end



c**********************************************************************
c**********************************************************************

      subroutine timemarching( ig,i,nlayer,chemthermod,n_comp_en_EQ,
     $     compmin,tmin,timefrac_sec, deltat,fmargin1 ) 

  
c Calculates the timestep of the model, deltat, and verifies if the species
c that can be in PE actually verify or not the PE condition
c
c     jul 2008      MA        Version en subrutina
c         2009      FGG       Adaptation to GCM
c**********************************************************************

      use iono_h
      use param_v4_h, only: tminco2,tmino2,tmino3p,tminco,tminh,tminoh,
     .  tminho2,tminh2,tminh2o,tmino1d,tminh2o2,tmino3,tminn,tminno,
     .  tminno2,tminn2,tminn2d,tminco2plus,tminoplus,tmino2plus,
     .  tmincoplus,tmincplus,tminnplus,tminnoplus,tminn2plus,
     .  tminhplus,tminhco2plus

      implicit none

c     arguments
c
      integer   i,ig,nlayer                              ! I. Layer  
      integer   chemthermod
      integer   n_comp_en_EQ(nlayer)         ! Number of species in PE
      integer   compmin(nlayer)              ! Species with minimum lifetime
      real*8    tmin(nlayer)                 ! Minimum lifetime
      real*8    timefrac_sec                   ! I. 
      real*8    deltat                         ! O. TimeMarching step

c     local variables
c
      integer   j
      real*8    tminaux

      real*8   fmargin1, fmargin2

ccccccccccccccc CODE STARTS 

!      fmargin1=1.
      fmargin2=5.

      !Internal timestep as the minimum of the external (GCM) timestep
      !and the minimum lifetime of the species not in PE divided by a 
      !given factor
      tminaux = min( timefrac_sec, tmin(i)/fmargin1 )
      
      !Given the internal timestep, we verify if the species that can be
      !in PE verify or not the PE condition (lifetime < deltat/fmargin2)

      ! 6 neutral species that can be in PE
      !      [ O1D , OH , HO2 , H , N2D , NO2 ]

      ! O1D
      if ( (o1d_eq(i).eq.'Y') .and. 
     &     (tmino1d(i).gt.tminaux/fmargin2) ) then 
         o1d_eq(i)='N'
         n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
      endif
      ! OH
      if ( (oh_eq(i).eq.'Y') .and. 
     &     (tminoh(i).gt.tminaux/fmargin2) ) then
         oh_eq(i)='N'
         n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
      endif
      ! HO2
      if ( (ho2_eq(i).eq.'Y') .and. 
     &     (tminho2(i).gt.tminaux/fmargin2) ) then
         ho2_eq(i)='N'
         n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
      endif
      ! H
      if ( (h_eq(i).eq.'Y') .and. 
     &     (tminh(i).gt.tminaux/fmargin2) ) then
         h_eq(i)='N'
         n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
      endif
      
      !Only if N or ion chemistry requested
      if(chemthermod.ge.2) then
      ! N2D
         if ( (n2d_eq(i).eq.'Y') .and. 
     &        (tminn2d(i).gt.tminaux/fmargin2) ) then
            n2d_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! NO2
         if ( (no2_eq(i).eq.'Y') .and. 
     &        (tminno2(i).gt.tminaux/fmargin2) ) then
            no2_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
         
      endif   !Of chemthermod.ge.2


      !
      ! 9 ions 
      ! 
      
      !Only if ion chemistry requested
      if(chemthermod.eq.3) then
      ! C+
         if ( (cplus_eq(i).eq.'Y') .and. 
     &        (tmincplus(i).gt.tminaux/fmargin2) ) then
            cplus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! CO+
         if ( (coplus_eq(i).eq.'Y') .and. 
     &        (tmincoplus(i).gt.tminaux/fmargin2) ) then
            coplus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! O+
         if ( (oplus_eq(i).eq.'Y') .and. 
     &        (tminoplus(i).gt.tminaux/fmargin2) ) then
            oplus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! N2+
         if ( (n2plus_eq(i).eq.'Y') .and. 
     &        (tminn2plus(i).gt.tminaux/fmargin2) ) then
            n2plus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! H+
         if ( (hplus_eq(i).eq.'Y') .and. 
     &        (tminhplus(i).gt.tminaux/fmargin2) ) then
            hplus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! CO2+
         if ( (co2plus_eq(i).eq.'Y') .and. 
     &        (tminco2plus(i).gt.tminaux/fmargin2) ) then
            co2plus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! O2+
         if ( (o2plus_eq(i).eq.'Y') .and. 
     &        (tmino2plus(i).gt.tminaux/fmargin2) ) then
            o2plus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! NO+
         if ( (noplus_eq(i).eq.'Y') .and. 
     &        (tminnoplus(i).gt.tminaux/fmargin2) ) then
            noplus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
      ! N+
         if ( (nplus_eq(i).eq.'Y') .and. 
     &        (tminnplus(i).gt.tminaux/fmargin2) ) then
            nplus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif
       ! HCO2+
         if ( (hco2plus_eq(i).eq.'Y') .and. 
     &        (tminhco2plus(i).gt.tminaux/fmargin2) ) then
            hco2plus_eq(i)='N'
            n_comp_en_EQ(i) = n_comp_en_EQ(i) - 1
         endif

      endif   !Of chemthermod.eq.3



         ! And we set the internal timestep
         !
      deltat = tminaux

      return
c END
      end



c**********************************************************************
c**********************************************************************

      subroutine prodsandlosses ( ig,i,nlayer,chemthermod,zenit,zx,
     &                 jdistot8, jdistot8_b, jion8,
     &                              co2xinput, o2xinput, o3pxinput,
     &                              coxinput,  h2xinput, o3xinput,
     &                              h2oxinput, nxinput,  noxinput, 
     &                              h2o2xinput, n2xinput, 
     &                           o1dxinput, ohxinput,  ho2xinput,
     &                           hxinput,   n2dxinput, no2xinput,
     &                 co2plusxinput,  o2plusxinput, coplusxinput,
     &                 oplusxinput,    cplusxinput,  noplusxinput,
     &                 n2plusxinput,   hplusxinput,  nplusxinput,
     &                 hco2plusxinput,electxinput )


  
c Calculates the productions and losses of each species in each reaction
c and each altitude
c
c     jul 2008      MA        Version en subrutina
c     apr 2009      FGG       Compact version to save CPU time, adaptation
c                             to GCM
c**********************************************************************

      use param_v4_h
      implicit none

c     arguments
c
      integer   ig,nlayer
      integer   i                                 ! I. Layer  
      integer   chemthermod
      real      zx(nlayer)
      real      zenit
      real*8    jdistot8(nabs,nlayer)
      real*8    jdistot8_b(nabs,nlayer)
      real*8    jion8(nabs,nlayer,4)
      real*8    co2xinput,o2xinput,o3pxinput,coxinput
      real*8    ho2xinput,h2xinput,hxinput,ohxinput
      real*8    h2o2xinput,o1dxinput,o3xinput,h2oxinput
      real*8    nxinput,noxinput,n2xinput,n2dxinput,no2xinput
      real*8    co2plusxinput,coplusxinput,oplusxinput,o2plusxinput
      real*8    cplusxinput,noplusxinput,n2plusxinput,hplusxinput
      real*8    electxinput,nplusxinput,hco2plusxinput

c     local variables
c
      integer   j 

      external  ionsec_nplus
      real*8    ionsec_nplus

      external  ionsec_n2plus
      real*8    ionsec_n2plus

      external  ionsec_oplus
      real*8    ionsec_oplus
     
      external  ionsec_coplus
      real*8    ionsec_coplus

      external  ionsec_co2plus
      real*8    ionsec_co2plus

      external  ionsec_o2plus
      real*8    ionsec_o2plus

      logical firstcall
      save firstcall
      data firstcall /.true./

ccccccccccccccc CODE STARTS 

      !Initialization'
!      if (firstcall) then
!        firstcall= .false.
        do j=1,nreact
           Lco2(i,j) = 0.d0
           Pco2(i,j) = 0.d0
           Lo2(i,j) = 0.d0
           Po2(i,j) = 0.d0
           Lo3p(i,j) = 0.d0
           Po3p(i,j) = 0.d0
           Lco(i,j) = 0.d0
           Pco(i,j) = 0.d0
           Ph(i,j) = 0.d0
           Lh(i,j) = 0.d0
           Poh(i,j) = 0.d0
           Loh(i,j) = 0.d0
           Pho2(i,j) = 0.d0
           Lho2(i,j) = 0.d0
           Ph2(i,j) = 0.d0
           Lh2(i,j) = 0.d0
           Ph2o(i,j) = 0.d0
           Lh2o(i,j) = 0.d0
           Po1d(i,j) = 0.d0
           Lo1d(i,j) = 0.d0
           Ph2o2(i,j) = 0.d0
           Lh2o2(i,j) = 0.d0
           Po3(i,j) = 0.d0
           Lo3(i,j) = 0.d0
           Pn(i,j) = 0.d0
           Ln(i,j) = 0.d0
           Pno(i,j) = 0.d0
           Lno(i,j) = 0.d0
           Pn2(i,j) = 0.d0
           Ln2(i,j) = 0.d0
           Pn2d(i,j) = 0.d0
           Ln2d(i,j) = 0.d0
           Pno2(i,j) = 0.d0
           Lno2(i,j) = 0.d0			   
           Lco2plus(i,j) = 0.d0
           Pco2plus(i,j) = 0.d0
           Loplus(i,j) = 0.d0
           Poplus(i,j) = 0.d0
           Lo2plus(i,j) = 0.d0
           Po2plus(i,j) = 0.d0
           Pcoplus(i,j) = 0.d0
           Lcoplus(i,j) = 0.d0
           Pcplus(i,j) = 0.d0
           Lcplus(i,j) = 0.d0
           Pnplus(i,j) = 0.d0
           Lnplus(i,j) = 0.d0
           Pnoplus(i,j) = 0.d0
           Lnoplus(i,j) = 0.d0
           Pn2plus(i,j) = 0.d0
           Ln2plus(i,j) = 0.d0
           Phplus(i,j) = 0.d0
           Lhplus(i,j) = 0.d0
           Phco2plus(i,j) = 0.d0
           Lhco2plus(i,j) = 0.d0
           Pelect(i,j) = 0.d0
           Lelect(i,j) = 0.d0
        end do
        Pco2tot(i) = 0.d0
        Lco2tot(i) = 0.d0
        Po2tot(i) = 0.d0
        Lo2tot(i) = 0.d0
        Po3ptot(i) = 0.d0
        Lo3ptot(i) = 0.d0
        Pcotot(i) = 0.d0
        Lcotot(i) = 0.d0
        Phtot(i) = 0.d0
        Lhtot(i) = 0.d0
        Pohtot(i) = 0.d0
        Lohtot(i) = 0.d0
        Pho2tot(i) = 0.d0
        Lho2tot(i) = 0.d0
        Ph2tot(i) = 0.d0
        Lh2tot(i) = 0.d0
        Ph2otot(i) = 0.d0
        Lh2otot(i) = 0.d0
        Po1dtot(i) = 0.d0
        Lo1dtot(i) = 0.d0
        Ph2o2tot(i) = 0.d0
        Lh2o2tot(i) = 0.d0
        Po3tot(i) = 0.d0
        Lo3tot(i) = 0.d0
        Pntot(i) = 0.d0
        Lntot(i) = 0.d0
        Pnotot(i) = 0.d0
        Lnotot(i) = 0.d0
        Pn2tot(i) = 0.d0
        Ln2tot(i) = 0.d0
        Pn2dtot(i) = 0.d0
        Ln2dtot(i) = 0.d0
        Pno2tot(i) = 0.d0
        Lno2tot(i) = 0.d0
                                !
        Pco2plustot(i) = 0.d0
        Lco2plustot(i) = 0.d0
        Loplustot(i) = 0.d0
        Poplustot(i) = 0.d0
        Lo2plustot(i) = 0.d0
        Po2plustot(i) = 0.d0
        Lcoplustot(i) = 0.d0
        Pcoplustot(i) = 0.d0
        Lcplustot(i) = 0.d0
        Pcplustot(i) = 0.d0
        Lnplustot(i) = 0.d0
        Pnplustot(i) = 0.d0
        Lnoplustot(i) = 0.d0
        Pnoplustot(i) = 0.d0
        Ln2plustot(i) = 0.d0
        Pn2plustot(i) = 0.d0
        Lhplustot(i) = 0.d0
        Phplustot(i) = 0.d0
        Lhco2plustot(i) = 0.d0
        Phco2plustot(i) = 0.d0
        Pelecttot(i) = 0.0d0
        Lelecttot(i) = 0.0d0
!      endif

      

      !!! Productions and losses reaction by reaction

c     R1: CO2 + hv -> CO + O

      Lco2(i,1) = jdistot8(1,i)
      Pco(i,1) = co2xinput * jdistot8(1,i)
      Po3p(i,1) = Pco(i,1)!co2xinput * jdistot8(1,i)


c     R16(1b): CO2 + hv -> CO + O1D
      
      Lco2(i,16) = jdistot8_b(1,i)
      Pco(i,16) = co2xinput * jdistot8_b(1,i)
      Po1d(i,16) = Pco(i,16)!co2xinput * jdistot8_b(1,i)


c     R2: H + O2 + CO2 -> HO2 + CO2

      Lh(i,2) = ch2 * o2xinput * co2xinput
      Lo2(i,2) = ch2 * hxinput * co2xinput
      Pho2(i,2) = ch2 * hxinput * o2xinput * co2xinput


c     R3: O + HO2 -> OH + O2

      Lo3p(i,3) = ch3 * ho2xinput 
      Lho2(i,3) = ch3 * o3pxinput
      Poh(i,3) = ch3 * o3pxinput * ho2xinput
      Po2(i,3) = Poh(i,3)!ch3 * o3pxinput * ho2xinput
      

c     R4: CO + OH -> CO2 + H
      
      Lco(i,4) = ch4 * ohxinput
      Loh(i,4) = ch4 * coxinput
      Pco2(i,4) = ch4 * coxinput * ohxinput
      Ph(i,4) = Pco2(i,4)!ch4 * coxinput * ohxinput


c     R5: 2HO2 -> H2O2 + O2

      Lho2(i,5) = 2.d0 * ch5 * ho2xinput 
      Po2(i,5) = ch5 * ho2xinput * ho2xinput 
      Ph2o2(i,5) = Po2(i,5)!ch5 * ho2xinput * ho2xinput


c     R6: H2O2 + hv -> 2OH

      Lh2o2(i,6) = jdistot8(6,i) 
      Poh(i,6) = 2.d0 * jdistot8(6,i) * h2o2xinput


c     R7: OH + HO2 -> H2O + O2

      Loh(i,7) = ch7 * ho2xinput
      Lho2(i,7) = ch7 * ohxinput
      Po2(i,7) = ch7 * ohxinput * ho2xinput
      Ph2o(i,7) = Po2(i,7)!ch7 * ohxinput * ho2xinput
      

c     R8: H20 + hv -> H + OH

      Lh2o(i,8) = jdistot8(4,i)
      Ph(i,8) = jdistot8(4,i) * h2oxinput
      Poh(i,8) = Ph(i,8)!jdistot8(4,i) * h2oxinput


c     R9: O(1D) + H2O -> 2OH

      Lo1d(i,9) = ch9 * h2oxinput
      Lh2o(i,9) = ch9 * o1dxinput
      Poh(i,9) = 2.d0 * ch9 * o1dxinput * h2oxinput 


c     R10: 2O + CO2 -> O2 + CO2

      Lo3p(i,10) = 2.d0 * ch10 * o3pxinput * co2xinput  
      Po2(i,10) = ch10 * o3pxinput * o3pxinput * co2xinput
      

c     R11: O + OH -> O2 + H

      Lo3p(i,11) = ch11 * ohxinput
      Loh(i,11) = ch11 * o3pxinput
      Po2(i,11) = ch11 * o3pxinput * ohxinput
      Ph(i,11) = Po2(i,11)!ch11 * o3pxinput * ohxinput


c     R12: O2 + hv -> 2O

      Lo2(i,12) = jdistot8(2,i)
      Po3p(i,12) = 2.d0 * jdistot8(2,i) * o2xinput


c     R17(12b): O2 + hv -> O + O1D

      Lo2(i,17) = jdistot8_b(2,i)
      Po3p(i,17) = jdistot8_b(2,i) * o2xinput
      Po1d(i,17) = Po3p(i,17)!jdistot8_b(2,i) * o2xinput


c     R13: H + HO2 -> H2 + O2

      Lh(i,13) = ch13 * ho2xinput
      Lho2(i,13) = ch13 * hxinput
      Ph2(i,13) = ch13 * hxinput * ho2xinput
      Po2(i,13) = Ph2(i,13)!ch13 * hxinput * ho2xinput


c     R14: O(1D) + H2 -> H + OH

      Lo1d(i,14) = ch14 * h2xinput
      Lh2(i,14) = ch14 * o1dxinput
      Ph(i,14) = ch14 * o1dxinput * h2xinput
      Poh(i,14) = Ph(i,14)!ch14 * o1dxinput * h2xinput


c     R15: OH + H2 -> H + H20

      Loh(i,15) = ch15 * h2xinput
      Lh2(i,15) = ch15 * ohxinput
      Ph(i,15) = ch15 * ohxinput * h2xinput
      Ph2o(i,15) = Ph(i,15)!ch15 * ohxinput * h2xinput


c     R18: OH + H2O2 -> H2O + HO2

      Loh(i,18) = ch18 * h2o2xinput
      Lh2o2(i,18) = ch18 * ohxinput
      Ph2o(i,18) = ch18 * ohxinput * h2o2xinput
      Pho2(i,18) = Ph2o(i,18)!ch18 * ohxinput * h2o2xinput


c     R19: O(1D) + CO2 -> O + CO2

      Lo1d(i,19) = ch19 * co2xinput
      Po3p(i,19) = ch19 * o1dxinput * co2xinput


c     R20: O(1D) + O2 -> O + O2

      Lo1d(i,20) = ch20 * o2xinput
      Po3p(i,20) = ch20 * o1dxinput * o2xinput


c     R21: O + O2 + CO2 -> O3 + CO2

      Lo3p(i,21) = ch21 * o2xinput * co2xinput
      Lo2(i,21) = ch21 * o3pxinput * co2xinput
      Po3(i,21) = ch21 * o3pxinput * o2xinput * co2xinput


      !Only if O3, N or ion chemistry requested
      if(chemthermod.ge.1) then

c     R22: O3 + H -> OH + O2

         Lo3(i,22) = ch22 * hxinput
         Lh(i,22) = ch22 * o3xinput
         Poh(i,22) = ch22 * o3xinput * hxinput
         Po2(i,22) = Poh(i,22)  !ch22 * o3xinput * hxinput


c     R23: O3 + OH -> HO2 + O2

         Lo3(i,23) = ch23 * ohxinput
         Loh(i,23) = ch23 * o3xinput
         Pho2(i,23) = ch23 * ohxinput * o3xinput
         Po2(i,23) = Pho2(i,23) !ch23 * ohxinput * o3xinput


c     R24: O3 + HO2 -> OH + 2O2

         Lo3(i,24) = ch24 * ho2xinput
         Lho2(i,24) = ch24 * o3xinput
         Poh(i,24) = ch24 * o3xinput * ho2xinput
         Po2(i,24) = Poh(i,24)  !2.d0 * ch24 * o3xinput * ho2xinput


c     R25: O3 + hv -> O2 + O3P

         Lo3(i,25) = jdistot8(7,i)
         Po2(i,25) = jdistot8(7,i) * o3xinput
         Po3p(i,25) = Po2(i,25) !jdistot8(7,i) * o3xinput


c     R26 (R25_b): O3 + hv -> O2 + O1D

         Lo3(i,26) = jdistot8_b(7,i)
         Po2(i,26) = jdistot8_b(7,i) * o3xinput
         Po1d(i,26) = Po2(i,26) !jdistot8_b(7,i) * o3xinput

      endif    !Of chemthermod.ge.1


c     R27: H2 + hv -> 2H

      Lh2(i,27) = jdistot8(5,i)
      Ph(i,27) = 2.d0 * jdistot8(5,i) * h2xinput
			

      !Only if N or ion chemistry requested
      if(chemthermod.ge.2) then

c     R28: N2 + hv -> N + N2D

         Ln2(i,28) = jdistot8(8,i)
         Pn(i,28) = jdistot8(8,i) * n2xinput
         Pn2d(i,28) = Pn(i,28)  !jdistot8(8,i) * n2xinput


c     R29: NO + hv -> N + O
            
         Lno(i,29) = jdistot8(10,i)
         Pn(i,29) = jdistot8(10,i) * noxinput
         Po3p(i,29) = Pn(i,29)  !jdistot8(10,i) * noxinput


c     R30: N + NO -> N2 + O

         Ln(i,30) = ch30 * noxinput
         Lno(i,30) = ch30 * nxinput
         Pn2(i,30) = ch30 * nxinput * noxinput
         Po3p(i,30) = Pn2(i,30) !ch30 * nxinput * noxinput


c     R31: N2 + O1D -> N2 + O

         Lo1d(i,31) = ch31 * n2xinput
         Po3p(i,31) = ch31 * n2xinput * o1dxinput


c     R32: N + O2 -> NO + O
            
         Ln(i,32) = ch32 * o2xinput
         Lo2(i,32) = ch32 * nxinput
         Pno(i,32) = ch32 * o2xinput * nxinput
         Po3p(i,32) = Pno(i,32) !ch32 * o2xinput * nxinput

      
c     R33: N + OH -> NO + H

         Ln(i,33) = ch33 * ohxinput
         Loh(i,33) = ch33 * nxinput
         Pno(i,33) = ch33 * nxinput * ohxinput
         Ph(i,33) = Pno(i,33)   !Pch33 * nxinput * ohxinput


c     R34: N + O3 -> NO + O2

         Ln(i,34) = ch34 * o3xinput
         Lo3(i,34) = ch34 * nxinput
         Pno(i,34) = ch34 * nxinput * o3xinput
         Po2(i,34) = Pno(i,34)  !ch34 * nxinput * o3xinput


c     R35: N + HO2 -> NO + OH

         Ln(i,35) = ch35 * ho2xinput
         Lho2(i,35) = ch35 * nxinput
         Pno(i,35) = ch35 * nxinput * ho2xinput
         Poh(i,35) = Pno(i,35)  !ch35 * nxinput * ho2xinput


c     R36: N2D + O -> N + O

         Ln2d(i,36) = ch36 * o3pxinput
         Pn(i,36) = ch36 * n2dxinput * o3pxinput
         

c     R37: N2D + N2 -> N + N2

         Ln2d(i,37) = ch37 * n2xinput
         Pn(i,37) = ch37 * n2dxinput * n2xinput


c     R38: N2D + CO2 -> NO + CO

         Ln2d(i,38) = ch38 * co2xinput
         Lco2(i,38) = ch38 * n2dxinput
         Pno(i,38) = ch38 * n2dxinput * co2xinput
         Pco(i,38) = Pno(i,38)  !ch38 * n2dxinput * co2xinput


c     R39: NO + HO2 -> NO2 + OH
            
         Lno(i,39) = ch39 * ho2xinput
         Lho2(i,39) = ch39 * noxinput
         Pno2(i,39) = ch39 * noxinput * ho2xinput
         Poh(i,39) = Pno2(i,39) !ch39 * noxinput * ho2xinput


c     R40: O + NO + CO2 -> NO2 + CO2

         Lo3p(i,40) = ch40 * noxinput * co2xinput
         Lno(i,40) = ch40 * o3pxinput * co2xinput
         Pno2(i,40) = ch40 * o3pxinput * noxinput * co2xinput
      

c     R41: O + NO2 -> NO + O2

         Lo3p(i,41) = ch41 * no2xinput
         Lno2(i,41) = ch41 * o3pxinput
         Pno(i,41) = ch41 * o3pxinput * no2xinput
         Po2(i,41) = Pno(i,41)  !ch41 * o3pxinput * no2xinput


c     R42: NO + O3 -> NO2 + O2

         Lno(i,42) = ch42 * o3xinput
         Lo3(i,42) = ch42 * noxinput
         Pno2(i,42) = ch42 * noxinput * o3xinput
         Po2(i,42) = Pno2(i,42) !ch42 * noxinput * o3xinput


c     R43: H + NO2 -> NO + OH

         Lh(i,43) = ch43 * no2xinput
         Lno2(i,43) = ch43 * hxinput
         Pno(i,43) = ch43 * no2xinput * hxinput
         Poh(i,43) = Pno(i,43)  !ch43 * no2xinput * hxinput


c     R44: NO2 + hv -> NO + O

         Lno2(i,44) = jdistot8(13,i)
         Pno(i,44) = jdistot8(13,i) * no2xinput
         Po3p(i,44) = Pno(i,44) !jdistot8(13,i) * no2xinput
      

c     R45: N + O -> NO

         Ln(i,45) = ch45 * o3pxinput
         Lo3p(i,45) = ch45 * nxinput
         Pno(i,45) = ch45 * o3pxinput * nxinput

      endif    !Of chemthermod.ge.2


      ! >>>> IONOSFERA 

      !Only if ion chemistry requested
      if(chemthermod.eq.3) then

c     R46: CO2+  + O2  ->  CO2 + O2+

         Lco2plus(i,46) = ch46 * o2xinput
         Lo2(i,46)      = ch46 * co2plusxinput
         Pco2(i,46)     = ch46 * co2plusxinput * o2xinput
         Po2plus(i,46)  = Pco2(i,46) !ch46 * co2plusxinput * o2xinput


c     R47: CO2+ + O ->  O+ + CO2

         Lco2plus(i,47) = ch47 * o3pxinput
         Lo3p(i,47)     = ch47 * co2plusxinput
         Pco2(i,47)     = ch47 * co2plusxinput * o3pxinput
         Poplus(i,47)   = Pco2(i,47) !ch47 * co2plusxinput * o3pxinput


c     R48:  CO2+  O  ->    O2+  +  CO

         Lco2plus(i,48) = ch48 * o3pxinput  
         Lo3p(i,48)     = ch48 * co2plusxinput
         Po2plus(i,48)  = ch48 * co2plusxinput * o3pxinput	
         Pco(i,48)      = Po2plus(i,48) !ch48 * co2plusxinput * o3pxinput


c     R49:  O2+  +  e -> O  + O

         Lo2plus(i,49) = ch49 * electxinput
         Lelect(i,49) = ch49 * o2plusxinput
         Po3p(i,49)    = ch49 * o2plusxinput * electxinput *2.d0

	    	
c     R50:   O+  + CO2   -> O2+  + CO
      
         Loplus(i,50) = ch50 * co2xinput
         Lco2(i,50)   = ch50 * oplusxinput
         Po2plus(i,50)= ch50 * oplusxinput * co2xinput
         Pco(i,50)    = Po2plus(i,50) !ch50 * oplusxinput * co2xinput


c     R51: CO2 + hv -> CO2+ + e

         Lco2(i,51)     = jion8(1,i,1)
         Pco2plus(i,51) = co2xinput * jion8(1,i,1)
         Pelect(i,51)   = Pco2plus(i,51) !co2xinput * jion8(1,i,1)


c     R52: CO2 + hv --> O+  + CO + e
      
         Lco2(i,52)   = jion8(1,i,2)
         Pco(i,52)    = co2xinput * jion8(1,i,2)
         Poplus(i,52) = Pco(i,52) !co2xinput * jion8(1,i,2)
         Pelect(i,52) = Pco(i,52) !co2xinput * jion8(1,i,2)


c     R53:  CO2 + hv -->  CO+  + O  + e

         Lco2(i,53)    = jion8(1,i,3)
         Pcoplus(i,53) = co2xinput * jion8(1,i,3)
         Po3p(i,53)    = Pcoplus(i,53) !co2xinput * jion8(1,i,3)
         Pelect(i,53)  = Pcoplus(i,53) !co2xinput * jion8(1,i,3)
      

c     R54:  CO2 + hv  -->   C+   +  O2 + e

         Lco2(i,54)   = jion8(1,i,4)
         Pcplus(i,54) = co2xinput * jion8(1,i,4)
         Po2(i,54)   = Pcplus(i,54) !co2xinput * jion8(1,i,4)
         Pelect(i,54)  = Pcplus(i,54) !co2xinput * jion8(1,i,4)

      
c     R55:   CO2+ +  e -->   CO  +   O

         Lco2plus(i,55) = ch55 * electxinput
         Lelect(i,55) = ch55 * co2plusxinput
         Pco(i,55)      = ch55 * co2plusxinput * electxinput
         Po3p(i,55)     = Pco(i,55) !ch55 * co2plusxinput * electxinput


c     R56:    O+ + CO2  -->  O2   +   CO+   Probablemente no se dá

         Lco2(i,56)   = ch56 * oplusxinput
         Loplus(i,56) = ch56 * co2xinput
         Po2(i,56)    = ch56 * co2xinput * oplusxinput
         Pcoplus(i,56)= Po2(i,56) !ch56 * co2xinput * oplusxinput


c     R57:    CO+ + CO2  --> CO2+ + CO

         Lco2(i,57)    = ch57 * coplusxinput
         Lcoplus(i,57) = ch57 * co2xinput
         Pco2plus(i,57)= ch57 *  coplusxinput * co2xinput
         Pco(i,57)     = Pco2plus(i,57) !ch57 *  coplusxinput * co2xinput


c     R58:   CO+ + O  -->   O+  + CO

         Lo3p(i,58)    = ch58 * coplusxinput
         Lcoplus(i,58) = ch58 * o3pxinput
         Poplus(i,58)  = ch58 * coplusxinput * o3pxinput
         Pco(i,58)     = Poplus(i,58) !ch58 * coplusxinput * o3pxinput

      
c     R59:  C+  +  CO2 -->  CO+  + CO 

         Lco2(i,59)   = ch59 * cplusxinput
         Lcplus(i,59) = ch59 * co2xinput
         Pcoplus(i,59)= ch59 *  cplusxinput * co2xinput
         Pco(i,59)    = Pcoplus(i,59) !ch59 *  cplusxinput * co2xinput


c     R60:   O2 + hv   -->   O2+   +  e

         Lo2(i,60)     = jion8(2,i,1)
         Po2plus(i,60) = o2xinput * jion8(2,i,1)
         Pelect(i,60)  = Po2plus(i,60) !o2xinput * jion8(2,i,1)


c     R61:   O + hv    -->    O+   +  e
         
         Lo3p(i,61)    = jion8(3,i,1)
         Poplus(i,61)  = o3pxinput * jion8(3,i,1)
         Pelect(i,61)  = Poplus(i,61) !o3pxinput * jion8(3,i,1)


c     R62 :   CO2+  +  NO   -->  NO+  +  CO2

         Lco2plus(i,62)    = ch62 * noxinput
         Lno(i,62)         = ch62 * co2plusxinput
         Pnoplus(i,62)     = ch62 *  co2plusxinput * noxinput
         Pco2(i,62)        = Pnoplus(i,62) !ch62 *  co2plusxinput * noxinput
         

c     R63 :   CO2+  +  N    -->  NO  + CO+

         Lco2plus(i,63)    = ch63 * nxinput
         Ln(i,63)          = ch63 * co2plusxinput
         Pcoplus(i,63)     = ch63 * co2plusxinput * nxinput
         Pno(i,63)         = Pcoplus(i,63) !ch63 * co2plusxinput * nxinput


c     R64 :   O2+  +  NO    -->  NO+  + O2

         Lo2plus(i,64)    = ch64 * noxinput
         Lno(i,64)        = ch64 * o2plusxinput
         Pnoplus(i,64)    = ch64 * o2plusxinput * noxinput
         Po2(i,64)        = Pnoplus(i,64) !ch64 * o2plusxinput * noxinput


c     R65 :   O2+  +  N2    -->  NO+  + NO

         Lo2plus(i,65)    = ch65 * n2xinput
         Ln2(i,65)        = ch65 * o2plusxinput
         Pnoplus(i,65)    = ch65 * o2plusxinput * n2xinput
         Pno(i,65)        = Pnoplus(i,65) !ch65 * o2plusxinput * n2xinput


c     R66:   O2+  +  N    -->  NO+  + O

         Lo2plus(i,66)    = ch66 * nxinput
         Ln(i,66)         = ch66 * o2plusxinput
         Pnoplus(i,66)    = ch66 * o2plusxinput * nxinput
         Po3p(i,66)       = Pnoplus(i,66) !ch66 * o2plusxinput * nxinput


c     R67 :   O+  +  N2    -->  NO+  + N

         Loplus(i,67)    = ch67 * n2xinput
         Ln2(i,67)       = ch67 * oplusxinput
         Pnoplus(i,67)   = ch67 * oplusxinput * n2xinput
         Pn(i,67)        = Pnoplus(i,67) !ch67 * oplusxinput * n2xinput


c     R68 :   N2+  +  CO2    -->  CO2+  + N2

         Ln2plus(i,68)    = ch68 * co2xinput
         Lco2(i,68)       = ch68 * n2plusxinput
         Pco2plus(i,68)   = ch68 * n2plusxinput * co2xinput
         Pn2(i,68)        = Pco2plus(i,68) !ch68 * n2plusxinput * co2xinput


c     R69 :   N2+  +  O3p    -->  NO+  + N

         Ln2plus(i,69)    = ch69 * o3pxinput
         Lo3p(i,69)       = ch69 * n2plusxinput
         Pnoplus(i,69)    = ch69 * n2plusxinput * o3pxinput
         Pn(i,69)         = Pnoplus(i,69) !ch69 * n2plusxinput * o3pxinput


c     R70 :   N2+  +  CO    -->  N2  + CO+

         Ln2plus(i,70)    = ch70 *  coxinput
         Lco(i,70)        = ch70 *  n2plusxinput
         Pcoplus(i,70)    = ch70 *  n2plusxinput * coxinput
         Pn2(i,70)        = Pcoplus(i,70) !ch70 *  n2plusxinput * coxinput


c     R71 :   N2+  +  e-    -->  N  + N

         Ln2plus(i,71)   = ch71 * electxinput
         Lelect(i,71)    = ch71 * n2plusxinput
         Pn(i,71)        = 2. * ch71 *  n2plusxinput * electxinput


c     R72 :   N2+  +  O3p    -->  O+  + N2

         Ln2plus(i,72) = ch72 * o3pxinput
         Lo3p(i,72)    = ch72 * n2plusxinput
         Poplus(i,72)  = ch72 *  n2plusxinput * o3pxinput
         Pn2(i,72)     = Poplus(i,72) !ch72 *  n2plusxinput * o3pxinput


c     R73 :   CO+  +  H    -->  H+  + CO
         
         Lcoplus(i,73) = ch73 * hxinput
         Lh(i,73)      = ch73 * coplusxinput
         Phplus(i,73)  = ch73 * coplusxinput * hxinput
         Pco(i,73)     = Phplus(i,73) !ch73 * coplusxinput * hxinput

c     R74 :   O+  +  H    -->  H+  + O

         Loplus(i,74) = ch74 * hxinput
         Lh(i,74)     = ch74 * oplusxinput
         Phplus(i,74) = ch74*  oplusxinput * hxinput
         Po3p(i,74)   = Phplus(i,74) !ch74 * oplusxinput * hxinput
		

c     R75:   NO+  +  e-  -->  N  +  O

         Lnoplus(i,75)    = ch75 * electxinput
         Lelect(i,75)    = ch75 * noplusxinput
         Pn(i,75)= ch75 *  noplusxinput * electxinput
         Po3p(i,75)= Pn(i,75)   !ch75 *  noplusxinput * electxinput


c     R76:   H+  +  O3p  -->  O+  +  H

         Lhplus(i,76)  = ch76 * o3pxinput
         Lo3p(i,76)    = ch76 * hplusxinput
         Poplus(i,76)  = ch76 *  hplusxinput * o3pxinput
         Ph(i,76)      = Poplus(i,76) !ch76 *  hplusxinput * o3pxinput


c     R77:   CO + hv   -->  CO+ + e-

         Lco(i,77)    = jion8(11,i,1)
         Pcoplus(i,77) = coxinput * jion8(11,i,1)
         Pelect(i,77)  = Pcoplus(i,77) !coxinput * jion8(11,i,1)


c     R78:   CO + hv   -->  C+   +   O    + e-

         Lco(i,78)     = jion8(11,i,2)
         Pcplus(i,78)  = coxinput * jion8(11,i,2)
         Po3p(i,78)    = Pcplus(i,78) !coxinput * jion8(11,i,2)
         Pelect(i,78)  = Pcplus(i,78) !coxinput * jion8(11,i,2)


c     R79:   CO + hv   -->  C  +   O+   + e-

!     Lco(i,79)    = jion8(11,i,3)
!     Pc(i,79) = coxinput * jion8(11,i,3)
!     Poplus(i,79) = Pc(i,79)!Pc(i,79)coxinput * jion8(11,i,3)
!     Pelect(i,79)  = Pc(i,79)!coxinput * jion8(11,i,3)


c     R80:   NO + hv   -->  NO+  +  e-

         Lno(i,80)    = jion8(10,i,1)
         Pnoplus(i,80) = noxinput * jion8(10,i,1)
         Pelect(i,80)  = Pnoplus(i,80) !noxinput * jion8(10,i,1)


c     R81:   N2 + hv   -->  N2+  +   e-

         Ln2(i,81)    = jion8(8,i,1)
         Pn2plus(i,81) = n2xinput * jion8(8,i,1)
         Pelect(i,81)  = Pn2plus(i,81) !n2xinput * jion8(8,i,1)


c     R82:   N2 + hv   -->  N+  +   N  +  e-

         Ln2(i,82)     = jion8(8,i,2)
         Pnplus(i,82)  = n2xinput * jion8(8,i,2)
         Pn(i,82)      = Pnplus(i,82) !n2xinput * jion8(8,i,2)
         Pelect(i,82)  = Pnplus(i,82) !n2xinput * jion8(8,i,2)


c     R83:   H + hv   -->   H+   +  e

         Lh(i,83)     = jion8(12,i,1)
         Phplus(i,83)     = hxinput * jion8(12,i,1)
         Pelect(i,83)     = Phplus(i,83) !hxinput * jion8(12,i,1)


c     R84:   N + hv   -->   N+   +  e

         Ln(i,84)     = jion8(9,i,1)
         Pnplus(i,84)     = nxinput * jion8(9,i,1)
         Pelect(i,84)     = Pnplus(i,84) !nxinput * jion8(9,i,1)


c     R85:   N+ + CO2   -->   CO2+  +   N

         Pn(i,85)       = ch85 * co2xinput * nplusxinput
         Pco2plus(i,85) = Pn(i,85) !ch85 * co2xinput * nplusxinput
         Lnplus(i,85)   = ch85 * co2xinput
         Lco2(i,85)     = ch85 * nplusxinput


c     R86:   H2 + CO2+   --> H + HCO2+
         Ph(i,86) = ch86 * co2plusxinput * h2xinput
         Lco2plus(i,86) = ch86 * h2xinput
         Lh2(i,86) = ch86 * co2plusxinput
                                !!! 

c     R87: HCO2+ + e --> H + CO2
         Ph(i,87) = ch87 * hco2plusxinput * electxinput
         Pco2(i,87) = Ph(i,87)
         Lhco2plus(i,87) = ch87 * electxinput
         Lelect(i,87) = ch87 * hco2plusxinput


c     R88:   N  +  e  -->  N+  +  e  +  e
      
         Pnplus(i,88)   = ionsec_nplus(zenit,zx(i))*Pnplus(i,84)
         Pelect(i,88)    = Pnplus(i,88)
         if(nxinput.ge.1.d-20) then
            Ln(i,88)    = Pnplus(i,88)/nxinput
         else
            Ln(i,88)    = 1.d-30
         endif


c     R89:   N2  +  e  -->  N2+  +  e  +  e

         Pn2plus(i,89)  = ionsec_n2plus(zenit,zx(i))*Pn2plus(i,81)
         Pelect(i,89)    = Pn2plus(i,89)
         if(n2xinput.ge.1.d-20) then
            Ln2(i,89)   = Pn2plus(i,89)/n2xinput
         else
            Ln2(i,89)   = 1.d-30
         endif


c     R90:  O  +  e  -->  O+  +  e  +  e

         Poplus(i,90)   = ionsec_oplus(zenit,zx(i))*Poplus(i,61)
         Pelect(i,90)    = Poplus(i,90)
         if(o3pxinput.ge.1.d-20) then
            Lo3p(i,90)    = Poplus(i,90)/o3pxinput
         else
            Lo3p(i,90)    = 1.d-30
         endif


c     R91:  CO  +  e  -->  CO+  +  e  +  e
      
         Pcoplus(i,91)  = ionsec_coplus(zenit,zx(i))*Pcoplus(i,77)
         Pelect(i,91)    = Pcoplus(i,91)
         if(coxinput.ge.1.d-20) then
            Lco(i,91)   = Pcoplus(i,91)/coxinput
         else
            Lco(i,91)   = 1.d-30
         endif


c     R92:  CO2  +  e  -->  CO2+  +  e  +  e
      
         Pco2plus(i,92) = ionsec_co2plus(zenit,zx(i))*
     $        Pco2plus(i,51)
         Pelect(i,92)    = Pco2plus(i,92)
         if(co2xinput.ge. 1.d-20) then
            Lco2(i,92)  = Pco2plus(i,92)/co2xinput
         else
            Lco2(i,92)  = 1.d-30
         endif


c     R93:  O2  +  e  -->  O2+  +  e
         Po2plus(i,93)  = ionsec_o2plus(zenit,zx(i))*Po2plus(i,60)
         Pelect(i,93)    = Po2plus(i,93)
         if(o2xinput.ge.1.d-20) then
            Lo2(i,93)   = Po2plus(i,93)/o2xinput
         else
            Lo2(i,93)   = 1.d-30
         endif
         
      endif   !Of chemthermod.eq.3



c     Total production and loss for each species. To save CPU time, we 
c     do not sum over all reactions, but only over those where the species
c     participates

c     CO2: 
c     Prod: R4, R46, R47, R62, R87
c     Loss: R1, R16, R38, R50, R51, R52, R53, R54, R57, R59, R68, R85, R92
      
      Pco2tot(i) = Pco2(i,4) + Pco2(i,46) + Pco2(i,47) + Pco2(i,62)
     $     + Pco2(i,87)
      Lco2tot(i) = Lco2(i,1) + Lco2(i,16) + Lco2(i,38) + Lco2(i,50) +
     $     Lco2(i,51) + Lco2(i,52) + Lco2(i,53) + Lco2(i,54) +
     $     Lco2(i,56) + Lco2(i,57) + Lco2(i,59) + Lco2(i,68) + 
     $     Lco2(i,85) + Lco2(i,92)

c     O2
c     Prod: R3, R5, R7, R10, R11, R13, R22, R23, R24, R25, R26, R34, R41, R42,
c           R54, R64
c     Loss: R2, R12, R17, R21, R32, R46, R60, R93

      Po2tot(i) = Po2(i,3) + Po2(i,5) + Po2(i,7) + Po2(i,10) + 
     $     Po2(i,11) + Po2(i,13) + Po2(i,22) + Po2(i,23) + Po2(i,24) +
     $     Po2(i,25) + Po2(i,26) + Po2(i,34) + Po2(i,41) + Po2(i,42) +
     $     Po2(i,54) + Po2(i,56) + Po2(i,64)
      Lo2tot(i) = Lo2(i,2) + Lo2(i,12) + Lo2(i,17) + Lo2(i,21) +
     $     Lo2(i,32) + Lo2(i,46) + Lo2(i,60) + Lo2(i,93)

      
c     O(3p)
c     Prod: R1, R12, R17, R19, R20, R25, R29, R30, R31, R32, R44, R49, R53,
c           R55, R66, R74, R75, R78
c     Loss: R3, R10, R11, R21, R40, R41, R45, R47, R48, R58, R61, R69, R72,
c           R76, R90

      Po3ptot(i) = Po3p(i,1) + Po3p(i,12) + Po3p(i,17) + Po3p(i,19) +
     $     Po3p(i,20) + Po3p(i,25) + Po3p(i,29) + Po3p(i,30) + 
     $     Po3p(i,31) + Po3p(i,32) + Po3p(i,44) + Po3p(i,49) +
     $     Po3p(i,53) + Po3p(i,55) + Po3p(i,66) + Po3p(i,74) +
     $     Po3p(i,75) + Po3p(i,78)
      Lo3ptot(i) = Lo3p(i,3) + Lo3p(i,10) + Lo3p(i,11) + Lo3p(i,21) +
     $     Lo3p(i,40) + Lo3p(i,41) + Lo3p(i,45) + Lo3p(i,47) +
     $     Lo3p(i,48) + Lo3p(i,58) + Lo3p(i,61) + Lo3p(i,69) +
     $     Lo3p(i,72) + Lo3p(i,76) + Lo3p(i,90)


c     CO
c     Prod: R1, R16, R38, R48, R50, R52, R55, R57, R58, R59, R73
c     Loss: R4, R70, R77, R78, R91
      
      Pcotot(i) = Pco(i,1) + Pco(i,16) + Pco(i,38) + Pco(i,48) +
     $     Pco(i,50) + Pco(i,52) + Pco(i,55) + Pco(i,57) + Pco(i,58) +
     $     Pco(i,59) + Pco(i,73)
      Lcotot(i) = Lco(i,4) + Lco(i,70) + Lco(i,77) + Lco(i,78) +
     $     Lco(i,91)


c     H
c     Prod: R4, R8, R11, R14, R15, R27, R33, R76, R86, R87
c     Loss: R2, R13, R22, R43, R73, R74, R83
      
      Phtot(i) = Ph(i,4) + Ph(i,8) + Ph(i,11) + Ph(i,14) + Ph(i,15) +
     $     Ph(i,27) + Ph(i,33) + Ph(i,76) + Ph(i,86) + Ph(i,87)
      Lhtot(i) = Lh(i,2) + Lh(i,13) + Lh(i,22) + Lh(i,43) + Lh(i,73) +
     $     Lh(i,74) + Lh(i,83)
      

c     OH
c     Prod: R3, R6, R8, R9, R14, R22, R24, R35, R39, R43
c     Loss: R4, R7, R11, R15, R18, R23, R33

      Pohtot(i) = Poh(i,3) + Poh(i,6) + Poh(i,8) + Poh(i,9) +
     $     Poh(i,14) + Poh(i,22) + Poh(i,24) + Poh(i,35) + Poh(i,39) +
     $     Poh(i,43)
      Lohtot(i) = Loh(i,4) + Loh(i,7) + Loh(i,11) + Loh(i,15) +
     $     Loh(i,18) + Loh(i,23) + Loh(i,33)


c     HO2
c     Prod: R2, R18, R23
c     Loss: R3, R5, R7, R13, R24, R35, R39

      Pho2tot(i) = Pho2(i,2) + Pho2(i,18) + Pho2(i,23)
      Lho2tot(i) = Lho2(i,3) + Lho2(i,5) + Lho2(i,7) + Lho2(i,13) +
     $     Lho2(i,24) + Lho2(i,35) + Lho2(i,39)


c     H2
c     Prod: R13
c     Loss: R14, R15, R27, R86

      Ph2tot(i) = Ph2(i,13)
      Lh2tot(i) = Lh2(i,14) + Lh2(i,15) + Lh2(i,27) + Lh2(i,86)
      

c     H2O
c     Prod: R7, R15, R18
c     Loss: R8, R9

      Ph2otot(i) = Ph2o(i,7) + Ph2o(i,15) + Ph2o(i,18)
      Lh2otot(i) = Lh2o(i,8) + Lh2o(i,9)
      

c     O(1d)
c     Prod: R16, R17, R26
c     Loss: R9, R14, R19, R20, R31

      Po1dtot(i) = Po1d(i,16) + Po1d(i,17) + Po1d(i,26)
      Lo1dtot(i) = Lo1d(i,9) + Lo1d(i,14) + Lo1d(i,19) + Lo1d(i,20) +
     $     Lo1d(i,31)

c     H2O2
c     Prod: R5
c     Loss: R6, R18

      Ph2o2tot(i) = Ph2o2(i,5)
      Lh2o2tot(i) = Lh2o2(i,6) + Lh2o2(i,18)
      

      !Only if O3, N or ion chemistry requested
      if(chemthermod.ge.1) then

c     O3
c     Prod: R21
c     Loss: R22, R23, R24, R25, R26, R34, R42

         Po3tot(i) = Po3(i,21)
         Lo3tot(i) = Lo3(i,22) + Lo3(i,23) + Lo3(i,24) + Lo3(i,25) +
     $        Lo3(i,26) + Lo3(i,34) + Lo3(i,42)

      endif


      !Only if N or ion chemistry requested
      if(chemthermod.ge.2) then

c     N
c     Prod: R28, R29, R36, R37, R67, R69, R71, R75, R82, R85
c     Loss: R30, R32, R33, R34, R35, R45, R63, R66, R84, R88
         
         Pntot(i) = Pn(i,28) + Pn(i,29) + Pn(i,36) + Pn(i,37) +Pn(i,67)+
     $        Pn(i,69) + Pn(i,71) + Pn(i,75) + Pn(i,82) + Pn(i,85)
         Lntot(i) = Ln(i,30) + Ln(i,32) + Ln(i,33) + Ln(i,34) +Ln(i,35)+
     $        Ln(i,45) + Ln(i,63) + Ln(i,66) + Ln(i,84) + Ln(i,88)
      

c     NO
c     Prod: R32, R33, R34, R35, R38, R41, R43, R44, R45, R63, R65 
c     Loss: R29, R30, R39, R40, R42, R62, R64, R80

         Pnotot(i) = Pno(i,32) + Pno(i,33) + Pno(i,34) + Pno(i,35) +
     $        Pno(i,38) + Pno(i,41) + Pno(i,43) + Pno(i,44) + Pno(i,45)+
     $        Pno(i,63) + Pno(i,65)
         Lnotot(i) = Lno(i,29) + Lno(i,30) + Lno(i,39) + Lno(i,40) +
     $        Lno(i,42) + Lno(i,62) + Lno(i,64) + Lno(i,80)
      

c     N2
c     Prod: R30, R68, R70, R72
c     Loss: R28, R65, R67, R81, R82, R89

         Pn2tot(i) = Pn2(i,30) + Pn2(i,68) + Pn2(i,70) + Pn2(i,72)
         Ln2tot(i) = Ln2(i,28) + Ln2(i,65) + Ln2(i,67) + Ln2(i,81) + 
     $        Ln2(i,82) + Ln2(i,89)
      

c     N(2d)
c     Prod: R28
c     Loss: R36, R37, R38

         Pn2dtot(i) = Pn2d(i,28)
         Ln2dtot(i) = Ln2d(i,36) + Ln2d(i,37) + Ln2d(i,38)
      

c     NO2
c     Prod: R39, R40, R42
c     Loss: R41, R43, R44

         Pno2tot(i) = Pno2(i,39) + Pno2(i,40) + Pno2(i,42)
         Lno2tot(i) = Lno2(i,41) + Lno2(i,43) + Lno2(i,44)
      
      endif    !Of chemthermod.ge.2
                                ! 
      
      !Only if ion chemistry requested
      if(chemthermod.eq.3) then

c     CO2+
c     Prod: R51,R57, R68, R85, R92
c     Loss: R46, R47, R48, R55, R62, R63, R86

         Pco2plustot(i) = Pco2plus(i,51) + Pco2plus(i,57) + 
     $        Pco2plus(i,68) + Pco2plus(i,85) + Pco2plus(i,92)
         Lco2plustot(i) = Lco2plus(i,46) + Lco2plus(i,47) + 
     $        Lco2plus(i,48) + Lco2plus(i,55) + Lco2plus(i,62) + 
     $        Lco2plus(i,63) + Lco2plus(i,86)
      

c     O+
c     Prod: R47, R52, R58, R61, R72, R76, R90
c     Loss: 50,67,74

         Poplustot(i) = Poplus(i,47) + Poplus(i,52) + Poplus(i,58) +
     $        Poplus(i,61) + Poplus(i,72) + Poplus(i,76) + Poplus(i,90)
         Loplustot(i) = Loplus(i,50) + Loplus(i,56) + Loplus(i,67) + 
     $        Loplus(i,74)

c     O2+
c     Prod: R46, R48, R50, R60, R93
c     Loss: R49, R64, R65, R66

         Po2plustot(i) = Po2plus(i,46) + Po2plus(i,48) + Po2plus(i,50) +
     $        Po2plus(i,60) + Po2plus(i,93)
         Lo2plustot(i) = Lo2plus(i,49) + Lo2plus(i,64) + Lo2plus(i,65) +
     $        Lo2plus(i,66)


c     CO+
c     Prod: R53, R59, R63, R70, R77, R91
c     Loss: R57, R58, R73
      
         Pcoplustot(i) = Pcoplus(i,53) + Pcoplus(i,56) + Pcoplus(i,59) + 
     $        Pcoplus(i,63) + Pcoplus(i,70) + Pcoplus(i,77) + 
     $        Pcoplus(i,91)
         Lcoplustot(i) = Lcoplus(i,57) + Lcoplus(i,58) + Lcoplus(i,73)
      

c     C+
c     Prod: R54, R78
c     Loss: R59

         Pcplustot(i) = Pcplus(i,54) + Pcplus(i,78)
         Lcplustot(i) = Lcplus(i,59)
      

c     N+
c     Prod: R82, R84, R88
c     Loss: R85

         Pnplustot(i) = Pnplus(i,82) + Pnplus(i,84) + Pnplus(i,88)
         Lnplustot(i) = Lnplus(i,85)
      

c     N2+
c     Prod: R81, R89
c     Loss: R68, R69, R70, R71, R72

         Pn2plustot(i) = Pn2plus(i,81) + Pn2plus(i,89)
         Ln2plustot(i) = Ln2plus(i,68) + Ln2plus(i,69) + Ln2plus(i,70) +
     $        Ln2plus(i,71) + Ln2plus(i,72)
      

c     NO+
c     Prod: R62, R64, R65, R66, R67, R69, R80
c     Loss: R75

         Pnoplustot(i) = Pnoplus(i,62) + Pnoplus(i,64) + Pnoplus(i,65) +
     $        Pnoplus(i,66) + Pnoplus(i,67) + Pnoplus(i,69) + 
     $        Pnoplus(i,80)
         Lnoplustot(i) = Lnoplus(i,75)
      

c     H+
c     Prod: R73, R74, R83
c     Loss: R76

         Phplustot(i) = Phplus(i,73) + Phplus(i,74) + Phplus(i,83)
         Lhplustot(i) = Lhplus(i,76)
      

c     HCO2+
c     Prod: R86
c     Loss: R87

         Phco2plustot(i) = Phco2plus(i,86)
         Lhco2plustot(i) = Lhco2plus(i,87)

c     e-
c     Prod: R51, R52, R53, R54, R60, R61, R77, R78, R80, R81, R82, R83, R84, 
c           R88, R89, R90, R91, R92, R93
c     Loss: R49, R55, R71, R75, R87

         Pelecttot(i) = Pelect(i,51) + Pelect(i,52) + Pelect(i,53) + 
     $        Pelect(i,54) + Pelect(i,60) + Pelect(i,61) + Pelect(i,77)+
     $        Pelect(i,78) + Pelect(i,80) + Pelect(i,81) + Pelect(i,82)+
     $        Pelect(i,83) + Pelect(i,84) + Pelect(i,88) + Pelect(i,89)+
     $        Pelect(i,90) + Pelect(i,91) + Pelect(i,92) + Pelect(i,93)
         Lelecttot(i) = Lelect(i,49) + Lelect(i,55) + Lelect(i,71) +
     $        Lelect(i,75) + Lelect(i,87)
      
      endif    !Of chemthermod.eq.3

      return
c END
      end




c**********************************************************************
c**********************************************************************

      subroutine EF_oscilacion
     &               (ig,i,nlayer,paso,chemthermod,zenit, zx,
     &                 jdistot8, jdistot8_b,jion8, 
     &                 tminaux, 
     &                           co2xoutput,     co2xinput,
     $                           o2xoutput,      o2xinput,
     $                           o3pxoutput,     o3pxinput,
     $                           coxoutput,      coxinput,
     $                           h2xoutput,      h2xinput,
     $                           h2oxoutput,     h2oxinput,
     $                           h2o2xoutput,    h2o2xinput,
     $                           o3xoutput,      o3xinput,
     $                           nxoutput,       nxinput,
     $                           noxoutput,      noxinput,
     $                           n2xoutput,      n2xinput,
     &                           o1dxoutput,     o1dxinput, 
     &                           ohxoutput,      ohxinput, 
     &                           ho2xoutput,     ho2xinput, 
     &                           hxoutput,       hxinput, 
     &                           n2dxoutput,     n2dxinput, 
     &                           no2xoutput,     no2xinput, 
     &                         co2plusxoutput, co2plusxinput, 
     &                         o2plusxoutput,  o2plusxinput, 
     &                         coplusxoutput,  coplusxinput, 
     &                         oplusxoutput,   oplusxinput, 
     &                         cplusxoutput,   cplusxinput, 
     &                         noplusxoutput,  noplusxinput, 
     &                         n2plusxoutput,  n2plusxinput, 
     &                         hplusxoutput,   hplusxinput, 
     &                         nplusxoutput,   nplusxinput, 
     $                         hco2plusxoutput,hco2plusxinput,
     &                         electxoutput,   electxinput, 
     &                           electxoutput_timemarching )

  
c Calculates the concentrations of the fast species in PE. Includes a
c procedure to avoid oscillations

c

c         2009      FGG       Adaptation to GCM
c     jul 2008      MA        Version en subrutina
c     nov 2007      MA        Original. Evita oscilacion en EF
c**********************************************************************

      use iono_h
      use param_v4_h, only: nabs,
     .  ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,ch13,ch14,ch15,ch18,
     .  ch19,ch20,ch21,ch22,ch23,ch24,ch30,ch31,ch32,ch33,ch34,
     .  ch35,ch36,ch37,ch38,ch39,ch40,ch41,ch42,ch43,ch45,
     .  ch46,ch47,ch48,ch49,ch50,ch55,ch56,ch57,ch58,ch59,ch62,
     .  ch63,ch64,ch65,ch66,ch67,ch68,ch69,ch70,ch71,
     .  ch72,ch73,ch74,ch75,ch76,ch85,ch86,ch87


      implicit none

c     arguments
      
      integer   ig,nlayer
      integer   i         ! I. Layer                              
      integer   paso      ! I. paso temporal del timemarching, 1,numpasos
      integer   chemthermod
      real*8    tminaux   ! I.
      real      zx(nlayer)
      real      zenit
      real*8    jdistot8(nabs,nlayer)                  ! I.
      real*8    jdistot8_b(nabs,nlayer)                ! I.
      real*8    jion8(nabs,nlayer,4)

      real*8    co2xoutput,o2xoutput,o3pxoutput,coxoutput,h2xoutput
      real*8    h2o2xoutput,o3xoutput,h2oxoutput
      real*8    nxoutput,noxoutput,n2xoutput
      real*8    ho2xoutput,      hxoutput,       ohxoutput    ! O.
      real*8    o1dxoutput,      n2dxoutput,     no2xoutput   ! O.

      real*8    co2xinput,o2xinput,o3pxinput,coxinput,h2xinput
      real*8    h2o2xinput,o3xinput,h2oxinput
      real*8    nxinput,noxinput,n2xinput
      real*8    ho2xinput,       hxinput,        ohxinput     ! I.
      real*8    o1dxinput,       n2dxinput,      no2xinput    ! I. 

      real*8    co2plusxoutput,  coplusxoutput,  oplusxoutput   ! O.
      real*8    o2plusxoutput,   cplusxoutput,   noplusxoutput  ! O.
      real*8    n2plusxoutput,   hplusxoutput                   ! O.
      real*8    electxoutput,    nplusxoutput, hco2plusxoutput  ! O.
      real*8    co2plusxinput,   coplusxinput,   oplusxinput    ! I.
      real*8    o2plusxinput,    cplusxinput,    noplusxinput   ! I.
      real*8    n2plusxinput,    hplusxinput                    ! I.
      real*8    electxinput,     nplusxinput, hco2plusxinput    ! I.

      real*8    electxoutput_timemarching           ! I. 

   

c     local variables

      integer   npasitos
      parameter (npasitos=20)
      real*8    electxoutput_neutr

      real*8    ho2xoutput2, hxoutput2, ohxoutput2
      real*8    o1dxoutput2, n2dxoutput2, no2xoutput2
      real*8    co2plusxoutput2, coplusxoutput2, oplusxoutput2
      real*8    o2plusxoutput2, hplusxoutput2, nplusxoutput2
      real*8    cplusxoutput2, noplusxoutput2, n2plusxoutput2
      real*8    hco2plusxoutput2
      !real*8    electxoutput2

      integer   correc_oscil,  ip 
      real*8    avg_pares, avg_impar, dispersion
      real*8    dif_impar, dif_pares , dif_pares_impar
      real*8    ho2xpares(3), hxpares(3), ohxpares(3)
      real*8    o1dxpares(3), n2dxpares(3), no2xpares(3)
      real*8    co2plusxpares(3), coplusxpares(3), oplusxpares(3)
      real*8    o2plusxpares(3), hplusxpares(3), nplusxpares(3)
      real*8    cplusxpares(3), noplusxpares(3), n2plusxpares(3)
      real*8    hco2plusxpares(3)
      real*8    ho2ximpar(3), hximpar(3), ohximpar(3)
      real*8    o1dximpar(3), n2dximpar(3), no2ximpar(3)
      real*8    co2plusximpar(3), coplusximpar(3), oplusximpar(3)
      real*8    o2plusximpar(3), hplusximpar(3), nplusximpar(3)
      real*8    cplusximpar(3), noplusximpar(3), n2plusximpar(3)
      real*8    hco2plusximpar(3)
      real*8    auxp,auxl,ca, cb,cc, cd1,cd2,ce
      real*8    IonMostAbundant

      real*8    cocimin
      parameter (cocimin=1.d-30)

      integer   cociopt
      parameter (cociopt=1)

c external functions
c
      external cociente
      real*8   cociente

      external ionsec_nplus
      real*8   ionsec_nplus

      external ionsec_n2plus
      real*8   ionsec_n2plus

      external ionsec_oplus
      real*8   ionsec_oplus

      external ionsec_coplus
      real*8   ionsec_coplus

      external ionsec_co2plus
      real*8   ionsec_co2plus

      external ionsec_o2plus
      real*8   ionsec_o2plus

      external avg
      real*8   avg

      external dif
      real*8   dif

      real*8 log1
      real*8 log2
      real*8 log3

ccccccccccccccc CODE STARTS 

           !!! Starts iteration to avoid oscilations

      correc_oscil=0

      do ip = 1, npasitos 

         if (ip.eq.1) then 
            if (o1d_eq(i).eq.'Y') o1dxoutput = o1dxinput 
            if (oh_eq(i).eq.'Y') ohxoutput = ohxinput 
            if (ho2_eq(i).eq.'Y') ho2xoutput = ho2xinput 
            if (h_eq(i).eq.'Y') hxoutput = hxinput 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dxoutput = n2dxinput 
               if (no2_eq(i).eq.'Y') no2xoutput = no2xinput
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.ge.3) then
               if (n2plus_eq(i).eq.'Y') n2plusxoutput=n2plusxinput 
               if (cplus_eq(i).eq.'Y') cplusxoutput=cplusxinput 
               if (coplus_eq(i).eq.'Y') coplusxoutput=coplusxinput 
               if (oplus_eq(i).eq.'Y') oplusxoutput=oplusxinput 
               if (hplus_eq(i).eq.'Y') hplusxoutput=hplusxinput 
               if (co2plus_eq(i).eq.'Y') co2plusxoutput=co2plusxinput 
               if (noplus_eq(i).eq.'Y') noplusxoutput=noplusxinput 
               if (o2plus_eq(i).eq.'Y') o2plusxoutput=o2plusxinput 
               if (nplus_eq(i).eq.'Y') nplusxoutput=nplusxinput 
               if (hco2plus_eq(i).eq.'Y') hco2plusxoutput=hco2plusxinput
               electxoutput = electxinput 
            endif
         endif

              
         ! 6 neutral , O1D, OH, HO2, H, N2D, NO2

         !O1D
         if (o1d_eq(i).eq.'Y') then
            auxp = jdistot8_b(1,i) * co2xoutput 
     &           + jdistot8_b(2,i) * o2xoutput
     &           + jdistot8_b(7,i) * o3xoutput
            auxl = ch9 * h2oxoutput
     &           + ch14 * h2xoutput
     &           + ch19 * co2xoutput
     &           + ch20 * o2xoutput
     &           + ch31 * n2xoutput
            o1dxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
         end if
           
           !OH
         if (oh_eq(i).eq.'Y') then
            auxp = ch3 * o3pxoutput * ho2xoutput
     &           + 2.d0 * jdistot8(6,i) * h2o2xoutput
     &           + jdistot8(4,i) * h2oxoutput
     &           + 2.d0 * ch9 * o1dxoutput * h2oxoutput 
     &           + ch14 * o1dxoutput * h2xoutput
     &           + ch22 * o3xoutput * hxoutput
     &           + ch24 * o3xoutput * ho2xoutput
     &           + ch35 * nxoutput * ho2xoutput
     &           + ch39 * noxoutput * ho2xoutput
     &           + ch43 * no2xoutput * hxoutput
            auxl = ch4 * coxoutput
     &           + ch7 * ho2xoutput
     &           + ch11 * o3pxoutput
     &           + ch15 * h2xoutput
     &           + ch18 * h2o2xoutput
     &           + ch23 * o3xoutput
     &           + ch33 * nxoutput
            ohxoutput2 = cociente(auxp,auxl,cocimin,cociopt)

         end if

			
         !HO2
         if (ho2_eq(i).eq.'Y') then
            auxp = ch2 * hxoutput * o2xoutput * co2xoutput
     &           + ch18 * ohxoutput * h2o2xoutput
     &           + ch23 * ohxoutput * o3xoutput
            auxl = ch3 * o3pxoutput
     &           + 2.d0 * ch5 * ho2xoutput 
     &           + ch7 * ohxoutput
     &           + ch13 * hxoutput
     &           + ch24 * o3xoutput
     &           + ch35 * nxoutput
     &           + ch39 * noxoutput
            ho2xoutput2 = cociente(auxp,auxl,cocimin,cociopt)
         end if

         !H
         if (h_eq(i).eq.'Y') then
            auxp = ch4 * coxoutput * ohxoutput
     &           + jdistot8(4,i) * h2oxoutput
     &           + ch11 * o3pxoutput * ohxoutput
     &           + ch14 * o1dxoutput * h2xoutput
     &           + ch15 * ohxoutput * h2xoutput
     &           + 2.d0 * jdistot8(5,i) * h2xoutput
     &           +  ch33 * nxoutput * ohxoutput
     &           + ch76 *  hplusxoutput * o3pxoutput
     &           + hxoutput * jion8(12,i,1) 
     $           + ch86 * h2xoutput * co2plusxoutput
     $           + ch87 * hco2plusxoutput * electxoutput
            auxl = ch2 * o2xoutput * co2xoutput
     &           + ch13 * ho2xoutput
     &           + ch22 * o3xoutput
     &           + ch43 * no2xoutput
     &           + ch73 * coplusxoutput
     &           + ch74 * oplusxoutput
     &           + jion8(12,i,1)
            hxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
         end if
			
         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
         !N2D
            if (n2d_eq(i).eq.'Y') then
               auxp = jdistot8(8,i) * n2xoutput
               auxl = ch36 * o3pxoutput
     &              + ch37 * n2xoutput
     &              + ch38 * co2xoutput
               n2dxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            end if

         !NO2
            if (no2_eq(i).eq.'Y') then
               auxp = ch39 * noxoutput * ho2xoutput
     &              + ch40 * o3pxoutput * noxoutput * co2xoutput
     &              + ch42 * noxoutput * o3xoutput
               auxl = ch41 * o3pxoutput
     &              + ch43 * hxoutput
     &              + jdistot8(13,i)
               no2xoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            end if
            
         endif    !Of chemthermod.ge.2
         
         ! 9 ions 

         !Only if ion chemistry requested
         if(chemthermod.eq.3) then
         ! N2+
            if (n2plus_eq(i).eq.'Y') then
               auxp = jion8(8,i,1)*n2xoutput*
     $              (1.d0+ionsec_n2plus(zenit,zx(i)))
               auxl = ch68*co2xoutput
     &              + ch69*o3pxoutput 
     &              + ch70*coxoutput
     $              + ch71*electxoutput
     &              + ch72*o3pxoutput
               n2plusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            endif

         ! C+
            if (cplus_eq(i).eq.'Y') then
               auxp = jion8(1,i,4)*co2xoutput 
     &              + jion8(11,i,2)*coxoutput
               auxl = ch59*co2xoutput
               cplusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            end if
            
         ! CO+
            if (coplus_eq(i).eq.'Y') then
               auxp = jion8(1,i,3)*co2xoutput 
     $              + ch59*cplusxoutput *co2xoutput 
     $              + ch63*co2plusxoutput*nxoutput 
     $              + ch70*n2plusxoutput*coxoutput 
     $              + jion8(11,i,1)*coxoutput*
     $              (1.d0+ionsec_coplus(zenit,zx(i))) 
               auxl = ch57*co2xoutput 
     &              + ch58*o3pxoutput 
     &              + ch73*hxoutput
               coplusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            end if

         ! CO2+
            if (co2plus_eq(i).eq.'Y') then
               auxp = jion8(1,i,1)*co2xoutput*
     $              (1.d0+ionsec_co2plus(zenit,zx(i))) 
     &              + ch57*coplusxoutput*co2xoutput
     $              + ch68*n2plusxoutput*co2xoutput 
     &              + ch85*nplusxoutput*co2xoutput
               auxl = ch46*o2xoutput 
     &              + ch47*o3pxoutput 
     &              + ch55*electxoutput
     $              + ch62*noxoutput 
     &              + ch63*nxoutput 
     &              + ch48*o3pxoutput
     $              + ch86*h2xoutput
               co2plusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            end if

         !
         !   [O+] [H+] are linked: 2 equations with 2 unknowns 
         !
         !      [H+] = (ca + cb [O+])/cd1
         !      [O+] = (cc + ce [H+])/cd2
         !        
            ca  = ch73 * coplusxoutput * hxoutput 
     &           + jion8(12,i,1)*hxoutput 
            cb  = ch74 * hxoutput
            cd1 = ch76 * o3pxoutput              
            cc  = ch47*co2plusxoutput*o3pxoutput 
     &           + jion8(1,i,2)*co2xoutput
     $           + jion8(3,i,1)*o3pxoutput*
     $           (1.d0+ionsec_oplus(zenit,zx(i))) 
     &           + ch58*coplusxoutput*o3pxoutput
     $           + ch72*n2plusxoutput*o3pxoutput 
            ce  = ch76 * o3pxoutput
            cd2 = ch50*co2xoutput + ch67*n2xoutput + ch74*hxoutput
            
         ! O+
            if (oplus_eq(i).eq.'Y') then
               auxp = cc*cd1 + ce*ca
               auxl = cd1*cd2 -ce*cb
               oplusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            end if
         
         ! H+
            if (hplus_eq(i).eq.'Y') then 
               auxp = ca + cb * oplusxoutput
               auxl = cd1
               hplusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            endif
           
         ! O2+
            if (o2plus_eq(i).eq.'Y') then
               auxp = ch46*co2plusxoutput*o2xoutput  
     $              + ch48*co2plusxoutput*o3pxoutput
     $              + ch50*oplusxoutput*co2xoutput 
     $              + jion8(2,i,1)*o2xoutput*
     $              (1.d0+ionsec_o2plus(zenit,zx(i)))
               auxl = ch49*electxoutput + ch64*noxoutput 
     $              + ch65*n2xoutput + ch66*nxoutput
               o2plusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            endif

         ! NO+
            if(noplus_eq(i).eq.'Y') then
               auxp = ch62*coplusxoutput*noxoutput 
     $              + ch64*o2plusxoutput*noxoutput
     $              + ch65*o2plusxoutput*n2xoutput 
     $              + ch66*o2plusxoutput*nxoutput
     $              + ch67*oplusxoutput*n2xoutput 
     $              + ch69*n2plusxoutput*o3pxoutput
     $              + jion8(10,i,1)*noxoutput 
               auxl = ch75*electxoutput
               noplusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            endif

         ! N+
            if(nplus_eq(i).eq.'Y') then
               auxp = jion8(8,i,2)*n2xoutput 
     &              + jion8(9,i,1)*nxoutput*
     $              (1.d0+ionsec_nplus(zenit,zx(i))) 
               auxl = ch85*co2xoutput
               nplusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            endif

         ! HCO2+
            if(hco2plus_eq(i).eq.'Y') then
               auxp = ch86*h2xoutput*co2plusxoutput
               auxl = ch87*electxoutput
               hco2plusxoutput2 = cociente(auxp,auxl,cocimin,cociopt)
            endif

         endif    !Of chemthermod.eq.3

         ! Detection of oscilations and elimination
           
         if (ip.eq.4 .or. ip.eq.14) then ! ***pares(1) 
            
            if (o1d_eq(i).eq.'Y') o1dxpares(1)=o1dxoutput2
            if (oh_eq(i).eq.'Y')  ohxpares(1)=ohxoutput2 
            if (ho2_eq(i).eq.'Y') ho2xpares(1)=ho2xoutput2 
            if (h_eq(i).eq.'Y')   hxpares(1)=hxoutput2 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dxpares(1)=n2dxoutput2 
               if (no2_eq(i).eq.'Y') no2xpares(1)=no2xoutput2
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (n2plus_eq(i).eq.'Y')  n2plusxpares(1)=n2plusxoutput2 
               if (cplus_eq(i).eq.'Y')   cplusxpares(1)=cplusxoutput2 
               if (coplus_eq(i).eq.'Y')  coplusxpares(1)=coplusxoutput2 
               if (oplus_eq(i).eq.'Y')   oplusxpares(1)=oplusxoutput2 
               if (hplus_eq(i).eq.'Y')   hplusxpares(1)=hplusxoutput2 
               if (co2plus_eq(i).eq.'Y')co2plusxpares(1)=co2plusxoutput2 
               if (noplus_eq(i).eq.'Y')  noplusxpares(1)=noplusxoutput2 
               if (o2plus_eq(i).eq.'Y')  o2plusxpares(1)=o2plusxoutput2 
               if (nplus_eq(i).eq.'Y')   nplusxpares(1)=nplusxoutput2
               if (hco2plus_eq(i).eq.'Y') 
     $              hco2plusxpares(1)=hco2plusxoutput2
            endif

         elseif (ip.eq.6 .or. ip.eq.16) then ! ***pares(2)

            if (o1d_eq(i).eq.'Y') o1dxpares(2)=o1dxoutput2
            if (oh_eq(i).eq.'Y')  ohxpares(2)=ohxoutput2 
            if (ho2_eq(i).eq.'Y') ho2xpares(2)=ho2xoutput2 
            if (h_eq(i).eq.'Y')   hxpares(2)=hxoutput2 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dxpares(2)=n2dxoutput2 
               if (no2_eq(i).eq.'Y') no2xpares(2)=no2xoutput2
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (n2plus_eq(i).eq.'Y')  n2plusxpares(2)=n2plusxoutput2 
               if (cplus_eq(i).eq.'Y')   cplusxpares(2)=cplusxoutput2 
               if (coplus_eq(i).eq.'Y')  coplusxpares(2)=coplusxoutput2 
               if (oplus_eq(i).eq.'Y')   oplusxpares(2)=oplusxoutput2 
               if (hplus_eq(i).eq.'Y')   hplusxpares(2)=hplusxoutput2 
               if (co2plus_eq(i).eq.'Y')co2plusxpares(2)=co2plusxoutput2 
               if (noplus_eq(i).eq.'Y')  noplusxpares(2)=noplusxoutput2 
               if (o2plus_eq(i).eq.'Y')  o2plusxpares(2)=o2plusxoutput2 
               if (nplus_eq(i).eq.'Y')   nplusxpares(2)=nplusxoutput2 
               if (hco2plus_eq(i).eq.'Y')
     $              hco2plusxpares(2)=hco2plusxoutput2
            endif
            
         elseif (ip.eq.8 .or. ip.eq.18) then ! ***pares(3)
            
            if (o1d_eq(i).eq.'Y') o1dxpares(3)=o1dxoutput2
            if (oh_eq(i).eq.'Y')  ohxpares(3)=ohxoutput2 
            if (ho2_eq(i).eq.'Y') ho2xpares(3)=ho2xoutput2 
            if (h_eq(i).eq.'Y')   hxpares(3)=hxoutput2 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dxpares(3)=n2dxoutput2 
               if (no2_eq(i).eq.'Y') no2xpares(3)=no2xoutput2
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (n2plus_eq(i).eq.'Y')  n2plusxpares(3)=n2plusxoutput2 
               if (cplus_eq(i).eq.'Y')   cplusxpares(3)=cplusxoutput2 
               if (coplus_eq(i).eq.'Y')  coplusxpares(3)=coplusxoutput2 
               if (oplus_eq(i).eq.'Y')   oplusxpares(3)=oplusxoutput2 
               if (hplus_eq(i).eq.'Y')   hplusxpares(3)=hplusxoutput2 
               if (co2plus_eq(i).eq.'Y')co2plusxpares(3)=co2plusxoutput2 
               if (noplus_eq(i).eq.'Y')  noplusxpares(3)=noplusxoutput2 
               if (o2plus_eq(i).eq.'Y')  o2plusxpares(3)=o2plusxoutput2 
               if (nplus_eq(i).eq.'Y')   nplusxpares(3)=nplusxoutput2 
               if (hco2plus_eq(i).eq.'Y')
     $              hco2plusxpares(3)=hco2plusxoutput2
            endif

         elseif (ip.eq.5 .or. ip.eq.15) then ! ***impar(1) 
            
            if (o1d_eq(i).eq.'Y') o1dximpar(1)=o1dxoutput2
            if (oh_eq(i).eq.'Y')  ohximpar(1)=ohxoutput2 
            if (ho2_eq(i).eq.'Y') ho2ximpar(1)=ho2xoutput2 
            if (h_eq(i).eq.'Y')   hximpar(1)=hxoutput2 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dximpar(1)=n2dxoutput2 
               if (no2_eq(i).eq.'Y') no2ximpar(1)=no2xoutput2
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (n2plus_eq(i).eq.'Y')  n2plusximpar(1)=n2plusxoutput2 
               if (cplus_eq(i).eq.'Y')   cplusximpar(1)=cplusxoutput2 
               if (coplus_eq(i).eq.'Y')  coplusximpar(1)=coplusxoutput2 
               if (oplus_eq(i).eq.'Y')   oplusximpar(1)=oplusxoutput2 
               if (hplus_eq(i).eq.'Y')   hplusximpar(1)=hplusxoutput2 
               if (co2plus_eq(i).eq.'Y')co2plusximpar(1)=co2plusxoutput2 
               if (noplus_eq(i).eq.'Y')  noplusximpar(1)=noplusxoutput2 
               if (o2plus_eq(i).eq.'Y')  o2plusximpar(1)=o2plusxoutput2 
               if (nplus_eq(i).eq.'Y')   nplusximpar(1)=nplusxoutput2 
               if (hco2plus_eq(i).eq.'Y')
     $              hco2plusximpar(1)=hco2plusxoutput2
            endif
            
         elseif (ip.eq.7 .or. ip.eq.17) then ! ***impar(2)
            
            if (o1d_eq(i).eq.'Y') o1dximpar(2)=o1dxoutput2
            if (oh_eq(i).eq.'Y')  ohximpar(2)=ohxoutput2 
            if (ho2_eq(i).eq.'Y') ho2ximpar(2)=ho2xoutput2 
            if (h_eq(i).eq.'Y')   hximpar(2)=hxoutput2 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dximpar(2)=n2dxoutput2 
               if (no2_eq(i).eq.'Y') no2ximpar(2)=no2xoutput2
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (n2plus_eq(i).eq.'Y')  n2plusximpar(2)=n2plusxoutput2 
               if (cplus_eq(i).eq.'Y')   cplusximpar(2)=cplusxoutput2 
               if (coplus_eq(i).eq.'Y')  coplusximpar(2)=coplusxoutput2 
               if (oplus_eq(i).eq.'Y')   oplusximpar(2)=oplusxoutput2 
               if (hplus_eq(i).eq.'Y')   hplusximpar(2)=hplusxoutput2 
               if (co2plus_eq(i).eq.'Y')co2plusximpar(2)=co2plusxoutput2 
               if (noplus_eq(i).eq.'Y')  noplusximpar(2)=noplusxoutput2 
               if (o2plus_eq(i).eq.'Y')  o2plusximpar(2)=o2plusxoutput2 
               if (nplus_eq(i).eq.'Y')   nplusximpar(2)=nplusxoutput2 
               if (hco2plus_eq(i).eq.'Y')
     $              hco2plusximpar(2)=hco2plusxoutput2
            endif
            
         elseif (ip.eq.9 .or. ip.eq.19) then ! ***impar(3)
            
            if (o1d_eq(i).eq.'Y') o1dximpar(3)=o1dxoutput2
            if (oh_eq(i).eq.'Y')  ohximpar(3)=ohxoutput2 
            if (ho2_eq(i).eq.'Y') ho2ximpar(3)=ho2xoutput2 
            if (h_eq(i).eq.'Y')   hximpar(3)=hxoutput2 
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') n2dximpar(3)=n2dxoutput2 
               if (no2_eq(i).eq.'Y') no2ximpar(3)=no2xoutput2
            endif
                                ! 
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (n2plus_eq(i).eq.'Y')  n2plusximpar(3)=n2plusxoutput2 
               if (cplus_eq(i).eq.'Y')   cplusximpar(3)=cplusxoutput2 
               if (coplus_eq(i).eq.'Y')  coplusximpar(3)=coplusxoutput2 
               if (oplus_eq(i).eq.'Y')   oplusximpar(3)=oplusxoutput2 
               if (hplus_eq(i).eq.'Y')   hplusximpar(3)=hplusxoutput2 
               if (co2plus_eq(i).eq.'Y')co2plusximpar(3)=co2plusxoutput2 
               if (noplus_eq(i).eq.'Y')  noplusximpar(3)=noplusxoutput2 
               if (o2plus_eq(i).eq.'Y')  o2plusximpar(3)=o2plusxoutput2 
               if (nplus_eq(i).eq.'Y')   nplusximpar(3)=nplusxoutput2 
               if (hco2plus_eq(i).eq.'Y')
     $              hco2plusximpar(3)=hco2plusxoutput2
            endif
            
            
            if (o1d_eq(i).eq.'Y') then
               log1 = log10(o1dxpares(1))
               log2 = log10(o1dxpares(2))
               log3 = log10(o1dxpares(3))
               avg_pares = avg(log1,log2,log3)
               dif_pares = dif(log1,log2,log3,avg_pares)
               log1 = log10(o1dximpar(1))
               log2 = log10(o1dximpar(2))
               log3 = log10(o1dximpar(3))
               avg_impar = avg(log1,log2,log3)
               dif_impar = dif(log1,log2,log3,avg_impar)
               dispersion = dif_pares + dif_impar
               dif_pares_impar = abs(avg_pares-avg_impar)
               if (dif_pares_impar .gt. 5.d0*dispersion) then
                  correc_oscil=correc_oscil+1
                  o1dxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
               endif
            endif
            
            if (oh_eq(i).eq.'Y') then
               log1 = log10(ohxpares(1))
               log2 = log10(ohxpares(2))
               log3 = log10(ohxpares(3))
               avg_pares = avg(log1,log2,log3)
               dif_pares = dif(log1,log2,log3,avg_pares)
               log1 = log10(ohximpar(1))
               log2 = log10(ohximpar(2))
               log3 = log10(ohximpar(3))
               avg_impar = avg(log1,log2,log3)
               dif_impar = dif(log1,log2,log3,avg_impar)
               dispersion = dif_pares + dif_impar
               dif_pares_impar = abs(avg_pares-avg_impar)
               if (dif_pares_impar .gt. 5.*dispersion) then
                  correc_oscil=correc_oscil+1
                  ohxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
               endif
               
            endif
            
            if (ho2_eq(i).eq.'Y') then
               log1 = log10(ho2xpares(1))
               log2 = log10(ho2xpares(2))
               log3 = log10(ho2xpares(3))
               avg_pares = avg(log1,log2,log3)
               dif_pares = dif(log1,log2,log3,avg_pares)
               log1 = log10(ho2ximpar(1))
               log2 = log10(ho2ximpar(2))
               log3 = log10(ho2ximpar(3))
               avg_impar = avg(log1,log2,log3)
               dif_impar = dif(log1,log2,log3,avg_impar)
               dispersion = dif_pares + dif_impar
               dif_pares_impar = abs(avg_pares-avg_impar)
               if (dif_pares_impar .gt. 5.*dispersion) then
                  correc_oscil=correc_oscil+1
                  ho2xoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
               endif
            endif
            
            if (h_eq(i).eq.'Y') then
               log1 = log10(hxpares(1))
               log2 = log10(hxpares(2))
               log3 = log10(hxpares(3))
               avg_pares = avg(log1,log2,log3)
               dif_pares = dif(log1,log2,log3,avg_pares)
               log1 = log10(hximpar(1))
               log2 = log10(hximpar(2))
               log3 = log10(hximpar(3))
               avg_impar = avg(log1,log2,log3)
               dif_impar = dif(log1,log2,log3,avg_impar)
               dispersion = dif_pares + dif_impar
               dif_pares_impar = abs(avg_pares-avg_impar)
               if (dif_pares_impar .gt. 5.*dispersion) then
                  correc_oscil=correc_oscil+1
                  hxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
               endif
            endif
            
            !Only if N or ion chemistry requested
            if(chemthermod.ge.2) then
               if (n2d_eq(i).eq.'Y') then
                  log1 = log10(n2dxpares(1))
                  log2 = log10(n2dxpares(2))
                  log3 = log10(n2dxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(n2dximpar(1))
                  log2 = log10(n2dximpar(2))
                  log3 = log10(n2dximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     n2dxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
            
               if (no2_eq(i).eq.'Y') then
                  log1 = log10(no2xpares(1))
                  log2 = log10(no2xpares(2))
                  log3 = log10(no2xpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(no2ximpar(1))
                  log2 = log10(no2ximpar(2))
                  log3 = log10(no2ximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     no2xoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
            endif    !Of chemthermod.ge.2

                                ! IONS
            !Only if ion chemistry requested
            if(chemthermod.eq.3) then
               if (cplus_eq(i).eq.'Y') then
                  log1 = log10(cplusxpares(1))
                  log2 = log10(cplusxpares(2))
                  log3 = log10(cplusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(cplusximpar(1))
                  log2 = log10(cplusximpar(2))
                  log3 = log10(cplusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     cplusxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (coplus_eq(i).eq.'Y') then
                  log1 = log10(coplusxpares(1))
                  log2 = log10(coplusxpares(2))
                  log3 = log10(coplusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(coplusximpar(1))
                  log2 = log10(coplusximpar(2))
                  log3 = log10(coplusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     coplusxoutput2=10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (oplus_eq(i).eq.'Y') then
                  log1 = log10(oplusxpares(1))
                  log2 = log10(oplusxpares(2))
                  log3 = log10(oplusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(oplusximpar(1))
                  log2 = log10(oplusximpar(2))
                  log3 = log10(oplusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     oplusxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (n2plus_eq(i).eq.'Y') then
                  log1 = log10(n2plusxpares(1))
                  log2 = log10(n2plusxpares(2))
                  log3 = log10(n2plusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(n2plusximpar(1))
                  log2 = log10(n2plusximpar(2))
                  log3 = log10(n2plusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     n2plusxoutput2 =10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (hplus_eq(i).eq.'Y') then
                  log1 = log10(hplusxpares(1))
                  log2 = log10(hplusxpares(2))
                  log3 = log10(hplusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(hplusximpar(1))
                  log2 = log10(hplusximpar(2))
                  log3 = log10(hplusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     hplusxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (co2plus_eq(i).eq.'Y') then
                  log1 = log10(co2plusxpares(1))
                  log2 = log10(co2plusxpares(2))
                  log3 = log10(co2plusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(co2plusximpar(1))
                  log2 = log10(co2plusximpar(2))
                  log3 = log10(co2plusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     co2plusxoutput2=10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (o2plus_eq(i).eq.'Y') then
                  log1 = log10(o2plusxpares(1))
                  log2 = log10(o2plusxpares(2))
                  log3 = log10(o2plusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(o2plusximpar(1))
                  log2 = log10(o2plusximpar(2))
                  log3 = log10(o2plusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     o2plusxoutput2 =10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (noplus_eq(i).eq.'Y') then
                  log1 = log10(noplusxpares(1))
                  log2 = log10(noplusxpares(2))
                  log3 = log10(noplusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(noplusximpar(1))
                  log2 = log10(noplusximpar(2))
                  log3 = log10(noplusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     noplusxoutput2 =10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif
               
               if (nplus_eq(i).eq.'Y') then
                  log1 = log10(nplusxpares(1))
                  log2 = log10(nplusxpares(2))
                  log3 = log10(nplusxpares(3))
                  avg_pares = avg(log1,log2,log3)
                  dif_pares = dif(log1,log2,log3,avg_pares)
                  log1 = log10(nplusximpar(1))
                  log2 = log10(nplusximpar(2))
                  log3 = log10(nplusximpar(3))
                  avg_impar = avg(log1,log2,log3)
                  dif_impar = dif(log1,log2,log3,avg_impar)
                  dispersion = dif_pares + dif_impar
                  dif_pares_impar = abs(avg_pares-avg_impar)
                  if (dif_pares_impar .gt. 5.*dispersion) then
                     correc_oscil=correc_oscil+1
                     nplusxoutput2 = 10.d0**(0.5*(avg_pares+avg_impar))
                  endif
               endif

               if (hco2plus_eq(i).eq.'Y') then
                 log1 = log10(hco2plusxpares(1))
                 log2 = log10(hco2plusxpares(2))
                 log3 = log10(hco2plusxpares(3))
                 avg_pares = avg(log1,log2,log3)
                 dif_pares = dif(log1,log2,log3,avg_pares)
                 log1 = log10(hco2plusximpar(1))
                 log2 = log10(hco2plusximpar(2))
                 log3 = log10(hco2plusximpar(3))
                 avg_impar = avg(log1,log2,log3)
                 dif_impar = dif(log1,log2,log3,avg_impar)
                 dispersion = dif_pares + dif_impar
                 dif_pares_impar = abs(avg_pares-avg_impar)
                 if (dif_pares_impar .gt. 5.*dispersion) then
                    correc_oscil=correc_oscil+1
                    hco2plusxoutput2=10.d0**(0.5*(avg_pares+avg_impar))
                 endif
              endif

            endif   !Of chemthermod.eq.3
            
         endif


          ! Preparation of next step
         
         if (o1d_eq(i).eq.'Y') o1dxoutput = o1dxoutput2
         if (oh_eq(i).eq.'Y') ohxoutput = ohxoutput2 
         if (ho2_eq(i).eq.'Y') ho2xoutput = ho2xoutput2 
         if (h_eq(i).eq.'Y') hxoutput = hxoutput2 
         !Only if N or ion chemistry requested
         if(chemthermod.ge.2) then
            if (n2d_eq(i).eq.'Y') n2dxoutput = n2dxoutput2 
            if (no2_eq(i).eq.'Y') no2xoutput = no2xoutput2
         endif
                                ! 
         !Only if ion chemistry requested
         if(chemthermod.eq.3) then
            if (n2plus_eq(i).eq.'Y') n2plusxoutput=n2plusxoutput2 
            if (cplus_eq(i).eq.'Y') cplusxoutput=cplusxoutput2 
            if (coplus_eq(i).eq.'Y') coplusxoutput=coplusxoutput2 
            if (oplus_eq(i).eq.'Y') oplusxoutput=oplusxoutput2 
            if (hplus_eq(i).eq.'Y') hplusxoutput=hplusxoutput2 
            if (co2plus_eq(i).eq.'Y') co2plusxoutput=co2plusxoutput2 
            if (noplus_eq(i).eq.'Y') noplusxoutput=noplusxoutput2 
            if (o2plus_eq(i).eq.'Y') o2plusxoutput=o2plusxoutput2 
            if (nplus_eq(i).eq.'Y') nplusxoutput=nplusxoutput2
            if (hco2plus_eq(i).eq.'Y') hco2plusxoutput=hco2plusxoutput2
            electxoutput = o2plusxoutput +
     @           co2plusxoutput +
     @           coplusxoutput +
     @           oplusxoutput +
     @           cplusxoutput +
     @           n2plusxoutput +
     @           nplusxoutput +
     @           noplusxoutput +
     @           hplusxoutput +
     $           hco2plusxoutput
            
            electxoutput_neutr = electxoutput 
         
            IonMostAbundant = o2plusxoutput
            IonMostAbundant = max( co2plusxoutput, IonMostAbundant) 
            IonMostAbundant = max( coplusxoutput, IonMostAbundant) 
            IonMostAbundant = max( oplusxoutput, IonMostAbundant) 
            IonMostAbundant = max( cplusxoutput, IonMostAbundant) 
            IonMostAbundant = max( n2plusxoutput, IonMostAbundant) 
            IonMostAbundant = max( noplusxoutput, IonMostAbundant) 
            IonMostAbundant = max( nplusxoutput, IonMostAbundant) 
            IonMostAbundant = max( hplusxoutput, IonMostAbundant) 
            IonMostAbundant = max( hco2plusxoutput, IonMostAbundant)
            IonMostAbundant = IonMostAbundant / electxoutput

         endif   !Of chemthermod.eq.3
         


      enddo
            !!! End of iteration


      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function avg(log1,log2,log3)
        implicit none
        real*8 avg
        real*8 log1,log2,log3
        avg = (log1+log2+log3)*0.333
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function dif(log1,log2,log3,avg)
        implicit none
        real*8 dif
        real*8 avg
        real*8 log1,log2,log3
        dif = (abs(log1-avg) +
     &                abs(log2-avg) +
     &                abs(log3-avg) ) * 0.333
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c**********************************************************************
c***********************************************************************
	function cociente (xnum, xdenom, minimo, option)                      

c	Returns the ratio between XNUM and XDENOM avoiding null values
c       according to the action given by OPTION. Checks that the are not
c       negative values

c       If XNUM = 0 -> cociente=0
c       If XDENOM < minimo, increases XDENOM in a factor=1e20 to avoid
c         overflow, and after the ratio XNUM/XDENOM the value is back to normal
c         if cociente < minimo -> cociente=minimo
c       If XDENOM = 0 then :  
c          option = 0  .... warning message and stop
c          option = 1  ....    "       "    and cociente=minimo

c       MALV    Jul-08    Original
c***********************************************************************
                                                
	implicit none             

! Arguments 
        real*8  cociente
        real*8  xnum
        real*8  xdenom
        real*8  minimo
        integer option

! Local variables 
        real*8  factor

!!!!!!! Program starts

	if (xnum .lt. 0.d0) then 
	   write (*,*) log10( xnum )
	   STOP ' ERROR. Negative productions. XNUM=0.'
	elseif (xdenom .lt. 0.d0) then 
	   STOP ' ERROR. Negative losses. XDENOM=0.'
	endif

        if (xnum .eq. 0.d0) then 
          cociente = minimo
          return
        endif

        if (xdenom .eq. 0.d0) then 
          if (option .eq. 0) then 
             STOP ' ERROR. xdenom=0. '
          elseif (option .eq. 1) then 
!             write (*,*) 'WARNING !!    xdenom=0  '
!             write (*,*) 'WARNING !!    option=2 => cociente=minimo',
!     $            xdenom
             cociente = minimo
             return
          else
             STOP ' ERROR. option undefined in call to cociente'
          endif
        endif
            
        if (xdenom .lt. minimo) then 
           factor = xdenom * 1.d20
           cociente = xnum / factor * 1.d20
        else
	   cociente = xnum / xdenom
	endif

	if (cociente .lt. minimo) cociente = minimo
                                                
	return                                         
c END
	end
