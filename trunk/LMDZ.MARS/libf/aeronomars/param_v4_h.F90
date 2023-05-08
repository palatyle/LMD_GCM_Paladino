MODULE param_v4_h

   IMPLICIT NONE

   integer, parameter :: ninter=36
   integer, parameter :: nabs=13
   integer, parameter :: nz2=253
   integer, parameter :: ninter2=16
   real*8, parameter :: kboltzman = 1.381e-16
   real*8, parameter :: n_avog = 6.023e23
   real*8, parameter :: gg = 6.67259e-8
   real*8, parameter :: masa = 6.4163e26
   real*8, parameter :: radio = 3390.
   integer, parameter :: nreact=93
   integer, parameter :: tapas=42

   real crscabsi2(nabs,16)    !cross section
   real c1_16(nz2,16)   !Col. int. 1 (o2+o+h2+n) (cm^-2)
   real c17_24(nz2)     !Col. int. 17-24 (co2+o2+n2+no+co+no2) (cm^-2)
   real c25_29(nz2)     !Col. int. 25-29 (co2+o2+h2o+h2o2+no+co+no2) (cm^-2)
   real c30_31(nz2)     !Col. int. 30-31 (co2+o2+h2o+h2o2+no+no2)
   real c32(nz2)        !col. int. 32 (co2+h2o2+o2+no+no2) (cm^-2)
   real c33(nz2)        !col. int. 33 (h2o2+o2+no2) (cm^-2)
   real c34(nz2)        !col. int. 34 (h2o2+o2+o3+no2) (cm^-2)
   real c35(nz2)        !col. int. 35 (h2o2+o3+no2) (cm^-2)
   real c36(nz2)        !col. int. 36 (o3+no2) (cm^-2)
   real co2crsc195(9)
   real co2crsc295(9)
   real t0(nz2)
   real fluxtop(ninter)
   real freccen(ninter)          !representative wavelenght
   real jabsifotsintpar(nz2,nabs,ninter)
   real e107,date_e107(669),e107_tab(669)
   real coefit0(ninter,nabs),coefit1(ninter,nabs)
   real coefit2(ninter,nabs)
   real coefit3(ninter,nabs),coefit4(ninter,nabs)

   !reaction rates
   real*8 ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,ch13,ch14,ch15,ch18
   real*8 ch19,ch20,ch21,ch22,ch23,ch24,ch30,ch31,ch32,ch33,ch34
   real*8 ch35,ch36,ch37,ch38,ch39,ch40,ch41,ch42,ch43,ch45
   real*8 ch46,ch47,ch48,ch49,ch50,ch55,ch56,ch57,ch58,ch59,ch62
   real*8 ch63,ch64,ch65,ch66,ch67,ch68,ch69,ch70,ch71
   real*8 ch72,ch73,ch74,ch75,ch76,ch85,ch86,ch87
   real*8 rcoef(61,3)

   real fluxtophr(ninter)
   real ct1(ninter),p1(ninter),ct2(ninter),p2(ninter)

   real efdisco2(ninter), efdiso2(ninter), efdish2o(ninter)
   real efdish2o2(ninter), efdish2(ninter), efdiso3(ninter)
   real efdiso(ninter), efdisn(ninter), efdish(ninter)
   real efdisno(ninter), efdisn2(ninter), efdisno2(ninter)
   real efdisco(ninter)
   real efionco2(ninter,4)
   real efiono2(ninter,2)
   real efionn2(ninter,2)
   real efionco(ninter,3)
   real efiono3p(ninter),efionn(ninter)
   real efionno(ninter),efionh(ninter)

   !photodissociation rates
   REAL,SAVE,ALLOCATABLE :: jfotsout(:,:,:)
   REAL,SAVE,ALLOCATABLE :: jdistot(:,:)
   REAL,SAVE,ALLOCATABLE :: jdistot_b(:,:)
   REAL,SAVE,ALLOCATABLE :: jion(:,:,:)

   REAL*8,SAVE,ALLOCATABLE :: Pco2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Po2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Po3p(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pco(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ph(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Poh(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pho2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ph2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ph2o(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Po1d(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ph2o2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Po3(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pn(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pno(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pno2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pn2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pn2d(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pco2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Poplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Po2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pelect(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pcoplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pcplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pnplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pnoplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pn2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Phplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Phco2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Pco2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Po2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Po3ptot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pcotot(:)
   REAL*8,SAVE,ALLOCATABLE :: Phtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pohtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pho2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Ph2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Ph2otot(:)
   REAL*8,SAVE,ALLOCATABLE :: Po1dtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Ph2o2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Po3tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pntot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pnotot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pno2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pn2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pn2dtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pco2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Poplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Po2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pelecttot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pcoplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pcplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pnplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pnoplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Pn2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Phplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Phco2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lco2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lo2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lo3p(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lco(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lh(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Loh(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lho2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lh2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lh2o(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lo1d(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lh2o2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lo3(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ln(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lno(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lno2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ln2(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ln2d(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lco2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Loplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lo2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lelect(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lcoplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lcplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lnplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lnoplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Ln2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lhplus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lhco2plus(:,:)
   REAL*8,SAVE,ALLOCATABLE :: Lco2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lo2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lo3ptot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lcotot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lhtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lohtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lho2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lh2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lh2otot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lo1dtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lh2o2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lo3tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lntot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lnotot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lno2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Ln2tot(:)
   REAL*8,SAVE,ALLOCATABLE :: Ln2dtot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lco2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Loplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lo2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lelecttot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lcoplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lcplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lnplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lnoplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Ln2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lhplustot(:)
   REAL*8,SAVE,ALLOCATABLE :: Lhco2plustot(:)
   REAL*8,SAVE,ALLOCATABLE :: tminco2(:)
   REAL*8,SAVE,ALLOCATABLE :: tmino2(:)
   REAL*8,SAVE,ALLOCATABLE :: tmino3p(:)
   REAL*8,SAVE,ALLOCATABLE :: tminco(:)
   REAL*8,SAVE,ALLOCATABLE :: tminh(:)
   REAL*8,SAVE,ALLOCATABLE :: tminoh(:)
   REAL*8,SAVE,ALLOCATABLE :: tminho2(:)
   REAL*8,SAVE,ALLOCATABLE :: tminh2(:)
   REAL*8,SAVE,ALLOCATABLE :: tminh2o(:)
   REAL*8,SAVE,ALLOCATABLE :: tmino1d(:)
   REAL*8,SAVE,ALLOCATABLE :: tminh2o2(:)
   REAL*8,SAVE,ALLOCATABLE :: tmino3(:)
   REAL*8,SAVE,ALLOCATABLE :: tminn(:)
   REAL*8,SAVE,ALLOCATABLE :: tminno(:)
   REAL*8,SAVE,ALLOCATABLE :: tminno2(:)
   REAL*8,SAVE,ALLOCATABLE :: tminn2(:)
   REAL*8,SAVE,ALLOCATABLE :: tminn2d(:)
   REAL*8,SAVE,ALLOCATABLE :: tminco2plus(:)
   REAL*8,SAVE,ALLOCATABLE :: tminoplus(:)
   REAL*8,SAVE,ALLOCATABLE :: tmino2plus(:)
   REAL*8,SAVE,ALLOCATABLE :: tmincoplus(:)
   REAL*8,SAVE,ALLOCATABLE :: tmincplus(:)
   REAL*8,SAVE,ALLOCATABLE :: tminnplus(:)
   REAL*8,SAVE,ALLOCATABLE :: tminnoplus(:)
   REAL*8,SAVE,ALLOCATABLE :: tminn2plus(:)
   REAL*8,SAVE,ALLOCATABLE :: tminhplus(:)
   REAL*8,SAVE,ALLOCATABLE :: tminhco2plus(:)

   CONTAINS

      SUBROUTINE fill_data_thermos

      IMPLICIT NONE

      ! data for the UV heating tabulation

      crscabsi2(1,1:16) = (/ 5.61031E-19,1.59677E-18,4.7072E-18,&
          1.48254e-17,2.07445e-17,2.573e-17,2.901e-17,3.083e-17,&
          3.217e-17,3.539e-17,3.658e-17,3.63e-17,3.41239e-17,&
          2.71019e-17,4.93677e-17,1.64e-17 /) 

      crscabsi2(2,1:16) = (/ 0.27250E-18,0.11650E-17,0.39250E-17,&
          0.10630E-16,0.15590E-16,0.17180E-16,0.19270E-16,0.22860E-16,&
          0.24270E-16,0.24440E-16,0.25020E-16,0.26600E-16,0.25400E-16,&
          0.35800E-16,0.25590E-16,0.16740E-16 /)

      crscabsi2(3,1:16) = (/ 0.2776E-18,0.9792E-18,0.3313E-17,&
          0.6621E-17,0.8481E-17,0.9146E-17,0.9414E-17,0.1039E-16,&
          0.1012E-16,0.1033E-16,0.1033E-16,0.1033E-16,0.8268E-17,&
          0.6563E-17,0.3506E-17,0.3470E-17 /) 

      crscabsi2(5,1:16) = (/ .5E-20,.1077607E-19,.5670491E-19,&
          .3322716E-18,.1054509E-17,.1700005E-17,.3171188E-17,&
          .4734241E-17,.5108741E-17,.6022236E-17,.6741537E-17,&
          .7277079E-17,.9070787E-17,.9708916E-17,.4026281E-17,0.0 /)

      crscabsi2(8,1:16) = (/ 0.0, 7.44175e-19, 2.23167e-18,&
          8.46200e-18,1.18275e-17,1.54900e-17,2.32475e-17,2.41373e-17,&
          2.55482e-17,2.38431e-17,2.28600e-17,2.35067e-17,2.56000e-17,&
          2.64636e-17,2.86260e-17,3.26561e-17 /)

      crscabsi2(9,1:16) = (/ 3.48182e-20,3.37038e-19,1.03077e-18,&
          4.01364e-18,6.45e-18,7.8e-18,1.0e-17,1.13500e-17,1.15500e-17,&
          1.18000e-17,1.17500e-17,1.16000e-17,1.28667e-17,1.18500e-17,&
          1.11000e-17,9.50000e-18 /)

      crscabsi2(10,1:16) = (/ 0.0,9.39833e-19,2.87714e-18,&
          9.66900e-18,1.37063e-17,1.61820e-17,2.30450e-17,2.63373e-17,&
          2.63773e-17,2.67677e-17,2.64100e-17,2.53000e-17,2.18100e-17,&
          2.04941e-17,2.28160e-17,2.93550e-17 /)

      crscabsi2(11,1:16) = (/ 0.0,9.58555e-19,2.52767e-18,&
          8.29700e-18,1.21850e-17,1.40500e-17,1.97025e-17,2.12018e-17,&
          2.14673e-17,2.20331e-17,2.21500e-17,2.21600e-17,2.33200e-17,&
          2.67800e-17,2.56400e-17,3.58561e-17 /)

      crscabsi2(12,1:16) = (/ 0.0,1.0e-20,2.5e-20,1.30400e-19,&
          2.93800e-19,4.36000e-19,8.49400e-19,1.29400e-18,1.40500e-18,&
          1.67600e-18,1.93400e-18,2.12200e-18,2.75800e-18,3.48400e-18,&
          4.17200e-18,5.26000e-18 /)

      crscabsi2(13,1:16) = (/ 0.0,1.60e-18,4.99111e-18,1.48496e-17,&
          2.17395e-17,2.55857e-17,2.87754e-17,3.65571e-17,3.85691e-17,&
          4.16286e-17,4.15117e-17,4.05901e-17,3.64000e-17,2.99670e-17,&
          2.46796e-17,2.51789e-17 /)

      freccen(1:ninter)=(/ 3.4,7.5,14.5,23.0,30.3,34.1,&
           49.6,50.5,52.5,56.0,&
           59.0,61.5,68.7,73.1,78.4,83.1,92.4,97.5,99.3,100.1,100.7,102.1,&
           104.5,116.8,121.3,127.0,130.6,153.7,162.8,171.4,&
           195.6,206.3,222.0,236.0,289.0,600. /)

      co2crsc195(1:9)=(/ 2.05864e-17,5.90557e-20,3.1027e-19,6.70653e-19,&
           4.55132e-19,8.87122e-20,1.32138e-20,7.22244e-23,2.88002e-26 /)

      co2crsc295(1:9)=(/2.05897e-17,6.71104e-20,3.45509e-19,7.45711e-19,&
           4.82752e-19,1.11594e-19,1.98308e-20,1.3853e-22,2.1414e-25 /)

      END SUBROUTINE fill_data_thermos

      SUBROUTINE allocate_param_thermos(nlayer)

        IMPLICIT NONE

        INTEGER :: nlayer
        allocate(jdistot(nabs,nlayer))    
        allocate(jdistot_b(nabs,nlayer))
        allocate(jion(nabs,nlayer,4))
        allocate(jfotsout(ninter,nabs,nlayer))
        allocate(Pco2(nlayer,nreact))
        allocate(Po2(nlayer,nreact))
        allocate(Po3p(nlayer,nreact))
        allocate(Pco(nlayer,nreact))
        allocate(Ph(nlayer,nreact))
        allocate(Poh(nlayer,nreact))
        allocate(Pho2(nlayer,nreact))
        allocate(Ph2(nlayer,nreact))
        allocate(Ph2o(nlayer,nreact))
        allocate(Po1d(nlayer,nreact))
        allocate(Ph2o2(nlayer,nreact))
        allocate(Po3(nlayer,nreact))
        allocate(Pn(nlayer,nreact))
        allocate(Pno(nlayer,nreact))
        allocate(Pno2(nlayer,nreact))
        allocate(Pn2(nlayer,nreact))
        allocate(Pn2d(nlayer,nreact))
        allocate(Pco2plus(nlayer,nreact))
        allocate(Poplus(nlayer,nreact))
        allocate(Po2plus(nlayer,nreact))
        allocate(Pelect(nlayer,nreact))
        allocate(Pcoplus(nlayer,nreact))
        allocate(Pcplus(nlayer,nreact))
        allocate(Pnplus(nlayer,nreact))
        allocate(Pnoplus(nlayer,nreact))
        allocate(Pn2plus(nlayer,nreact))
        allocate(Phplus(nlayer,nreact))
        allocate(Phco2plus(nlayer,nreact))
        allocate(Pco2tot(nlayer))
        allocate(Po2tot(nlayer))
        allocate(Po3ptot(nlayer))
        allocate(Pcotot(nlayer))
        allocate(Phtot(nlayer))
        allocate(Pohtot(nlayer))
        allocate(Pho2tot(nlayer))
        allocate(Ph2tot(nlayer))
        allocate(Ph2otot(nlayer))
        allocate(Po1dtot(nlayer))
        allocate(Ph2o2tot(nlayer))
        allocate(Po3tot(nlayer))
        allocate(Pntot(nlayer))
        allocate(Pnotot(nlayer))
        allocate(Pno2tot(nlayer))
        allocate(Pn2tot(nlayer))
        allocate(Pn2dtot(nlayer))
        allocate(Pco2plustot(nlayer))
        allocate(Poplustot(nlayer))
        allocate(Po2plustot(nlayer))
        allocate(Pelecttot(nlayer))
        allocate(Pcoplustot(nlayer))
        allocate(Pcplustot(nlayer))
        allocate(Pnplustot(nlayer))
        allocate(Pnoplustot(nlayer))
        allocate(Pn2plustot(nlayer))
        allocate(Phplustot(nlayer))
        allocate(Phco2plustot(nlayer))
        allocate(Lco2(nlayer,nreact))
        allocate(Lo2(nlayer,nreact))
        allocate(Lo3p(nlayer,nreact))
        allocate(Lco(nlayer,nreact))
        allocate(Lh(nlayer,nreact))
        allocate(Loh(nlayer,nreact))
        allocate(Lho2(nlayer,nreact))
        allocate(Lh2(nlayer,nreact))
        allocate(Lh2o(nlayer,nreact))
        allocate(Lo1d(nlayer,nreact))
        allocate(Lh2o2(nlayer,nreact))
        allocate(Lo3(nlayer,nreact))
        allocate(Ln(nlayer,nreact))
        allocate(Lno(nlayer,nreact))
        allocate(Lno2(nlayer,nreact))
        allocate(Ln2(nlayer,nreact))
        allocate(Ln2d(nlayer,nreact))
        allocate(Lco2plus(nlayer,nreact))
        allocate(Loplus(nlayer,nreact))
        allocate(Lo2plus(nlayer,nreact))
        allocate(Lelect(nlayer,nreact))
        allocate(Lcoplus(nlayer,nreact))
        allocate(Lcplus(nlayer,nreact))
        allocate(Lnplus(nlayer,nreact))
        allocate(Lnoplus(nlayer,nreact))
        allocate(Ln2plus(nlayer,nreact))
        allocate(Lhplus(nlayer,nreact))
        allocate(Lhco2plus(nlayer,nreact))
        allocate(Lco2tot(nlayer))
        allocate(Lo2tot(nlayer))
        allocate(Lo3ptot(nlayer))
        allocate(Lcotot(nlayer))
        allocate(Lhtot(nlayer))
        allocate(Lohtot(nlayer))
        allocate(Lho2tot(nlayer))
        allocate(Lh2tot(nlayer))
        allocate(Lh2otot(nlayer))
        allocate(Lo1dtot(nlayer))
        allocate(Lh2o2tot(nlayer))
        allocate(Lo3tot(nlayer))
        allocate(Lntot(nlayer))
        allocate(Lnotot(nlayer))
        allocate(Lno2tot(nlayer))
        allocate(Ln2tot(nlayer))
        allocate(Ln2dtot(nlayer))
        allocate(Lco2plustot(nlayer))
        allocate(Loplustot(nlayer))
        allocate(Lo2plustot(nlayer))
        allocate(Lelecttot(nlayer))
        allocate(Lcoplustot(nlayer))
        allocate(Lcplustot(nlayer))
        allocate(Lnplustot(nlayer))
        allocate(Lnoplustot(nlayer))
        allocate(Ln2plustot(nlayer))
        allocate(Lhplustot(nlayer))
        allocate(Lhco2plustot(nlayer))
        allocate(tminco2(nlayer))
        allocate(tmino2(nlayer))
        allocate(tmino3p(nlayer))
        allocate(tminco(nlayer))
        allocate(tminh(nlayer))
        allocate(tminoh(nlayer))
        allocate(tminho2(nlayer))
        allocate(tminh2(nlayer))
        allocate(tminh2o(nlayer))
        allocate(tmino1d(nlayer))
        allocate(tminh2o2(nlayer))
        allocate(tmino3(nlayer))
        allocate(tminn(nlayer))
        allocate(tminno(nlayer))
        allocate(tminno2(nlayer))
        allocate(tminn2(nlayer))
        allocate(tminn2d(nlayer))
        allocate(tminco2plus(nlayer))
        allocate(tminoplus(nlayer))
        allocate(tmino2plus(nlayer))
        allocate(tmincoplus(nlayer))
        allocate(tmincplus(nlayer))
        allocate(tminnplus(nlayer))
        allocate(tminnoplus(nlayer))
        allocate(tminn2plus(nlayer))
        allocate(tminhplus(nlayer))
        allocate(tminhco2plus(nlayer))

      END SUBROUTINE allocate_param_thermos

END MODULE param_v4_h
