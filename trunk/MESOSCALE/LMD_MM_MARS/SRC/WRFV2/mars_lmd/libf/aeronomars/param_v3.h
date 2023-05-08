c**********************************************************************

c	param.cmn
	
c	fgg     mar 02          first version
c**********************************************************************

      common/uv/ anchint,crscabsi2,freccen,c123,c1,c2,c12,c23,ch2o2
     $ ,jabsifotsint,jfotsout,co2crsc195,co2crsc295,fluxtop,t0
     
     


      real       anchint(ninter)          !widht of each interval
      real       crscabsi2(nabs,ninter)    !cross section
      real       freccen(ninter)          !representative wavelenght
      real       c123(ninter2,nz2)        !column densities (cm^-2)
                                          !(CO2+O2+O3P)
      real       c1(nz2)                  !idem (CO2) (cm^-2)
      real       c2(nz2)                  !idem (O2) (cm^-2)
      real       c12(nz2)                 !idem (CO2+O2) (cm^-2)
      real       c23(nz2)                 !idem (O2+O3P) (cm^-2)
      real       ch2o2(nz2)               !idem (H2O2) (cm^-2)
      real       jabsifotsint(ninter,nabs,nzmax) !photoabsorption
                                                 !coefficients (fot/s)
      real       jfotsout(ninter,nabs,nzmax)
      real       co2crsc195(9)
      real       co2crsc295(9)
      real       fluxtop(ninter)
      real       t0(nzmax)
      


      common/fotoquimica/  ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,
     $                     ch13,ch14,ch15,ch18,ch19,ch20,ch21
     $			   ,jdistot,jdistot_b

      real*8 ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,ch13,ch14,ch15,ch18
      real*8 ch19,ch20,ch21
						  !reaction rates
      real jdistot(nabs,nzmax)    !photodissociation rates
      real jdistot_b(2,nzmax)



      common/prod/ Pco2,Po2,Po3p,Pco,Ph,Poh,Pho2,Ph2,Ph2o,Po1d,Ph2o2,
     $		   Pco2tot,Po2tot,Po3ptot,Pcotot,Phtot,Pohtot,
     $             Pho2tot,Ph2tot,Ph2otot,Po1dtot,Ph2o2tot
	
      real*8 Pco2(nzmax,nreact),Po2(nzmax,nreact),Po3p(nzmax,nreact)
      real*8 Pco(nzmax,nreact),Ph(nzmax,nreact),Poh(nzmax,nreact)
      real*8 Pho2(nzmax,nreact),Ph2(nzmax,nreact),Ph2o(nzmax,nreact)
      real*8 Po1d(nzmax,nreact),Ph2o2(nzmax,nreact)
      real*8 Pco2tot(nzmax),Po2tot(nzmax),Po3ptot(nzmax),Pcotot(nzmax)
      real*8 Phtot(nzmax)
      real*8 Pohtot(nzmax),Pho2tot(nzmax),Ph2tot(nzmax),Ph2otot(nzmax)
      real*8 Po1dtot(nzmax),Ph2o2tot(nzmax)


      common/loss/ Lco2,Lo2,Lo3p,Lco,Lh,Loh,Lho2,Lh2,Lh2o,Lo1d,Lh2o2,
     $		   Lco2tot,Lo2tot,Lo3ptot,Lcotot,Lhtot,Lohtot,
     $             Lho2tot,Lh2tot,Lh2otot,Lo1dtot,Lh2o2tot
	
      real*8 Lco2(nzmax,nreact),Lo2(nzmax,nreact),Lo3p(nzmax,nreact)
      real*8 Lco(nzmax,nreact),Lh(nzmax,nreact),Loh(nzmax,nreact)
      real*8 Lho2(nzmax,nreact),Lh2(nzmax,nreact),Lh2o(nzmax,nreact)
      real*8 Lo1d(nzmax,nreact),Lh2o2(nzmax,nreact)
      real*8 Lco2tot(nzmax),Lo2tot(nzmax),Lo3ptot(nzmax),Lcotot(nzmax)
      real*8 Lhtot(nzmax),Lohtot(nzmax),Lho2tot(nzmax),Lh2tot(nzmax)
      real*8 Lh2otot(nzmax),Lo1dtot(nzmax),Lh2o2tot(nzmax)


      common/vida/ tminco2,tmino2,tmino3p,tminco,tminh,tminoh,tminho2
     $               ,tminh2,tminh2o,tmino1d,tminh2o2

    	real*8 tminco2(nzmax),tmino2(nzmax),tmino3p(nzmax),tminco(nzmax)
    	real*8 tminh(nzmax),tminoh(nzmax),tminho2(nzmax),tminh2(nzmax)
    	real*8 tminh2o(nzmax),tmino1d(nzmax),tminh2o2(nzmax)


      common/flujoctes/ fluxtophr,ct1,p1,ct2,p2   
 
      real fluxtophr(ninter)
      real ct1(ninter),p1(ninter),ct2(ninter),p2(ninter)       

      common/phdisef/ efdisco2, efdiso2, efdish2, efdish2o, efdish2o2

      real efdisco2(ninter), efdiso2(ninter), efdish2o(ninter)
      real efdish2o2(ninter), efdish2(ninter)






