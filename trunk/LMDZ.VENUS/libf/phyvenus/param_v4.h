!**********************************************************************

!	param.cmn
	
!	fgg     mar 02          first version
!**********************************************************************



      common/uv/ co2crsc195,co2crsc295,crscabsi2,                                          &
     &     e107,date_e107,e107_tab,freccen               

      real co2crsc195(9)
      real co2crsc295(9),crscabsi2(13,16)
      real freccen(36)          !representative wavelenght
      real e107,date_e107(669),e107_tab(669)



      common/fotoquimica/  ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,                &
     &                     ch13,ch14,ch15,ch18,ch19,ch20,ch21,ch22,              &
     &                     ch23,ch24,ch30,ch31,ch32,ch33,ch34,ch35,              &
     &                     ch36,ch37,ch38,ch39,ch40,ch41,ch42,ch43,              &
     &                     ch45,ch46,ch47,ch48,ch49,ch50,ch55,ch56,              &
     &                     ch57,ch58,ch59,ch62,ch63,ch64,ch65,                   &
     &                     ch66,ch67,ch68,ch69,ch70,ch71,ch72,                   &
     &                     ch73,ch74,ch75,ch76,ch85,ch86,ch87                   
!     &                     jdistot,jdistot_b,jion

      real*8 ch2, ch3, ch4, ch5, ch7,ch9,ch10,ch11,ch13,ch14,ch15,ch18
      real*8 ch19,ch20,ch21,ch22,ch23,ch24,ch30,ch31,ch32,ch33,ch34
      real*8 ch35,ch36,ch37,ch38,ch39,ch40,ch41,ch42,ch43,ch45
      real*8 ch46,ch47,ch48,ch49,ch50,ch55,ch56,ch57,ch58,ch59,ch62 
      real*8 ch63,ch64,ch65,ch66,ch67,ch68,ch69,ch70,ch71
      real*8 ch72,ch73,ch74,ch75,ch76,ch85,ch86,ch87

						  !reaction rates
!      real jdistot(nabs,nlayermx)    !photodissociation rates
!      real jdistot_b(nabs,nlayermx)
     
!      real jion(nabs,nlayermx,4)



!      common/prod/ Pco2,Po2,Po3p,Pco,Ph,Poh,Pho2,Ph2,Ph2o,Po1d,Ph2o2,          &
!     &             Po3,Pn,Pno,Pno2,Pn2,Pn2d,                                   &
!     &             Pco2plus,Poplus,Po2plus,Pelect,Pcoplus,Pcplus,              &
!     &             Pnplus, Pnoplus, Pn2plus, Phplus, Phco2plus,                &
!     &  	   Pco2tot,Po2tot,Po3ptot,Pcotot,Phtot,Pohtot,                 &
!     &             Pho2tot,Ph2tot,Ph2otot,Po1dtot,Ph2o2tot,                    &
!     &             Po3tot,Pntot,Pnotot,Pno2tot,Pn2tot,Pn2dtot,                 &
!     &             Pcoplustot,Pcplustot,                                       &
!     &             Pco2plustot,Poplustot,Po2plustot,Pelecttot,                 &
!     &             Pnplustot,Pnoplustot,Pn2plustot,Phplustot,                  &
!     &             Phco2plustot  

	
!      real*8 Pco2(nlayermx,nreact),Po2(nlayermx,nreact)
!      real*8 Po3p(nlayermx,nreact)
!      real*8 Pco(nlayermx,nreact),Ph(nlayermx,nreact)
!      real*8 Poh(nlayermx,nreact)
!      real*8 Pho2(nlayermx,nreact),Ph2(nlayermx,nreact)
!      real*8 Ph2o(nlayermx,nreact)
!      real*8 Po1d(nlayermx,nreact),Ph2o2(nlayermx,nreact)
!      real*8 Po3(nlayermx,nreact),Pn(nlayermx,nreact)
!      real*8 Pno(nlayermx,nreact)
!      real*8 Pno2(nlayermx,nreact),Pn2(nlayermx,nreact)
!      real*8 Pn2d(nlayermx,nreact)
!      real*8 Pco2plus(nlayermx,nreact),Poplus(nlayermx,nreact)
!      real*8 Po2plus(nlayermx,nreact), Pelect(nlayermx,nreact)
!      real*8 Pcoplus(nlayermx,nreact),Pcplus(nlayermx,nreact)
!      real*8 Pnplus(nlayermx,nreact),Pnoplus(nlayermx,nreact)
!      real*8 Pn2plus(nlayermx,nreact),Phplus(nlayermx,nreact)
!      real*8 Phco2plus(nlayermx,nreact)
!      real*8 Pco2tot(nlayermx),Po2tot(nlayermx)
!      real*8 Po3ptot(nlayermx),Pcotot(nlayermx)
!      real*8 Phtot(nlayermx)
!      real*8 Pohtot(nlayermx),Pho2tot(nlayermx)
!      real*8 Ph2tot(nlayermx),Ph2otot(nlayermx)
!      real*8 Po1dtot(nlayermx),Ph2o2tot(nlayermx)
!      real*8 Po3tot(nlayermx),Pntot(nlayermx),Pnotot(nlayermx)
!      real*8 Pno2tot(nlayermx),Pn2tot(nlayermx),Pn2dtot(nlayermx)
!      real*8 Pco2plustot(nlayermx), Poplustot(nlayermx)
!      real*8 Po2plustot(nlayermx),Pelecttot(nlayermx)
!      real*8 Pcoplustot(nlayermx), Pcplustot(nlayermx)
!      real*8 Pnplustot(nlayermx), Pnoplustot(nlayermx)
!      real*8 Pn2plustot(nlayermx),Phplustot(nlayermx)
!      real*8 Phco2plustot(nlayermx)



!      common/loss/ Lco2,Lo2,Lo3p,Lco,Lh,Loh,Lho2,Lh2,Lh2o,Lo1d,Lh2o2,         &
!     &             Lo3,Ln,Lno,Lno2,Ln2,Ln2d,Lelect,                           &
!     &             Lco2plus,Lo2plus,Loplus, Lcoplus, Lcplus,                  &
!     &             Lnplus,Lnoplus,Ln2plus, Lhplus, Lhco2plus,                 &
!     &             Lco2tot,Lo2tot,Lo3ptot,Lcotot,Lhtot,Lohtot,                &
!     &             Lho2tot,Lh2tot,Lh2otot,Lo1dtot,Lh2o2tot,                   &
!     &             Lo3tot,Lntot,Lnotot,Lno2tot,Ln2tot,Ln2dtot,                &
!     &             Lcoplustot,Lcplustot,Lco2plustot,Loplustot,                &
!     &             Lo2plustot,Lelecttot,Lnplustot,Lnoplustot,                 &
!     &             Ln2plustot,Lhplustot,Lhco2plustot

	
!      real*8 Lco2(nlayermx,nreact),Lo2(nlayermx,nreact)
!      real*8 Lo3p(nlayermx,nreact)
!      real*8 Lco(nlayermx,nreact),Lh(nlayermx,nreact)
!      real*8 Loh(nlayermx,nreact)
!      real*8 Lho2(nlayermx,nreact),Lh2(nlayermx,nreact)
!      real*8 Lh2o(nlayermx,nreact)
!      real*8 Lo1d(nlayermx,nreact),Lh2o2(nlayermx,nreact)
!      real*8 Lo3(nlayermx,nreact),Ln(nlayermx,nreact)
!      real*8 Lno(nlayermx,nreact)
!      real*8 Lno2(nlayermx,nreact),Ln2(nlayermx,nreact)
!      real*8 Ln2d(nlayermx,nreact)
!      real*8 Lco2plus(nlayermx,nreact), Loplus(nlayermx,nreact)
!      real*8 Lo2plus(nlayermx,nreact), Lelect(nlayermx,nreact)
!      real*8 Lcoplus(nlayermx,nreact), Lcplus(nlayermx,nreact)
!      real*8 Lnplus(nlayermx,nreact), Lnoplus(nlayermx,nreact)
!      real*8 Ln2plus(nlayermx,nreact), Lhplus(nlayermx,nreact)
!      real*8 Lhco2plus(nlayermx,nreact)
!      real*8 Lco2tot(nlayermx),Lo2tot(nlayermx)
!      real*8 Lo3ptot(nlayermx),Lcotot(nlayermx)
!      real*8 Lhtot(nlayermx),Lohtot(nlayermx)
!      real*8 Lho2tot(nlayermx),Lh2tot(nlayermx)
!      real*8 Lh2otot(nlayermx),Lo1dtot(nlayermx)
!      real*8 Lh2o2tot(nlayermx)
!      real*8 Lo3tot(nlayermx),Lntot(nlayermx),Lnotot(nlayermx)
!      real*8 Lno2tot(nlayermx),Ln2tot(nlayermx)
!      real*8 Ln2dtot(nlayermx)
!      real*8 Lco2plustot(nlayermx), Loplustot(nlayermx)
!      real*8 Lo2plustot(nlayermx), Lelecttot(nlayermx)
!      real*8 Lcoplustot(nlayermx), Lcplustot(nlayermx)
!      real*8 Lnplustot(nlayermx), Lnoplustot(nlayermx)
!      real*8 Ln2plustot(nlayermx), Lhplustot(nlayermx)
!      real*8 Lhco2plustot(nlayermx)



!      common/vida/ tminco2,tmino2,tmino3p,tminco,tminh,tminoh,                &
!     &      tminho2,tminh2,tminh2o,tmino1d,tminh2o2,tmino3,tminn,             &
!     &      tminno,tminno2,tminn2,tminn2d,                                    &
!     &      tminco2plus,tminoplus,tmino2plus,tmincoplus,                      &
!     &      tmincplus,tminnplus,tminnoplus,                                   &
!     &      tminn2plus,tminhplus,tminhco2plus



!    	real*8 tminco2(nlayermx),tmino2(nlayermx)
!        real*8 tmino3p(nlayermx),tminco(nlayermx)
!    	real*8 tminh(nlayermx),tminoh(nlayermx)
!        real*8 tminho2(nlayermx),tminh2(nlayermx)
!    	real*8 tminh2o(nlayermx),tmino1d(nlayermx),tminh2o2(nlayermx)
!        real*8 tmino3(nlayermx),tminn(nlayermx),tminno(nlayermx)
!        real*8 tminno2(nlayermx),tminn2(nlayermx),tminn2d(nlayermx)
!        real*8 tminco2plus(nlayermx),tminoplus(nlayermx)
!        real*8 tmino2plus(nlayermx)
!        real*8 tmincoplus(nlayermx)
!        real*8 tmincplus(nlayermx)
!        real*8 tminnplus(nlayermx),tminnoplus(nlayermx)
!        real*8 tminn2plus(nlayermx)
!        real*8 tminhplus(nlayermx)
!        real*8 tminhco2plus(nlayermx)


!      common/flujoctes/ fluxtophr,ct1,p1,ct2,p2   
 
!        real fluxtophr(ninter)
!        real ct1(ninter),p1(ninter),ct2(ninter),p2(ninter)       

!      common/phdisef/ efdisco2, efdiso2, efdish2, efdish2o,                   &
!     & efdish2o2,efdiso3,efdiso,efdisn,efdish,efdisno,efdisn2,                &
!     & efdisno2,efdisco,efionco2,efionn2,efionco,                             &
!     & efiono3p,efionn,efionno,efionh
  

!        real efdisco2(ninter), efdiso2(ninter), efdish2o(ninter)
!        real efdish2o2(ninter), efdish2(ninter), efdiso3(ninter)
!        real efdiso(ninter), efdisn(ninter), efdish(ninter)
!        real efdisno(ninter), efdisn2(ninter), efdisno2(ninter)
!        real efdisco(ninter)
!        real efionco2(ninter,4)
!        real efionn2(ninter,2)
!        real efionco(ninter,3)
!        real efiono3p(ninter),efionn(ninter)
!        real efionno(ninter),efionh(ninter)






