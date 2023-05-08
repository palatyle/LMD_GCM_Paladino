      SUBROUTINE uniplot(im,jm,champ)
      INTEGER isys,system
      EXTERNAL system
      INTEGER imx,jmx
      PARAMETER (imx=400,jmx=200)
      INTEGER im,jm,i,j,nreg
      REAL champ(im,jm)
      REAL xinf,xsup,yinf,ysup
      REAL xorig,yorig,xlength,ylength,xfac,yfac
      REAL xreg(1),yreg(1)
      REAL x(imx),y(jmx)
      LOGICAL francais,flreg,flgri,fliso
      REAL txtheight,linewidth,fdismin,flonmin
      INTEGER ishade,inter_exter
      CHARACTER*200 txtstr(5)

      xorig=20.
      yorig=20.
      xlength=120.
      ylength=120.
      xfac=1.
      yfac=1.
      xinf=1.
      xsup=im
      yinf=1.
      ysup=jm
      nreg=1
      DO i=1,im
         x(i)=i
      ENDDO
      DO j=1,jm
         y(j)=j
      ENDDO

      DO i=1,5
       txtstr(i)=" "
      ENDDO

      txtstr(2)='i'
      txtstr(3)='j'

      OPEN(10,file='/tmp2/planetes/tmp.graph',
     s  status='unknown',form='formatted')

      WRITE(10,*) im,jm,nreg
      WRITE(10,*) xorig,yorig,xlength,ylength,xfac,yfac,
     s     xinf,xsup,yinf,ysup
      WRITE(10,*) (x(i),i=1,im),(y(j),j=1,jm),
     s     ((champ(i,j),i=1,im),j=1,jm),
     s     (xreg(i),i=1,nreg),(yreg(j),j=1,nreg)

      francais=.true.
      flreg=.false.
      flgri=.false.
      fliso=.true.
      txtheight=3.7
      linewidth=.2
      fdismin=0.
      flonmin=0.
      ishade=3
      inter_exter=0
      WRITE(10,*) francais,flreg,flgri,fliso
      WRITE(10,*) txtheight,linewidth,fdismin,flonmin
      WRITE(10,*) ishade,inter_exter
      DO i=1,5
      WRITE(10,'(a200)') txtstr(i)
      ENDDO
      CLOSE(10)

      isys=system(
     s  'rsh ella " xterm -display $DISPLAY -n controle_graph
     s   -e /u2/planetes/graph/lance" ' )

      END
