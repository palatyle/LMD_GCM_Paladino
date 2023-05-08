      subroutine blackl(blalong,blat,blae)
c.......................................................................
      implicit double precision (a-h,o-z)
c.......................................................................
c.....donnees de fond
      sigma=5.6693d-08
      pi=datan(1.d0)*4.d0
      c0=2.9979d+08
      h=6.6262d-34
      cbol=1.3806d-23
      rind=1.d0
      c=c0/rind
      c1=h*(c**2)
      c2=h*c/cbol
c.....fin de donnees de fond
c.......................................................................
      blae=2.d0*pi*c1/blalong**5/(dexp(c2/blalong/blat)-1.d0)
c.......................................................................
      return
      end
