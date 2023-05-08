      integer nfmx,imx,jmx,lmx,nf,nvarmx
      parameter(nfmx=10,imx=200,jmx=150,lmx=20,nvarmx=1000)

      real xd(imx,nfmx),yd(jmx,nfmx),zd(lmx,nfmx),dtime(nfmx)

      integer imd(imx),jmd(jmx),lmd(lmx)
      integer iid(imx),jid(jmx)
      integer ifd(imx),jfd(jmx)
      integer unit(nfmx),irec(nfmx),itime(nfmx),nld(nvarmx,nfmx)

      integer nvar(nfmx),ivar(nfmx)
      logical firsttime(nfmx)

      character*10 var(nvarmx,nfmx),fichier(nfmx)
      character*40 title(nfmx),tvar(nvarmx,nfmx)

      common/gradsdef/xd,yd,zd,dtime,
     s   imd,jmd,lmd,iid,jid,ifd,jfd,
     s   unit,irec,nvar,ivar,nf,itime,nld,firsttime,
     s   var,fichier,title,tvar
