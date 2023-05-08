*CDK comgeom2
      COMMON/comgeom/
     1 cu(iip1,jjp1),cv(iip1,jjm),unscu2(iip1,jjp1),unscv2(iip1,jjm)  ,
     2 aire(iip1,jjp1),airesurg(iip1,jjp1),aireu(iip1,jjp1)           ,
     3 airev(iip1,jjm),unsaire(iip1,jjp1),apoln,apols                 ,  
     4 unsairez(iip1,jjm),airuscv2(iip1,jjm),airvscu2(iip1,jjm)       ,
     5 aireij1(iip1,jjp1),aireij2(iip1,jjp1),aireij3(iip1,jjp1)       ,
     6 aireij4(iip1,jjp1),alpha1(iip1,jjp1),alpha2(iip1,jjp1)         ,
     7 alpha3(iip1,jjp1),alpha4(iip1,jjp1),alpha1p2(iip1,jjp1)        ,
     8 alpha1p4(iip1,jjp1),alpha2p3(iip1,jjp1),alpha3p4(iip1,jjp1)    ,
     9 fext(iip1,jjm),constang(iip1,jjp1), rlatu(jjp1),rlatv(jjm),
     1 rlonu(iip1),rlonv(iip1),cuvsurcv(iip1,jjm),cvsurcuv(iip1,jjm)  ,
     1 cvusurcu(iip1,jjp1),cusurcvu(iip1,jjp1)                        ,
     2 cuvscvgam1(iip1,jjm),cuvscvgam2(iip1,jjm),cvuscugam1(iip1,jjp1),
     3 cvuscugam2(iip1,jjp1),cvscuvgam(iip1,jjm),cuscvugam(iip1,jjp1) ,
     4 unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2                ,
     5 unsair_gam1(iip1,jjp1),unsair_gam2(iip1,jjp1)                  ,
     6 unsairz_gam(iip1,jjm),aivscu2gam(iip1,jjm),aiuscv2gam(iip1,jjm)
     7 , xprimu(iip1),xprimv(iip1)

c
      REAL 
     1 cu,cv,unscu2,unscv2,aire,airesurg,aireu,airev,apoln,apols,unsaire
     2 ,unsairez,airuscv2,airvscu2,aireij1,aireij2,aireij3,aireij4     ,
     3 alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 ,
     4 fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     ,
     5 cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1           , 
     6 unsapolnga2,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2     ,
     7 unsairz_gam,aivscu2gam,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu    ,
     8 cusurcvu,xprimu,xprimv
c
