c
c $Header$
c
c calculs statistiques distribution nuage ftion du regime dynamique 
c
c Ce calcul doit etre fait a partir de valeurs mensuelles ??
      CALL histo_o500_pctau(nbregdyn,pct_ocean,o500,fq_isccp,
     &histoW,nhistoW)
c
c nhistoWt = somme de toutes les nhistoW
      DO nreg=1, nbregdyn
       DO k = 1, kmaxm1
        DO l = 1, lmaxm1
         DO iw = 1, iwmax
          nhistoWt(k,l,iw,nreg)=nhistoWt(k,l,iw,nreg)+
     &    nhistoW(k,l,iw,nreg)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
c
cIM 190504 END
