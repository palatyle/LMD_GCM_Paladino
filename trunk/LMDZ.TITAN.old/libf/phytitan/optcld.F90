!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      MODULE optcld :
!
!      regroupe les variables communes au routines/fonctions pour l'extrapolation/interpolation
!      des proprietes optiques des gouttes.     
!
!      Contient les routines :
!        INITIALISATION  - lecture de la look-up table et initialisation des variables communes.
!        LOCATE          - Recherche de valeur dans une table.
!        INTERPOLEMOI    - interpolation d'une valeur a partir d'un jeu de donnees.
!        EXTRAPOLEMOI    - extrapolation des valeurs en dehors de la look-up table.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       MODULE optcld
         IMPLICIT NONE      
         REAL,SAVE             :: A_ex,A_ab,A_gg, B_ex,B_ab,B_gg
         REAL,SAVE             :: fmin_ex,fmin_ab,fmin_gg
         REAL,SAVE             :: r0cld
         REAL,SAVE             :: lmin,lmax
         INTEGER,SAVE          :: npts
         REAL,ALLOCATABLE,SAVE :: tq_ex(:),tq_ab(:),tq_gg(:),tq_wln(:)
         REAL,ALLOCATABLE,SAVE :: ltq_ex(:),ltq_ab(:),ltq_gg(:), &
                                  ltq_wln(:)
         REAL,SAVE             :: frac_c(3),fhvi_c,fhir_c,lseuil_c
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         CONTAINS 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        SUBROUTINE iniqcld :
!
!        Initialisation des variables commune et lecture de la look-up table.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           SUBROUTINE iniqcld() 
             IMPLICIT NONE
! ---------- LOCAL
             INTEGER             :: i,iun
             LOGICAL             :: ok
             CHARACTER*100       :: tmp
             iun=-1
!            rechercher une unite logique libre sur les 200 premieres.
             DO i=1,200
               INQUIRE (UNIT=i, OPENED=ok)
               IF (.not.ok) THEN
                 iun=i
                 exit
               ENDIF
             ENDDO
!            Une belle condition d'arret  ==> aucun unite dispo.
             IF (iun.eq.-1) THEN
               PRINT*,"CATASTROPHE !"
               PRINT*,"Impossible de trouver une unite logique libre", &
                      " sur les 200 premieres."
               PRINT*,"Je ne peux pas lire optcld.table."
               STOP "Je m'arrete (et salement en plus)..."
             ENDIF
!            Lecture de la table  !!!
             OPEN(iun,file='optcld.table')
             DO i=1,4
               READ(iun,*) tmp 
             ENDDO
             READ(iun,*) npts
!            Petite pause dans la lecture pour allouer les tableaux
             ALLOCATE(tq_ex(npts))
             ALLOCATE(tq_ab(npts))
             ALLOCATE(tq_gg(npts))
             ALLOCATE(tq_wln(npts))
             ALLOCATE(ltq_ex(npts))
             ALLOCATE(ltq_ab(npts))
             ALLOCATE(ltq_gg(npts))
             ALLOCATE(ltq_wln(npts))
!            Reprise de la              
             READ(iun,*) tmp
             READ(iun,'(ES14.7)') r0cld
             READ(iun,*) tmp 
             READ(iun,'(2(ES14.7,2X))') lmin,lmax 
             DO i=1,3
               READ(iun,*) tmp 
             ENDDO
             READ(iun,'(2(ES14.7,2X))') A_ex,B_ex 
             READ(iun,*) tmp 
             READ(iun,'(2(ES14.7,2X))') A_ab,B_ab 
             READ(iun,*) tmp 
             READ(iun,'(2(ES14.7,2X))') A_gg,B_gg
             DO i=1,3
               READ(iun,*) tmp 
             ENDDO
             READ(iun,'(ES14.7)') fmin_ex
             READ(iun,*) tmp 
             READ(iun,'(ES14.7)') fmin_ab
             READ(iun,*) tmp 
             READ(iun,'(ES14.7)') fmin_gg
             DO i=1,3
               READ(iun,*) tmp
             ENDDO
             READ(iun,*) (frac_c(i),i=1,3)
             DO i=1,2
               READ(iun,*) tmp
             ENDDO
             READ(iun,*) fhvi_c,fhir_c,lseuil_c
             READ(iun,*) tmp 

             DO i=1,npts
               READ(iun,'(4(ES23.15,1X))') &
               tq_wln(i),tq_ex(i),tq_ab(i),tq_gg(i)
! ------------ on passe tout en log pour les interpolations
               ltq_wln(i) = alog(tq_wln(i))
                ltq_ex(i) = alog(tq_ex(i))
                ltq_ab(i) = alog(tq_ab(i))
                ltq_gg(i) = alog(tq_gg(i))
             ENDDO
             CLOSE(iun)

             WRITE(*,*) &
             "LECTURE LOOK-UP optcld.table... TERMINEE :)"

           END SUBROUTINE iniqcld

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        SUBROUTINE locate(xx,n,x,j) :
!
!        Recherche de l'indice j de la valeur la plus proche par defaut de x dans le tableau xx de 
!        dimension n
!
!        ARGUMENTS D'ENTREE :
!           xx : tableau dans lequel rechercher l'indice.
!            x : valeur recherchee
!            n : dimension de xx
!
!        ARGUMENT DE SORTIE :
!           j : indice de la valeur la plus proche PAR DEFAUT de x 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           SUBROUTINE locate(xx,n,x,j)
             IMPLICIT NONE
! ---------- INPUT
             INTEGER,INTENT(in)  :: n
             REAL   ,INTENT(in)  :: x,xx(n)
             INTEGER,INTENT(out) :: j
! ---------- LOCAL
             INTEGER jl,jm,ju
             jl=0
             ju=n+1
             DO WHILE (ju-jl.gt.1)
               jm=(ju+jl)/2
               IF (jm.eq.0) STOP "ALERTE jm=0 !!"
               IF((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) THEN
                 jl=jm
               ELSE
                 ju=jm
               ENDIF
             ENDDO
             IF (x.eq.xx(1))THEN
               j=1
             ELSE IF (x.eq.xx(n)) THEN
               j=n-1
             ELSE
               j=jl
             ENDIF
             RETURN
           END SUBROUTINE locate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        SUBROUTINE interpolemoi(ii,x,xin,yin,npts,yout,ver,loga) :
!
!        Interpolation d'une valeur dans un jeu de donnees xin(npts),yin(npts).
!
!        ARGUMENTS D'ENTREE :
!                ii : ind dans xin de la valeur la plus proche de x par defaut (locate est ton ami)
!                 x : abscisse de la valeur a interpoler.
!           xin,yin : jeu de valeurs pour l'interpolation.
!              npts : nombre de points de xin,yin.
!              ver  : type d'interpolation (0 = lineaire / 1 = quadratique)
!              loga : interpolation en espace log
!
!        ARGUMENT DE SORTIE :
!           yout : valeur interpolee en x.                            !
!
!        NOTES : 
!           - Si loga est utilisee alors xin et yin doivent alors representer les logarithmes du
!             jeu de donnees. 
!           - Quelque soit la valeur de loga, yout est la valeur interpolee dans l'espace normal 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           SUBROUTINE interpolemoi(ii,x,xin,yin,npts,yout,ver,loga)
             IMPLICIT NONE
! ---------- INPUT
             INTEGER,INTENT(in)    :: npts
             INTEGER,INTENT(inout) :: ii,ver  ! ces 2 variables sont suceptible de changer
             LOGICAL,INTENT(in)    :: loga
             REAL   ,INTENT(in)    :: x,xin(npts),yin(npts)
! ---------- OUTPUT
             REAL   ,INTENT(out)   :: yout
! ---------- LOCAL
             REAL                  :: myx,ytmp,denom

! ---------- on pourrait ameliorer la condition ici et tester le cas
!            x plus proche de x(i+1) et utiliser dans ce cas l'interpolation
!            quadratique... mais est ce vraiment nécessaire ???
!            Si on est sur les 2 premiers ou les 2 derniers points :
!                  ===> interpolation lineaire
             IF (ii.eq.1.or.ii.eq.npts-1) ver = 0

!            Interpolation lineaire
             IF (ver.eq.0) THEN
               myx = x
               IF (loga) myx=alog(myx)
               denom = (xin(ii+1)-xin(ii))
               ytmp=((yin(ii+1)-yin(ii))*(myx-xin(ii)))/denom+yin(ii)
               IF (loga) THEN
                 yout=exp(ytmp)
               ELSE
                 yout=ytmp
               ENDIF
             ELSE
! ------------ Recherche de l'indice le plus proche de la valeur.
!              Permet de choisir si l'on interpole avec :
!              i-1;i;i+1   OU  i;i+1;i+2
               IF (x-xin(ii).gt.xin(ii+1)-x) ii = ii+1
               myx=x
               IF (loga) myx=alog(myx)
               ytmp = (myx-xin(ii))*(myx-xin(ii+1))   /         &
                      ((xin(ii-1)-xin(ii))                *     &
                      (xin(ii-1)-xin(ii+1)))              *     &
                      yin(ii-1)                                 &
                      +                                         &
                      (myx-xin(ii-1))*(myx-xin(ii+1)) /         &
                      ((xin(ii)-xin(ii-1))                *     &
                      (xin(ii)-xin(ii+1)))                *     &
                      yin(ii)                                   &
                      +                                         &
                      (myx-xin(ii-1))*(myx-xin(ii))   /         &
                      ((xin(ii+1)-xin(ii-1))              *     &
                      (xin(ii+1)-xin(ii)))                *     &
                      yin(ii)
               IF (loga) THEN
                 yout=exp(ytmp)
               ELSE
                 yout=ytmp
               ENDIF
             ENDIF
             RETURN
           END SUBROUTINE interpolemoi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        SUBROUTINE extrapolemoi(x,A,B,yout,loga) :
!
!        Extrapolation lineaire d'une valeur a partir de coeff A et B issus d'un jeu de donnees 
!        xin(npts), yin(npts). [ f(x) = A*x + B  OU exp(A*log(x)+B) ]
!
!        ARGUMENTS D'ENTREE :
!              x : abscisse de la valeur a extrapoler.
!              A : coeff directeur de la droite.
!              B : ordonnee a l'origine de la droite.
!           loga : interpolation en espace log
!
!        ARGUMENT DE SORTIE :
!           yout : valeur extrapolee en x.
!
!        NOTES : 
!           - Si loga est utilisee alors A et B n'ont pas la meme signification (voir forme f(x))
!           - Quelque soit la valeur de loga, yout est la valeur extrapolee dans l'espace normal 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           SUBROUTINE extrapolemoi(x,A,B,yout,loga)
             IMPLICIT NONE
! ---------- INPUT
             REAL   ,INTENT(in)  :: x,A,B
             LOGICAL,INTENT(in)  :: loga 
! ---------- OUTPUT
             REAL   ,INTENT(out) :: yout
             
             IF (loga) THEN
               yout = exp(A*alog(x)+B)
             ELSE
               yout = A*x+B
             ENDIF

           END SUBROUTINE extrapolemoi

       END MODULE optcld 
