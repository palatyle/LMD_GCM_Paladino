!
! $Header$
!
      if (ok_regdyn) then
      
      if (is_sequential) then


      ndex3d = 0
      itau_w = itau_phy + itap
c
       CALL histwrite(nid_regdyn,"hw1",itau_w,histoW(:,:,:,1),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nh1",itau_w,nhistoW(:,:,:,1),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nht1",itau_w,nhistoWt(:,:,:,1),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"hw2",itau_w,histoW(:,:,:,2),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nh2",itau_w,nhistoW(:,:,:,2),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nht2",itau_w,nhistoWt(:,:,:,2),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"hw3",itau_w,histoW(:,:,:,3),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nh3",itau_w,nhistoW(:,:,:,3),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nht3",itau_w,nhistoWt(:,:,:,3),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"hw4",itau_w,histoW(:,:,:,4),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nh4",itau_w,nhistoW(:,:,:,4),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nht4",itau_w,nhistoWt(:,:,:,4),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"hw5",itau_w,histoW(:,:,:,5),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nh5",itau_w,nhistoW(:,:,:,5),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)
c
       CALL histwrite(nid_regdyn,"nht5",itau_w,nhistoWt(:,:,:,5),
     &               kmaxm1*lmaxm1*iwmax,ndex3d)

      if (ok_sync) then
        call histsync(nid_regdyn)
      endif

      endif ! is_sequential

      endif
