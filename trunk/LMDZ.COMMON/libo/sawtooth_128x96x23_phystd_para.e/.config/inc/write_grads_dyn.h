!
! $Header$
!
      if (callinigrads) then

         string10='dyn'
         call inigrads(1,iip1
     s  ,rlonv,180./pi,-180.,180.,jjp1,rlatu,-90.,90.,180./pi
     s  ,llm,presnivs,1.
     s  ,dtvr*iperiod,string10,'dyn_zon ')

        callinigrads=.false.


      endif

      string10='ps'
      CALL wrgrads(1,1,ps,string10,string10)

      string10='u'
      CALL wrgrads(1,llm,unat,string10,string10)
      string10='v'
      CALL wrgrads(1,llm,vnat,string10,string10)
      string10='teta'
      CALL wrgrads(1,llm,teta,string10,string10)
      do iq=1,nqtot
         string10='q'
         write(string10(2:2),'(i1)') iq
         CALL wrgrads(1,llm,q(:,:,iq),string10,string10)
      enddo

