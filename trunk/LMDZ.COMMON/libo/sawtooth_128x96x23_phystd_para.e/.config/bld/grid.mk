# Automatic Make rule for grid

SRCDIR0__grid = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/libf/grid

dimensions.h: \
          $(SRCDIR0__grid)/dimensions.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

dimensions.h.idone: \
          $(SRCDIR0__grid)/dimensions.h
	touch $(FCM_DONEDIR)/$@

fxyprim.h: \
          $(SRCDIR0__grid)/fxyprim.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

fxyprim.h.idone: \
          $(SRCDIR0__grid)/fxyprim.h
	touch $(FCM_DONEDIR)/$@

fxy_new.h: \
          $(SRCDIR0__grid)/fxy_new.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

fxy_new.h.idone: \
          $(SRCDIR0__grid)/fxy_new.h
	touch $(FCM_DONEDIR)/$@

fxy_sin.h: \
          $(SRCDIR0__grid)/fxy_sin.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

fxy_sin.h.idone: \
          $(SRCDIR0__grid)/fxy_sin.h
	touch $(FCM_DONEDIR)/$@

fxy_reg.h: \
          $(SRCDIR0__grid)/fxy_reg.h
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

fxy_reg.h.idone: \
          $(SRCDIR0__grid)/fxy_reg.h
	touch $(FCM_DONEDIR)/$@

