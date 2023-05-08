#-
#- $Id: AA_make.ldef 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Local parameters related to the "tools" of IOIPSL (rebuild, config)
#---------------------------------------------------------------------
SHELL = /bin/ksh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
P_P =
MODEL_LIB = $(LIBDIR)/libioipsl.a
#-
#- $Id: AA_make.gdef 4976 2020-02-25 10:06:58Z aclsce $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../../util
#-
#-Q-
######-Q- ada      F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- curie  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- cur_mono  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
#- Global definitions for gfortran, generic system
M_K = make
P_C = cpp
FCM_ARCH = gfortran 
P_O = -fpreprocessed -P -C -traditional $(P_P)
F_C = gfortran -c -cpp
F_D =
F_P = -fdefault-real-8
w_w = -O3 -funroll-all-loops $(F_D) $(F_P) -I$(MODDIR)
F_O = $(w_w) -J$(MODDIR)
F_L = gfortran
M_M = 0
L_X = 0
L_O =
A_C = ar -rs
A_G = ar -x
C_C = cc -c
C_O =
C_L = cc
#-
NCDF_INC = /apps/local/easybuild/software/netCDF-Fortran/4.4.2-gmvolf-5.5.4/include
NCDF_LIB = -L/apps/local/easybuild/software/netCDF-Fortran/4.4.2-gmvolf-5.5.4/lib -lnetcdf
#-
####-Q- lxiv8    F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- occigen  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
#- $Id: AA_make 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Create "rebuild" related to IOIPSL
#---------------------------------------------------------------------

SXMODEL_LIB = $(MODEL_LIB)
ifeq ($(L_X),1)
SXMODEL_LIB = $(LIBDIR)/libsxioipsl.a
endif

rbld: libioipsl $(BINDIR)/rebuild $(BINDIR)/flio_rbld
	@echo rebuild is OK

$(BINDIR)/rebuild: ./rebuild
	cp -p ./rebuild $(BINDIR)/rebuild
	chmod a+rx $(BINDIR)/rebuild

$(BINDIR)/flio_rbld: $(SXMODEL_LIB) flio_rbld.f90
	$(F_C) $(F_O) -I$(NCDF_INC) flio_rbld.f90
	$(F_L) flio_rbld.o -o $(BINDIR)/flio_rbld \
               $(L_O) $(SXMODEL_LIB) $(NCDF_LIB)
	chmod a+rx $(BINDIR)/flio_rbld
	$(RM) flio_rbld.o

libioipsl:
	(cd ../src; $(M_K) -f Makefile)

clean:
	$(RM) *.*~ Makefile~ core *.o *.mod i.*.L *.L i.*.f90 *.M
	$(RM) $(BINDIR)/rebuild $(BINDIR)/flio_rbld

