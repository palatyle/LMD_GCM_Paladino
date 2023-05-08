#-
#- $Id: AA_make.ldef 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Local parameters for "libioipsl"
#---------------------------------------------------------------------
SHELL = /bin/bash
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
#- $Id: AA_make 105 2007-07-26 09:46:30Z mafoipsl $
#-
MAKE_NAM = $(MAKE)
ifneq ($(MAKE_NAM),$(M_K))
 @$(error You must invoke this Makefile with the $(M_K) command)
endif
USER_DIR = $(shell /bin/pwd)
MAKE_DIR = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON/ioipsl/modipsl/modeles/IOIPSL/src
ifneq ($(USER_DIR),$(MAKE_DIR))
 @$(error You must invoke this Makefile from its directory)
endif
#-
#---------------------------------------------------------------------
#- Create "libioipsl"
#---------------------------------------------------------------------

.PRECIOUS : $(MODEL_LIB)
SXMODEL_LIB = $(MODEL_LIB)
ifeq ($(L_X),1)
SXMODEL_LIB = $(LIBDIR)/libsxioipsl.a
.PRECIOUS : $(SXMODEL_LIB)
endif

MODS1 = defprec.f90   \
        errioipsl.f90 \
        stringop.f90  \
        mathelp.f90   \
        getincom.f90  \
        calendar.f90  \
        fliocom.f90   \
        flincom.f90   \
        histcom.f90   \
        restcom.f90   \
        ioipsl.f90
OBJSMODS1 = $(MODS1:.f90=.o)

all: $(MODEL_LIB)($(OBJSMODS1))
	@echo IOIPSL is OK

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
ifeq ($(L_X),1)
	$(A_X) $(SXMODEL_LIB) $*.o
endif
ifeq ($(M_M),1)
	mv $*.mod $(MODDIR)
endif
	$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)
ifeq ($(L_X),1)
	$(RM) $(SXMODEL_LIB)
endif
	$(RM) *.*~ Makefile~ core *.o *.mod i.*.L *.L i.*.f90
	$(RM) $(MODDIR)/*.mod $(MODDIR)/*.M *.M

#- Specific dependencies

$(MODEL_LIB)(errioipsl.o):   \
  $(MODEL_LIB)(defprec.o)

$(MODEL_LIB)(stringop.o):    \
  $(MODEL_LIB)(defprec.o)

$(MODEL_LIB)(mathelp.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(getincom.o):    \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(calendar.o):    \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(fliocom.o):     \
  $(MODEL_LIB)(defprec.o)    \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(calendar.o)

$(MODEL_LIB)(flincom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(calendar.o)

$(MODEL_LIB)(histcom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(mathelp.o)    \
  $(MODEL_LIB)(calendar.o)   \
  $(MODEL_LIB)(fliocom.o)

$(MODEL_LIB)(restcom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(mathelp.o)    \
  $(MODEL_LIB)(calendar.o)   \
  $(MODEL_LIB)(fliocom.o)

$(MODEL_LIB)(ioipsl.o):      \
  $(MODEL_LIB)(getincom.o)   \
  $(MODEL_LIB)(fliocom.o)    \
  $(MODEL_LIB)(flincom.o)    \
  $(MODEL_LIB)(histcom.o)    \
  $(MODEL_LIB)(restcom.o)

