! Copyright 2013-2015,2017 Université de Reims Champagne-Ardenne 
! Contributor: J. Burgalat (GSMA, URCA)
! email of the author : jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute
! microphysics processes using a two-moments scheme.
! 
! This library is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.

!! file: mm_mprec.F90
!! summary: Library floating point precision module.
!! author: J. Burgalat
!! date: 2013-2015,2017

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef PREC
#define PREC 64 
#elif (PREC != 32 && PREC != 64 && PREC != 80)
#undef PREC
#define PREC 64
#endif

MODULE MM_MPREC
  !! Library floating point computations precision module.
  !!
  !! This module only defines a single variable [[mm_mprec(module):mm_wp(variable)]] which sets
  !! the kind of floating point value used in all other part of the library (REAL(kind=mm_wp) 
  !! declaration statement).
  IMPLICIT NONE

#if (PREC == 32)
  !> Size of floating point variables in the library (single).
  INTEGER, PUBLIC, PARAMETER :: mm_wp = SELECTED_REAL_KIND(p=6)  ! 32 bits
#elif (PREC == 64)
  !> Size of floating point variables in the library (double).
  INTEGER, PUBLIC, PARAMETER :: mm_wp = SELECTED_REAL_KIND(p=15) ! 64 bits 
#elif (PREC == 80)
  !> Size of floating point variables in the library (extended-double).
  INTEGER, PUBLIC, PARAMETER :: mm_wp = SELECTED_REAL_KIND(p=18) ! 80 bits
#endif
END MODULE MM_MPREC

