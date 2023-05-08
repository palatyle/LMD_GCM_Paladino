/* Copyright Jérémie Burgalat (2010-2015,2017)
 * 
 * jeremie.burgalat@univ-reims.fr
 * 
 * This software is a computer program whose purpose is to provide configuration 
 * file and command line arguments parsing features to Fortran programs.
 * 
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 * 
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-B license and that you accept its terms.
 */

 /** 
  * @file defined.h
  * @brief CPP macro definitions files
  * @details This header defines few CPP symbols and macros that are used 
  * in the library source code.
  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/** @def ASSIGN_DTSTR(in,out)
 *  Performs string assignment
 *
 *  This macro definition depends on compiler's support for allocatable string
 *  in derived type:
 *    - If it actually supports this feature, the macro defines an allocation
 *      statement.
 *    - Otherwise it defines a simple assignment statement.
 */
#if ! HAVE_FTNDTSTR
#define ASSIGN_DTSTR(in,out) out = in
#else
#define ASSIGN_DTSTR(in,out) ALLOCATE(out,source=in)
#endif

/** @def OBJECT(name)
 *  Derived type declaration
 *
 *  This macro definition depends on compiler's support for Bounded procedures
 *  in derived type (more precisely, Fortran 2003 PROCEDURE keyword support): 
 *    - If it actually supports this feature, the macro defines derived type
 *      declaration as dummy argument of subroutine/function using CLASS keyword.
 *    - Otherwise, derived type dummy argument are declared using TYPE keyword.
 */
#if ! HAVE_FTNPROC
#define OBJECT(name) TYPE(name)
#else
#define OBJECT(name) CLASS(name)
#endif

/* Defines SSLEN if needed */
#ifndef SSLEN
#define SSLEN 250
#endif

/* Defines SLLEN if needed */
#ifndef SLLEN
#define SLLEN 2500
#endif
