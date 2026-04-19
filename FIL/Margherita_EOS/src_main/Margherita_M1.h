/*
/* =====================================================================================
 *
 *       Filename:  Margherita_M1.h
 *
 *    Description:  Header file that includes all Leakage related header files
 *
 *        Version:  1.0
 *        Created:  17/05/2017 17:46:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#ifndef _MARGHERITA_M1_
#define _MARGHERITA_M1_

enum M1_interactions { BETA=0,
		       PLASMON_DECAY,
		       BREMS,
		       PAIR,
		       NINTERACTIONS };
			   
#include "margherita.hh"
#include "M1/fugacities.hh"
#include "M1/M1_constants.hh"
#include "M1/M1.hh"
#include "M1/M1_find_betaeq.hh"

#endif
