/*
 * operationalParams.hpp
 *
 *  Created on: Feb 11, 2015
 *      Author: entropy
 */

#ifndef OPERATIONALPARAMS_HPP_
#define OPERATIONALPARAMS_HPP_

#include <cmath>
#include <time.h>

/**
* This data structure holds the parameters for operational constraints and variables
**/
class OperationalParams{
public:
	// These parameters depend on the logistic and operational constraints
	int recAvail;	    		/// Available receiver
	int orgNRL;					/// Original number of receiver line
	int orgNSL;					/// Original number of source line
	int NRL;					/// Number of receiver line for SA
	int NSL;					/// Number of source line for SA
	int NR;						/// Number of receivers
	int NS;						/// Number of sources
	float RI;					/// Receiver interval
	float SI;					/// Source interval
	float initSlOffset;			/// Initial source line offset
	float RLI;					/// Receiver line interval for SA
	float SLI;					/// Source line interval for SA
	float orgRLI;					/// Original Receiver line interval
	float orgSLI;					/// Original Source line interval
	int decNRL;					/// Decimated receiver line
	int decNSL; 				/// Decimated source line
	int decNS; 					/// Decimated sources in one line
	float maxRLI;				/// maximum receiver line interval
	float maxSLI;				/// maximum source line interval
	float minRLI;				/// maximum receiver line interval
	float minSLI;				/// maximum source line interval


	/******** clear configuration variables **********/
	void reset()
	{

	}
};

#endif /* OPERATIONALPARAMS_HPP_ */
