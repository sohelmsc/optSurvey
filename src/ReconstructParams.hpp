/*
 * reconstructParams.hpp
 *
 *  Created on: Feb 11, 2015
 *      Author: Sohel Bhuiyan
 */

#ifndef RECONSTRUCTPARAMS_HPP_
#define RECONSTRUCTPARAMS_HPP_

#include <cmath>
#include <time.h>
#include <vector>

/**
 * This data structure holds the parameters for reconstructional parameters and variables
 **/

class ReconstructParams{
public:
	// These parameters depend on the seismic data and velocity model
	float *mu;	    				/// Mutual coherency for all patches
	float *initMu;	    			/// Mutual coherency for all patches
	float avgMu;					/// Average mutual coherency
	std::vector<float> itrAvgMu;	/// Avergae mu for different outer iterations

    float simCriticalMu; 			/// simulated value of fractional critical mutual coherency
    float desCriticalMu; 			/// Desired value of fractional critical mutual coherency

	float desiredGridDensity;	    /// Desired grid density
	float *simGridDensity;	    	/// Simulated grid density
	float *initSimGridDensity; 		/// Initial Simulated Grid density
	float avgGridDensity; 			/// Average grid density

	float desiredCriticGridDensity;	/// Desired critical fraction of grid density
	float simCriticGridDensity;		/// Simulated critical fraction of grid density

	float desiredGridEfficiency;	/// Desired grid efficiency
	float simGridEfficiency;	    /// Simulated grid efficiency

	/******** clear configuration variables **********/
	void reset()
	{
		delete [] mu;
		delete [] simGridDensity;
	}
};
#endif /* RECONSTRUCTPARAMS_HPP_ */
