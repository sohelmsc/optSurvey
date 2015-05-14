/*
 * geophyParams.hpp
 *
 *  Created on: Feb 11, 2015
 *      Author: entropy
 */

#ifndef GEOPHYPARAMS_HPP_
#define GEOPHYPARAMS_HPP_

#include <cmath>
#include <time.h>

/**
 * This data structure holds the parameters for Geophysical parameters and variables
 **/

class GeophyParams{
public:
	// These parameters depend on the geophysical parameters
	float binX;	    			/// Bin size in crossline direction (m)
	float binY;	    			/// Bin size in inline direction (m)
	float simXmax;	    		/// Simulated Maximum offset (m)
	float simXmaxX;	    		/// Simulated Maximum offsetx (m)
	float simXmaxY;	    		/// Simulated Maximum offsety (m)
	float simXmin;	    		/// Simulated Minimum offset (m)
	float desXmax;	    		/// Desired Minimum offset (m)
	float desXmin;	    		/// Desired Minimum offset (m)

	/******** clear configuration variables **********/
	void reset()
	{
	}
};
#endif /* GEOPHYPARAMS_HPP_ */
