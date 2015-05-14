/*
 * Survey.hpp
 *
 *  Created on: Feb 11, 2015
 *      Author: entropy
 */

#ifndef SURVEY_HPP_
#define SURVEY_HPP_

#include "Params.hpp"
#include "GeophyParams.hpp"
#include "OperationalParams.hpp"
#include "ReconstructParams.hpp"
#include <algorithm>
#include <fstream>
using std::ofstream;
using std::ios;
#include <iostream>


class Survey{
	/**
	 * Public methods:
	 * Indeed, I compute the optimum survey design in this public function. SA with Augmented Lagrangian method
	 * has been implemented in this function.
	 **/
	public:
			void findOptSurvey(Params* params, GeophyParams* gpp, ReconstructParams* rpp, OperationalParams* opp);
			int ** alloc2DPointer(int n1, int n2, int **data);
			float **  alloc2DPointerFloat(int n1, int n2, float **data);
			void write2Ddata(std::string fileName, int n1, int n2, int **data);
			void write1Ddata(std::string fileName, int n1, float *data);
			void writeDdata(std::string fileName, Params* params, GeophyParams* gpp, ReconstructParams* rpp , int flag );

			// Write binary file for the optimised sources' and receivers' locations ****** //

	/**
	 * Access modifier of these methods are made private to prohibit the
	 * execution of these functions from anywhere except the object of this class
	 **/
	private:

			// Allocating 4D data pointer *************** //
			int **** alloc4DPointer(int n1, int n2, int n3, int n4, int ****data);

			// Computing standardization of the data *************** //
			float standardization(float *data, float avg, int totPatch);

			// Implementing SA ************************** //
			void SA(Params* params, GeophyParams* gpp, ReconstructParams* rpp, OperationalParams* opp, int outLoopCount);

			// Find the location of station lines ************************** //
			inline int* findLinLoc(int **Lineloc, int length);

			// Find next feasible locations of sources and receivers *** //
			int findLoc(int j, int *indx, int allowedMovement, int allowedMovementY, float inter, int deciStationLine, int preLoc, int totalStationLine, int totStation, int **currStationLoc);

			// Find source locations for a source line *** //
			int * permSrcLoc(int allowedMovementSrc, int inter, int srcNum, int totalSrc);

			// Compute CMP-Offset domain **************************//
			void computeCMPOffset(Params* params, OperationalParams* opp, GeophyParams* gpp, ReconstructParams* rpp, int flag, int backFlag, int currLoc,int totRecInLine, int totSrcInLine, int currLinNum);

			// Compute 4D patches ************************** //
			void computePatches(Params* params, GeophyParams* gpp, ReconstructParams* rpp,int flag);

			// Compute binning ************************** //
			void computeBinning(Params* params, GeophyParams* gpp, ReconstructParams* rpp, int allocFlag);

			void computePatch(Params* params, GeophyParams* gpp, ReconstructParams* rpp, int flag);

			// Compute mutual coherency for all patches ****************** //
			long computeMutualCoh(std::vector<std::vector<std::vector<std::vector<int> > > >&patch4D, float *mu, int patchNum);

			// Find source locations index in source line ****************** //
			int* findSrcLoc(int **Lineloc, int length, int lineNum);

			// Deallocate memeory for 4D pointer ****************** //
			void deallocPointer(int n1, int n2, int n3, int ****data);

			// Evaluate fitness function ********************************  //
			float evalFitness(Params* params, GeophyParams* gpp, ReconstructParams* rpp, int maxRLI, float RLI, int maxSLI, float SLI, int outLoopCount);

};
#endif /* SURVEY_HPP_ */
