/*  =================================================================================
 *  Name        : Survey.cpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Execute methods to find optimum survey
 *  Date        : Feb 12, 2015
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 *  =================================================================================
 */

#include <cmath>
#include <algorithm>
using std::sort;
using std::random_shuffle;
#include <stdlib.h>
#include <vector>
#include <iostream>
using std::cout;
using std::ofstream;
using std::endl;
#include <fstream>
#include <string>
#include <cstdlib>      // std::rand, std::srand
#include <complex.h>
#include "fftw3.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "Params.hpp"
#include "GeophyParams.hpp"
#include "OperationalParams.hpp"
#include "ReconstructParams.hpp"
#include "Survey.hpp"

#define PI 3.1415926535897932

std::ofstream logOut ("logFile.csv");

/***** Core method to optimise 3D seismic survey *****************/
void Survey::findOptSurvey(Params* params, GeophyParams* gpp, ReconstructParams* rpp, OperationalParams* opp)
{
	int j=0, totSrcInLine, totRecInLine, CMPOffsetFlag=1;
	int count=0;
	std::vector<int> indx;

	float alpha = 4.0, beta = 10.0, epsilon = 0.01, tempMax1, tempMax2, eqConstraint[3], uneqConstraint[3], kNew, k = 100000.0;
	clock_t start, finish;
	params->firstLocRecTem = alloc2DPointer(params->totSrc, 2, params->firstLocRecTem);
	totSrcInLine = params->totSrc/opp->decNSL;
	totRecInLine = params->totRec/opp->decNRL;
	computeCMPOffset(params, opp, gpp, rpp, CMPOffsetFlag, -1, -1, totRecInLine, totSrcInLine, 10000);
	computeBinning(params, gpp, rpp, 0);
	computePatch(params, gpp, rpp, 0);
	params->fitness = evalFitness(params, gpp, rpp, opp->maxRLI, opp->RLI, opp->maxSLI, opp->SLI, count);
	params->finalCons.push_back(k);
	rpp->itrAvgMu.push_back(rpp->avgMu);

	while(1){
		if(count>=18)
			break;
		indx.clear();
		start = clock();
		if(count > 0)
			params->initTemp = params->initTemp * powf(0.99, (1 + params->numCycle));

		SA(params, gpp, rpp, opp, j);
		finish = clock();

		logOut << "Time: " << (finish-start)/double(CLOCKS_PER_SEC) << " Seconds ";
		eqConstraint[0] = rpp->simCriticalMu;
		eqConstraint[1] = rpp->simCriticGridDensity;
		eqConstraint[2] = params->normGridEfficiency;
		tempMax1 = *std::max_element(eqConstraint, eqConstraint + 3);
		uneqConstraint[0] =(params->tempNormMin > -1.0*params->theta[j][2])?fabs(params->tempNormMin):fabs(params->theta[j][2]);
		uneqConstraint[1] =(params->tempNormMax > -1.0*params->theta[j][3])?fabs(params->tempNormMax) : fabs(params->theta[j][3]);
		uneqConstraint[2] =(params->tempGridDensity > -1.0*params->theta[j][5])?fabs(params->tempGridDensity):fabs(params->theta[j][5]);

		tempMax2 = *std::max_element(uneqConstraint, uneqConstraint + 3);
		kNew = (tempMax1 > tempMax2) ? tempMax1 : tempMax2;

		if(uneqConstraint[0] > k/alpha)
			indx.push_back(2);
		if(uneqConstraint[1] > k/alpha)
			indx.push_back(3);
		if(uneqConstraint[2] > k/alpha)
			indx.push_back(5);
		if(eqConstraint[0] > k/alpha)
			indx.push_back(0);
		if(eqConstraint[1] > k/alpha)
			indx.push_back(1);
		if(eqConstraint[2] > k/alpha)
			indx.push_back(4);

		if(kNew <= epsilon)
			break;
		else
		{
			if (kNew >= k)
			{
				for (unsigned int i = 0; i < indx.size(); i++)
				{
					params->gamma[indx[i]] = params->gamma[indx[i]]*beta;
					params->theta[j][indx[i]] = params->theta[j][indx[i]]/beta;
				}
			}
			else
			{
				for (unsigned int i = 0; i < indx.size(); i++)
				{
					if (indx[i] == 2)
						params->theta[j+1][indx[i]] = params->theta[j][indx[i]]+(params->tempNormMin> -1.0*params->theta[j][indx[i]])?(params->tempNormMin):(-params->theta[j][indx[i]]);
					else if (indx[i] == 3)
						params->theta[j+1][indx[i]] = params->theta[j][indx[i]]+ (params->tempNormMax> -1.0*params->theta[j][3])?(params->tempNormMax):(-params->theta[j][3]);
					else if (indx[i] == 5)
						params->theta[j+1][indx[i]] = params->theta[j][indx[i]]+ (params->tempGridDensity	> -1.0*params->theta[j][5]) ?(params->tempGridDensity) :(-params->theta[j][5]);
					else if (indx[i] == 0)
						params->theta[j+1][indx[i]] = params->theta[j][indx[i]]+ rpp->simCriticalMu;
					else if (indx[i] == 1)
						params->theta[j+1][indx[i]] = params->theta[j][indx[i]]+ rpp->simCriticGridDensity;
					else
						params->theta[j+1][indx[i]] = params->theta[j][indx[i]]+ params->normGridEfficiency;
				}

				if(kNew <= (k / alpha))
				{
					k = kNew;
					j++;
				}
				else
				{
					for (unsigned int i = 0; i < indx.size(); i++)
					{
						params->gamma[indx[i]] = params->gamma[indx[i]]*beta;
						params->theta[j + 1][indx[i]] = params->theta[j+1][indx[i]]/beta;
					}
					k = kNew;
					j++;
				}
			}
		}
		params->fitness = evalFitness(params, gpp, rpp, opp->minRLI, opp->RLI, opp->minSLI, opp->SLI, count);
		logOut<< k ;
		logOut <<"\n";
		params->finalCons.push_back(k);
		rpp->itrAvgMu.push_back(rpp->avgMu);

		logOut << "Outer Loop: " << j << std::endl;
		logOut << "Count: " << count << std::endl;
		count++;

		if(j>=18)
			break;
	}

	// To deallocate memory space ***************************//
	deallocPointer(params->binDimX, params->binDimY, params->offsetxDim, params->binning);
	params->binning = NULL;
}

/***** Implementing SA ***********************************************************/
void Survey::SA(Params* params, GeophyParams* gpp, ReconstructParams* rpp,	OperationalParams* opp, int outLoopCount)
{
	float currTemp, r, diffFitness, fitness;
	int numCycle = 1, *indxSrc, *indxRec, preLoc, currLoc, totSrcInLine,totRecInLine, *indxSrcLoc, *preIndxSrcLoc;
	int CMPOffsetFlag = 1;

	//** Initial source and receiver locations ***********************************/
	indxRec = findLinLoc(params->rec,opp->NRL);
	indxSrc = findLinLoc(params->src,opp->NSL);
	totSrcInLine = params->totSrc/opp->decNSL;
	totRecInLine = params->totRec/opp->decNRL;

	//End of Data Initialisation *************************************************/
	params->convergenceSAconstraint[0] = params->fitness; // cost_function;

	while(numCycle <= params->numCycle)
	{
		logOut << numCycle << std::endl;
		currTemp = params->initTemp * powf(0.99, (1 + numCycle));

		// Looping over number of source line ************************************/
		for (int j=0; j<opp->decNSL; j++)
		{
			if(params->srcLineDecimationRate < 1.0 && params->srcDecimationRate < 1.0)
			{
				// Find the previous source line location ******/
				preLoc = indxSrc[j];

				// Find the current source position in a particular source line before changing position ***/
				preIndxSrcLoc = findSrcLoc(params->src,ceil(params->surveySizeY/opp->SI),preLoc);

				// Rearrange the source position in a particular source line ******************************/
				indxSrcLoc = permSrcLoc(params->allowedMovementSrcX, params->intervalSrc, (totSrcInLine - 1), ceil(params->surveySizeY/opp->SI));

				// Find the new source line location ******/
				currLoc = findLoc(j, indxSrc, params->allowedMovementSrcLine, params->allowedMovementSrcY, params->intervalSrcLine, (opp->decNSL - 1), preLoc, opp->NSL, totSrcInLine, params->currSrcLoc);

				// IF current location is changed *********************************/
				indxSrc[j] = currLoc;
				if (currLoc != preLoc)
				{
					if (j > 0)
						opp->maxSLI = ((indxSrc[j] - indxSrc[j-1]) > opp->maxSLI) ?(indxSrc[j] - indxSrc[j - 1]) : opp->maxSLI;
						//opp->minSLI = ((indxRec[j] - indxRec[j-1]) < opp->minSLI) ? (indxRec[j] - indxRec[j-1]) : opp->minSLI;
					for (int k = 0; k < totSrcInLine; k++)
					{
						params->src[preLoc][preIndxSrcLoc[k]] = 0;
						params->src[currLoc][indxSrcLoc[k]] = 1;
						params->srcLoc[j * totSrcInLine +k][0] = currLoc;
						params->srcLoc[j * totSrcInLine +k][1] = indxSrcLoc[k];
					}
				}
				else
				{
					for (int k = 0; k < totSrcInLine; k++)
					{
						params->src[currLoc][preIndxSrcLoc[k]] = 0;
						params->src[currLoc][indxSrcLoc[k]] = 1;
						params->srcLoc[j * totSrcInLine+k][1] = indxSrcLoc[k];
					}
				}
				CMPOffsetFlag = 0;
				computeCMPOffset(params, opp, gpp, rpp,CMPOffsetFlag, -1, currLoc,	totRecInLine, totSrcInLine, j);

				//fitness = evalFitness(params, gpp, rpp,opp->maxRLI, opp->RLI, opp->maxSLI, opp->SLI, outLoopCount);
				fitness = evalFitness(params, gpp, rpp,opp->minRLI, opp->RLI, opp->minSLI, opp->SLI, outLoopCount);
				diffFitness = fitness - params->fitness;

				// Metropolis criteria  ******************//
				r = ((double) rand()/(RAND_MAX));
				if (diffFitness < 0 || (r < (expf(-diffFitness / currTemp))))
					params->fitness = fitness;
				else {
					indxSrc[j] = preLoc;
					for (int k = 0; k < totSrcInLine; k++)
					{
						params->src[preLoc][preIndxSrcLoc[k]] = 1;
						params->src[currLoc][indxSrcLoc[k]] = 0;
						params->srcLoc[j * totSrcInLine +k][0] = preLoc;
						params->srcLoc[j * totSrcInLine +k][1] = preIndxSrcLoc[k];
					}
				}
			}
			else if(params->srcDecimationRate < 1.0)
			{
				// Find the current source position in a particular source line before changing position ***/
				preIndxSrcLoc = findSrcLoc(params->src,ceil(params->surveySizeY/opp->SI),indxSrc[j]);

				// Rearrange the source position in a particular source line ******/
				indxSrcLoc = permSrcLoc(params->allowedMovementSrcX, params->intervalSrc, (totSrcInLine - 1), ceil(params->surveySizeY/opp->SI));

				for (int k=0; k<totSrcInLine;k++)
				{
					params->src[j][preIndxSrcLoc[k]] = 0;
					params->src[j][indxSrcLoc[k]] = 1;
					params->srcLoc[j * totSrcInLine +k][1] = indxSrcLoc[k];
				}

				CMPOffsetFlag = 0;
				computeCMPOffset(params, opp, gpp, rpp,CMPOffsetFlag, -1,currLoc,	totRecInLine, totSrcInLine, j);

				//fitness = evalFitness(params, gpp, rpp,opp->maxRLI, opp->RLI, opp->maxSLI, opp->SLI, outLoopCount);
				fitness = evalFitness(params, gpp, rpp,opp->minRLI, opp->RLI, opp->minSLI, opp->SLI, outLoopCount);
				diffFitness = fitness - params->fitness;

				// Metropolis criteria  ******************//
				r = ((double) rand() / (RAND_MAX));
				if (diffFitness < 0 || (r < (expf(-diffFitness / currTemp))))
					params->fitness = fitness;
				else
				{
					indxSrc[j] = preLoc;
					for (int k=0; k<totSrcInLine; k++)
					{
						params->src[j][preIndxSrcLoc[k]] = 1;
						params->src[j][indxSrcLoc[k]] = 0;
						params->srcLoc[j*totSrcInLine+k][1] = preIndxSrcLoc[k];
					}
				}
			}
			else
			{
				preLoc = indxSrc[j];
				currLoc = findLoc(j, indxSrc, params->allowedMovementSrcLine, params->allowedMovementSrcY, params->intervalSrcLine, (opp->decNSL - 1), preLoc, opp->NSL, totSrcInLine, params->currSrcLoc);

				// IF current location is changed *************************************/
				indxSrc[j] = currLoc;
				if (currLoc != preLoc)
				{
					if(j > 0)
						opp->maxSLI = ((indxSrc[j] - indxSrc[j-1]) > opp->maxSLI) ?(indxSrc[j] - indxSrc[j - 1]) : opp->maxSLI;
						//opp->minSLI = ((indxRec[j] - indxRec[j-1]) < opp->minSLI) ? (indxRec[j] - indxRec[j-1]) : opp->minSLI;
						if(j == 0)
						{
							for (int k = 0; k < totSrcInLine; k++)
							{
								params->src[preLoc][k] = 0;
								params->src[currLoc][k] = 1;
								params->srcLoc[j * totSrcInLine + k][0] = currLoc;
							}
						}
						else
						{
							for (int k = 0; k < totSrcInLine; k++)
							{
								params->src[preLoc+params->preSrcLoc[j][k]][k] = 0;
								params->src[currLoc+params->currSrcLoc[j][k]][k] = 1;
								params->srcLoc[j * totSrcInLine + k][0] = currLoc + params->currSrcLoc[j][k];
							}
						}
					}
					CMPOffsetFlag = 0;
					computeCMPOffset(params, opp, gpp, rpp,CMPOffsetFlag, -1,currLoc,	totRecInLine, totSrcInLine, j);
					//fitness = evalFitness(params, gpp, rpp,opp->maxRLI, opp->RLI, opp->maxSLI, opp->SLI, outLoopCount);
					fitness = evalFitness(params, gpp, rpp,opp->minRLI, opp->RLI, opp->minSLI, opp->SLI, outLoopCount);
					diffFitness = fitness - params->fitness;

					// Metropolis criteria  ******************//
					r = ((double) rand() / (RAND_MAX));
					if (diffFitness < 0 || (r < (expf(-diffFitness / currTemp))))
					{
						params->fitness = fitness;
						if(j > 0)
							for(int k=0; k<totSrcInLine; k++)
								params->preSrcLoc[j][k] = params->currSrcLoc[j][k];
					}
					else
					{
						indxSrc[j] = preLoc;
						if(j > 0)
						{
							for (int k = 0; k < totSrcInLine; k++)
							{
								params->src[preLoc+params->preSrcLoc[j][k]][k] = 1;
								params->src[currLoc+params->currSrcLoc[j][k]][k] = 0;
								params->srcLoc[j * totSrcInLine + k][0] = preLoc+params->preSrcLoc[j][k];
							}
						}
						else
						{
							for (int k = 0; k < totSrcInLine; k++)
							{
								params->src[preLoc][k] = 1;
								params->src[currLoc][k] = 0;
								params->srcLoc[j * totSrcInLine + k][0] = preLoc;
							}
						}
						CMPOffsetFlag = 0;
						computeCMPOffset(params, opp, gpp, rpp, CMPOffsetFlag, 3, currLoc,	totRecInLine, totSrcInLine, j);
					}
				} // End of else
    		} // End of for loop

		if(params->recDecimationRate < 1.0)
		{
			// Looping over number of receiver line *********//
			for(int j=0; j<opp->decNRL; j++)
			{
				preLoc = indxRec[j];
				currLoc = findLoc(j, indxRec, params->allowedMovementRec, 0, params->intervalRec, (opp->decNRL - 1), preLoc, opp->NRL, totRecInLine, params->currRecLoc);

				// IF current location is changed ************//
				indxRec[j] = currLoc;
				if (currLoc != preLoc) {
					if(j>0)
						//opp->maxRLI = ((indxRec[j] - indxRec[j-1]) > opp->maxRLI) ? (indxRec[j] - indxRec[j-1]) : opp->maxRLI;
						opp->minRLI = ((indxRec[j] - indxRec[j-1]) < opp->minRLI) ? (indxRec[j] - indxRec[j-1]) : opp->minRLI;
					if(j>0)
					{
						for (int k=0; k < totRecInLine; k++)
						{
							params->rec[preLoc+params->preRecLoc[j][k]][k] = 0;
							params->rec[currLoc+params->currRecLoc[j][k]][k] = 1;
							params->recLoc[j * totRecInLine + k][0] = currLoc+params->currRecLoc[j][k];
						}
					}
					else
					{
						for (int k=0; k < totRecInLine; k++)
						{
							params->rec[preLoc][k] = 0;
							params->rec[currLoc][k] = 1;
							params->recLoc[j * totRecInLine + k][0] = currLoc;
						}
					}
					CMPOffsetFlag = 2;
					computeCMPOffset(params, opp, gpp, rpp, CMPOffsetFlag, -1, currLoc,	totRecInLine, totSrcInLine, j);
					//fitness = evalFitness(params, gpp, rpp, opp->maxRLI, opp->RLI, opp->maxSLI, opp->SLI, outLoopCount);
					fitness = evalFitness(params, gpp, rpp,opp->minRLI, opp->RLI, opp->minSLI, opp->SLI, outLoopCount);
					diffFitness = fitness - params->fitness;

					// Metropolis criteria *****************************//
					r = ((double) rand() / (RAND_MAX));
					if (diffFitness < 0 || (r < (expf(-diffFitness / currTemp))))
					{
						params->fitness = fitness;
						for(int k=0; k<totRecInLine; k++)
							params->preRecLoc[j][k] = params->currRecLoc[j][k];
					}
					else
					{
						indxRec[j] = preLoc;
						if(j>0)
							for (int k = 0; k < totRecInLine; k++)
							{
								// Backward the cmp and offset, and patch *******//
								params->rec[preLoc+params->preRecLoc[j][k]][k] = 1;
								params->rec[currLoc+params->currRecLoc[j][k]][k] = 0;
								params->recLoc[j * totRecInLine + k][0] = preLoc+params->preRecLoc[j][k];
							}
						else
							for (int k = 0; k < totRecInLine; k++)
							{
								// Backward the cmp and offset, and patch *******//
								params->rec[preLoc][k] = 1;
								params->rec[currLoc][k] = 0;
								params->recLoc[j * totRecInLine + k][0] = preLoc;
							}

						CMPOffsetFlag = 2;
						computeCMPOffset(params, opp, gpp, rpp, CMPOffsetFlag, 3, currLoc,	totRecInLine, totSrcInLine, j);
					}
				}
			} // End of for loop
		} // End of IF condition
		params->convergenceSAconstraint[numCycle] = params->fitness;
		numCycle=numCycle+1;
	}
	params->convergenceSAunconstraint =	params->convergenceSAconstraint[numCycle-1];
}

/**** Find the next location of sources and receivers  */
int Survey::findLoc(int j, int *indx, int allowedMovement, int allowedMovementY, float inter, int deciStationLine, int preLoc, int totalStationLine, int totStation, int **currStationLoc)
{
	// Defining whether the perturbation will be decimal or integer number ***********//
	int numDigit = 0, count = 0, currLoc, flag, loc;
	float perturb, r;
	r = ((float) rand() / (RAND_MAX));

	// Compute the current receiver location ****************************//
	if (j == 0) {
		perturb = round((r * allowedMovement / 2 * pow(10, numDigit)))/ pow(10, numDigit);
		currLoc = inter * (j) + perturb;
	} else if (j > 0) {
		perturb = round((r * allowedMovement * pow(10, numDigit)))/ pow(10, numDigit);
		currLoc = floor(inter * (j)) + perturb - floor(allowedMovement / 2);
	} else if (j == deciStationLine) {
		perturb = round((r * allowedMovement * pow(10, numDigit))) / pow(10, numDigit);
		currLoc = floor(inter * (j)) + perturb - floor(allowedMovement / 2);
	}

	for(int i=0; i<totStation; i++)
	{
		r = ((float) rand() / (RAND_MAX));
		perturb = round((r * allowedMovementY * pow(10, numDigit)))/ pow(10, numDigit);
		r = ((float) rand() / (RAND_MAX));
		if(r<=0.5)
			perturb *=-1;
		currStationLoc[j][i] = perturb;
	}


	// Infinity loop until a new location is found ************************
	while (1) {
		flag = 0;
		for (int k = 0; k <= deciStationLine; k++)
			if (indx[k] == currLoc) {
				flag = 1;
				break;
			}

		if (flag == 0 && currLoc < totalStationLine)
			break;
		else {
			r = ((float) rand() / (RAND_MAX));
			if (j == 0) {
				perturb = round((r * allowedMovement / 2 * pow(10, numDigit)))/ pow(10, numDigit);
				currLoc = inter * (j) + perturb;
			} else if (j > 0) {
				perturb = round((r * allowedMovement * pow(10, numDigit)))/ pow(10, numDigit);
				currLoc = floor(inter * (j)) + perturb - floor(allowedMovement / 2);
			} else if (j == sizeof(indx) / sizeof(int) - 1) {
				perturb = round((r * allowedMovement * pow(10, numDigit)))	/ pow(10, numDigit);
				currLoc = floor(inter * (j)) + perturb - floor(allowedMovement / 2);
			}
			count = count + 1;
			if (count > 20) {
				currLoc = preLoc;
				break;
			}
		}
	}
	return currLoc;
}

int * Survey::permSrcLoc(int allowedMovementSrc, int inter, int srcNum , int totalSrc)
{
	// Defining whether the perturbation will be decimal or integer number ***********//
	int numDigit = 0, currLoc, *p;
	float perturb, r;

	p = (int *) malloc(sizeof(int) * (srcNum+1));

	// Compute the current source location ****************************//
	for(int j=0; j<=srcNum; j++)
	{
		r = ((float) rand() / (RAND_MAX));
		if (j == 0) {
			perturb = round((r * allowedMovementSrc / 2 * pow(10, numDigit)))/ pow(10, numDigit);
			currLoc = inter * (j) + perturb;
		} else if (j > 0) {
			perturb = round((r * allowedMovementSrc * pow(10, numDigit)))/ pow(10, numDigit);
			currLoc = inter * (j) + perturb - floor(allowedMovementSrc / 2);
		} else if (j == srcNum)
		{
			perturb = round((r * allowedMovementSrc * pow(10, numDigit)))/ pow(10, numDigit);
			currLoc = inter * (j) + perturb - floor(allowedMovementSrc / 2);
		}
		if(currLoc >= totalSrc)
			p[j] = totalSrc-1;
		else
			p[j] = currLoc;
	}

	return p;
}

//** Find critical portion of data set *********************************************/
float Survey::standardization(float *data, float avg, int totPatch)
{
	int i, count = 0;
	float sum = 0.0, criticalData, stdData, temp;

	for(i = 0; i < totPatch; i++)
		sum += powf((data[i] - avg), 2);

	stdData = sqrtf((sum / totPatch));

	for (i = 0; i < totPatch; i++) {
		temp = ((data[i] - avg) / stdData) - 1;
		if (temp > 0.0)
			count++;
	}
	criticalData = (float) (count/totPatch);

	return criticalData;
}

//** Find the location of source or receiver lines *********************************************/
inline int* Survey::findLinLoc(int **Lineloc, int length) {
	int *p, i, count = 0;
	p = (int *) malloc(sizeof(int) * length);
	for (i = 0; i < length; i++)
		if (Lineloc[i][0] == 1) {
			p[count] = i;
			count++;
		}
	return p;
}

//** Find sources in one line *************************************************************/
inline int* Survey::findSrcLoc(int **Lineloc, int length, int lineNum) {

	int *p, i, count = 0;
	p = (int *) malloc(sizeof(int) * length);
	for (i = 0; i < length; i++)
		if (Lineloc[lineNum][i] == 1) {
			p[count] = i;
			count++;
		}
	return p;
}

//*** Compute CMPOffset domain ************************************************************/
void Survey::computeCMPOffset(Params* params, OperationalParams* opp, GeophyParams* gpp, ReconstructParams* rpp, int flag, int backFlag, int currLoc,int totRecInLine, int totSrcInLine, int currLinNum)
{
	int i, preRecLoc, count = 0, countSalvo, chkSalvoY, chkSalvoX,countSrcLine = 0, maxPatchNum=-10;
	int temMovementX, tempLoc, indx, imx, imy, ihx, ihy, endWinX, endWinY, endOffsetX,endOffsetY, minPatchNum = 100000;
	float offsetx, offsety, radtodeg = 180/M_PI, azimuth, sumMu = 0.0, sumGridDensity = 0.0;
	long int countPopBin = 0;

	//** For all sources in the acquisition design. This is done only once at the beginning.
	if(flag == 1)
	{
		for(i=0; i<params->totSrc; i++)
		{
			params->cmpX.push_back(std::vector<std::vector<float> >());
			count = 0;

			//** Find the starting location of receiver template for the first source of a source line **//
			if (params->templateExist == 1)
			{
				if (i % (totSrcInLine) == 0)
				{
					countSalvo = 0;
					if(countSrcLine == 0)
					{
						params->firstLocRecTem[i][0] = params->recLoc[0][0]; /// Store the receiver line
						params->firstLocRecTem[i][1] = params->recLoc[0][1]; /// Store the receiver number
					}
					//*** shifting receiver template towards right direction ***************************
					else if ((params->firstLocRecTem[i-1][1] * opp->RI+ params->tempLenX / 2)	< params->srcLoc[i][0] * opp->SLI)
					{
						temMovementX = ceil((params->srcLoc[i][0] * opp->SLI- (params->firstLocRecTem[i-1][1]	* opp->RI + params->tempLenX / 2))/ opp->RI) + 1;
						chkSalvoX = (params->firstLocRecTem[i - 1][1] + temMovementX+ params->tempNR - totRecInLine);
						params->firstLocRecTem[i][1] =(params->firstLocRecTem[i-1][1] + temMovementX)- (chkSalvoX > 0 ? chkSalvoX : 0);/// Store the receiver number

						// For up going source direction **************************
						if (!(countSrcLine%2))
							params->firstLocRecTem[i][0] = params->rec[0][0];			/// Store the receiver line number

						// For down going source direction ************************
						else
							params->firstLocRecTem[i][0] = opp->decNRL - 1;				/// Store the receiver line number

					} else {
						params->firstLocRecTem[i][1] = params->firstLocRecTem[i - 1][1];/// Store the receiver number

						// For up going source direction **************************
						if (!(countSrcLine % 2))
							params->firstLocRecTem[i][0] = params->rec[0][0];			/// Store the receiver line number

						// For down going source direction ************************
						else
							params->firstLocRecTem[i][0] = opp->decNRL-1;	 			/// Store the receiver line number
					}
					countSrcLine++;
				}

				//** Find the starting location of receiver template  ********
				else if ((countSalvo + 1) >= params->salvo)
				{
					// For up going source direction *************************
					if (countSrcLine % 2)
					{
						chkSalvoY = (params->firstLocRecTem[i-1][0] +1+ params->tempNRL - opp->decNRL);
						if (chkSalvoY <= 0)
							params->firstLocRecTem[i][0] = params->firstLocRecTem[i-1][0] + 1;/// Store the receiver number line
						else
							params->firstLocRecTem[i][0] = params->firstLocRecTem[i-1][0];	/// Store the receiver number line
					}

					// For down going source direction **************************
					else {
						chkSalvoY = (params->firstLocRecTem[i - 1][0]- params->tempNRL + 1);
						if (chkSalvoY > 0)
							params->firstLocRecTem[i][0] = params->firstLocRecTem[i- 1][0] - 1;/// Store the receiver number line
						else
							params->firstLocRecTem[i][0] = params->firstLocRecTem[i- 1][0];	/// Store the receiver number line
					}
					params->firstLocRecTem[i][1] = params->firstLocRecTem[i-1][1];/// Store the receiver number
					countSalvo = 0;
				}
				else
				{
					params->firstLocRecTem[i][0] = params->firstLocRecTem[i-1][0];/// Store the receiver number
					params->firstLocRecTem[i][1] = params->firstLocRecTem[i-1][1];/// Store the receiver number
					countSalvo++;
				}
			}
			else  // If every source is listened by every receiver
			{
				params->firstLocRecTem[i][0] = params->recLoc[0][0];	/// Store the receiver line
				params->firstLocRecTem[i][1] = params->recLoc[0][1];	/// Store the receiver number
			}

			// Traversing the receiver template from lower to higher receivers ******************************//
			for (int j = 0; j < params->tempNRL; j++)
			{
				preRecLoc = j * totRecInLine;

				if(params->templateExist == 1)
				{
					// For up going source direction **************************
					if (countSrcLine % 2)
						tempLoc = (params->firstLocRecTem[i][0] + j);

					// For down going source direction **************************
					else
						tempLoc = (params->firstLocRecTem[i][0] - j);
				}
				else
					 tempLoc = params->recLoc[preRecLoc][0];

				for (int k = 0; k < params->tempNR; k++)
				{
					params->cmpX[i].push_back(std::vector<float>()); // i is the index for sources ***************//

					params->cmpX[i][count].push_back(((params->srcLoc[i][0]*opp->SLI + opp->initSlOffset) + (params->recLoc[preRecLoc+k][1])*opp->RI)/2); // CMPX

					params->cmpX[i][count].push_back(params->recLoc[preRecLoc+k][0]);       // RLN

					params->cmpX[i][count].push_back(params->recLoc[preRecLoc+k][1]);	  	// RN

					params->cmpXVal.push_back(((params->srcLoc[i][0] * opp->SLI + opp->initSlOffset) + (params->recLoc[preRecLoc+k][1])*opp->RI)/2);  // CMP-X

					params->cmpY.push_back((params->srcLoc[i][1] * opp->SI+ (params->recLoc[preRecLoc + k][0]) * opp->RLI)/2);  // CMP-Y

					offsetx = fabs((params->recLoc[preRecLoc + k][1]) * opp->RI - (params->srcLoc[i][0] * opp->SLI+ opp->initSlOffset));

					offsety = fabs((params->recLoc[preRecLoc + k][0]) * opp->RLI - params->srcLoc[i][1] * opp->SI);

                    params->offSetX.push_back(offsetx);  // hx = Rx - Sx

                    params->offSetY.push_back(offsety);  // hy =  Ry - Sy

					params->offSet.push_back(sqrtf(powf(offsetx, 2) + powf(offsety, 2)));	//offset

					azimuth = (atan2(offsety, offsetx) * radtodeg);

					params->aziMuth.push_back(azimuth < 0 ? azimuth + 360 : azimuth); // Azimuth

					count++;
				}
			}
		}
	}

	//** For the perturbation of one source line in the acquisition design *********************//
	else if(flag == 0)
	{
		if(currLinNum > 0)
		{
			if(params->templateExist == 1)
			{
				if((params->firstLocRecTem[currLinNum * totSrcInLine - 1][1]*opp->RI + params->tempLenX / 2)< params->srcLoc[currLinNum * totSrcInLine][0] * opp->SLI)
				{
					temMovementX = ceil((params->srcLoc[currLinNum * totSrcInLine][0] * opp->SLI- (params->firstLocRecTem[currLinNum* totSrcInLine - 1][1] * opp->RI+ params->tempLenX / 2)) / opp->RI) + 1;
					chkSalvoX = (params->firstLocRecTem[currLinNum * totSrcInLine- 1][1] + temMovementX + params->tempNR - totRecInLine);
					params->firstLocRecTem[currLinNum * totSrcInLine][1] = (params->firstLocRecTem[currLinNum * totSrcInLine - 1][1]+ temMovementX)- (chkSalvoX > 0 ? chkSalvoX : 0);/// Store the receiver number
				}
			}
		}

		//*** Traversing all sources in that source line *******************************************//
		for(int m = 0; m < totSrcInLine; m++)
		{
			count = 0;
			//*** Traversing all receiver line for one source **************************************//
			for (int j = 0; j < params->tempNRL; j++)
			{
				if(params->templateExist == 1)
				{
					// For up going template movement *****************************************//
					if (!(currLinNum % 2))
						tempLoc =(params->firstLocRecTem[currLinNum * totSrcInLine][0]+ j);    // Current receiver line number

					// For down going template movement ********************************//
					else
						tempLoc = (params->firstLocRecTem[currLinNum * totSrcInLine][0]- j);   // Current receiver line number
				}

				preRecLoc = j * totRecInLine;

				//***** Traversing every receiver in that receiver line ***************//
				for(int k = 0; k<params->tempNR; k++)
				{
					ihx = floor(params->offSetX[currLinNum * totSrcInLine * params->tempNRL* params->tempNR + m * (params->tempNRL * params->tempNR) + count] / params->dhx);
					ihy = floor(params->offSetY[currLinNum * totSrcInLine * params->tempNRL* params->tempNR + m * (params->tempNRL * params->tempNR) + count] / params->dhy);
					imx = floor(params->cmpXVal[currLinNum * totSrcInLine * params->tempNRL* params->tempNR	+ m * (params->tempNRL * params->tempNR) + count] / gpp->binX);
					imy = floor(params->cmpY[currLinNum * totSrcInLine * params->tempNRL* params->tempNR	+ m * (params->tempNRL * params->tempNR) + count] / gpp->binY);

					//iaz = floor(params->aziMuth[currLinNum * totSrcInLine * params->tempNRL* params->tempNR	+ m * (params->tempNRL * params->tempNR) + count] / params->daz);
					//ih = (ih>=params->offsetDim)?params->offsetDim-1:ih;

					ihx = (ihx>=params->offsetxDim)?params->offsetxDim-1:ihx;
					ihy = (ihy>=params->offsetyDim)?params->offsetyDim-1:ihy;

//					** Turn off the previous trace in that bin **********************/
					if (params->binning[imx][imy][ihx][ihy] == 1)
					{
						params->binning[imx][imy][ihx][ihy] = 0;
						params->cmpOffsetBinPopulate--;
					}

					params->cmpX[currLinNum * totSrcInLine + m][count][0] = (((params->srcLoc[currLinNum * totSrcInLine + m][0]* opp->SLI + opp->initSlOffset) + (params->recLoc[preRecLoc+k][1]) * opp->RI) / 2); // CMPX

					params->cmpX[currLinNum * totSrcInLine + m][count][2] =	(params->recLoc[preRecLoc+k][1]);	  // RN

					params->cmpXVal[currLinNum * totSrcInLine * params->tempNRL* params->tempNR	+ m * (params->tempNRL * params->tempNR) + count] = ((params->srcLoc[currLinNum * totSrcInLine + m][0] * opp->SLI + opp->initSlOffset) + (params->recLoc[preRecLoc + k][1])* opp->RI)/2;  // CMP-Y

					params->cmpY[currLinNum * totSrcInLine * params->tempNRL* params->tempNR+ m * (params->tempNRL * params->tempNR) + count] =((params->srcLoc[currLinNum * totSrcInLine + m][1]* opp->SI + (params->recLoc[preRecLoc + k][0])* opp->RLI) / 2);  // CMP-Y

					offsetx = fabs((params->recLoc[preRecLoc + k][1]) * opp->RI - (params->srcLoc[currLinNum * totSrcInLine + m][0]* opp->SLI + opp->initSlOffset)); // hx

					offsety =  fabs((params->recLoc[preRecLoc + k][0]) * opp->RLI - params->srcLoc[currLinNum * totSrcInLine + m][1]* opp->SI); // hy

					params->offSetX[currLinNum * totSrcInLine * params->tempNRL* params->tempNR+ m * (params->tempNRL * params->tempNR) + count] =offsetx; //offsetx

					params->offSetY[currLinNum * totSrcInLine * params->tempNRL* params->tempNR+ m * (params->tempNRL * params->tempNR) + count] =offsety; //offsety

					params->offSet[currLinNum * totSrcInLine * params->tempNRL* params->tempNR+ m * (params->tempNRL * params->tempNR) + count] =(sqrtf(powf(offsetx, 2) + powf(offsety, 2)));//offset

					azimuth = (atan2(offsety, offsetx) * radtodeg);

					params->aziMuth[currLinNum * totSrcInLine * params->tempNRL* params->tempNR+ m * (params->tempNRL * params->tempNR) + count] =(azimuth < 0 ? azimuth + 360 : azimuth); // Azimuth

					count++;
				} // End of receiver
			}// End of receiver line
			params->firstLocRecTem[currLinNum * totSrcInLine + m][1] =params->firstLocRecTem[currLinNum * totSrcInLine][1];
		}// End of source

		if(backFlag != 3)
		{
			computeBinning(params, gpp, rpp, 1);
			computePatch(params, gpp, rpp, 1);
		}
		else if(backFlag == 3)
		{
			computeBinning(params, gpp, rpp, 1);
		}
	}

	//** For the perturbation of one receiver line in the acquisition design *****//
	else
	{
		count = 0;
		for (i=0; i<params->totSrc; i++)
		{
			for(int j=0; j<totRecInLine; j++)
			{
				//***Turn off the previous trace in that bin **************/
				imx = floor(params->cmpXVal[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] / gpp->binX);
				imy = floor(params->cmpY[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] / gpp->binY);
				ihx = floor(params->offSetX[i * params->tempNRL* params->tempNR  +  currLinNum * (params->tempNR) + j] / params->dhx);
				ihy = floor(params->offSetY[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] / params->dhy);

				ihx = (ihx>=params->offsetxDim)?params->offsetxDim-1:ihx;
				ihy = (ihy>=params->offsetyDim)?params->offsetyDim-1:ihy;

				if (params->binning[imx][imy][ihx][ihy] == 1)
				{
					params->binning[imx][imy][ihx][ihy] = 0;
					params->cmpOffsetBinPopulate--;
				}

				// Compute new offset, and CMP bins *****************//
				params->cmpX[i][currLinNum*params->tempNR+j][0] = (((params->srcLoc[i][0]*opp->SLI + opp->initSlOffset) + (params->recLoc[currLinNum * totRecInLine + j][1]) * opp->RI) / 2); // CMPX

				params->cmpXVal[i * params->tempNRL * params->tempNR +  currLinNum * params->tempNR + j] = (((params->srcLoc[i][0] * opp->SLI + opp->initSlOffset) + (params->recLoc[currLinNum * totRecInLine+j][1]) * opp->RI) / 2);  //CMP-X

				params->cmpY[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] = ((params->srcLoc[i][1] * opp->SI + (params->recLoc[currLinNum * totRecInLine+j][0])* opp->RLI) / 2);  // CMP-Y

				offsetx = fabs((params->recLoc[currLinNum * totRecInLine + j ][1])	* opp->RI - (params->srcLoc[i][0] * opp->SLI + opp->initSlOffset));  // hx

				offsety = fabs((params->recLoc[currLinNum * totRecInLine+j][0])* opp->RLI - params->srcLoc[i][1] * opp->SI);  // hy

				params->offSetX[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] = offsetx;

				params->offSetY[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] = offsety;

				params->offSet[i * params->tempNRL* params->tempNR + currLinNum * (params->tempNR) + j] = (sqrtf(powf(offsetx,2)+powf(offsety,2)));		//offset

				azimuth = (atan2(offsety, offsetx) * radtodeg);

				params->aziMuth[i * params->tempNRL* params->tempNR +  currLinNum * (params->tempNR) + j] = (azimuth < 0 ? azimuth + 360 : azimuth); // Azimuth
			}
		}
		if(backFlag != 3)
		{
			computeBinning(params, gpp, rpp,1);
			computePatch(params, gpp, rpp, 1);
		}
		else if(backFlag == 3)
		{
			computeBinning(params, gpp, rpp, 1);
		}
	}
}

//*** Compute 4D binning based on the CMPOffset domain *********/
void Survey::computeBinning(Params* params, GeophyParams* gpp, ReconstructParams* rpp, int allocFlag)
{
	int ihx, ihy, imx, imy, count = 0;
	params->binDimX = ceil((params->surveySizeX)/gpp->binX);
	params->binDimY = ceil((params->surveySizeY)/gpp->binY);
	params->azDim = 360/params->daz;
	gpp->simXmax = *std::max_element(params->offSet.begin(),params->offSet.end());

	gpp->simXmaxX = *std::max_element(params->offSetX.begin(),params->offSetX.end());
	gpp->simXmaxY = *std::max_element(params->offSetY.begin(),params->offSetY.end());

	//*** 0.1 is added to solve the indexing problem *****************************/
	params->offsetDim = ceil((gpp->simXmax + 0.001) / params->dh);
	params->offsetxDim = ceil((gpp->simXmaxX + 0.001) / params->dhx);
	params->offsetyDim = ceil((gpp->simXmaxY + 0.001) / params->dhy);
	params->cmpOffsetBinSize = params->binDimX * params->binDimY * params->offsetxDim * params->offsetyDim;

	// Allocating memory for binning *********************************************/
	if(allocFlag == 0)
		params->binning = alloc4DPointer(params->binDimX, params->binDimY, params->offsetxDim, params->offsetyDim, params->binning);

	// Binning the cmp-offset domain**********************************************/
	for (unsigned int i = 0; i < params->offSetX.size(); i++)
	{
		ihx = floor(params->offSetX[i] / params->dhx);
		ihy = floor(params->offSetY[i] / params->dhy);
		imx = floor(params->cmpXVal[i] / gpp->binX);
		imy = floor(params->cmpY[i] / gpp->binY);
		//iaz = floor(params->aziMuth[i] / params->daz);

		// Checking the bin is populated earlier ********************************/
		if(params->binning[imx][imy][ihx][ihy] == 0) {
			params->binning[imx][imy][ihx][ihy] = 1;
			count++;
		}
	}
	if(allocFlag == 0)
		params->cmpOffsetBinPopulate  = count;
	else
		params->cmpOffsetBinPopulate += count;
}

//*** Compute 5D patch **********************************************************/
void Survey::computePatch(Params* params, GeophyParams* gpp, ReconstructParams* rpp, int flag)
{
	int spatialWindowNum[4];
	int j, k, l, m, countPatch, initCmpX, initCmpY, initOffsetY, initOffsetX, initAzimuth, n, p, q, r;
	int endCmpX, endCmpY, endOffsetX, endOffsetY, countPopBin, countBinx,countBiny, countHx;
	float sumMu = 0.0, sumGridDensity = 0.0;

	params->patch4D = std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > >();

	// Compute the number of patches in 4 spatial window directions *************/
	spatialWindowNum[0] = floor((params->binDimX - params->olapX) / (params->winX - params->olapX)); 		// CMPX
	spatialWindowNum[1] = floor((params->binDimY - params->olapY) / (params->winY - params->olapY)); 		// CMPY
	spatialWindowNum[2] = floor((params->offsetxDim - params->olapHx) / (params->winHx - params->olapHx));  // OFFSETx
	spatialWindowNum[3] = floor((params->offsetyDim - params->olapHy) / (params->winHy - params->olapHy));  // OFFSETy
	//spatialWindowNum[3] = floor((params->azDim - params->olapAz) / (params->winAz - params->olapAz));		// AZIMUTH

	params->totPatches = spatialWindowNum[0] * spatialWindowNum[1] * spatialWindowNum[2] * spatialWindowNum[3];
	countPatch = 0;

	rpp->mu = NULL;
	rpp->simGridDensity = NULL;

	rpp->mu = new float[params->totPatches];
	rpp->simGridDensity = new float[params->totPatches];
	params->patch4D.reserve(params->totPatches);

	//** Computing 4D patches from binning *************************************/
		for (j = 0; j < spatialWindowNum[0]; j++) {
			for (k = 0; k < spatialWindowNum[1]; k++) {
				for (l = 0; l < spatialWindowNum[2]; l++) {
					for (m = 0; m < spatialWindowNum[3]; m++) {
						countPopBin = 0;
						params->patch4D.push_back(std::vector<std::vector<std::vector<std::vector<int> > > >());
						initCmpX = j * (params->winX - params->olapX);
						if (spatialWindowNum[0] == j + 1)
							endCmpX = params->binDimX - 1;
						else
							endCmpX = initCmpX + params->winX - 1;

						initCmpY = k * (params->winY - params->olapY);

						if (spatialWindowNum[1] == k + 1)
							endCmpY = params->binDimY - 1;
						else
							endCmpY = initCmpY + params->winY - 1;

						initOffsetX = l * (params->winHx - params->olapHx);
						if (spatialWindowNum[2] == l + 1)
							endOffsetX = params->offsetxDim - 1;
						else
							endOffsetX = initOffsetX + params->winHx - 1;

						initOffsetY = m * (params->winHy - params->olapHy);

						if (spatialWindowNum[3] == m + 1)
							endOffsetY = params->offsetyDim - 1;
						else
							endOffsetY = initOffsetY + params->winHy - 1;

						params->patch4D[countPatch].reserve(endCmpX - initCmpX+1);

						countBinx = -1;
						for (n = initCmpX; n <= endCmpX; n++) {
							countBinx++;
							countBiny = -1;
							params->patch4D[countPatch].push_back(std::vector<std::vector<std::vector<int> > >());
							params->patch4D[countPatch][countBinx].reserve(endCmpY - initCmpY+1);

							//tempPatchData.push_back(std::vector<std::vector<std::vector<int> > >());
							for (p = initCmpY; p <= endCmpY; p++) {
								countBiny++;
								countHx = -1;
								params->patch4D[countPatch][countBinx].push_back(std::vector<std::vector<int> >());
								params->patch4D[countPatch][countBinx][countBiny].reserve(endOffsetX - initOffsetX+1);
								//tempPatchData[countBinx].push_back(std::vector<std::vector<int> >());
								for (q = initOffsetX; q <= endOffsetX; q++) {
									countHx++;
									params->patch4D[countPatch][countBinx][countBiny].push_back(std::vector<int>());
									params->patch4D[countPatch][countBinx][countBiny][countHx].reserve(endOffsetY - initOffsetY+1);
									//tempPatchData[countBinx][countBiny].push_back(std::vector<int>());
									for (r = initOffsetY; r <= endOffsetY; r++) {
										params->patch4D[countPatch][countBinx][countBiny][countHx].push_back(params->binning[n][p][q][r]);
										//tempPatchData[countBinx][countBiny][countHx].push_back(params->binning[n][p][q][r]);
									}
								}
							}
						}
						countPopBin = computeMutualCoh(params->patch4D[countPatch], rpp->mu, countPatch);
						sumMu += rpp->mu[countPatch];

						// Computing grid density for every patch *******************//
						rpp->simGridDensity[countPatch] = ((float)countPopBin)/((float) (params->patch4D[countPatch].size()*params->patch4D[countPatch][0].size()*params->patch4D[countPatch][0][0].size()*params->patch4D[countPatch][0][0][0].size()));
						//std::cout<<  rpp->simGridDensity[countPatch] << ","<<rpp->mu[countPatch]<<std::endl;
						sumGridDensity += rpp->simGridDensity[countPatch];
						countPatch++;
					}
				}
			}
		} // End of making 4D patch

		if(flag == 0)
		{
			rpp->initMu = new float[params->totPatches];
			rpp->initSimGridDensity = new float[params->totPatches];
			for(int i=0; i<params->totPatches; i++)
			{
				rpp->initMu[i] = rpp->mu[i];
				rpp->initSimGridDensity[i] = rpp->simGridDensity[i];
			}
		}

		// Computing average mutual coherency *******************/
		rpp->avgMu = (sumMu / params->totPatches);

		// Computing average grid density ***********************/
		rpp->avgGridDensity = sumGridDensity / params->totPatches;
	}

	//******* Compute mutual Coherency for every patch in cmp-offset grid ***************************/
	long Survey::computeMutualCoh(std::vector<std::vector<std::vector<std::vector<int> > > >&patch, float *mu, int patchNum)
	{
		int n[4], j, k, l;
		long int lexDim , popBinCount = 0;
		float maxCoeff, muCoherency=-10000.0, dataOut;
		fftw_plan p1;
		fftw_complex *data;

		n[0] = patch.size();
		n[1] = patch[0].size();
		n[2] = patch[0][0].size();
		n[3] = patch[0][0][0].size();

		// Dimensions of the 4D spatial cases //
		data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n[0] * n[1] * n[2] * n[3]);

		for (j = 0; j < n[0]; j++)
			for (k = 0; k < n[1]; k++)
				for (l = 0; l < n[2]; l++)
					for (int m = 0; m < n[3]; m++)
					{
						if(patch[j][k][l][m] == 1)
							popBinCount++;
						data[j * n[1] * n[2] * n[3] + k * n[2] * n[3] + l * n[3] + m] = (double) patch[j][k][l][m]+ 0.0*I;
					}

		p1 = fftw_plan_dft(4, n, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(p1);
		maxCoeff = sqrtf(powf(crealf(data[0]), 2.0) + powf(cimagf(data[0]), 2.0));
		lexDim = n[0] * n[1] * n[2] * n[3];

		for(j = 1; j < lexDim; j++)
		{
			dataOut = (sqrtf(powf(crealf(data[j]), 2) + powf(cimagf(data[j]), 2))/maxCoeff);
			if (muCoherency < dataOut)
				muCoherency = dataOut;
		}

		mu[patchNum] = (maxCoeff<=0.0)?1.0:muCoherency;
		fftw_free(data);
		data = NULL;
		fftw_destroy_plan(p1);

		return popBinCount;
	}

	//**** Evaluate the fitness value of the cost function ************************/
	float Survey::evalFitness(Params* params, GeophyParams* gpp,ReconstructParams* rpp, int minRLI, float RLI, int minSLI, float SLI, int outLoopCount)
	{
		float fitness = 0.0, normXmin, normXmax, normGridDensity;
		int densityCount = 0;
		float normGridEfficiency;

		//gpp->simXmin = sqrtf(powf(maxRLI * RLI, 2) + powf(maxSLI * SLI, 2));

		gpp->simXmin = sqrtf(powf(minRLI * RLI, 2) + powf(minSLI * SLI, 2));

		// Computing critical mutual coherency **********************************/
		rpp->simCriticalMu = standardization(rpp->mu, rpp->avgMu, params->totPatches);

		// Computing critical grid density **************************************/
		for (int i = 0; i < params->totPatches; i++) {
			if (rpp->simGridDensity[i] < rpp->desiredGridDensity)
				densityCount++;
		}

		rpp->simCriticGridDensity = ((float) (densityCount) / (float) params->totPatches);

		// Computing normalised simulated tempGridDensity ***********************/
		params->tempGridDensity = ((rpp->desiredGridDensity - rpp->avgGridDensity)/ rpp->desiredGridDensity);

		// Computing normalised simulated GridDensity ***************************/
		normGridDensity =powf((((params->tempGridDensity + params->theta[outLoopCount][5])> 0.0) ? params->tempGridDensity : 0.0), 2);

		// Computing simulated grid efficiency **********************************/
		rpp->simGridEfficiency = ((float) params->cmpOffsetBinPopulate)/((float) params->offSet.size());

		// Computing normalised simulated grid Efficiency ***********************/
		params->normGridEfficiency = ((rpp->desiredGridEfficiency - rpp->simGridEfficiency) / rpp->desiredGridEfficiency);

		normGridEfficiency = powf((params->normGridEfficiency + params->theta[outLoopCount][4]), 2);

		// Computing normalised simulated Xmax **********************************/
		params->tempNormMax = ((gpp->desXmax - gpp->simXmax) / gpp->desXmax);
		params->tempNormMax = (fabs(params->tempNormMax)>1.0)?(params->tempNormMax/fabs(params->tempNormMax)):params->tempNormMax;

		normXmax = powf((((params->tempNormMax + params->theta[outLoopCount][3]) > 0.0)?params->tempNormMax : 0.0), 2);

		// Computing normalised simulated Xmin **********************************/
		params->tempNormMin = (((gpp->desXmin - gpp->simXmin) / gpp->desXmin));
		params->tempNormMin = (fabs(params->tempNormMin)>1.0)?(params->tempNormMin/fabs(params->tempNormMin)):params->tempNormMin;

		normXmin = powf((((params->tempNormMin + params->theta[outLoopCount][2]) < 0.0)?fabs(params->tempNormMax) : 0.0), 2);

		// Evaluate the fitness function ****************************************/
		fitness =rpp->avgMu + 0.5*(params->gamma[0]* powf((rpp->simCriticalMu+ params->theta[outLoopCount][0]),2)+
							 params->gamma[1]* powf((rpp->simCriticGridDensity+ params->theta[outLoopCount][1]),2)+
							 params->gamma[2] * normXmin+params->gamma[3] * normXmax+params->gamma[4]*normGridEfficiency+
							 params->gamma[5] * normGridDensity);

		return fitness;
	}

	/***** Allocate 2D pointer ***************************************/
	int ** Survey::alloc2DPointer(int n1, int n2, int **data) {
		if (!(data = new int *[n1])) {
			std::cout << " Error out of memory " << std::endl;
			exit(1);
		} else {
			for (int i = 0; i < n1; i++)
				data[i] = new int[n2];
		}
		return data;
	}

	/***** Allocate 2D pointer ***************************************/
	float ** Survey::alloc2DPointerFloat(int n1, int n2, float **data) {
		if (!(data = new float *[n1])) {
			std::cout << " Error out of memory " << std::endl;
			exit(1);
		} else {
			for (int i = 0; i < n1; i++)
				data[i] = new float[n2];
		}
		return data;
	}

	/***** Allocate 4D pointer ***************************************/
	int **** Survey::alloc4DPointer(int n1, int n2, int n3, int n4, int ****data) {
		int i, j, k, l;

		data = new int ***[n1];
		if (!(data = new int ***[n1])) {
			std::cout << " Error: out of memory " << std::endl;
			exit(1);
		} else {
			for (i = 0; i < n1; i++) {
				data[i] = new int**[n2];
				for (j = 0; j < n2; j++) {
					data[i][j] = new int*[n3];
					for (k = 0; k < n3; k++) {
						data[i][j][k] = new int[n4];
						for (l = 0; l < n4; l++) {
							data[i][j][k][l] = 0;
						}
					}
				}
			}
		}
		return data;
	}

	/***** Allocate 4D pointer ***************************************/
	void Survey::deallocPointer(int n1, int n2, int n3, int ****data)
	{
		int i, j, k;
		for (i = 0; i < n1; i++)
			for (j = 0; j < n2; j++)
				for (k = 0; k < n3; k++)
				{
					delete [] data[i][j][k];
					data[i][j][k] = NULL;
				}
		delete [] data;
		data = NULL;
	}

	// Write binary file for the optimised sources' and receivers' locations ****** //
	void Survey::write2Ddata(std::string fileName, int n1, int n2, int **data)
	{
		int i, j;
		std::ofstream output;
		std::cout << fileName <<std::endl;
		output.open(fileName.c_str(), ios::out);

		try
		{
			if(n2==2)
			{
				output<<"Source Number";
				output<<",";
				output<<"Receiver line Number";
				output<<",";
				output<<"Receiver Number";
				output<<"\n";
			}

			for(i=0; i<n1; i++)
			{
				if(n2==2)
				{
					output<<i+1;
					output<<",";
					output<<data[i][0]+1;
					output<<",";
					output<<data[i][1]+1;
				}
				else
					for(j=0; j<n2; j++)
					{
						output<<data[i][j];
						if(j<n2-1)
							output<<",";
					}
				output<<"\n";
			}
			output.close();
		}
		catch (int n)
		{
			std::cout <<"ERROR: " << n <<std::endl;
		}
	}

	// Write binary file for the optimised sources' and receivers' locations ****** //
	void Survey::write1Ddata(std::string fileName, int n1, float *data)
	{
		int i;
		std::ofstream output1 (fileName.c_str());
		std::cout << fileName <<std::endl;
		try
		{
			for(i=0; i<n1; i++)
			{
				output1<<data[i];
				output1<<"\n";
			}
		}
		catch (int n)
		{
			std::cout <<"ERROR: " << n <<std::endl;
		}
	}

	// Write binary file for the optimised sources' and receivers' locations ****** //
	void Survey::writeDdata(std::string fileName, Params* params, GeophyParams* gpp, ReconstructParams* rpp, int flag )
	{
		std::ofstream output1 (fileName.c_str());

		std::cout << fileName <<std::endl;

		if(flag == 1)
		{
			output1<<"Number of receiver in Template: " << params->tempNR*params->tempNRL;
			output1<<"\n";
			output1<<"Number of receiver line in Template: " << params->tempNRL;
			output1<<"\n";
			output1<<"Simulated Maximum offset: " << gpp->simXmax;
			output1<<"\n";
			output1<<"Simulated Largest minimum offset: " << gpp->simXmin;
			output1<<"\n";
			output1<<"Simulated grid efficiency: " << rpp->simGridEfficiency;
			output1<<"\n";
		}
		else if(flag == 2)
		{
			for(unsigned int i=0; i<params->finalCons.size(); i++)
			{
				output1<<params->finalCons[i];
				output1<<"\n";
			}
		}
		else
		{
			for(unsigned int i=0; i<rpp->itrAvgMu.size(); i++)
			{
				output1<<rpp->itrAvgMu[i];
				output1<<"\n";
			}
		}
	}
