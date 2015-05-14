/*  =================================================================================
 *  Name        : RefVelParams.hpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Store reference velocities information for PSPI
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

/**
 * This data structure holds the parameters for S-R domain
 */

#include <iostream>
using std::cout;
using std::ostream;
using std::ofstream;
using std::endl;
#include <fstream>
#include <sstream>
#include <vector>

/**
 * This data structure holds the dimension and grid size of model, and seismic trace constant parameters
 */

enum surveyType{ORTHO, STAGG}; /// Different migration types

class Params{
public:
	// These parameters depend on the seismic data and velocity model *******************/
	int ****binning;					/// 4D binning
	std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > >patch4D; /// 4D patch
	std::vector<std::vector<std::vector<float> > >cmpX; /// CMP in X direction plus
	std::vector<float>cmpXVal;			/// CMP in X direction
	std::vector<float>cmpY;				/// CMP in Y direction
	std::vector<float>offSet;			/// offSet
	std::vector<float>offSetX;			/// offSetX
	std::vector<float>offSetY;			/// offSetY
	std::vector<float>aziMuth;			/// Azimuth
	std::vector<float>finalCons;		/// Final constraint values
	std::vector<int>chkPatch;			/// Check modified patch number
	std::vector<std::vector<int> > patchTrack; // Track the bin index in every patch.
	int **firstLocRecTem;				/// Store the first location of receiver template for every source
	int **recLoc;          				/// Receiver location
	int **srcLoc;        				/// Source location
	int **src;          				/// Source geometry
	int **rec;        					/// Receiver Geometry
	float **theta;						/// Optimisation parameter1
	float *gamma;						/// Optimisation parameter2
	float *convergenceSAconstraint;		/// Convergence SA constrained
	int **currSrcLoc;					/// Current perturbed locations of Sources (src)
	int **currRecLoc;					/// Current perturbed locations of Receivers (Rec)
	int **preSrcLoc;					/// Current perturbed locations of Sources (src)
	int **preRecLoc;					/// Current perturbed locations of Receivers (Rec)
	float convergenceSAunconstraint;	/// Convergence SA Unconstrained
	bool templateExist; 				/// Template Exist or not
	float fitness;						/// Fitness value of the cost function
	int totPatches;						/// Number of total patch
	long int cmpOffsetBinSize; 			/// Dimension of the CMPOffset bin size
	long int cmpOffsetBinPopulate; 		/// Populated number of CMPOffset bins
	float tempGridDensity;				/// Normalised grid density
	float tempNormMax;					/// Normalised Xmax
	float tempNormMin;					/// Normalised Xmin
	float normGridEfficiency;			/// Normalised grid efficiency
	float fMax;	        				/// Maximum frequency used in shot (Hz)
	float fDom;	        				/// Dominant frequency used in shot (Hz)
	float vMin;	        				/// Minimum velocity in the subsurface layer (m/s)
	float vMax;          				/// Maximum velocity in the subsurface layer (m/s)
	float vAve;          				/// Average velocity in the subsurface layer (m/s)
	float srcLineDecimationRate;       /// Source Line decimation rate (e.g., 50% ----> 0.5)
	float srcDecimationRate;       		/// Source decimation rate (e.g., 50% ----> 0.5)
	float recDecimationRate;          	/// Receiver decimation rate (e.g., 50% ----> 0.5)
	float depthSurvey;          		/// Depth of survey (m)
	float targetSize;			 		/// Size of the smallest target (m)
	float depthShHorizon;			 	/// Depth of shallow horizon
	float surveySizeX;					/// Survey size in X-direction
	float surveySizeY;					/// Survey size in Y-direction
	int perturbSrc;						/// Perturbation of Source line
	int perturbRec;                     /// Perturbation of Receiver line
	int perturbSrcNode;					/// Perturbation of Source node
	int perturbRecNode;					/// Perturbation of Receiver node
	int allowedMovementSrcLine;			/// Allowed movement of source line
	int allowedMovementSrcX;			/// Allowed movement of source in X direction
	int allowedMovementSrcY;			/// Allowed movement of source in Y direction
	int allowedMovementRec;				/// Allowed movement for receiver line
	float intervalSrcLine;				/// Source line interval
	int intervalSrc;					/// Source interval
	float intervalRec;					/// Receiver interval
	int totSrc;							/// Total Source
	int totRec;							/// Total Receiver
	int salvo;							/// Number of shots before rolling the receiver patch
	int winX;							/// Window in the inline direction
	int winY;							/// Window in the crossline direction
	int winHx;							/// Window in the offsetX direction
	int winHy;							/// Window in the offsetY direction
	int winAz;							/// Window in the offsetY direction
	int olapX;							/// Window in the overlapX direction
	int olapY;							/// Window in the overlapY direction
	int olapHx;							/// Window in the overlapHx direction
	int olapHy;							/// Window in the overlapHy direction
	int olapAz;							/// Window in the overlapHx direction
	int binDimX;						/// Dimension of cmp in inline direction
	int binDimY;						/// DImension of cmp in crossline direction
	int	offsetDim;						/// Dimension of bin in offset direction
	int	offsetxDim;						/// Dimension of bin in offsetx direction
	int	offsetyDim;						/// Dimension of bin in offsety direction
	int azDim;							/// Dimension of bin in azimuth direction
	float dh;							/// Offset bin length
	float daz;							/// Azimuth bin length
	float dhx;							/// Offsetx bin length
	float dhy;							/// Offsety bin length
	float tempLenX;						/// Template length in X direction
	float tempLenY;						/// Template length in Y direction
	int tempNR;							/// Number of receivers in X direction
	int tempNRL;						/// Number of receiver line in Y direction
	int constraintNum;					/// Number of constraint
	float initTemp;						/// Temperature of SA
	int numCycle;						/// Number of cycles in SA
	int numthreads;						/// Number of threads
	int testFlag; 						/// Flag for testing the code; 1==> do test, 0==> don't test
	surveyType surType;        			/// Assign initial survey type
	std::string swathShootTrac, srcGeo, recGeo, mutualCoherency, convergence1, convergence2, initMutualCoherency;
	std::string otherInfo, simGridDensity, initSimGridDensity, finalFitness, finalAvgMu, line, token, logFile;
	std::ofstream logOut;

	/// Clear configuration variables ///
	void reset()
	{
		delete [] src; src = NULL;
		delete [] rec; rec = NULL;
		delete [] theta; theta = NULL;
		delete [] gamma; gamma = NULL;
		delete [] convergenceSAconstraint; convergenceSAconstraint = NULL;
		delete [] binning; binning = NULL;
		cmpX = std::vector<std::vector<std::vector<float> > >(); /// CMP in X direction plus
		cmpXVal = std::vector<float>();			/// CMP in X direction
		cmpY = std::vector<float>();			/// CMP in Y direction
		offSet = std::vector<float>();			/// offSet
		offSetX = std::vector<float>();			/// offSetx
		offSetY = std::vector<float>();			/// offSety
		aziMuth = std::vector<float>();			/// Azimuth
		chkPatch = std::vector<int>();			/// Check modified patch number
		delete [] firstLocRecTem; firstLocRecTem = NULL;
		delete [] recLoc; recLoc = NULL;
		delete [] srcLoc; srcLoc = NULL;
	}
};

#endif /* PARAMS_HPP_ */
