/*  =================================================================================
 *  Name        : Optsurvey.cpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 1.0
 *  Purpose     : Optimum survey design. We have decimated receiver and source lines.
 *  			: Here we have moved the receiver and source lines in inline and X-line
 *  			: directions.
 *  Date        : February 11, 2015
 *  Affiliation : University of Alberta, Physics department, (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 *  =================================================================================
 */

#include "Params.hpp"
#include "OperationalParams.hpp"
#include "GeophyParams.hpp"
#include "ReconstructParams.hpp"
#include "Survey.hpp"

int main()
{
   clock_t start, finish;
   float binSize[3];
   int count=0, val, lineCounter = 0;
   std::vector<int> rndSrcLoc; 		  // For random source locations
   std::vector<int> rndRecLoc; 		  // For random receiver locations

   /** Creating the objects of the dependence classes ************************************/
   OperationalParams opp;
   GeophyParams gpp;
   ReconstructParams rpp;
   Params params;
   Survey sur;

   /** Parameters for simulation and initial ACQ design **********************************/
   params.fMax = 25;   			// Hz
   params.fDom = 18;    		// Hz
   params.vAve = 2150; 			// (m/s)
   params.vMax = 2800; 			// (m/s)
   params.vMin = 1500; 			// (m/s)
   params.depthSurvey = 4000;	// (m)
   params.targetSize = 100;		// (m)
   params.depthShHorizon = 80;	// (m)
   params.surveySizeX = 1475;   // (m)
   params.surveySizeY = 1475;   // (m)
   params.constraintNum = 6;    // Number of constraint

   /** Compute bin size parameter ********************************************************/
   binSize[0] = floor(params.targetSize/3);
   binSize[1] = floor(params.vAve/(4*params.fDom));
   binSize[2] = floor(params.vMin/(4*params.fMax));

   /** Initialisation of Geophysical parameters. *****************************************/
   gpp.desXmax = 1.4*params.depthSurvey;
   gpp.desXmin = 1.4*params.depthShHorizon;
   gpp.binX = *std::min_element(binSize, binSize+3);
   gpp.binY = gpp.binX;

   /** Initialisation of Operational parameters. *****************************************/
   opp.RI = 2*gpp.binX;				  					 	 // Receiver interval (m)
   opp.SI = 2*gpp.binY;				  					 	 // Source interval (m)
   opp.orgRLI = round((1/sqrtf(5.0))*gpp.desXmin);	  	 	 // Receiver line interval (m)
   opp.orgSLI = 2*opp.orgRLI;				  				 // Source line interval (m)
   opp.orgNRL = floor(params.surveySizeY/opp.orgRLI)+1;  	 // Number of receiver line
   opp.orgNSL = floor(params.surveySizeX/opp.orgSLI)+1;   	 // Number of source line
   opp.SLI = 25;
   opp.RLI = 25;
   opp.NRL = floor(params.surveySizeY/opp.RLI)+1;
   opp.NSL = floor(params.surveySizeX/opp.SLI)+1;
   opp.NR = (floor(params.surveySizeX/opp.RI)+1)*opp.orgNRL;   // Total number of receivers in the whole survey
   opp.NS = (floor(params.surveySizeY/opp.SI)+1)*opp.orgNSL;   // Total number of sources in the whole survey
   opp.initSlOffset = opp.RI/2;							  	    //  Source line offset from the first receiver

   /** Initialisation of Reconstructional parameters. *********************************/
   rpp.desiredGridDensity = 0.15;                		  // grid density = (populated bin)/(total number of bins)
   rpp.desiredCriticGridDensity = 0.0;
   rpp.desCriticalMu = 0.0;
   rpp.desiredGridEfficiency = 1.0;						  // Grid efficiency = (total number of traces in CMP-Offset domain) / (total number of traces in S-R domain)
   opp.recAvail = 6400;

   /** Initialising Lagrangian multipliers ********************************************/
   params.gamma = new float [params.constraintNum];
   params.theta = sur.alloc2DPointerFloat(30, params.constraintNum, params.theta);

   for(int i=0;i<params.constraintNum;i++)
   	   params.gamma[i] = 0.1;
   for(int j=0;j<30;j++)
	   for(int i=0;i<params.constraintNum;i++)
		   params.theta[j][i] = 0.1;

   // Initialisation of windowing parameters *****************************************/
   params.winX = 20;
   params.winY = 25;
   params.winHx = 13;
   params.winHy = 13;
   params.winAz = 4;
   params.olapHx = (int) floor(params.winHx*0.25);
   params.olapHy = (int) floor(params.winHy*0.25);
   params.olapX = (int) floor(params.winX*0.25);
   params.olapY = (int) floor(params.winY*0.25);
   params.olapAz = (int) floor(params.winAz*0.25);
   // Initialisation of windowing parameters *****************************************/

   /** Initialisation of template design parameters. *********************************/
   params.templateExist = 0;
   params.testFlag = 0;
   params.tempLenX = 2*gpp.desXmax/sqrt(5);				      // Template length in X direction
   params.tempLenY = gpp.desXmax/sqrt(5);				      // Template length in Y direction
   params.tempNR = ceil(params.tempLenX/(2*gpp.binX))+1;     // Number of receivers in one template line
   params.tempNRL = ceil(params.tempLenY/(opp.orgRLI))+1;    // Number of receiver line in the template
   params.srcLineDecimationRate = 0.67;					     // Rate of available in source line
   params.srcDecimationRate = 1.0;						     // Rate of decimation in source
   params.recDecimationRate = 0.67;						     // Rate of available in receiver line
   params.allowedMovementSrcLine = 2;                        // Predefined based on decimation and source Line perturbation
   params.allowedMovementSrcX = 1;                            // Predefined based on decimation and source perturbation
   params.allowedMovementSrcY = 1;                            // Predefined based on decimation and source perturbation
   /** End of initialisation of template design parameters. *************************/

   // Receiver and source distribution after decimation *****************************/
   opp.decNRL = round(opp.orgNRL * params.recDecimationRate);    				  // Decimated number of receiver line
   params.intervalRec = floor((float)(opp.NRL-1)/(float)(opp.decNRL-1))*1.0;    // Receiver line interval
   params.allowedMovementRec = 1;                              					  // Predefined based on decimation and receiver perturbation
   opp.decNSL = round(opp.orgNSL * params.srcLineDecimationRate);    			  		   // Decimated number of source line
   opp.decNS = round((floor(params.surveySizeY/opp.SI)+1) * params.srcDecimationRate);    // Decimated number of sources in one line
   params.intervalSrcLine = floor((float) (opp.NSL-1)/(float)(opp.decNSL - 1));     	   // Source line interval
   params.intervalSrc = round((floor(params.surveySizeX/opp.SI))/(opp.decNS - 1));   	   // Sources interval in one line
   // Receiver and source distribution after decimation *****************************/

   // If every source is listened by every receiver
   if(params.templateExist != 1 )
   {
   	  params.tempLenX = params.surveySizeX;											  // Template length in X direction
   	  params.tempLenY = params.surveySizeY;											  // Template length in Y direction
   	  params.tempNR = floor(params.surveySizeX/opp.RI)+1; 							  // Number of receivers in one template line
   	  params.tempNRL = opp.decNRL;   												  // Number of receiver line in the template
   }

   // Allocating memory for 2D pointer for src and receiver locations.
   // The first column of array is for the line number of a source or a receiver.
   // The second column of array is for the position of a source or a receiver in that line ***//
   params.recLoc = sur.alloc2DPointer(opp.decNRL*(floor(params.surveySizeX/opp.RI)+1), 2, params.recLoc);
   params.srcLoc = sur.alloc2DPointer(opp.decNSL*opp.decNS, 2, params.srcLoc);

   // Memory allocation of the sources and receivers distribution as 2D array *********//
   params.rec = sur.alloc2DPointer(opp.NRL, floor(params.surveySizeX/opp.RI)+1, params.rec);
   params.src = sur.alloc2DPointer(opp.NSL, floor(params.surveySizeY/opp.SI)+1, params.src);

   // Initialisation of the shots and receivers locations in a 2D array ***************//
   for(int i=0; i<floor(params.surveySizeX/opp.RLI)+1; i++)
	   for(int j=0; j<floor(params.surveySizeX/opp.RI)+1; j++)
		   params.rec[i][j] = 0;

   for(int i=0; i<floor(params.surveySizeY/opp.SLI)+1;i++)
   	   for(int j=0; j<floor(params.surveySizeY/opp.SI)+1;j++)
   		   params.src[i][j] = 0;

   opp.maxRLI = params.intervalRec;
   opp.maxSLI = params.intervalSrc;

   opp.minRLI = params.intervalRec;
   opp.minSLI = params.intervalSrc;

   //** Random interval of decimated receiver lines  ****************************//
   if(params.testFlag == 0)
   {
	   if(params.recDecimationRate>0.9)
	   {
		  for(int i=0;i<opp.NRL-1; i++)
			 rndRecLoc.push_back(i+1);
		  std::random_shuffle(rndRecLoc.begin(),rndRecLoc.end());
		  std::sort(rndRecLoc.begin(), rndRecLoc.begin()+opp.decNRL-1);

		  int i = 0;
		  while(i<opp.decNRL)
		  {
			  if(i>opp.NRL)
				  i = opp.NRL;
			  for(int j=0; j<floor(params.surveySizeX/opp.RI)+1; j+=params.intervalRec)
			  {
				  if(i == 0)
				  {
					   params.recLoc[count][0]=i;   		 	// Receiver line
					   params.rec[i][j] = 1;
				  }
				  else
				  {
					   params.recLoc[count][0]=rndRecLoc[i-1];  // Receiver line
					   params.rec[rndRecLoc[i-1]][j] = 1;
				  }
				  params.recLoc[count][1]=j;   					// Receiver number in receiver line
				  count++;
			  }
			  i++;
		  }
	  }
	  else //** Regular interval of decimated receiver lines  ****************************//
	  {
		for(int i=0; i<opp.NRL; i+=params.intervalRec)
		{
		   if(i>opp.NRL)
			   i = opp.NRL;
		   for(int j=0; j<floor(params.surveySizeX/opp.RI)+1; j++)
		   {
			   params.recLoc[count][0]=i;   // Receiver line
			   params.recLoc[count][1]=j;	// Receiver location in that line
			   params.rec[i][j] = 1;
			   count++;
		   }
		 }
	  }
   }
   else
   {
	   // Read the receiverGeoFile.csv file to compute the grid density and grid efficiency
	   std::ifstream recFileInfo("/home/entropy/workspace/optSurvey/src/receiverGeoFile.csv");
	   while(!recFileInfo.eof())
	   {
	   	    int j = 0;
			getline(recFileInfo, params.line);
			std::istringstream info_line (params.line);
			while(getline(info_line,params.token,','))
			{
				std::istringstream temp(params.token);
				temp>>val;
				if(val == 1)
				{
					params.recLoc[count][0] = lineCounter;   // Receiver line
					params.recLoc[count][1] = j;   	   	 	 // Receiver number in receiver line
					params.rec[lineCounter][j] = 1;
					count++;
				}
				j++;
			}
			lineCounter++;
	   }
   }

  params.totRec = count;
  count = 0;
  lineCounter = 0;

   //** Random interval of decimated source lines  ****************************//
   if(params.testFlag == 0)
   {
	   if(params.srcLineDecimationRate > 0.9)
	   {
		  for(int i=1;i<opp.NSL; i++)
			  rndSrcLoc.push_back(i);
		  std::random_shuffle(rndSrcLoc.begin(),rndSrcLoc.end());
		  std::sort(rndSrcLoc.begin(), rndSrcLoc.begin()+opp.decNSL-1);
		  int i = 0;

		  while(i<opp.decNSL)
		  {
			if(i>opp.NSL)
				 i = opp.NSL;
			for(int j=0; j<floor(params.surveySizeY/opp.SI)+1; j+=params.intervalSrc)
			{
				if(i == 0)
				{
					params.srcLoc[count][0]=i;   		 // Source line
					params.src[i][j] = 1;
				}
				else
				{
					params.srcLoc[count][0]=rndSrcLoc[i-1];  // Source line
					params.src[rndSrcLoc[i-1]][j] = 1;
				}
			    params.srcLoc[count][1]=j;   // Source number in source line
			    count++;
			}
			i++;
		  }
	   }
	   else //** Regular interval of decimated source lines  *********************//
	   {
		   for(int i=0; i<opp.NSL; i+=params.intervalSrcLine)
		   {
			   if(i>opp.NSL)
					 i = opp.NSL;
			   for(int j=0; j<floor(params.surveySizeY/opp.SI)+1; j+=params.intervalSrc)
			   {
				   params.srcLoc[count][0]=i;   // Source line
				   params.srcLoc[count][1]=j;   // Source number in source line
				   params.src[i][j] = 1;
				   count++;
			   }
		   }
	   }
   }
   else
   {
	   std::ifstream srcFileInfo("/home/entropy/workspace/optSurvey/src/sourceGeoFile.csv");
	   int lineCounter = 0;
	   while(!srcFileInfo.eof())
	   {
		    int j = 0;
			getline(srcFileInfo, params.line);
			std::istringstream info_line (params.line);
			while(getline(info_line,params.token,','))
			{
				std::istringstream temp(params.token);
				temp>>val;
				if(val == 1)
				{
					params.srcLoc[count][0] = lineCounter;   // Source line
					params.srcLoc[count][1] = j;   	   	 	 // Source number in source line
					params.src[lineCounter][j] = 1;
					count++;
				}
				j++;
			}
			lineCounter++;
	   }
   }

   params.totSrc = count;
   params.perturbRecNode = 1;
   params.perturbSrcNode = 1;
   params.salvo = (ceil) ((params.intervalRec*opp.orgRLI)/opp.SI)+1;	// Number of shot before the movement of the template
   gpp.binX = 29;														// Size of CMPx bin
   gpp.binY = gpp.binX;													// Size of CMPy bin
   params.dh = 40;		  												// Size of offset bin (meter)
   params.dhx = 30;		  												// Size of offsetX bin (meter)
   params.dhy = 30;		  												// Size of offsetY bin (meter)
   params.daz = 45;       												// Size of azimuth bin (degree)


   params.currSrcLoc = new int *[opp.decNSL];
   params.preSrcLoc = new int *[opp.decNSL];
   params.currRecLoc = new int *[opp.decNRL];
   params.preRecLoc = new int *[opp.decNRL];


   for (int i = 0; i <opp.decNRL; i++)
   {
	   params.currRecLoc[i] = new int[params.totRec/opp.decNRL];
	   params.preRecLoc[i] = new int[params.totRec/opp.decNRL];
	   for(int j=0; j<floor(params.surveySizeX/opp.RI); j++)
	   	   	   params.preRecLoc[i][j] = 0;
   }

   for (int i = 0; i <opp.decNSL; i++)
   {
	   params.currSrcLoc[i] = new int[params.totSrc/opp.decNSL];
	   params.preSrcLoc[i] = new int[params.totSrc/opp.decNSL];
	   for(int j=0; j<floor(params.surveySizeY/opp.SI); j++)
	   	       params.preSrcLoc[i][j] = 0;
   }

   params.numthreads = 1;

   // Display message ********************************* //
   std::cout <<"Acquisition parameters for ideal case: " <<std::endl;
   std::cout <<"--------------------------------------" <<std::endl;
   std::cout <<"Total number of required receivers: " << opp.NR <<std::endl;
   std::cout <<"Total number of required sources: " << opp.NS <<std::endl;
   std::cout <<"Receiver interval: " << opp.RI <<std::endl;
   std::cout <<"Receiver line interval: " << opp.orgRLI <<std::endl;
   std::cout <<"Source interval: " << opp.SI <<std::endl;
   std::cout <<"Source line interval: " << opp.orgSLI <<std::endl;
   std::cout <<"Required least number of available receivers (based on maximum offset and target depth) or template size: " << params.tempNR * params.tempNRL <<std::endl;
   std::cout <<"Swath shooting path will be provided after the simulation: " << params.tempNR * params.tempNRL <<std::endl;
   std::cout <<"Simulation started ... ..." <<std::endl;

   // SA parameters ***********************************************/
   params.numCycle = 100;
   params.initTemp = 100.0;

   // Survey design parameters
   params.surType = ORTHO;

   // File name to save *******************************************/
   params.swathShootTrac = "swathShootingFile.csv";
   params.srcGeo = "sourceGeoFile.csv";
   params.recGeo = "receiverGeoFile.csv";
   params.convergence1 = "saConvergence.csv";
   params.convergence2 = "aulConvergence.csv";
   params.mutualCoherency = "mutualCoherency.csv";
   params.simGridDensity = "simGridDensity.csv";
   params.initMutualCoherency = "initMutualCoherency.csv";
   params.initSimGridDensity = "initSimGridDensity.csv";
   params.otherInfo = "otherInfo.csv";
   params.finalFitness = "finalFitness.csv";
   params.finalAvgMu = "finalAvgMu.csv";
   params.logFile = "logFile.csv";

   params.convergenceSAconstraint = NULL;
   params.convergenceSAconstraint = new float[params.numCycle+1];

   // Executing ACQ design optimisation algorithm ******************/
   gpp.desXmax = 2100;
   start = clock();
   sur.findOptSurvey(&params, &gpp, &rpp, &opp);
   finish = clock();
   std::cout << "Time: " << (finish-start)/double(CLOCKS_PER_SEC) << " Seconds "<<std::endl;

   // Writing the files *******************************************************************/
   sur.write2Ddata(params.swathShootTrac, params.totSrc, 2, params.firstLocRecTem);
   sur.write2Ddata(params.recGeo, opp.NRL, floor(params.surveySizeX/opp.RI)+1, params.rec);
   sur.write2Ddata(params.srcGeo, opp.NSL, floor(params.surveySizeY/opp.SI)+1, params.src);
   sur.write1Ddata(params.mutualCoherency, params.totPatches, rpp.mu);
   sur.write1Ddata(params.simGridDensity, params.totPatches, rpp.simGridDensity);
   sur.write1Ddata(params.initMutualCoherency, params.totPatches, rpp.initMu);
   sur.write1Ddata(params.initSimGridDensity, params.totPatches, rpp.initSimGridDensity);
   sur.write1Ddata(params.convergence1, params.numCycle, params.convergenceSAconstraint);
   sur.writeDdata(params.otherInfo, &params, &gpp, &rpp, 1);
   sur.writeDdata(params.finalFitness, &params, &gpp, &rpp, 2);
   sur.writeDdata(params.finalAvgMu, &params, &gpp, &rpp, 0);

   /** Clearing the memory space **********************************************************/
   params.reset();
   opp.reset();
   gpp.reset();
   rpp.reset();
   std::cout << "Cleared the memory space." <<std::endl;
 return 0;
}
