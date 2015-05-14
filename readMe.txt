# OptSurvey
Software package for optimizing 3D seismic survey 
Developed by: Mafijul Bhuiyan

About
=====

optSurvey is a C++ program which can deduce optimal 3D seismic survey design applying global optimization method named "Simulated Annealing  (SA) with augmented lagrangian method and mutual coherency criteria. Mutual coherency is computed using 5D fast fourier transformation (FFT). 


Installing
==========

Required packages: 
 1. FFTW 3.3.0 or higher
 2. g++
 
Installation of optSurvey is a straightforward process. You just need to copy the whole software package in your machine
and run make all command in your terminal. 

Make sure the path of FFTW plug-in is correct. 


Running
=======
Before running, check the input parameters in input.txt file. Input.txt file contains the dimensions of 3D curvelet  
transformation. 

1. cd Debug
2. make clean 
3. make all 
4. ./fdct3d 


That's all. 

Please contact: mbhuiyan@ualberta.ca
