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
From optSurvey directory in terminal run the following commands: 

1. cd Debug
2. make clean 
3. make all 
4. ./optSurvey 

Output
======
The program will generate and write several output file in Debug directory. Please check the files for the final results 
and convergence of the algorithm.


That's all. 

Please contact: mbhuiyan@ualberta.ca for further inquiries. 
