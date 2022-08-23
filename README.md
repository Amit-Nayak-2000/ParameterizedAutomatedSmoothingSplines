# Parameterized Automated Smoothing Splines

By: Amit Nayak & Dr. David Kopriva

This software takes in noisy data and outputs smoothed data to within a user defined length scale (minimum radius of curvature) and tolerance. 
The parameterized smoothing splines that are used are based on the SMOOTH algorithm in A Practical Guide to Splines by Carl de Boor [1]. 
A C++ compiler is required to build the software. GNU and CLANG compilers were tested on MacOS.
The user must provide the noisy data in a .csv format and input their parameters into the control file.
A template of the control file can be seen in userdef.control.
The software will output the smoothed data in a .csv file format. 

The commands to compile and run the software are as follows:
g++ -std=c++11 main.cpp
./a.out yourcontrolfile.control


If you have any questions, please email me at anaya085@uottawa.ca


References:
[1] C. D. Boor, A Practical Guide to Splines. New York: Springer, 2001. 
