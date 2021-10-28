# SWE_1D
This is a Shallow Water Equations solver (in 1D), which employs the HLL approximate Riemann solver.

The code is based on the HLL solver, as presented in Toro, E.F. (2001) <em>Shock-Capturing Methods for Free-Surface Shallow Flows</em>. Great Britain: Wiley. Bed friction is included via a semi-implicit method. The bedslope term is incorporated using the method described by Valiani and Begnudelli (<em>J. Hydraul. Eng.</em>, 2006, 132(7): 652-665). 

Matlab scripts are provided for setting up and visualising simulations.

[IMAGES HERE]

# Installation

The following installation instructions are aimed at Civil Encgineering undergrad students.

1. Make sure you have a <b>gfortran compiler</b>, which is part of GCC (the GNU Compiler Collection). This is a free compiler for which this code has been tested (latest version tested is gcc 11.1.0, both for MacOS and Windows). A simple way of knowing wether you have this compiler is by openning a Terminal (MacOS) or Command Prompt (Windows) and typing `gfortran`. If you have the compiler, you will get a message like: "gfortran: fatal error: no input files" (i.e. the compiler exists but cannot find a file to compile). Otherwise you will get a message along the lines of "command not found/recognised". You may also try the command `gfortran -v` which will give you information of the version you have installed (if you have one). For Windows, I recommend the compiler by [equation.com](http://www.equation.com/servlet/equation.cmd?fa=fortran). 
