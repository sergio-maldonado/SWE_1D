# SWE_1D
This is a Shallow Water Equations solver (in 1D), which employs the HLL approximate Riemann solver.

The code is based on the HLL solver, as presented in Toro, E.F. (2001) <em>Shock-Capturing Methods for Free-Surface Shallow Flows</em>. Great Britain: Wiley. Bed friction is included via a semi-implicit method. The bedslope term is incorporated using the method described by Valiani and Begnudelli <em>J. Hydraul. Eng.</em>, 2006, 132(7): 652-665. 

Matlab scripts are provided for setting up and visualising simulations.

[IMAGES HERE]

# Installation

The following installation instructions are aimed at Civil Engineering undergrad students.

1. Make sure you have a <b>gfortran compiler</b>, which is part of GCC (the GNU Compiler Collection). This is a free compiler for which this code has been tested (latest version tested is gcc 11.1.0, for both MacOS and Windows). A simple way of knowing wether you have this compiler is by openning a Terminal (MacOS) or Command Prompt (Windows) and typing `gfortran`. If you have the compiler, you will get a message like: "gfortran: fatal error: no input files" (i.e. the compiler exists but cannot find a file to compile). Otherwise you will get a message along the lines of "command not found/recognised". You may also try the command `gfortran -v` which will give you information on the version you have installed (if you have one). If you do not have the compiler installed, you will need to download and install it. For Windows, I recommend the compiler by [equation.com](http://www.equation.com/servlet/equation.cmd?fa=fortran).

2. Once you have downloaded this code, navigate to the <b>src</b> folder. You can do this from the Terminal/Command Prompt by typing `cd ` followed by the path; for example: `cd C:\Users\John\Documents\SWE_1D-main\src`. From here:
    - <u>If you are using <b>MacOS</b>:</u> just type `make` and that should be it (ignore warnings). An executable called <b>SWE1D_HLLC</b> should have been created in the same folder. Go straight to point 4 below.
    - <u>If you are using <b>Windows</b>:</u> you will need to make a few simple changes to the <b>Makefile</b> as described below.

3. <b>(for Windows only)</b> In the same <b>src</b> folder there is a file called <b>Makefile</b>. Open this file with any text editor (e.g. Notepad). Uncomment (remove the #) Line 21 and comment (add a #) Line 22. Then in Line 49 delete the following: <b>$(FLAGS)</b>. Save and close. Back in the Command Prompt type `make` (ignore warnings). If successful, a file called <b>SWE1D_HLLC.exe</b> should have been created in the same folder.

4. In the <b>src</b> folder there is a text file called <b>working_directory.inp</b>. Open this file with any text editor (e.g. Notepad or TextEdit) and modify the path to whichever path you plan to work in; for example: `"C:\Users\John\Documents\SWE_1D-main\runs"`.

5. Open the Matlab file (in the folder <b>scripts</b>) called <b>Tests_generate_IC.m</b>. What to do next should be self-explanatory (make sure you read the comments carefully). Choose any test and run the script. <b>NOTE:</b> If you activate the option to run automatically the SWE solver (from Matlab), go to the end of the script and comment/uncomment depending on whether you are using MacOS or Windows.

6. Visualise the results with the code <b>view_results.m</b> (make sure you change the path accordingly).
