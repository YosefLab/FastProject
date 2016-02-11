FastProject
===========
Here we present FastProject, a software package for two-dimensional visualization of single cell data, which utilizes a plethora of projection methods and provides a way to systematically investigate the biological relevance of these low dimensional representations by incorporating domain knowledge.  

FastProject analyzes a gene expression matrix and produces a dynamic output report in which two-dimensional projections of the data can be explored.  Annotated gene sets (referred to as 'signatures') are incorporated so that features in the projections can be understood in relation to the biological processes they might represent.  FastProject provides a novel method of scoring each cell against a gene signature so as to minimize the effect of missed transcripts as well as a method to rank signature-projection pairings so that meaningful associations can be quickly identified. Additionally, FastProject is written with a modular architecture and designed to serve as a platform for incorporating and comparing new projection methods and gene selection algorithms.


Installing FastProject
-----------------
FastProject is written in Python 2.7 and has a few dependencies.  To simplify the install process, we have created an install file for each platform that will take care of everything.

*Download the file, and move into the folder where you'd like to install the application before running*

- [Linux](https://rawgit.com/YosefLab/FastProject/master/FP_Linux_Install.sh) (Install using `bash FP_Linux_Install.sh`)
- [OSX](https://rawgit.com/YosefLab/FastProject/master/FP_OSX_Install.sh) (Install using `bash FP_OSX_Install.sh`)
- [Windows](https://rawgit.com/YosefLab/FastProject/master/FP_Windows_Install.ps1) (Right click and select "Run with PowerShell")

This installs FastProject and all dependencies as an isolated Python environment.  See [this page](https://github.com/YosefLab/FastProject/wiki/Install-Instructions) in the wiki for instructions on how to install into an existing Python environment for development.
 
How to Run FastProject
----------------------

Sample FastProject command:
 
    fastproject <data file> -s <signature file(s)> -p <precomputed signature file(s)> -k <housekeeping file> -o <output directory>
 
The -k and -o flags are optional.  Only need one of -s or -p though can include both
For signature files and precomputed signature files, multiple files can be specified, separated by a space.
 
To display help and additional options: 
    fastproject -h
 
Sample Output
-------------

Additional Info
---------------
All additional documentation is in the [FastProject Wiki](https://github.com/YosefLab/FastProject/wiki)
