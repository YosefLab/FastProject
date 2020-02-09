FastProject
===========

**Note: This repository is no longer undergoing active maintainence.**

**Instead, see [VISION](https://github.com/yoseflab/vision) for the successor to this project released in 2019.**
  - [https://github.com/yoseflab/VISION](https://github.com/yoseflab/vision)

Here we present FastProject, a software package for two-dimensional visualization of single cell data, which utilizes a plethora of projection methods and provides a way to systematically investigate the biological relevance of these low dimensional representations by incorporating domain knowledge.  

FastProject analyzes a gene expression matrix and produces a dynamic output report in which two-dimensional projections of the data can be explored.  Annotated gene sets (referred to as 'signatures') are incorporated so that features in the projections can be understood in relation to the biological processes they might represent.  FastProject provides a novel method of scoring each cell against a gene signature so as to minimize the effect of missed transcripts as well as a method to rank signature-projection pairings so that meaningful associations can be quickly identified. Additionally, FastProject is written with a modular architecture and designed to serve as a platform for incorporating and comparing new projection methods and gene selection algorithms.


Installing FastProject
-----------------
FastProject is written in Python (2.7 or 3.3+) and has a few dependencies.

To simplify the install process, we have created an install file for each platform that will take care of everything.  This method is recommended for users that do not regularly use Python.

*Download the file, and move into the folder where you'd like to install the application before running*

- [Linux](https://rawgit.com/YosefLab/FastProject/master/FP_Linux_Install.sh) (Install using `bash FP_Linux_Install.sh`)
- [OSX](https://rawgit.com/YosefLab/FastProject/master/FP_OSX_Install.sh) (Install using `bash FP_OSX_Install.sh`)
- [Windows](https://rawgit.com/YosefLab/FastProject/master/FP_Windows_Install.ps1) (Right click and select "Run with PowerShell")

This installs FastProject and all dependencies as an isolated Python environment.  See [this page](https://github.com/YosefLab/FastProject/wiki/Install-Instructions) in the wiki for instructions on how to install into an existing Python environment for development.

### Installing into an existing Python environment
A stable version of FastProject is maintained on the Python Package Index [https://pypi.python.org/pypi/FastProject/](https://pypi.python.org/pypi/FastProject/) and can be installed with:

    pip install FastProject

 
How to Run FastProject
----------------------

Sample FastProject command:
 
    fastproject <data file> -s <signature file(s)> -p <precomputed signature file(s)>
 
FastProject must be run with either the -s input or the -p input, though both can be used together.
For signature files and precomputed signature files, multiple files can be specified, separated by a space.
 
To display help and additional options: 
    fastproject -h
 
Sample Output
-------------
![FastProject Output Sample Image](/SampleOutput.png?raw=true)


Additional Info
---------------
All additional documentation is in the [FastProject Wiki](https://github.com/YosefLab/FastProject/wiki)
