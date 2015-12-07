FastProject
===========

Described here:

* Installing Python (required to run FastProject)
* Installing FastProject
* Running FastProject
  * Proper format for input files
  * Description of outputs


Installing Python
-----------------
### Instructions for Mac/Unix
Most Mac computers and Unix distributions come with python already installed.
 
FastProject requires Python 2.7.X   (where X can be any number)
 
To test if you already have Python installed, and test which version, open a terminal and enter this command:
 
    python --version    
 
If you do not have python installed, or have an old version, install the newest 2.7 version here:

* https://www.python.org/downloads/
 

### Instructions for Windows
FastProject requires Python 2.7.X   (where X can be any number)
 
To test if you already have Python installed, and test which version, open a command prompt and enter this command:
 
    python --version
 
If you do not have python installed, I recommend using the Anaconda Python distribution which makes the process of installing python dependencies much easier.

* https://www.continuum.io/downloads
 
Setting up FastProject
----------------------
### Installing from a tar or zip archive
 
If you are using the Anaconda distribution, run this command first to install dependencies:
 
    conda install numpy scipy matplotlib scikit-learn pandas
 
Then run:
 
    pip install FastProject-0.7.5-py2.7.tar.gz
 
If you are using OSX or Unix, just run
 
    pip install FastProject-0.7.5-py2.7.tar.gz 

How to Run FastProject
----------------------

Sample FastProject command:
 
    FastProject <data file> -s <signature file(s)> -p <precomputed signature file(s)> -k <housekeeping file> -o <output directory>
 
The -k and -o flags are optional.  Only need one of -s or -p though can include both
For signature files and precomputed signature files, multiple files can be specified, separated by a space.
 
To display help and additional options: 
    FastProject -h
 
*Note:  Depending on how your system's path is setup, the FastProject script might not be accessible.  If after installing FastProject you try to run "FastProject -h" and the system cannot find "FastProject", try running "python -m FastProject -h".  Then running FastProject will use this format instead:*
 
    python -m FastProject <data file> -s <signature file> -p <precomputed signature file> -k <housekeeping file> -o <output directory>

Input File Formats
------------------
### Data File
*Input data matrix of gene expression data for each sample.*
 
Text, Tab-delimited matrix with row and column labels
Rows are genes, columns are samples (cells)

First row in file is tab-separated list of sample labels
Then, each following row has:

    <Gene Name> TAB <Expression Value 1> TAB <Expression Value 2> TAB … (etc)
 
### Signature File (file extension matters, see below)
*List of gene signatures*
 
Text, tab-delimited
Genes should match the row labels in the input data matrix (matching is case insensitive)
 
#### Format A (File ends in .txt extension)
One gene per line (each signature has many lines)
 
    <Signature Name> TAB <Signature Sign> TAB <Gene Name>
 
    <Signature Sign> can be either "plus", "minus", or "both" if it's unsigned.

Alternately, you can just omit the second column and all genes will be treated as unsigned.
 
Example:

<table>
<tbody>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>plus</td><td>TEK</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>plus</td><td>AGPAT5</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>plus</td><td>HSPA4L</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>plus</td><td>FAM126A</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>minus</td><td>UBR1</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>minus</td><td>CYSLTR2</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>minus</td><td>ZNF205</td></tr>
<tr><td>Lin-neg_cell_vs_NKT_cell</td><td>minus</td><td>UBXN11</td></tr>
<tr><td>B_cell_vs_CD4T_cell</td><td>plus</td><td>RAB39</td></tr>
<tr><td>B_cell_vs_CD4T_cell</td><td>plus</td><td>EVI2B</td></tr>
<tr><td>B_cell_vs_CD4T_cell</td><td>plus</td><td>GALNS</td></tr>
<tr><td>B_cell_vs_CD4T_cell</td><td>plus</td><td>CRIP3</td></tr>
<tr><td>B_cell_vs_CD4T_cell</td><td>plus</td><td>HES6</td></tr>
<tr><td>B_cell_vs_CD4T_cell</td><td>plus</td><td>HMGXB4</td></tr>
<tr><td>...</td><td>...</td><td>...</td></tr>
</tbody>
</table>

 
#### Format B (File ends in .gmt)
One signature per line, (however "plus" and "minus" genes are split into two lines)
Each line: 

    <Signature Name> TAB <Signature Description> TAB <Gene1> TAB <Gene2> … (etc) 
 
The signature description is ignored by FastProject: that column just exists so the format conforms to the standard .gmt format.
 
To denote signature signs, use two lines to show the signature, with the "plus" genes in one and the "minus" genes in the other.  Add "_plus" to the signature name on the line with the "plus" genes and "_minus" to the signature name on the line with the "minus" genes.  

Example:

<tbody>
<table>
<tr><td>MEMORY_VS_NAIVE_CD8_TCELL_plus</td><td>GSE16522P</td><td>RHOC</td><td>OFD1</td><td>MLF1</td><td>...</td></tr>
<tr><td>MEMORY_VS_NAIVE_CD8_TCELL_minus</td><td>GSE16522</td><td>PTPRK</td><td>S100A5</td><td>IL1A</td><td>...</td></tr>
<tr><td>BCELL_VS_LUPUS_BCELL_plus</td><td>GSE10325</td><td>OXT</td><td>KCNH2</td><td>BTBD7</td><td>...</td></tr>
<tr><td>BCELL_VS_LUPUS_BCELL_minus</td><td>GSE10325</td><td>VAMP5</td><td>WSB2</td><td>CCR2</td><td>...</td></tr>
</tbody>
</table>

### Precomputed Signature File
*Covariates for each sample to also evaluate against projections.  This is where you can input information such as time or stimulus or batch*
 
Text, Tab-delimited

First line is a tab-separated list of sample labels.  Labels here should match sample labels in the expression data matrix file.
Each subsequent line contains:
 
    <Name> TAB <Signature Type> TAB <Value for sample 1> TAB <Value for sample 2> TAB … etc
 
<Signature Type> should be either “numerical” or “factor”
Numerical signature values must be parse-able to a float.  Factors can be anything (including numbers), and are treated differently with regards to calculating significance of signature/projection pairings.

Example:

<table>
<thead>
<tr><th></th><th></th><th>S1</th><th>S2</th><th>S3</th><th>S4</th><th>...</th></tr>
</thead>
<tbody>
<tr><td>Type</td><td>Factor</td><td>Wild</td><td>Wild</td><td>KO</td><td>Wild</td><td>...</td></tr>
<tr><td>Time</td><td>Numerical</td><td>0</td><td>4</td><td>4</td><td>6</td><td>...</td></tr>
<tr><td>Cluster</td><td>Factor</td><td>11</td><td>7</td><td>1</td><td>8</td><td>...</td></tr>
</tbody>
</table>
 
### Housekeeping Genes File
*List of genes to use as housekeeping genes for the estimation of False-Negative probabilities*

Text file containing one gene name per line
If omitted, FastProject uses the defaults stored in:

    <FastProject Root>/FastProject/Housekeeping Genes/

Outputs
-------

Outputs are located either in the directory specified by the -o option, or by default, a FastProject_Output directory is created.
The root of the directory contains Results.html.  This is the main output report. Running the file should launch your default web browser.

*Note: The “_resources” folder contains files needed by Results.html to render properly.  If needed, raw outputs (coordinates, signature scores, cluster assignments) can be found in the other sub-directories*
