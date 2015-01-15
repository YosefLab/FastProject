#! /usr/bin/env python
"""
Python wrapper to execute c++ tSNE implementation
for more information on tSNE, go to :
http://ticc.uvt.nl/~lvdrmaaten/Laurens_van_der_Maaten/t-SNE.html

HOW TO USE
Just call the method calc_tsne(dataMatrix)

Created by Philippe Hamel
hamelphi@iro.umontreal.ca
October 24th 2008

Modified by David DeTomaso
1/1/2015
To execute bh-tSNE


"""



from struct import *
import sys
import os
from numpy import *

def calc_tsne(dataMatrix,NO_DIMS=2,PERPLEX=30,INITIAL_DIMS=-1,THETA=1, DISTANCE_TYPE=1):
    """
    This is the main function.
    dataMatrix is a 2D numpy array containing your data (each row is a data point)
	THETA is a parameter for the gradient accuracy in bh-tSNE
	DISTANCE_TYPE chooses between probabilistic and euclidian distance
		0 = Euclidian
		1 = Probabilistic 
    """
    
    if DISTANCE_TYPE == 0 and INITIAL_DIMS != -1:
        dataMatrix=PCA(dataMatrix,INITIAL_DIMS)
    
    writeDat(dataMatrix,THETA,PERPLEX,DISTANCE_TYPE)
    tSNE()
    Xmat,LM,costs=readResult()
    #bh-tsne doesn't use landmarks, so no re-ording is necessary
    #clearData()
    #X=reOrder(Xmat,LM)  
    return Xmat


def PCA(dataMatrix, INITIAL_DIMS) :
    """
    Performs PCA on data.
    Reduces the dimensionality to INITIAL_DIMS
    """
    print 'Performing PCA'

    dataMatrix= dataMatrix-dataMatrix.mean(axis=0)

    if dataMatrix.shape[1]<INITIAL_DIMS:
        INITIAL_DIMS=dataMatrix.shape[1]

    (eigValues,eigVectors)=linalg.eig(cov(dataMatrix.T))
    perm=argsort(-eigValues)
    eigVectors=eigVectors[:,perm[0:INITIAL_DIMS]]
    dataMatrix=dot(dataMatrix,eigVectors)
    return dataMatrix

def readbin(type,file) :
    """
    used to read binary data from a file
    """
    return unpack(type,file.read(calcsize(type)))

def writeDat(dataMatrix,THETA,PERPLEX,DISTANCE_TYPE):
    """
    Generates data.dat
    """
    print 'Writing data.dat'
    print 'Theta : %f \nPerplexity : %f \nDistance Type : %i'%(THETA,PERPLEX,DISTANCE_TYPE)
    n,d = dataMatrix.shape
    f = open('data.dat', 'wb')
    f.write(pack('=iidd',n,d,THETA,PERPLEX))
    f.write(pack('=i',DISTANCE_TYPE))
    for inst in dataMatrix :
        for el in inst :
            f.write(pack('=d',el))
    f.close()


def tSNE():
    """
    Calls the tsne c++ implementation depending on the platform
    """
    platform=sys.platform
    print'Platform detected : %s'%platform
    if platform in ['mac', 'darwin'] :
        cmd='./tSNE_maci'
    elif platform == 'win32' :
        cmd= os.path.dirname(os.path.abspath(__file__)) + '\\bh_tsne_win64.exe'
    elif platform == 'linux2' :
        cmd='./tSNE_linux'
    else :
        print 'Not sure about the platform, we will try linux version...'
        cmd='./tSNE_linux'
    print 'Calling executable "%s"'%cmd
    os.system(cmd)
    

def readResult():
    """
    Reads result from result.dat
    """      
    print 'Reading result.dat'
    f=open('result.dat','rb')
    n,ND=readbin('ii',f)
    Xmat=empty((n,ND))
    for i in range(n):
        for j in range(ND):
            Xmat[i,j]=readbin('d',f)[0]
    LM=readbin('%ii'%n,f)
    costs=readbin('%id'%n,f)
    f.close()
    return (Xmat,LM,costs)

def reOrder(Xmat, LM):
    """
    Re-order the data in the original order
    Call only if LANDMARKS==1
    """
    print 'Reordering results'
    X=zeros(Xmat.shape)
    for i,lm in enumerate(LM):
        X[lm]=Xmat[i]
    return X

def clearData():
    """
    Clears files data.dat and result.dat
    """
    print 'Clearing data.dat and result.dat'
    platform=sys.platform
    print'Platform detected : %s'%platform
    if platform in ['mac', 'darwin'] :
        os.system('rm data.dat')
        os.system('rm result.dat')
    elif platform == 'win32' :
        os.system('del data.dat')
        os.system('del result.dat')
    elif platform == 'linux2' :
        os.system('rm data.dat')
        os.system('rm result.dat')
    else :
        os.system('rm data.dat')
        os.system('rm result.dat')
    
