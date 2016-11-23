# -*- coding: utf-8 -*-
"""Hartigans dip test

Python implementation of the old FORTRAN code
HDT_Sig_batch has been augmented to calculate
significance more quickly when running on many rows

"""
from __future__ import absolute_import, print_function, division;

import numpy as np;
from . import ProgressBar;
from ..Global import RANDOM_SEED;
import logging;

def HDT_Sig_batch(xpdf_matrix, nboot, progressbar=True):
    """
    Saves time if calculating dip test on many sample-s of same size.
    Only generates the background distribution for significance once.

    :param xpdf_matrix: 2-dimensional numpy ndarray
    :param nboot: number of times to generate a random test statistic for p-value calculation
    :return: dips: dip statistic for each row in xpdf_matrix
                ps: p-value for each row in xpdf_matrix
                xlows: xlow for each row in xpdf_matrix
                xups: xup for each row in xpdf_matrix
    """

    np.random.seed(RANDOM_SEED);
    dips = np.zeros(xpdf_matrix.shape[0]);
    xlows = np.zeros(xpdf_matrix.shape[0]);
    xups = np.zeros(xpdf_matrix.shape[0]);

    if(progressbar): pbar = ProgressBar(xpdf_matrix.shape[0] + nboot);

    for i in range(xpdf_matrix.shape[0]):
        (dip, xlow, xup, ifault, gcm, lcm, mn, mj) = DipTest(xpdf_matrix[i,:]);
        dips[i] = dip;
        xlows[i] = xlow;
        xups[i] = xup;
        if(progressbar): pbar.update();


    bootDip=np.zeros(nboot);
    for i in np.arange(nboot):
        unifpdf=np.sort(np.random.rand(xpdf_matrix.shape[1]))
        bootDip[i] = DipTest(unifpdf)[0];
        if(progressbar): pbar.update();

    dips = np.expand_dims(dips, axis=1);        #Make dips Nx1
    bootDip = np.expand_dims(bootDip, axis=0);  #Make bootDip 1xnboot

    ps = np.sum(dips < bootDip, axis=1) / float(nboot);

    if(progressbar): pbar.complete();

    return(dips, ps, xlows, xups)

def HDT_Sig(xpdf,nboot):
    """Dip test with significance
    """

    (dip,xlow,xup,ifault,gcm,lcm,mn,mj)=DipTest(xpdf)

    bootDip=np.zeros(nboot);
    for i in np.arange(nboot):
        unifpdf=np.sort(np.random.rand(xpdf.shape[0]))
        bootDip[i] = DipTest(unifpdf)[0];
    
    p=np.sum(np.less(dip,bootDip))/float(nboot)
    
    return (dip,p,xlow,xup)

def DipTest(xpdf):
    """Hartigan's dip test. 

    This is a copy 
    """
    x=np.sort(xpdf)
    N=x.shape[0]
    mn=np.zeros(x.shape, dtype='int')
    mj=np.zeros(x.shape, dtype='int')
    lcm=np.zeros(x.shape, dtype='int')
    gcm=np.zeros(x.shape, dtype='int')
    ifault=False
    
    #Check that N is positive
    if N<=0:
        ifault=1
        print('\nHartigansDipTest.    InputError :  ifault=', ifault)
        return (0.0,0.0,0.0,ifault,gcm,lcm,mn,mj)
    
    #check if N is one
    if N==1:
        xl=x[0]
        xu=x[0]
        dip=0.0
        ifault=2
        print('\nHartigansDipTest.    InputError :  ifault=', ifault)
        return (dip,xl,xu,ifault,gcm,lcm,mn,mj)
    
    #check for case 1<N<4 or all identical values
    if N<=4 or x[N-1]==x[0]:
        xl=x[0]
        xu=x[0]
        dip=0.0
        ifault=4
        print('\nHartigansDipTest.    InputError :  ifault=', ifault)
        return (dip,xl,xu,ifault,gcm,lcm,mn,mj)
    
    #check if x is perfectly unimodal
    xsign=-np.sign(np.diff(np.diff(xpdf)))
    posi=np.greater(xsign,0.0)
    negi=np.less(xsign,0.0)
    if np.sum(posi)==0 or np.sum(negi)==0 or np.sum(np.less(posi,np.min(negi)))==N:
        #A unimodal function is its own best unimodal approximation, 
        #with a zero corresponding dip
        xl=x[0]
        xu=x[N-1]
        dip=0.0
        ifault=5
        print('\nHartigansDipTest.    InputError :  ifault=', ifault)
        return (dip,xl,xu,ifault,gcm,lcm,mn,mj)

    # LOW  contains the index of the current estimate of the lower 
    # end of the modal interval
    # HIGH contains the index of the current estimate of the upper 
    # end of the modal interval
    fn=N
    low=1
    high=N    
    dip=1./fn
    xl=x[low -1]
    xu=x[high -1]

    # establish the indices over which combination is necessary 
    # for the convex minorant fit
    mn[0]=1
    for j in range(2,N+1):
        mn[j-1]=j-1
        
        mnj=mn[j-1]
        mnmnj=mn[mnj-1]
        a=mnj-mnmnj
        b=j-mnj
        while not ((mnj==1) or (x[j-1]-x[mnj-1])*a < (x[mnj-1]-x[mnmnj-1])*b):
            mn[j-1]=mnmnj
            mnj=mn[j-1]
            mnmnj=mn[mnj-1]
            a=mnj-mnmnj
            b=j-mnj

    # establish the indices over which combination is necessary 
    # for the concave majorant fit
    mj[N-1]=N
    na=N-1
    for jk in range(1,na+1):
        k=N-jk
        mj[k-1]=k+1
        
        mjk=mj[k-1]
        mjmjk=mj[mjk-1]
        a=mjk-mjmjk
        b=k-mjk
        while not ( (mjk==N) or (x[k-1]-x[mjk-1])*a<(x[mjk-1]-x[mjmjk-1])*b):
            mj[k-1]=mjmjk
            mjk=mj[k-1]
            mjmjk=mj[mjk-1]
            a=mjk-mjmjk
            b=k-mjk

    itarate_flag=True

    iter_num = 0;
    MAX_ITER = 100;
    while itarate_flag:

        # Hartigans Dip Test has an issue where it may freeze
        # This establishes a max number of iterations and outputs 
        # An insignificant statistic so that the gene is rejected
        iter_num += 1;
        if(iter_num == MAX_ITER):
            dip = 1e99;
            logmessage = ', '.join([str(xx) for xx in xpdf]); 
            logging.info(logmessage);
            return (dip,xl,xu,ifault,gcm,lcm,mn,mj)

        ic=1
        gcm[0]=high
        igcm1=gcm[ic-1]
        ic+=1
        gcm[ic-1]=mn[igcm1-1]
        
        while gcm[ic-1]>low:
            igcm1=gcm[ic-1]
            ic+=1
            gcm[ic-1]=mn[igcm1-1]

        icx=ic

        # collect the change points for the LCM from LOW to HIGH
        ic=1
        lcm[0]=low
        lcm1=lcm[ic-1]
        ic+=1
        lcm[ic-1]=mj[lcm1-1]
        while lcm[ic-1]<high:
            lcm1=lcm[ic-1]
            ic+=1
            lcm[ic-1]=mj[lcm1-1]
            
        icv=ic

        # ICX, IX, IG are counters for the convex minorant
        # ICV, IV, IH are counters for the concave majorant
        ig=icx
        ih=icv
        
        # find the largest distance greater than 'DIP' 
        # between the GCM and the LCM from low to high
        
        ix=icx-1
        iv=2
        d=0.0
        
        if not (icx!=2 or icv!=2):
            d=1./fn
        else:
            iterate_BP50=True
            
            while iterate_BP50:
                igcmx=gcm[ix-1]
                lcmiv=lcm[iv-1]
                if not (igcmx > lcmiv):
                    # if the next point of either the GCM or 
                    # LCM is from the LCM then calculate distance 
                    #here OTHERWISE, GOTO BREAK POINT 55
                    
                    lcmiv1=lcm[iv-1-1]
                    a=lcmiv-lcmiv1
                    b=igcmx-lcmiv1-1
                    dx=(x[igcmx-1]-x[lcmiv1-1])*a/(fn*(x[lcmiv-1]-x[lcmiv1-1]))-b/fn
                    ix-=1
                    if dx<d:
                        goto60=True
                    else:
                        d=dx
                        ig=ix+1
                        ih=iv
                        goto60=True
                else:
                    # if the next point of either the GCM or 
                    # LCM is from the GCM then calculate distance 
                    # here CODE BREAK POINT 55
                    lcmiv=lcm[iv-1]
                    igcm=gcm[ix-1]
                    igcm1=gcm[ix+1-1]
                    a=lcmiv-igcm1+1
                    b=igcm-igcm1
                    dx=a/fn - ((x[lcmiv-1]-x[igcm1-1])*b)/(fn*(x[igcm-1]-x[igcm1-1]))
                    iv+=1
                
                    if not dx<d:
                        d=dx
                        ig=ix+1
                        ih=iv-1
                        
                    goto60=True

                if goto60:
                    if ix<1 : ix=1
                    if iv>icv : iv=icv
                    iterate_BP50 = gcm[ix-1] != lcm[iv-1]

        itarate_flag= not d<dip
        if itarate_flag:
            # if itarate_flag is true, then continue 
            # calculations and the great iteration cycle
            #if itarate_flag is NOT true, then stop 
            # calculations here, and break out of 
            #great iteration cycle to BREAK POINT 100
   
            # calculate the DIPs for the current LOW and HIGH

            #the DIP for the convex minorant
            dl=0.
            if ig !=icx:
                icxa=icx-1
                for j in range(ig,icxa+1):
                    temp=1./fn
                    jb=gcm[j+1-1]
                    je=gcm[j-1]
                    if not (je-jb<=1):
                        if not (x[je-1]==x[jb-1]):
                            a=je-jb
                            const=a/(fn*(x[je-1]-x[jb-1]))
                            for jr in range(int(jb),int(je+1)):
                                b=jr-jb+1
                                t=b/fn-(x[jr-1]-x[jb-1])*const
                                if (t>temp): temp=t
                    if dl<temp: dl=temp

            du=0.
            if not(ih==icv):
                icva=icv-1
                for k in range(ih,icva+1):
                    temp=1./fn
                    kb=lcm[k-1]
                    ke=lcm[k+1-1]
                    if not (ke-kb<=1):
                        if not (x[ke-1]==x[kb-1]):
                            a=ke-kb
                            const=a/(fn*(x[ke-1]-x[kb-1]))
                            for kr in range(int(kb),int(ke+1)):
                                b=kr-kb-1
                                t=(x[kr-1]-x[kb-1])*const-b/fn
                                if t>temp: temp=t
                    if du<temp: du=temp

            dipnew=dl
            if du>dl: dipnew=du
            if dip<dipnew: dip=dipnew
            low=gcm[ig-1]
            high=lcm[ih-1]
        #end if itarate_flag

    dip*=0.5
    xl=x[low-1]
    xu=x[high-1]
    
    return (dip,xl,xu,ifault,gcm,lcm,mn,mj)


if __name__=="__main__":
    """
    xpdf=np.array([-0.316502,
           -0.4760215,  
           -0.2745295,  
           -0.043012,  
           -0.439517,  
           0.033979,  
           0.1215105,  
           0.0189845,  
           -0.246039,  
           0.1554785,  
           -0.333038,  
           -0.3330585,  
           0.095484,  
           0.1445005,  
           -0.2930415,  
           -0.192532,  
           -0.011058,  
           -0.035037,  
           -0.3665275,  
           0.084498,  
           0.012455,  
           -0.3340235,  
           -0.2225445,  
           -0.1030285,  
           -0.161511,  
           -0.257043,  
           -0.257026,  
           -0.251022,  
           -0.2730175,  
           -0.017697,  
           0.00498,  
           -0.252074,  
           -0.1900265,  
           -0.288551,  
           -0.1440105,  
           -0.273011,  
           -0.30153,  
           -0.3870235,  
           -0.302027,  
           -0.2765615,  
           -0.3289815,  
           0.0689785,  
           -0.02304,  
           -0.2210335,  
           -0.175028,  
           -0.165512,  
           -0.077005,  
           -0.270514,  
           -0.3465375,  
           -0.242042,  
           0.298475,  
           -0.2780885,  
           0.063467,  
           -0.34901  ])
           """
    xpdf=np.random.randn(100)
    xpdf+=np.random.randn(100)+5
    print(HDT_Sig(xpdf,1000))
    #print DipTest(xpdf)
