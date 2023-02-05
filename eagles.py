# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 10:04:07 2022

@authors: Rob Jeffries, Richard Jackson

EAGLES - Estimating AGes from Lithium Equivalent widthS

A code that uses a model of equivalent width of lithium, as a function of age and Teff in
order to estimate the age of a star, or group of stars, with measured lithum and Teff.

See Jeffries et al. 2023, 
The Gaia-ESO Survey: Empirical estimates of stellar ages from lithium equivalent widths (EAGLES)
MNRAS submitted

"""

import sys
import getopt
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# This will get the name of this file
script_name = os.path.basename(__file__)

##
# @brief Help document for this script.  See the main function below.
#
help = f'''
    {script_name} input_file output_file  [-c] [-s] [-m] [-p prior] [-z constant] [--lagesmin min_age] [--lagesmax max_age] [--lapkmin min_peak_age] [--nage age_steps] [-h]"
    
    Estimates the ages for a set of stars or a cluster based on their lithium EWs and Teff.
     Results are output for individual stars in a csv table. The posterior probability distribution
     for an individual star or a cluster are also written to a csv file. Summary plots can be saved
     as pdf files.
    
    input_file is an ascii file containing >=1 row and 4 or 5 columns (with no header). .
        Col 1 is a string identifier with no spaces,
        Col 2 is an effective temperature in K,
        Col 3 is an *optional* 1-sigma error in Teff in K (see paper)
        Col 4 is Li equivalent width in milli-angstroms,
        Col 5 is the error in equivalent width in milli-angstroms. 
        
    output_file is the stem for the output files, with no suffix.
        <output_file>.csv is a csv containing the individual results (median, best-fit etc.) 
        for each star in the input_file. If -c is set then the last line will be the combined
        results for the cluster.
        
        The first columns are those of the input file. The last 6 columns are:
        lApk is the more probable age (peak of the posterior probability distribution)
        siglo is the lower bound of the asymmetric 68% confidence interval
        sighi is the upper bound of the asymmetric 68% confidence interval
        limup is the 95% age upper limit
        limlo is the 95% lower limit
        lmed is the median age
        A value of -1 means the value cannot be or is not calculated.
        
        <output_file>_pos.csv is a csv file containing P(log age) of the (combined) posterior, which
         is only produced if the input is a cluster (-c) or a single star.
                
    examples:
        Read data for a single star and estimate its age, (input_file would have one row),  
         based on a prior age probability that is flat in age, saving the output plots
        
            {script_name} input_file output_file -s -p 1
    
        Read data for a cluster of stars and calculate their combined posterior age
         probability distribution using a prior that is flat in log age. Saves a combined posterior
         distribution and a combined probability plot and one of these for each star in the input file.
        
            {script_name} input_file output_file -c -s -e
    
    -c 
    --cluster
        Indicates that the list of targets belong to a coeval cluster and their log 
         likelihoods will be summed to produce the final result. NB: output_file.csv 
         will still contain the results for each star individually.
    -s 
    --save
        Indicates that the plots should be saved: 
         <output_file>_prob.pdf shows the P(log age) plot if the input is a cluster or a single star; 
         <output_file>_iso.pdf also shows the EWLi vs Teff isochrone if there is >1 star 
         and the -c flag is set.    
        
    -e
    --every
        Indicates that a posterior age probability file should be saved for each star and,
         if -s is set, will save a plot of the probability distribution for each star.
         Files will be <output_file>_id_pos.csv and <outputfile>_id_prob.pdf, where id is the
         identifier of the star in the input file.
    -p 
    --prior
        Sets the prior age probability: 
         0 = flat in log age; 1 = flat in age.
        
        The default prior is 0 .
        
    -z 
        A likelihood regularisation constant to avoid outliers skewing the results
         badly in clusters. Only used if -c is set.
        
        The default value of z is 1.0e-12 .
        
    --lagesmin 
        Is the minimum log age/yr for the prior 
        
        The default value of lagesmin is 6.

    --lagesmax 
        Is the maximum log age/yr for the prior 
        
        The default value of lagesmax is 10.1.
        
    --lapkmin 
        Is the minimum log age/yr at which the likelihood is calculated.
         Likelihoods at lower ages are set to the value at this age.
        
        The default value of lapkmin is 6.699 .
        
    --nage 
        The number of log ages at which the posterior probability is calculated 
         between lagemin and lagemax.
        
        The default value of nage is 820, giving a log age bin size of 0.005 dex 
         for the defualt values of --lagesmin and --lagesmax

    -h 
    --help
        prints this message.
'''

def print_help():
    print(help, file=sys.stderr)


def AT2EWm(lTeff, lAge):

# The model for the Li EWs  (called EW_m in the paper)    
# Input arguments are log10(Teff/K) and log10(Age/yr)

    # constants defining model fit (Table 3)

    lTc =  3.524
    CT0 = -6.44
    CT1 =  4.06
    AAc =  291.3
    AA1 =  20699.0
    BBc =  0.111
    BBt = -164.8
    CCc =  7.131

    # calculate parameters

    AM    = AAc - AA1*(lTeff -lTc)*(lTeff -lTc)/(2*lTc)
    BM    = BBc - BBt*(lTeff -lTc)*(lTeff -lTc)/(2*lTc)
    CM    = np.zeros(np.shape(lTeff)[0])

    index = np.nonzero(lTeff <= lTc)
    if np.shape(index)[1] > 0 :
        CM[index]=CCc + CT0*(lTeff[index]-lTc)

    index = np.nonzero(lTeff > lTc)   
    if np.shape(index)[1] > 0 :
        CM[index]=CCc + CT1*(lTeff[index]-lTc)

    # return model Lithium EW

    return AM*(1-np.sinh((lAge-CM)/BM)/np.cosh((lAge-CM)/BM)) 



def eAT2EWm(lTeff, lAge):

# This is the model for the intrinsic dispersion (labelled rho_m in the paper)

    # constants defining model fit (Table 3)

    lTc =  3.524
    CT0 = -6.44
    CT1 =  4.06
    AAc =  291.3
    AA1 =  20699.0
    BBc =  0.111
    BBt = -164.8
    CCc =  7.131
    
    EE0 = 84.3
    EE1 = 1.8
    EE2 = 0.47
    FF0 = 0.079
    FF1 = 0.219
        
    # calculate parameters

    AM    = AAc - AA1*(lTeff -lTc)*(lTeff -lTc)/(2*lTc)
    BM    = BBc - BBt*(lTeff -lTc)*(lTeff -lTc)/(2*lTc)
    CM    = np.zeros(np.shape(lTeff)[0])

    
    index = np.nonzero(lTeff <= lTc)
    if np.shape(index)[1] > 0 :
        CM[index]=CCc + CT0*(lTeff[index]-lTc)

    index = np.nonzero(lTeff > lTc)   
    if np.shape(index)[1] > 0 :
        CM[index]=CCc + CT1*(lTeff[index]-lTc)

    # model dispersion (see eqn. 4)
    eEEM   = EE0/np.exp(EE2*(lAge-6)) + EE1

    eFFM = FF0*np.ones(np.shape(lTeff)[0])

    index = np.nonzero((lTeff < 3.716) & (lTeff > 3.62325))
    if np.shape(index)[1] > 0 :
      eFFM[index] = eFFM[index] + FF1/np.exp( (lAge - 8.0)**2/0.08)  
    
    eFFM = eFFM*(AM/BM)/np.cosh((lAge-CM)/BM)**2
    
    return np.sqrt(eEEM**2 + eFFM**2)
     


  
def bounds(lAge, like, prior=None, lApkmin=None, pcuthi=None, limcut=None):
    
    # Determines the parameters of the posterior distribution after
    # applying the prior to the likelihood distribution.

    #INPUT
    # lAge    - array of uniformly spaced log ages
    # like    - likelihood at each value of lage
    #OUTPUT (Returns)
    # prob    - likelihood weighted by prior
    # p  - a numpy array containing the following elements
    # p[0] lApk    - the log age of the peak of the posterior
    # p[1] siglo   - offset in logage of 68% lower bound
    # p[2] sighi   - offset in logage of 68% upper bound
    # p[3] limup   - 95% upper limit where no clear peak
    # p[4] limlo   - 95% lower limit where no clear peak
    # p[5] lmed    - the median log Age from the probability distribution
    # values set to -1 if  parameter  undefined
    #OPTIONS
    # prior   - constant defining weighting of age scale
    #         - 0 for uniform with log10(age) 1 for uniform with age
    # lApkmin - lowest age for valid model of EW_Li
    #         - Likelihood is constant below this age
    # pcuthi  - level defining resolvablelog age decrease at limit
    # limcut  - fraction defining upper/lower limit    
    # median - if set to true calculates median log age, else returns most probable log age

    # default limits
    if prior is None:
        prior = 0
    if lApkmin is None:
        lApkmin = np.log10(5)+6
    if pcuthi is None:
        pcuthi = np.exp(-0.5)
    if limcut is None:
        limcut = 0.95

    # default values of output array = -1 until set
    p = np.full(6, -1.0)

 
    # fill low end of likelihood curve below lApkmin
    Nmax = np.shape(lAge)[0]
    index = np.nonzero(lAge < lApkmin)
    Nlo = np.shape(index)[1]
    if Nlo > 0 :
        like[0:Nlo] = like[Nlo-1]

    # scale likelihood by the prior
    prob = like*10**(prior*lAge) 

    # find the median log age and its index 
    cdf = np.cumsum(prob)/np.sum(prob)
    index = np.nonzero(cdf  <= 0.5)
    Npk = np.shape(index)[1]
    p[5] = lAge[Npk]
    
    # now find most probable age
    Npk = np.argmax(prob)
    p[0] = lAge[Npk]
    
    # handle case of very low lApk
    if Nlo > Npk :
        Npk = Nlo 
            
    # define max values
    
    Pmax = prob[Npk]
    lAstep =  lAge[Nmax-1]-lAge[Nmax-2]
#    print('lApk, lMed, Npk, Nlo', p[0], p[5], Npk, Nlo)

    # case with resolvable peak
    if prob[Nlo]/Pmax < pcuthi and prob[Nmax-1]/Pmax < pcuthi :
        
        # get lower bound
        plo = prob[0:Npk][::-1]  # reverses a slice of prob
        cdflo = np.cumsum(plo)/np.sum(plo)
        index = np.nonzero(cdflo > 0.68)
        indlo = np.shape(index)[1]
        p[1] = lAge[Npk] - lAge[indlo] + lAstep

        # get upper bound
        phi = prob[Npk:Nmax]
        cdfhi = np.cumsum(phi)/np.sum(phi)
        index = np.nonzero(cdfhi <= 0.68)
        indhi = np.shape(index)[1]
        p[2] = lAge[indhi + Npk] - lAge[Npk] + lAstep
      
    # case of upper limit
    if prob[Nlo]/Pmax > pcuthi and prob[Nmax-1]/Pmax < pcuthi :
        cdf = np.cumsum(prob)/np.sum(prob)     
        index = np.nonzero(cdf <= limcut)
        Ncut = np.shape(index)[1]
        p[3] = lAge[Ncut] + lAstep
        p[0] = -1
        
    # case of lower limit
    if prob[Nlo]/Pmax < pcuthi and prob[Nmax-1]/Pmax > pcuthi :
        cdf = np.cumsum(prob)/np.sum(prob)
        index = np.nonzero(cdf <= 1.0-limcut)
        Ncut = np.shape(index)[1]
        p[4] = lAge[Ncut] - lAstep
        p[0] = -1
        
    return prob, p   
  
def age_fit(lTeff, LiEW, eLiEW, lAges, z, lApkmin, nTeff, prior=None, elTeff=None):

    # generates the likelihood distribution
    # and determines the parameters of the posterior age distribution    
    
    # Outputs
    # llike  - the log likelihood distribution (normalised to a peak of zero)
    # lprob  - the log of the posterior age probability (normalised to a peak of zero)
    # p  - a numpy array containing the following elements
    # p[0] lApk    - the log age of the peak of the posterior
    # p[1] siglo   - offset in logage of 68% lower bound
    # p[2] sighi   - offset in logage of 68% upper bound
    # p[3] limup   - 95% upper limit where no clear peak
    # p[4] limlo   - 95% lower limit where no clear peak
    # p[5] lmed    - the median log Age from the probability distribution
    # chisq - a chi-squared for the combined fit of an isochrone to cluster data
    
    # default value of prior
    if prior is None:
        prior = 0      
    nStar = np.shape(LiEW)[0]
    nAges = np.shape(lAges)[0] 
    pfit = np.zeros(nAges) 
    dpfit = np.zeros(nAges)

    if elTeff is None:
        
        # get raw log likelihood at each log age 
        for k in range(0, nAges):
            EWm = AT2EWm(lTeff, lAges[k])
            xLiEW = eAT2EWm(lTeff, lAges[k])
            sLiEW = np.sqrt(xLiEW**2 + eLiEW**2)
            pstar = 1.0/np.sqrt(2.0*np.pi)/sLiEW* \
                np.exp(-(LiEW-EWm)*(LiEW-EWm)/(2.0*sLiEW*sLiEW)) + z
            pfit[k] = np.sum(np.log(pstar))
            
    else:
        
        # need to loop over a set of temperature AND over a set of ages
        dlTeff = 4.0*elTeff/nTeff
        
        for j in range(0, nTeff):
            # newTeff is the new (log) Teff to evaluate the likelihood
            newTeff = lTeff + (j - (nTeff-1)/2)*dlTeff
            # factor is the number by which the log likelihood is decreased
            # at that Teff. Assumes a Gaussian distribution.
            factor = -1.0*(lTeff-newTeff)**2/(2.0*elTeff**2)
            
            for k in range(0, nAges):
                EWm = AT2EWm(newTeff, lAges[k] )
                xLiEW = eAT2EWm(newTeff, lAges[k])
                sLiEW = np.sqrt(xLiEW**2 + eLiEW**2)
                pstar = 1.0/np.sqrt(2.0*np.pi)/sLiEW* \
                    np.exp(-(LiEW-EWm)*(LiEW-EWm)/(2.0*sLiEW*sLiEW)) + z
                dpfit[k] = np.sum(np.log(pstar*np.exp(factor)))
            
            pfit = pfit + np.exp(dpfit)
            
        pfit = np.log(pfit)
            
            
    # handle low likelihoods (set ln L to -70)
    llike = pfit - np.amax(pfit)
    index = np.nonzero(llike < -70)
    if np.shape(index)[1] > 0:
        llike[index] = -70

    # handle low ages - set likelhiood to have its value at lApkmin
    like = np.exp(llike - np.amax(llike))
    index = np.nonzero(lAges < lApkmin)
    Nlo = np.shape(index)[1]
    if Nlo > 0:
        like[0:Nlo-1] = like[Nlo-1]
        
    # normalise    
    like = like/np.sum(like)
    
    # get posterior, peak age, median, bounds and limits. NB This also applies the prior.
    prob, p = bounds(lAges, like, prior=prior)
    
    # log probability
    lprob = np.log(prob) - np.amax(np.log(prob))
    
    # get reduced chisqr if a cluster
    if p[0] > -1 and nStar > 1 :
        xLiEW = eAT2EWm(lTeff, p[0])
        sLiEW = np.sqrt(xLiEW**2 + eLiEW**2)
        dLiEW = (AT2EWm(lTeff, p[0]) - LiEW)
        chisq = np.sum(dLiEW**2/sLiEW**2)/nStar
    else:
        chisq = -1
    
    # normalise peak of log likelihood to zero
    llike = np.log(like) - np.amax(np.log(like))  
    
    return llike, lprob, p, chisq


def get_li_age(LiEW, eLiEW, Teff, eTeff=None, lagesmin=6.0, lagesmax=10.1, 
               lApkmin=6.699, nAge=820, z=1.0e-12, nTeff=21, prior=None) :
    
 
    # setup some default values
    
    # default value of prior
    if prior is None:
        prior = 0     
    
    # bounds for Teff
    Teffmin = 3000.0
    Teffmax = 6500.0
    
    lagesmin = lagesmin # fitted function is in terms of log Age/yr
    lagesmax = lagesmax
    lApkmin = lApkmin
    nStar = np.shape(LiEW)[0] # number of stars in input file

    if (np.amax(Teff) > Teffmax or np.amin(Teffmin) < Teffmin):
        raise RuntimeError("All temperatures must be between 3000K and 6500K")
    
    if (np.amax(LiEW-eLiEW) > 800 or np.amin(LiEW+eLiEW) < -200):
        raise RuntimeError("LiEW outside a sensible range")

# set array of ages and temperatures
    lAges = lagesmin+(np.arange(nAge)+0.5)*(lagesmax-lagesmin)/float(nAge)
    lTeff = np.log10(Teff)    
    
    # setup an array of log Teff around the central value if necessary
    if eTeff is not None:
        elTeff = 0.5*(np.log10(Teff+eTeff)-np.log10(Teff-eTeff))
    else:
        elTeff=None
  
# get the likelihood, posterior, ages and limits, and chi-squared (for a cluster)
    llike, lprob, p, chisq = \
      age_fit(lTeff, LiEW, eLiEW, lAges, z, lApkmin, nTeff, prior=prior, elTeff=elTeff)


    # outputs
    # lAges - the array of logarithmic ages used
    # llike - the ln likelihood distribution
    # lprob - the ln posterior age probability distribution
    # p  - a numpy array containing the following elements
    # p[0] lApk    - the log age of the peak of the posterior
    # p[1] siglo   - offset in logage of 68% lower bound
    # p[2] sighi   - offset in logage of 68% upper bound
    # p[3] limup   - 95% upper limit where no clear peak
    # p[4] limlo   - 95% lower limit where no clear peak
    # p[5] lmed    - the median log Age from the probability distribution
    # chisq - a chi-squared for the combined fit of an isochrone to cluster data
    
    return lAges, llike, lprob, p, chisq


def make_plots(lAges, lprob, p, chisq, lagesmin, lagesmax, \
               ID, LiEW, eLiEW, Teff, filename, is_cluster=False, savefig=False):


    
    # produce a probability plot, either for the cluster or for an individual star
    nStar = np.shape(LiEW)[0]
    plt.xlabel('log Age/yr')
    plt.ylabel('P(log Age)')

    if is_cluster and nStar > 1 : # case of a cluster of >1 stars
        index = np.nonzero(lprob > -20)
        nlo = np.amin(index)
        nhi = np.amax(index)
        plt.xlim([lAges[nlo], lAges[nhi]])
        label = 'Cluster'
    else: # case of a single star
        plt.xlim([lagesmin,lagesmax])
        label = ID
    
    # plot the probability distribution
    ax = plt.gca()
    plt.ylim([-0.05, 1.2])
    plt.plot(lAges, np.exp(lprob), c='b')
    ax.text(0.02,0.94, label, transform=ax.transAxes)

    # create shaded areas for the limits
    # handle cases where a peak is found and also limits
    if p[0] != -1: # a peak is found
        plt.fill_between(lAges, np.exp(lprob), \
            where = (lAges > p[0]-p[1]) & (lAges < p[0]+p[2]), color='0.7' )
        
        ax.text(0.02,0.87, "Age: {:.1f} +{:.1f}/-{:.1f} Myr".format(10**(p[0]-6), \
                10**(p[0]+p[2]-6)-10**(p[0]-6), 10**(p[0]-6)-10**(p[0]-p[1]-6)), \
                transform=ax.transAxes)  
    
    if p[3] != -1: # an upper limit
        plt.fill_between(lAges, np.exp(lprob), \
            where = (lAges < p[3]), color='0.7' )
        ax.text(0.02,0.87, "Age < {:.1f} Myr".format(10**(p[3]-6)), transform=ax.transAxes)

    if p[4] != -1: # a lower limit
        plt.fill_between(lAges, np.exp(lprob), \
            where = (lAges > p[4]), color='0.7' )
        ax.text(0.02,0.87, "Age > {:.1f} Myr".format(10**(p[4]-6)), transform=ax.transAxes)

    # find and plot the median age    
    index = np.nonzero(lAges <= p[5])
    nmed = np.shape(index)[1]
    plt.plot([p[5],p[5]], [0.0,np.exp(lprob[nmed])], 'k--')  
    
    if savefig:    
        if is_cluster :
            plt.savefig(filename+'_prob.pdf')
        else:
            plt.savefig(filename+'_'+ID+'_prob.pdf')
    plt.show()
    
    # produce a plot of LiEW vs Teff with the best-fitting isochrone also showing the 
    # intrinsic dispersion
    # NB only makes sense for a cluster, not for a file of unrelated stars or single star
    if (is_cluster and nStar > 1):    
        
        # set the Teff array and limits at which to plot isochrones
        lTeff = np.arange(3.4771, 3.8130, 0.001) #between 3000K and 6500K
        plt.xlim(6600, 2900)
        
       # if nStar == 1 :
       #     lTeff = np.arange(max(np.log10(Teff*0.9),3.4771), \
       #                       min(np.log10(Teff*1.1),3.8130), 0.001)   
       #     plt.xlim(min(6600, Teff*1.1), max(2900, Teff*0.9))
       # else :
       #     lTeff = np.arange(3.4771, 3.8130, 0.001) #between 3000K and 6500K
       #     plt.xlim(6600, 2900)
        
        ax = plt.gca()

        # plot an isochrone and the dispersion around the isochrone
        # handle cases where age is found or just limits
        if p[0] != -1:
            ax.text(0.02,0.94, "Most probable Age: {:.1f} +{:.1f}/-{:.1f} Myr".format(10**(p[0]-6), \
                    10**(p[0]+p[2]-6)-10**(p[0]-6), 10**(p[0]-6)-10**(p[0]-p[1]-6)), \
                    transform=ax.transAxes)
            ax.text(0.02,0.87, "Chisqr: {:.2f}, shaded = dispersion at best-fit age".format(chisq), \
                    transform=ax.transAxes)
            # calculate isochrone and dispersion around the isochrone    
            EWm = AT2EWm(lTeff, p[0])
            EWm_hi = EWm + eAT2EWm(lTeff, p[0])
            EWm_lo = EWm - eAT2EWm(lTeff, p[0])
            
        if p[3] != -1:
            ax.text(0.02,0.94, "Age < {:.1f} Myr (95%)".format(10**(p[3]-6)), transform=ax.transAxes)
            ax.text(0.02,0.87, "shaded = dispersion at upper limit", \
                    transform=ax.transAxes)
            # calculate isochrone and dispersion at upper limit
            EWm = AT2EWm(lTeff, p[3])
            EWm_hi = EWm + eAT2EWm(lTeff, p[3])
            EWm_lo = EWm - eAT2EWm(lTeff, p[3])

        if p[4] != -1:
            ax.text(0.02,0.94, "Age > {:.1f} Myr (95%)".format(10**(p[4]-6)), transform=ax.transAxes)
            ax.text(0.02,0.87, "shaded = dispersion at lower limit", \
                    transform=ax.transAxes)
            # calculate isochrone and dispersion at lower limit    
            EWm = AT2EWm(lTeff, p[4])
            EWm_hi = EWm + eAT2EWm(lTeff, p[4])
            EWm_lo = EWm - eAT2EWm(lTeff, p[4])

        plt.xlabel('Teff (K)')
        plt.ylabel('LiEW (mA)')
        plt.ylim(min(np.min(EWm_lo),-10, np.min(LiEW-eLiEW)), 1.2*max(np.max(EWm_hi), np.max(LiEW+eLiEW)))
        # plot the dispersion, the data and the isochrone
        plt.fill_between(10**lTeff, EWm_lo, EWm_hi, color='0.8')
        plt.errorbar(Teff, LiEW, yerr=eLiEW, color='b', fmt='.')
        plt.plot(10**lTeff, EWm, color='k')

        # plot dashed isochones at the upper and lower error bar/limit
        if p[0] != -1:
            plt.plot(10**lTeff, AT2EWm(lTeff, p[0]+p[2]), 'k--')
            plt.plot(10**lTeff, AT2EWm(lTeff, p[0]-p[1]), 'k--')

        if p[3] != -1:
            plt.plot(10**lTeff, AT2EWm(lTeff, lagesmin), 'k--' )

        if p[4] != -1:
            plt.plot(10**lTeff, AT2EWm(lTeff, lagesmax), 'k--' )

        if savefig:
            plt.savefig(filename+"_iso.pdf")
                       
        plt.show()

# Main code

def main(argv):
        
    # default values - see help file
    prior=0
    lagesmin = 6.0
    lagesmax = 10.1
    lApkmin = np.log10(5)+6    
    nAge = 820
    nTeff = 21 # not currently an optional parameter - number of integration steps
               # if marginalising over supplied Teff uncertainties
    z = 1.0e-12
    eTeff = None # over-ridden if input has 5 columns
    is_cluster = False # over-ridden by the -c flag
    savefig = False # over-ridden by the -s flag
    everystar = False # over-ridden by the -e flag
    
    # check there are at least two arguments
    if (len(argv) < 2 or '-h' in argv):
        print_help()
        sys.exit(())

    # open the input file and read the data
    inputfile = argv[0]
    if os.path.isfile(inputfile): # check file exists
                    
        with open(inputfile) as f:
            line = f.readline()
            ncol = len(line.split())
            f.close()
            if ncol == 4 or ncol == 5: # check if there are 4 or 5 columns
                if ncol == 4: # no Teff errors
                    cols = ['ID', 'Teff', 'LiEW', 'eLiEW']  
                            
                if ncol == 5: # additional Teff errors
                    cols = ['ID', 'Teff', 'eTeff', 'LiEW', 'eLiEW']
                       
                df = pd.read_csv(inputfile, header=None, delim_whitespace=True, names=cols) 
                IDraw = df[['ID']].to_numpy()
                ID = np.reshape(IDraw, len(IDraw))
                Teffraw = df[['Teff']].to_numpy()
                Teff = np.reshape(Teffraw, len(Teffraw))
                LiEWraw = df[['LiEW']].to_numpy()
                LiEW = np.reshape(LiEWraw, len(LiEWraw))
                eLiEWraw = df[['eLiEW']].to_numpy()
                eLiEW = np.reshape(eLiEWraw, len(eLiEWraw))
                nStar = np.shape(LiEW)[0]
                if ncol == 5:
                    eTeffraw = df[['eTeff']].to_numpy()
                    eTeff = np.reshape(eTeffraw, len(eTeffraw))
                        
            else:
                print("Wrong number of columns in input file!\n")
                sys.exit()
                    
    else:
        raise Exception("Input file does not exist") 

    # get the output filename stem; NB no check for valid path
    filename = argv[1]

    try:
        options, args = getopt.getopt(argv[2:], "p:z:cseh", \
                    ["help", "cluster", "save", "every", "prior=", "lagesmin=", \
                     "lagesmax=", "lapkmin=", "nage="])

    except:
        print_help()
    
    try:
        for name, value in options:
        
            if name in ['-c', '--cluster']:
                is_cluster = True
                if nStar < 2:
                    print("Warning: a cluster should have more than one star!")
       
            if name in ['-s', '--save']:
                savefig = True
                           
            if name in ['-e', '--every']:
                everystar = True
                print('Warning: -e can produce lots of files!')
                
            if name in ['-p', '--prior']:
                prior = int(value)
                if prior != 0 and prior != 1:
                    print("Prior must be either 0 or 1\n")
                    sys.exit()
                    
            if name in ['-h', '--help']:
                print_help()
                sys.exit()
    
            if name in ['-z']:
                z = float(value)
                print("Note that z will be reset to zero if the input file contains only 1 row\n")
    
            if name in ['--nage']:
                nAge = int(value)
                
            if name in ['--lagesmin']:
                lagesmin = float(value)
                if lagesmin < 6.0 :
                    print("Warning: lagesmin is very young!\n")
                
            if name in ['--lagesmax']:
                lagesmax = float(value)
                if lagesmax > 10.1 :
                    print("Warning: lagesmax is older than the Galaxy!\n")
            
            if lagesmax < lagesmin :
                print("lagesmax cannot be less than lagesmin!\n")
                sys.exit()
                
            if name in ['lapkmin']:
                lApkmin = float(value)
                if lApkmin < 6.3 :
                    print("Warning: lApkmin is very low.\n \
                          Likelihood is completely unconstrained below 2 Myr!\n")
  
    except (IndexError, ValueError):
        print_help()
     
    
    
    # setup the output lists
    l_lApk = []
    l_siglo = []
    l_sighi = []
    l_limup = []
    l_limlo = []
    l_lmed = []
    
    
    # Do the calculations in two passes
    # First treat each star as an individual and write the results into
    # lists that will be output to filename.csv
    
    for i in range(0, nStar):
        if ncol == 4:
            lAges, llike, lprob, p, chisq = \
            get_li_age(LiEW[i:i+1], eLiEW[i:i+1], Teff[i:i+1], lagesmax=lagesmax, lagesmin=lagesmin, \
                       lApkmin=lApkmin, z=0.0, nAge=nAge, prior=prior)
        
        if ncol == 5:
            lAges, llike, lprob, p, chisq = \
            get_li_age(LiEW[i:i+1], eLiEW[i:i+1], Teff[i:i+1], lagesmax=lagesmax, lagesmin=lagesmin, \
                       lApkmin=lApkmin, z=0.0, nAge=nAge, prior=prior, eTeff=eTeff[i:i+1])
 
        # make plots for the individual stars and save them if required
        if everystar :
            make_plots(lAges, lprob, p, chisq, \
                       lagesmin, lagesmax, ID[i], LiEW, eLiEW, Teff, filename, is_cluster=False, savefig=savefig)  

            # write out a posterior probability distribution for each star
            # format is log Age/yr probability (normalised to 1 at its peak)
            data = np.column_stack([lAges, np.exp(lprob)])
            np.savetxt(filename+"_"+ID[i]+"_pos.csv", data, fmt=['%7.4f', '%7.3e'], 
                           delimiter = ',', header= 'log (Age/yr)  Probability')

        l_lApk.append(p[0])
        l_siglo.append(p[1])
        l_sighi.append(p[2])
        l_limup.append(p[3])
        l_limlo.append(p[4])
        l_lmed.append(p[5])
 
    
 
    # Append new columns to the input dataframe
    # this is done whether they are in a cluster or not
    
    df['lApk'] = l_lApk
    df['siglo'] = l_siglo
    df['sighi'] = l_sighi
    df['limup'] = l_limup
    df['limlo'] = l_limlo
    df['lMed'] = l_lmed
    
    #Then if it is a cluster run again with the full list combined and
    # print out all the combined numerical results to the screen
 
    print("*********************************************")
    print("****************EAGLES V1.0******************")
    print("*********************************************\n")
    if is_cluster and nStar > 1:
        
        if ncol == 4:
            lAges, llike, lprob, p, chisq = \
                get_li_age(LiEW, eLiEW, Teff, lagesmax=lagesmax, lagesmin=lagesmin, \
                       lApkmin=lApkmin, z=z, nAge=nAge, prior=prior)
            # append the combined results for the cluster as the last tow in the dataframe
            df.loc[len(df)] = ['Cluster', 0, 0, 0, p[0], p[1], p[2], p[3], p[4], p[5]]
       
        if ncol == 5:
            lAges, llike, lprob, p, chisq = \
                get_li_age(LiEW, eLiEW, Teff, lagesmax=lagesmax, lagesmin=lagesmin, \
                       lApkmin=lApkmin, z=z, nAge=nAge, prior=prior, eTeff=eTeff)
            # append the combined results for the cluster as the last tow in the dataframe
            df.loc[len(df)] = ['Cluster', 0, 0, 0, 0, p[0], p[1], p[2], p[3], p[4], p[5]]
        


        # make the plots for a cluster and save if required
        make_plots(lAges, lprob, p, chisq, \
                   lagesmin, lagesmax, ID, LiEW, eLiEW, Teff, filename, is_cluster, savefig)

        print('Cluster of %i stars' % nStar)
        print('chi-squared of fit = %6.2f' % chisq)
        
    # print out the results for a cluster (or if there's only one star)
    if is_cluster or nStar == 1 :
       
        if nStar == 1 :
            print('LiEW = (%5.1f+/-%5.1f) mA Teff = %6.1f K' % (LiEW[0], eLiEW[0], Teff[0]))
        
        if p[0] != -1 :
            print('most probable log (Age/yr) = %5.3f +%5.3f/-%5.3f' % (p[0], p[2], p[1]))
            print('most probable Age (Myr) = %6.1f +%6.1f/-%6.1f' % (10**(p[0]-6), 10**(p[0]+p[2]-6)-10**(p[0]-6), \
                                10**(p[0]-6)-10**(p[0]-p[1]-6)))          

        if p[3] != -1 :
            print('log (Age/yr) < %4.2f' % (p[3]))
            print('Age (Myr) < %6.1f (95 per cent limit)' % 10**(p[3]-6))

        if p[4] != -1 :
            print('log (Age/yr) > %4.2f' % (p[4]))
            print('Age (Myr) > %6.1f (95 per cent limit)' % 10**(p[4]-6))
        
        print('median log (Age/yr) = %5.3f' % (p[5]))
        print('median Age (Myr) = %6.1f' % (10**(p[5]-6)))


        # write out the combined posterior probability to a csv file
        # format is log Age/yr probability (normalised to 1 at its peak)
        if (is_cluster) or nStar == 1 :
        
            data = np.column_stack([lAges, np.exp(lprob)])
            np.savetxt(filename+"_pos.csv", data, fmt=['%7.4f', '%7.3e'], 
                       delimiter = ',', header= 'log (Age/yr)  Probability')


    print("*********************************************\n")
    

    
    df.to_csv(filename+'.csv', float_format="%7.3f", index=False) 
    
    
###THE END###


if  __name__ == "__main__":   
    main(sys.argv[1:])



