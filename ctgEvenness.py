#!/home/gbonilla/.pyenv/shims/python
#/usr/bin/env python2.7

############################################################################
#################### Equitability of Mapped Reads to a contig ##############
############################################################################
# This module calculates the Shannon's evenness of the per-base coverage from
# a "samtools DEPTH" file. Then it calculates its Pielou's equitability as
# the ratio of the observed Shannon's evenness to the maximum evenness possible
# for a contig of the same length with the same number of total reads mapped.

# It is based on Pielou's equitability index:
# E = Ho / Hmax = -sum(pi*ln*pi) / lnL 
# where
# pi = numer of reads that cover contig site i
# L = contig length (number of sites in the contig)

import sys
import os
#print(sys.version)
import numpy as np

import scipy.stats
from scipy import stats
from scipy.stats import kurtosis
from scipy.stats import variation
from scipy.stats import skew
from scipy.stats import scoreatpercentile
import pandas
from pandas import read_csv

import math
from math import log
from math import pow
from math import sqrt
from math import ceil

from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

depthFile = str(sys.argv[1])
#outFile = str(sys.argv[2])

# Read the samtools depth file 
depthTab = pandas.read_table(depthFile, sep="\t", engine='python', header = None, names=["contig","site","cov"])    
ctgNames = depthTab.loc[1:1,"contig"]
ctgLengt = len(depthTab['cov'])
totReads = sum(depthTab['cov'])
covDensity = ( totReads/ctgLengt ) / 1000 # number of reads mapped per kilobase
covKurto = kurtosis(depthTab['cov'],fisher=True)
covMedia = np.median(depthTab['cov'])
covCVar = variation(depthTab['cov'])
covSkew = skew(depthTab['cov'])
covQ1 = stats.scoreatpercentile(depthTab['cov'],25)
covQ3 = stats.scoreatpercentile(depthTab['cov'],75)
# Quartile coefficient of dispersion: The closer the quartiles are (less disperesed), -> 0... the more away they are (more dispersed) -> 1
covQCD = ((covQ3 -covQ1)/2) / ((covQ1 + covQ3)/2)

obsSim = 0
obsEve = 0
coverdBases = 0

for index, row in depthTab.iterrows():
    cvg = float(row['cov'])
    if cvg != 0:
        coverdBases += 1
        pix = cvg / totReads
        obsEve += -math.log(pix,2)*pix
        #obsEve += -row['cov']*math.log(row['cov'],2)
        #obsSim += row['cov']**2
        obsSim += pix**2       

avgCover = math.ceil(coverdBases / totReads)
gnomeCov = float("{0:.10f}".format(coverdBases / ctgLengt))
#gnomeCov = round((coverdBases / ctgLengt),4)
#print("covered bases: {0} / contig length {1} = percentage covered genome {2}".format(coverdBases,ctgLengt,percGenomeCov))
ln2 = math.log(ctgLengt,2)
eqit = -obsEve/ln2
#maxEve = -ctgLengt*( (avgCover/totReads)*math.log((avgCover/totReads),2))
#maxSimp = ctgLengt*((avgCover/totReads)**2)
#equitSimp = obsSim/maxSimp

#print("contig \t ctgLen \t totReads \t CoverageDensity \t coverdBases \t percGenomeCov \t covMedian \t obsEve \t equitability \t covCVar \t covQCD \t covKurto \t covSkew" )
print("{0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7} \t {8} \t {9} \t {10} \t {11} \t {12}".format(ctgNames[1], ctgLengt, totReads, covDensity, coverdBases, gnomeCov, covMedia, obsEve, eqit, covCVar, covQCD, covKurto,covSkew ))
