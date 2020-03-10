#!/home/gbonilla/.pyenv/shims/python
#/usr/bin/env python3
#/home/gbonilla/Software/Komodo-IDE-10.2.2-89895-linux-x86_64/INSTALLDIR/lib/python/bin
############################################################################
        #################### ANI/AAI file parser ##################
############################################################################
# This script was develped for analysing the ANI/AAI values of viral clusters obtained
# from metagenomic datasets through vContACT. ANI/AAI values are expected in Rodriguez-Konstantiinidis format of enveomics.
# Usage: ANIparser.py [full path to ANIvalues.out] [AAI/ANI]

import sys
#print(sys.version)
import os
import subprocess
aniFile = str(sys.argv[1])
aaiFile = str(sys.argv[2])
#targetPath = os.listdir(sys.argv[2])

import pandas
from pandas import read_table

# Read the file with all Matrix Values (ANI enveomics)
dfAni = pandas.read_csv(aniFile, sep="\t", engine='python')
dfAAI = pandas.read_csv(aaiFile, sep="\t", engine='python')    

## First for ANI
# Count how many pairwise comparisons are below 85% (frx<85% over the shorter genome) for calculating vOTUs; or ANI below 80% which is its lower limit
totCtgs  = len(set(dfAni['SeqA']))
numFRX85 = len(dfAni[dfAni['Frx'] >= 85])
numvOTUs = len(dfAni[(dfAni['Frx'] >= 85) & (dfAni['ANI'] >= 95)])
numani80 = len(dfAni[dfAni['ANI'] >= 80])

# Calculate average nucleotide identity in the cluster
frxANI = dfAni['Frx'].mean()
avgANI = dfAni['ANI'].mean()
stdANI = dfAni['ANI'].std()
meanANI75 = dfAni[(dfAni.Frx >= 75)].mean()
stdvANI75 = dfAni[(dfAni.Frx >= 75)].std()
meanANI85 = dfAni[(dfAni.Frx >= 85)].mean()
stdvANI85 = dfAni[(dfAni.Frx >= 85)].std()

## AAI
# Get the total number of proteins for each contig
clusterID = aniFile.split("/")[8]
protPerc = []
totNProt = []
for index, row in dfAAI.iterrows():
    ctgName = (row['SeqA']).split(".")[0]
    comm = ("grep -c '>' ~/Documents/phages/spades/NEW/08_clusters/" + clusterID + "/" + ctgName + ".faa")
    totProt = (subprocess.Popen(comm, shell=True, stdout=subprocess.PIPE, universal_newlines=True).communicate()[0]).rstrip()
    totNProt.append(float(totProt))
    protPerc.append(float(row['N'])/float(totProt))
    ctgName = None
    comm = None

dfAAI['protPerc'] = protPerc
dfAAI['totNProt'] = totNProt
#print(dfAAI)

#Count how many genome pairwise comparisons are from the same genus (40%) or subfamily (20%)
aaiGEN = len(dfAAI[dfAAI['protPerc'] >= 0.4])
aaiFAM = len(dfAAI[dfAAI['protPerc'] >= 0.2])
aaiTOT = len(dfAAI)

# Calculate average aminoacid identity in the cluster
totAAI = dfAAI['totNProt'].mean()
frxAAI = dfAAI['protPerc'].mean()
avgAAI = dfAAI['AAI'].mean()
stdAAI = dfAAI['AAI'].std()

print("VCluster \t totCtg \t MeanANI \t StDevANI \t meanFrx \t avgANI75 \t stdANI75 \t avgFRX75 \t avgANI85 \t stdANI85 \t avgFRX85 \t aaiTOT \t aaiGEN \t aaiFAM \t frxAAI \t avgAAI \t stdAAI \t totAAI")
print("{0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7} \t {8} \t {9} \t {10} \t {11} \t {12} \t {13} \t {14} \t {15} \t {16} \t {17}".format(aniFile.split("/")[8],totCtgs,avgANI,stdANI,frxANI,meanANI75['ANI'],stdvANI75['ANI'],meanANI75['Frx'],meanANI85['ANI'],stdvANI85['ANI'],meanANI85['Frx'],aaiTOT,aaiGEN,aaiFAM,frxAAI,avgAAI,stdAAI,totAAI))


