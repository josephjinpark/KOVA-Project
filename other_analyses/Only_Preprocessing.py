#!/home/lab/bin/python2.7

import pysam
import sys, os
import re
import math
import random
import subprocess as sp
import numpy as np
from copy import deepcopy
from scipy.stats import binom, fisher_exact
from collections import Counter


lArgv = sys.argv
gsChrID = lArgv[1]
gsBamPathT, gsBamPathN = lArgv[2], lArgv[3]
gsOutDir = lArgv[4]
gsViewRange = lArgv[5] # e.g 1-1000 | 1001-2000
gsVCFFilePath = lArgv[6]
gbVCF, gbPLP = lArgv[7], lArgv[8]
gsOut_bam_path = lArgv[10]
gsPLPFilePath = gsVCFFilePath.replace('vcf', 'plp')
gfLODT, gfLODN, gfLODN4DBSNP = 1.0, 2.2, 5.5             # check this if contamination sample
gfLODTStrandBias = 2.0
gfAlphaT, gfAlphaN = float(lArgv[9]), 0 #default 0,0     #0.02, 0
giIndelWindowSize = 5
gfSensitivity = 0.9
gfClusteredPos1, gfClusteredPos2 = 10.0, 3.0

giMQ = 0
giBQS = -1 #Base Quality Score threshold
giQSMismatch = 100 # Sum of the Quality Scores of the Mismatches

##########new################
iShortage_blank_cutoff    = 0       # If The left window is over 0. That happen to make under N window.
iMax_window_range         = 12
fIndel_normal_cutoff      = 0.2
fRepeat_proportion_cutoff = 0.5
fTriallelic_pvalue_cutoff = 0.05
fStrandBias_pvalue_cutoff = 0.01

#### for filter valiable
gbDBSNP = True #1
gbProximalGap = True #2
gbStrandBias = True #3
gbClusteredPos = True#4
gbObservedInControl = True #5

gdResult = {'P':'PASS',
            'R':'REJECT',
            'p':'ProximalGap',
            'n':'ProximalGap_Normal',
            'r': 'Repeat',
            't':'Triallelic',
            's':'StrandBias',
            'c':'ClusteredPos',
            'd':'dbSNP',
            'o':'ObservedInControl',
            'x':'possible contamination'}
            # new condition : pn

samtools = 'samtools1.1'
gsRefBasePath = '/extdata5/Jinman/reference_genome/hg19/chromosomes/'
gsDBSNPPath = '/extdata5/Jinman/reference_genome/hg19/annotation/dbsnp/chromosomes/dbsnp_132.hg19.'+gsChrID
gsDBSNP = ''

giPQSType = giSanger = 33
gsCigarID = 'MIDNSHP=X'
giT, giN = 0, 1

gdSamFlag = {
    'p':0x0001, # the read is paired in sequencing
    'P':0x0002, # the read is mapped in a proper pair
    'u':0x0004, # the query sequence itself is unmapped
    'U':0x0008, # the mate is unmapped
    'r':0x0010, # strand of the query (1 for reverse)
    'R':0x0020, # strand of the mate
    '1':0x0040, # the read is the first read in a pair
    '2':0x0080, # the read is the second read in a pair
    's':0x0100, # the alignment is not primary
    'f':0x0200, # the read fails platform/vendor quality checks
    'd':0x0400  # the read is either a PCR or an optical duplicate
    }

class cPileup:
    def __init__(self):
        self.sTumReadSeq   = ''
        self.sTumQualScore = ''
        self.fTumLODT      = 0.0
        self.lTumContext   = []
        self.lTumDist5     = []
        self.lTumDist3     = []

        self.sNorReadSeq   = ''
        self.sNorQualScore = ''
        self.fNorLODN      = 0.0
        self.lNorContext   = []

        self.sReference                = '' # sRef
        self.sAltAllele                = ''
        self.sPassReject               = 'P'
        self.SensStrBias               = ''
        self.fMedianABSDev             = ''
        self.lStrandBiasCount          = [0,0,0,0]
        self.fStrandBias_PValue        = 0.0
        self.iRepeat                   = 0
        self.fHomoLODT                 = 0.0
        self.fHeteroLODT               = 0.0
        self.sGenotype                 = ''
    #end1: def __init__
#end: class cPileup


def fSetError(q): return math.pow(10, -1*(q/10.0)) # error rate of the position = 10^(-q/10)
def iDecodePQS(c): return ord(c) - giPQSType
def sEncodePQS(i): return chr(i + giPQSType)
def lMapCigar(c): return re.findall(r'(\d+)(\w)',c)

def fMedian(l):
    even = (0 if len(l) % 2 else 1) + 1
    half = (len(l) - 1) / 2
    return sum(sorted(l)[half:half + even]) / float(even)

def sSetFlagStrandBias(i, b): # int(0x10) = 16 (SEQ being reverse complemented)
    if i & 16: return (b.lower(), ',')[len(b)==0] # reverse
    else: return (b.upper(), '.')[len(b)==0] # forward

def sSetStrand(i): return '-' if i & 16 else '+'

def SamtoolsView(sCmd):
    sCmd = '%s view %s' %(gsSamtools, sCmd)
    sv = sp.Popen(sCmd, shell=True, stdout=sp.PIPE, bufsize=1)
    return sv.stdout

def SetSamByChr(sBamPath, iType): #../07_recal_bam/out/050005T_1T.recal.bam
    if iType==giT: # Tumor # PT1. XT != M
       # sCmd = '%s chr%s:%s | grep -v XT:A:M' %(sBamPath, gsChrID, gsViewRange) ###########################chr remove###########
        sCmd = '%s %s:%s | grep -v XT:A:M' %(sBamPath, gsChrID, gsViewRange)

    else: # Normal

    # sCmd = '%s chr%s:%s' %(sBamPath, gsChrID, gsViewRange)                  ###########################chr remvoe###########
        sCmd = '%s %s:%s' %(sBamPath, gsChrID, gsViewRange)

    return SamtoolsView(sCmd)

def tPreprocessingRead(iPos, sCigar, sSeq, sBQ,sRefBase):
    iPos -= 1
    sPrepSeq = ''
    sPrepBQ = ''
    lInsertion = []
    i = 0 # index for Insertion
    for mc in lMapCigar(sCigar): # (iCigarValue, sCigarKey)
        sCigarKey = mc[1]        # ex) M, M, M, ...
        iCigarValue = int(mc[0]) # ex) 101, 89, 101, 101, 95, ...

        if sCigarKey=='M': # M/0, I/1, D/2, N/3, S/4, H/5, P/6, =/7, X/8
            sPrepSeq += sSeq[:iCigarValue]
            sPrepBQ += sBQ[:iCigarValue]
            sSeq = sSeq[iCigarValue:]
            sBQ = sBQ[iCigarValue:]
            i += iCigarValue
        elif sCigarKey=='I': # skip
            lInsertion.append(i) #[i, sSeq[:iCigarValue]  #sSeq[:iCigarValue] -> insert sequence
            sSeq = sSeq[iCigarValue:]

            sBQ = sBQ[iCigarValue:]
        elif sCigarKey=='D':
            sPrepSeq += '*'*iCigarValue
            sPrepBQ += '!'*iCigarValue
            i += iCigarValue
        elif sCigarKey=='S': # skip
            sSeq = sSeq[iCigarValue:]
            sBQ = sBQ[iCigarValue:]


    sRefSeq = sRefBase[iPos:iPos+len(sPrepSeq)].upper()
    j = len(sRefSeq)
    return (sRefSeq, sPrepSeq[:j], sPrepBQ[:j], lInsertion)
#end tPreprocessingRead

def sPreprocessingBaseQualityScore(sPrepSeq, sPrepBQ): #### PT4 / PN1
    sSeq = ''
    for z in zip(sPrepSeq, sPrepBQ):
        if z[0] == '*':
            sSeq += '*'
        else:
            sSeq += ('E', z[0])[iDecodePQS(z[1]) >= giBQS] # E(Error)

    return sSeq

def bPreprocessingSoftClipped(sCigar, iPos, iType):
    if iType==giT: # Tumor
        iSoftClipped = iTotal = 0
        for t in lMapCigar(sCigar):
            if t[1]=='S': iSoftClipped += int(t[0])
            iTotal += int(t[0])
        return iSoftClipped / float(iTotal) * 100 >= 30
    else: return False

def bPreprocessingMismatches(sRefSeq, sPrepSeq, sPrepBQ, iType, iXM, sQName):
    
    iTotalBQ = 0
    #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': 
    #    print 'Ref', sRefSeq
    #    print 'Seq', sPrepSeq
    #    print 'bq', sPrepBQ

    if iType == giT and iXM > 2: # Tumor
        iTotalBQ = 0
        for z in zip(sRefSeq, sPrepSeq, sPrepBQ):

            if z[1] not in '*E,.' and z[0].upper() != z[1].upper(): #if z[0] != z[1]:
                
                iTotalBQ += iDecodePQS(z[2])
                
                #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print 'bqsum', z, iTotalBQ
                
                if iTotalBQ > giQSMismatch: break
            #end if
        #end for
    #end if


    return iTotalBQ > giQSMismatch # >100

def sPreprocessingRead2PLP(sRefSeq, sPrepSeq, iFlag):
    sPLPSeq = ''
    for z in zip(sRefSeq, sPrepSeq):
        if z[1] in '*NE': # *(Deletion), E(Error)
            sPLPSeq += z[1]
        else: # ACGT
            sPLPSeq += (sSetFlagStrandBias(iFlag, z[1]), sSetFlagStrandBias(iFlag, ''))[z[0]==z[1]]
    ### Example of sPLPSeq ###########################################################################################
    #  .........................................................C...........................................  (101)  #
    #  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  (101)  #
    ##################################################################################################################
    return sPLPSeq

def sInitOverlappingRead(sSeq): return sSeq.replace(',', '.').upper()

def lPreprocessingOverlappingRead(l1, l2, iType, iForward, iReverse, dPileup):

    iPos1S, sSeq1, sBQ1, sRef1, sStrand1 = l1 # S:Start, Seq = '.,ATGCatgc'
    iPos2S, sSeq2, sBQ2, sRef2, sStrand2 = l2

    i = 0
    j = 0

    if iPos1S > iPos2S:
        iPos1S, sSeq1, sBQ1, sRef1, sStrand1 = l2
        iPos2S, sSeq2, sBQ2, sRef2, sStrand2 = l1

    iSeq1Size = len(sSeq1)
    iPos1E    = iPos1S + iSeq1Size - 1 # E:End

    if iPos1E < iPos2S: # non overlap
        return False

    else: # overlap
        iSeq2Size = len(sSeq2)
        iPos2E = iPos2S + iSeq2Size - 1

        if iPos2E > iPos1E: # intersection
            i = iPos1E - iPos2S + 1
            sSeq1OL = sSeq1[-i:] # OL:Overlap
            sSeq2OL = sSeq2[:i]
            sBQ1OL = sBQ1[-i:]
            sBQ2OL = sBQ2[:i]

            sSeq1RD2 = '' # RD: Residue
            sBQ1RD2  = ''
            sSeq2RD1 = ''
            sBQ2RD1  = ''

            if sStrand1=='-':
                #sSeq1RD1 = 'E'*len(sSeq1[:-i])
                #sSeq2RD2 = 'E'*len(sSeq2[i:])

                sSeq1RD1 = sSeq1[:-i]
                sSeq2RD2 = sSeq2[i:]

            else:
                sSeq1RD1 = sSeq1[:-i]
                sBQ1RD1  = sBQ1[:-i]

                sSeq1RD1 = sPreprocessingBaseQualityScore(sSeq1RD1, sBQ1RD1)     # OL overlap, RD residue, each case BQ check.

                sSeq2RD2 = sSeq2[i:]
                sBQ2RD2  = sBQ2[i:]

                sSeq2RD2 = sPreprocessingBaseQualityScore(sSeq2RD2, sBQ2RD2)
        #end: if iPos2E

        else: # subset
            i = iPos2S - iPos1S
            j = i + iSeq2Size
            sSeq1OL = sSeq1[i:j]
            sSeq2OL = sSeq2
            sBQ1OL = sBQ1[i:j]
            sBQ2OL = sBQ2
            sSeq2RD1 = ''
            sBQ2RD1  = ''
            sSeq2RD2 = ''
            sBQ2RD2  = ''

            if sStrand1=='-':
                #sSeq1RD1 = 'E'*len(sSeq1[:i])
                sSeq1RD1 = sSeq1[:i]
                sBQ1RD1  = sBQ1[:i]

                sSeq1RD2 = sSeq1[j:]
                sBQ1RD2  = sBQ1[j:]

                sSeq1RD2 = sPreprocessingBaseQualityScore(sSeq1RD2, sBQ1RD2)

            else:
                sSeq1RD1 = sSeq1[:i]
                sBQ1RD1  = sBQ1[:i]

                sSeq1RD1 = sPreprocessingBaseQualityScore(sSeq1RD1, sBQ1RD1)

                sSeq1RD2 = sSeq1[j:]
                sBQ1RD1  = sBQ1[j:]
                #sSeq1RD2 = 'E'*len(sSeq1[j:])
        #end: else

        iPosC = iPos1S + len(sSeq1RD1) - 1 #iPosCurrent

        for z in zip(sSeq1OL, sSeq2OL, sBQ1OL, sBQ2OL):
            iPosC += 1

            if sInitOverlappingRead(z[0]) == sInitOverlappingRead(z[1]): # overlap and match

                sError = 'E'

                if z[0] in 'ACGTacgt':
                    sError = 'e'                                              # e is mismatch error, E is ref error

                if z[2]>=z[3]: s1, s2, = z[0], sError  # need to dig deeper

                #elif z[2]==z[3]: s1, s2 = z[0], 'E'
                elif z[2]<z[3]: s1, s2 = sError, z[1]

                if iDecodePQS(z[2]) < giBQS and iDecodePQS(z[3]) < giBQS:
                    s1, s2 = 'E', 'E'

                sSeq1RD1 += s1
                sSeq2RD1 += s2

            #end: if sInitOverlappingRead
            else: # overlap and mismatch

                """
                if z.count('*')==1:
                    if z[0]=='*':
                        sSeq1RD1 += 'E'
                        sSeq2RD1 += z[1]
                    else:
                        sSeq1RD1 += z[0]
                        sSeq2RD1 += 'E'
                    continue
                """
                if iType==giT: # Tumor

                    if z[0] in ',.':
                        sSeq1RD1 += 'E'
                    else: sSeq1RD1 += 'e'

                    if z[1] in ',.':
                        sSeq2RD1 += 'E'
                    else: sSeq2RD1 += 'e'

                else: # Normal
                    try:
                        sAlt = dPileup[iPosC].sAltAllele

                        if z[0].upper()==sAlt:

                            if iDecodePQS(z[2]) < giBQS:
                                sSeq1RD1 += 'E'
                                sSeq2RD1 += 'E'
                            else: sSeq1RD1 += z[0]; sSeq2RD1 += 'E'


                        elif z[1].upper()==sAlt:

                            if iDecodePQS(z[3]) < giBQS:
                                sSeq1RD1 += 'E'
                                sSeq2RD1 += 'E'
                            else: sSeq1RD1 += 'E'; sSeq2RD1 += z[1]

                        else:

                            if z[2]>=z[3]: s1, s2 = z[0], 'E' # > # need to dig deeper
                            elif z[2]<z[3]: s1, s2 = 'E', z[1]

                            if iDecodePQS(z[2]) < giBQS and iDecodePQS(z[3]) < giBQS:
                                s1, s2 = 'E', 'E'

                            sSeq1RD1 += s1 #z[0] #'E'
                            sSeq2RD1 += s2 #z[1] #'E'

                    #end: try
                    except:
                        sSeq1RD1 += 'E'
                        sSeq2RD1 += 'E'
                #end: else
        #end: for z

        lBQ1 = list(sBQ1)
        lBQ2 = list(sBQ2)
        lOverlap_BQ1 = []
        lOverlap_BQ2 = []

        if iType == 0:  #only tumor change the quality                      # if overlap, "forward, reverse base" change high score.
            if j == 0:  # intersection

                for sBQ1base, sBQ2base in zip(lBQ1[-i:], lBQ2[:i]):

                    if sBQ1base != sBQ2base:                                # Diffrent bases overlap are discarded both seq, so don't change higer quality.
                        lOverlap_BQ1 += sBQ1base
                        lOverlap_BQ2 += sBQ2base
                    #end: if sBQ1base
                    else:
                        if sBQ1base >= sBQ2base:                            # If same bases is change higher qaulity.
                            lOverlap_BQ1 += sBQ1base
                            lOverlap_BQ2 += sBQ1base
                        #end: if sBQ1base
                        else: lOverlap_BQ1 += sBQ2base; lOverlap_BQ2 += sBQ2base
                    #end: else
                #end for sSeq1base
                else: lBQ1[-i:] = lOverlap_BQ1; lBQ2[:i] = lOverlap_BQ2     # Chage the higher quality in overlap bases.

            else: # subset
                for sBQ1base, sBQ2base in zip(lBQ1[i:j], lBQ2):
                    if sBQ1base != sBQ2base:
                        lOverlap_BQ1 += sBQ1base
                        lOverlap_BQ2 += sBQ2base
                    #end: if sBQ1base
                    else:
                        if sBQ1base >= sBQ2base:
                            lOverlap_BQ1 += sBQ1base
                            lOverlap_BQ2 += sBQ1base
                        #end: if sBQ1base
                        else: lOverlap_BQ1 += sBQ2base; lOverlap_BQ2 += sBQ2base
                    #end: else
                #end: for sSeq1base
                else: lBQ1[i:j] = lOverlap_BQ1; lBQ2 = lOverlap_BQ2
            #end: j
        #end: itype

        sBQ1 = ''.join(lBQ1)
        sBQ2 = ''.join(lBQ2)

        return [sSeq1RD1 + sSeq1RD2, sSeq2RD1 + sSeq2RD2, sBQ1, sBQ2]
#end lPreprocessingOverlappingRead

def tRemoveErrorAlt(sXT, sQT, lD5, lD3, sAlt): # Remove except '.,'+sAlt
    x = q = ''
    l5 = []
    l3 = []
    for z in zip(sXT, sQT, lD5, lD3):
        if z[0].upper() in '.,'+sAlt:
            x += z[0]
            q += z[1]
            l5.append(z[2])
            l3.append(z[3])

    return(x, q, l5, l3)

def tRemoveErrorT(sX, sQ):
    x = q = ''
    for i, b in enumerate(sX):
        if b not in '*nNE':
            x += b
            q += sQ[i]
    return(x, q)

def tRemoveErrorN(sX, sQ, sAlt):
    x = q = ''
    for i, b in enumerate(sX):
        if b in '.,':
            x += b
            q += sQ[i]
        elif b.upper()==sAlt:
            x += b
            q += sQ[i]
    return(x, q)

def bCountBaseType(sX):
    for x in sX.upper():
        if x in 'ATGC': return True

    return False

def iCountBaseType(sX):
    sAltSeq = ''
    for x in sX.upper():
        if x in 'ATGC': sAltSeq += x

    return len(set(sAltSeq))


def SoVar_Processing():
    """
    PT1. XT != M
    PT2. MQS > 0
    PT3. BQS >= 5
    PT4. Soft clipped < 30%
    PT5. Mismatches <= 100
    PT6. Overlapping
    PN1. BQS >= 5
    PN2. Overlapping
    """
    #lFlag_check = [65,73,89,113,129,137,153,177] # set 1 reject
    #lFlag_check = [65,73,89,113,129,137,153,177,81,145] # set 2 reject
    #lFlag_check = [83,99,147,163] #set 3 proper_read
    #if not iFlag in lFlag_check:
    #    continue

    fp = open(gsRefBasePath+gsChrID+'.fa.base')
    sRefBase = fp.read()
    fp.close()
    dPileup = {}

    samfile_input = pysam.Samfile(gsBamPathT, 'rb')
    print gsOut_bam_path
    samfile_make  = pysam.Samfile(gsOut_bam_path, 'wb', template=samfile_input)
    sRange = gsViewRange.split('-')
    iStart = int(sRange[0])
    iEnd   = int(sRange[1])-101

    #i = int(lArgv[11]) ########################################################## 0 tumor, 1 normal ############################################
    dPrepRead = {}
    i = 0
    #for line in SetSamByChr(sBamPath, i): #### PT1

    for line in samfile_input.fetch(gsChrID, iStart-1, iEnd):
        sQName = line.qname
        iFlag = int(line.flag)

        iPos = int(line.pos+1) # 1-based leftmost mapping Position

        iMQ = int(line.mapq)
        
        if i == 0 and iMQ <= giMQ: continue #### PT2
        
        #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '1', iPos

        sCigar  = line.cigarstring
        sSeq    = line.query
        sBQ     = line.qual

        sXT_tag = ''

        for tag in line.tags:
            if tag[0] == 'XT' and tag[1] == 'M':
                continue
        
        #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '2'

        sRefSeq, sPrepSeq, sPrepBQ, lInsertion = tPreprocessingRead(iPos, sCigar, sSeq, sBQ,sRefBase)

        if bPreprocessingSoftClipped(sCigar, iPos, i): continue #### PT4
        
        iXM = 0

        for tag in line.tags:
            if tag[0] == 'XM': # XM = Mismatch
                iXM = int(tag[1])
        
        #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '3'

        sPrepSeq = sPreprocessingRead2PLP(sRefSeq, sPrepSeq, iFlag)   # make the seq to '...,A..,..'
        sStrand = sSetStrand(iFlag) # +|-

        if dPrepRead.has_key(sQName): # sQName horizontal direction read
            dPrepRead[sQName].append([iPos, sPrepSeq, sPrepBQ, sRefSeq, sStrand, lInsertion, iXM, i, iFlag, sXT_tag, line])

        else: dPrepRead[sQName] = [[iPos, sPrepSeq, sPrepBQ, sRefSeq, sStrand, lInsertion, iXM, i, iFlag, sXT_tag, line]]
        #                          [0]    [1]       [2]      [3]      [4]      [5]         [6] [7] [8]    [9]
        ### Examples 'dPrepRead'
        # {'D0ENMACXX111207:2:2201:2129:87437': [[9413215, '.....................................................................................................', '?@=BDFCFEGEFAFFEFEEFJGEHFFFFEFGGGFIGFFFFEFFHGGFIGFFFGHAICDG@FCDGAFFFFFFDFHGEFFGGEFFGEFFIDFHFGGFCE>=AB', 'CAATAGGCAGTATCAATTTAGGTCTATTTTCCATGAATATTTTCTCAGCAACTGTGGTGTTATGATATATATTGGTTTTCATCCACAGTTCCTGGCTTATA', '+', [], 0, 0, 99, 'R'], [9413479, ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,', 'C@>@DFEFCFFEEFGFFFGGDFGCCDGFFGGFFEEEFFFDHGFIHFCEGGFFFFF<GFEGGGIFFIEFFGFGGFFAEEHEGFEFFFEFEDFGDFCCCA??=', 'GAGGGAGGAATGCTGATATCATGAAACTTCCATAAAAATCCAGGAGGACAGGGTTCAGTGAGCTTCTGGGTAGTTGAACACATGGATGTTCCTGTAGGGTG', '-', [], 0, 0, 147, 'UorM']],
        #  sQName,                              iPosition, sPrepSeq,                                                                                              sPrepBQ,                                                                                                 sRefSeq,                                                                                               sStrand, lInsertion, iXM, i, iFlag, sXT_tag
        ###
    #end: for line
    #print 'pre', dPrepRead
    dProcessed_Read = {}
    for sQName, llPrepRead in dPrepRead.items():

        iForward_pos    = llPrepRead[0][0]
        sForward_seq    = llPrepRead[0][1]
        sForward_BQ     = llPrepRead[0][2]
        sForward_refseq = llPrepRead[0][3]
        sForward_Strand = llPrepRead[0][4]
        lForward_Insert = llPrepRead[0][5]
        sForward_XM     = llPrepRead[0][6]
        iForward_Type   = llPrepRead[0][7]
        sForward_XT_tag = llPrepRead[0][9]
        sForward_Samformat = llPrepRead[0][10]

        if len(llPrepRead) == 1:
            sForward_seq = sPreprocessingBaseQualityScore(sForward_seq, sForward_BQ) #### PT3 / PN1
            
            #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '1_4'
            if not bPreprocessingMismatches(sForward_refseq, sForward_seq, sForward_BQ, iForward_Type, sForward_XM, sQName):
                #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '1_5'
                dProcessed_Read[sQName] = [[iForward_pos, sForward_seq, sForward_BQ, sForward_refseq, sForward_Strand, lForward_Insert, sForward_XT_tag, sForward_Samformat]]

        #end: if len(llPrepRead)

        else:

            iReverse_pos    = llPrepRead[1][0]
            sReverse_seq    = llPrepRead[1][1]
            sReverse_BQ     = llPrepRead[1][2]
            sReverse_refSeq = llPrepRead[1][3]
            sReverse_Strand = llPrepRead[1][4]
            lReverse_Insert = llPrepRead[1][5]
            iReverse_XM     = llPrepRead[1][6]
            iReverse_Type   = llPrepRead[1][7]
            sReverse_XT_tag = llPrepRead[1][9]
            sReverse_Samformat = llPrepRead[1][10]

            rval = False # return value

            l1 = [iForward_pos, sForward_seq, sForward_BQ, sForward_refseq, sForward_Strand]
            l2 = [iReverse_pos, sReverse_seq, sReverse_BQ, sReverse_refSeq, sReverse_Strand]
            rval = lPreprocessingOverlappingRead(l1, l2, iForward_Type, iForward_pos, iReverse_pos, dPileup) #### PT5 / PN2  # i == iType

            #print rval
            
            if rval: # overlap

                if not bPreprocessingMismatches(sForward_refseq, rval[0], rval[2], iForward_Type, sForward_XM, sQName):
                    dProcessed_Read[sQName] = [[iForward_pos, rval[0], rval[2], sForward_refseq, sForward_Strand, lForward_Insert, sForward_XT_tag, sForward_Samformat]]
                #end: if not bPreprocessingMismatches

                if not bPreprocessingMismatches(sReverse_refSeq, rval[1], rval[3], iReverse_Type, iReverse_XM, sQName):

                    if dProcessed_Read.has_key(sQName):
                        dProcessed_Read[sQName].append([iReverse_pos, rval[1], rval[3], sReverse_refSeq, sReverse_Strand, lReverse_Insert, sReverse_XT_tag, sReverse_Samformat])
                    else: dProcessed_Read[sQName] = [[iReverse_pos, rval[1], rval[3], sReverse_refSeq, sReverse_Strand, lReverse_Insert, sReverse_XT_tag, sReverse_Samformat]]
                #end: if not bPreprocessingMismatches
            #end: if rval

            else: # non overlap

                sForward_seq = sPreprocessingBaseQualityScore(sForward_seq, sForward_BQ)
                sReverse_seq = sPreprocessingBaseQualityScore(sReverse_seq, sReverse_BQ)
                
                #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '4'
                if not bPreprocessingMismatches(sForward_refseq, sForward_seq, sForward_BQ, iForward_Type, sForward_XM, sQName):
                    #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '5'
                    dProcessed_Read[sQName] = [[iForward_pos, sForward_seq, sForward_BQ, sForward_refseq, sForward_Strand, lForward_Insert, sForward_XT_tag, sForward_Samformat]]

                if not bPreprocessingMismatches(sReverse_refSeq, sReverse_seq, sReverse_BQ, iReverse_Type, iReverse_XM, sQName):
                    #if sQName == 'D0EN0ACXX111207:4:1203:1634:151662': print '6'
                    if dProcessed_Read.has_key(sQName):
                        dProcessed_Read[sQName].append([iReverse_pos, sReverse_seq, sReverse_BQ, sReverse_refSeq, sReverse_Strand, lReverse_Insert, sReverse_XT_tag, sReverse_Samformat])

                    else:
                        dProcessed_Read[sQName] = [[iReverse_pos, sReverse_seq, sReverse_BQ, sReverse_refSeq, sReverse_Strand, lReverse_Insert, sReverse_XT_tag, sReverse_Samformat]]
                #end: if not bPreprocessingMismatches
            
            #end: else
        #end: else
    #end: for sQName, llPrepRead

    for sExtract_Samformat in dProcessed_Read.values():
        #print 'ex', sExtract_Samformat
        
        Samformat_forward = sExtract_Samformat[0][-1]       
        samfile_make.write(Samformat_forward)
        #print 'final', Samformat_forward
        
        try:
            Samformat_reverse = sExtract_Samformat[1][-1]
            samfile_make.write(Samformat_reverse)
            #print Samformat_reverse
        
        except IndexError:
            pass 

    samfile_input.close()
    samfile_make.close()


def Main():
    SoVar_Processing()
#end: def Main

import timeit
print 'start'
start = timeit.default_timer()

Main()

print 'Preprocessing_done'
print 'end'
stop = timeit.default_timer()
print 'Runtime: %.3f sec' %(stop - start)
