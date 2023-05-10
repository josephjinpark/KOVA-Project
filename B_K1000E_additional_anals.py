#!/extdata6/Jinman/opt/python3/bin/python3

bDEBUGMODE = 0

import os, sys, pickle, time, re, copy, subprocess, uuid, math, struct, scipy, itertools
from scipy import stats
from operator import itemgetter
import numpy as np


## region Globals

sAVINPUT          = '/extdata6/Jinman/util/annovar/convert2annovar.pl'
sANNOVAR          = '/extdata6/Jinman/util/annovar/table_annovar.pl'
sBASE_DIR         = '/extdata6/Kyuhong/05_K1000E_vcf_comparison/04_ANNOVAR'
sREF_DIR          = '/extdata6/Kyuhong/05_K1000E_vcf_comparison/04_ANNOVAR/ref'
fPSEUDO_COUNT     = float('1.0e-300')
sTIME_STAMP       = '%s' % (time.ctime().replace(' ', '-').replace(':', '_') )


#V-S Check Thresholds
nCOSMIC_FILE_COLUMN_SIZE    = 34
nGWAS_FILE_COLUMN_SIZE      = 34

## endregion


## region Common Functions

def copy_temp_core_script (sWorkDir):

    os.makedirs('%s/temp' % sWorkDir, exist_ok=True)
    os.system('cp %s/B_K1000E_additional_anals.py %s/temp/tmp_script_%s.py'
              % (sWorkDir, sWorkDir, sTIME_STAMP))

    return '%s/temp/tmp_script_%s.py' % (sWorkDir, sTIME_STAMP)


def mean (list_fValues):
    return sum(list_fValues)/len(list_fValues)


def median(values):
    '''median'''
    values = sorted(values)
    n      = len(values)
    q, r = divmod(n, 2)
    if r == 1: ret =  values[q]
    else:      ret = (values[q-1] + values[q]) / 2
    return ret
#def END: median

def stdev(list_fValues, fOption):
    if len(list_fValues) < 2:
        return None

    fSD   = 0.0
    fSum  = 0.0
    fMean = mean(list_fValues)

    for i in range(len(list_fValues)):
        fDiff = list_fValues[i] - fMean
        fSum += fDiff * fDiff

    fSD = math.sqrt(fSum / (len(list_fValues) - fOption))
    return fSD
#def END: stdev


def stderror(list_fValues):
    try:
        return stdev(list_fValues, 0.0)/math.sqrt(len(list_fValues))*1.96
    except TypeError:
        return 'NULL'
#def END: stderror


def print_start_msg (sTaskName, sFileDir=''):
    print('START...%s %s %s' % (sTaskName, sFileDir, time.ctime()))
#def END: print_start_msg


def print_done_msg (sTaskName, sFileDir=''):
    print('DONE...%s %s %s' % (sTaskName, sFileDir, time.ctime()))
#def END: print_done_msg


def annotate_with_database (list_cDatabase, list_cAnno, sDatabase, sDisease):

    print_start_msg('Annotating with %s database-%s' % (sDatabase, sDisease))

    if sDatabase == 'COSMIC':
        # By Position
        dict_cData       = {int(cCos.sPos) : cCos  for cCos in list_cDatabase if '-' not in cCos.sPos }
        dict_cAnno       = {cAnno.nStartPos: cAnno for cAnno in list_cAnno}

        list_sKeys_data  = list(dict_cData.keys())
        list_sKeys_anno  = list(dict_cAnno.keys())

    elif sDatabase == 'GWAS':
        # By Position
        dict_cData       = {cGwas.sSNP_ID  : cGwas for cGwas in list_cDatabase if cGwas.sSNP_ID != 'NA'}
        dict_cAnno       = {cAnno.sDBSNP_ID: cAnno for cAnno in list_cAnno if cAnno.sDBSNP_ID != 'NA'}

        list_sKeys_data  = list(dict_cData.keys())
        list_sKeys_anno  = list(dict_cAnno.keys())

    elif sDatabase == 'ExAC':
        # By Position
        dict_cData       = {cExAC.nPos : cExAC  for cExAC in list_cDatabase}
        dict_cAnno       = {cAnno.nStartPos: cAnno for cAnno in list_cAnno}

        list_sKeys_data  = list(dict_cData.keys())
        list_sKeys_anno  = list(dict_cAnno.keys())

    elif sDatabase == 'TCGA':
        # By Position
        dict_cData = {cVCF.nPos: cVCF for cVCF in list_cDatabase}
        dict_cAnno = {cAnno.nStartPos: cAnno for cAnno in list_cAnno}

        list_sKeys_data = list(dict_cData.keys())
        list_sKeys_anno = list(dict_cAnno.keys())
    #if END:

    print('%s-PosCount' % sDatabase, len(list_sKeys_data))
    print('%s-PosCount' % sDisease,  len(list_sKeys_anno))

    list_sIntersection = list(set(list_sKeys_data) & set(list_sKeys_anno))
    print('IntersectionCount', len(list_sIntersection))

    list_cAnno         = [dict_cAnno[sKey] for sKey in list_sIntersection if dict_cAnno[sKey].sGeneFunc != 'intergenic']

    #for cAnno in list_cAnno:
        #print(cAnno.sNMID, cAnno.sGeneSym, cAnno.nStartPos, dict_cData[cAnno.nStartPos].sDelete, dict_cData[cAnno.nStartPos].sMutaStatus)

    print('Final cAnno Cnt', len(list_cAnno))
    return list_cAnno
#def END: annotate_with_database


def output_for_maf_distribution (sOutDir, list_fMAF, sDisease, sSource):

    ## MAF Histogram ##
    dict_fMAF = {}
    fBinRange = max(list_fMAF) / 50 # 50 Bins
    for nGroup, lCnt in itertools.groupby(list_fMAF, key=lambda n: n//fBinRange):

        if nGroup not in dict_fMAF:
            dict_fMAF[nGroup] = []

        dict_fMAF[nGroup] = dict_fMAF[nGroup] + list(lCnt)
    #loop END: nGroup, lCnt

    nSumCheck = 0
    OutFile     = open('%s/MAF_Dist_%s_%s_%s.txt' % (sOutDir, sDisease, sSource, fBinRange), 'w')
    sHeader     = 'Bin-%s\t%s\t%s\t%s\n' % (fBinRange, 'Count', 'Mean', 'SEM')
    OutFile.write(sHeader)

    for nGroup in dict_fMAF:

        nSumCheck += len(dict_fMAF[nGroup])
        sOut       = '%s\t%s\t%0.5f\t%0.5f\n' % \
                     (nGroup, len(dict_fMAF[nGroup]), mean(dict_fMAF[nGroup]),
                      stderror(dict_fMAF[nGroup]) if stderror(dict_fMAF[nGroup])  != 'NULL' else 0)
        OutFile.write(sOut)
    #loop END: nGroup

    ## Subset of MAF distribution in first bin ##
    fBinRange       = float('%0.4f' % (max(dict_fMAF[0]) / 20)) # 20 Bins
    dict_fMAF_sub   = {}
    sHeader         = 'Bin-%s\t%s\t%s\t%s\n' % (fBinRange, 'Count', 'Mean', 'SEM')
    OutFile.write(sHeader)
    for nGroup, lCnt in itertools.groupby(dict_fMAF[0], key=lambda n: n//fBinRange):

        if nGroup not in dict_fMAF_sub:
            dict_fMAF_sub[nGroup] = []

        dict_fMAF_sub[nGroup] = dict_fMAF_sub[nGroup] + list(lCnt)
    #loop END: nGroup, lCnt

    for nGroup in dict_fMAF_sub:

        sOut       = '%s\t%s\t%0.5f\t%0.5f\n' % \
                     (nGroup, len(dict_fMAF_sub[nGroup]), mean(dict_fMAF_sub[nGroup]),
                      stderror(dict_fMAF_sub[nGroup]) if stderror(dict_fMAF_sub[nGroup])  != 'NULL' else 0)
        OutFile.write(sOut)
    #loop END: nGroup
    OutFile.close()

    #V-S Check
    if len(list_fMAF) != nSumCheck:
        sys.exit('Invalid Dictionary : dict_fMAF Size= %d : SumCheck= %d' % (len(list_fMAF), nSumCheck))
#def END: output_for_maf_distribution


def output_for_maf_diff (sOutDir, list_cAnno, sDisease):

    nCntCutoff     = 2 # Mutations per Gene

    ## Calculate MAF Difference
    list_fMAF_fc   = []
    list_fMAF_diff = []
    dict_sGene     = {}
    list_sOutput   = []
    for cAnno in list_cAnno:
        fDiff = cAnno.fMAF_PNorm - cAnno.fMAF_HNorm
        fFC   = math.log(cAnno.fMAF_PNorm/(cAnno.fMAF_HNorm + fPSEUDO_COUNT))

        if cAnno.sGeneSym not in dict_sGene:
            dict_sGene[cAnno.sGeneSym] = 0
        dict_sGene[cAnno.sGeneSym] += 1

        list_fMAF_diff.append(fDiff)
        list_fMAF_fc.append(fFC)
    #loop END: cAnno

    for sGene in dict_sGene:
        list_sOutput.append([sGene, dict_sGene[sGene]])

    list_sOutput = sorted(list_sOutput, key=lambda e:e[1], reverse=True)
    nGeneCnt     = 0
    OutFile      = open('%s/%s_MAF_gene_varcnt%s.txt' % (sOutDir, sDisease, nCntCutoff), 'w')
    sHeader      = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                   % ('VarCnt','GeneSym','MeanMAF_PN','MeanMAF_HN','SEM_PN','SEM_HN','Pvalue','Significant')
    OutFile.write(sHeader)
    for sGene, nMutaCnt in list_sOutput:
        list_fMAF_p = [cAnno.fMAF_PNorm for cAnno in list_cAnno if cAnno.sGeneSym == sGene]
        list_fMAF_h = [cAnno.fMAF_HNorm for cAnno in list_cAnno if cAnno.sGeneSym == sGene]
        fAvgMAF_P   = mean(list_fMAF_p)
        fAvgMAF_H   = mean(list_fMAF_h)
        fSEM_P      = stderror(list_fMAF_p) if stderror(list_fMAF_p) != 'NULL' else 0
        fSEM_H      = stderror(list_fMAF_h) if stderror(list_fMAF_h) != 'NULL' else 0
        Z, fPvalue  = stats.ranksums(list_fMAF_p, list_fMAF_h)
        sSig        = '*' if fPvalue < 0.05 else ''

        sOut = '%s\t%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%s\n' \
               % (nMutaCnt, sGene, fAvgMAF_P, fAvgMAF_H, fSEM_P, fSEM_H, fPvalue, sSig)

        ## Pass ones where healthy MAF is greater
        #if fAvgMAF_H > fAvgMAF_P: continue

        if nMutaCnt >= nCntCutoff:
            #print(sOut[:-1])
            nGeneCnt += 1
            OutFile.write(sOut)
    #loop END: sGene, sMutaCnt
    OutFile.close()

    print('TotalGeneCnt', len(dict_sGene))
    print('FinalGeneCnt', nGeneCnt)
#def END: output_for_maf_fc


def get_maf_from_freq_file (sInFile):

    dict_sOutput = {}
    print('Loading %s' % sInFile)
    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:
        '''
        CHROM   POS     N_ALLELES       N_CHR   {ALLELE:FREQ}
        chr1    14717   2       388     G:0.989691      A:0.0103093
        '''
        if sReadLine.startswith('CHROM'): continue
        list_sColumn = sReadLine.strip('\n').split('\t')

        sChrID       = list_sColumn[0]
        nPos         = int(list_sColumn[1])
        nAlleleCnt   = int(list_sColumn[2])
        nChrCnt      = int(list_sColumn[3])
        sMajorAllele = list_sColumn[4].split(':')[0]
        if sMajorAllele not in ['A','C','G','T']:continue

        fMajorAF     = float(list_sColumn[4].split(':')[1])
        sMinorAllele = list_sColumn[5].split(':')[0]

        if sMinorAllele not in ['A','C','G','T']:continue

        fMinorAF     = float(list_sColumn[5].split(':')[1])

        if nPos not in dict_sOutput:
            dict_sOutput[nPos] = ''
        dict_sOutput[nPos] = fMinorAF
    #loop END: sReadLine
    InFile.close()

    #V-S Check:
    if not dict_sOutput:
        sys.exit('Invalid Dictionary : get_maf_from_freq_file : dict_sOutput size= %d' % (len(dict_sOutput)))

    return dict_sOutput
#def END: dict_fMAF_tcga


def tcga_1000gp_validation (sTag, sCancer, sWorkDir, list_cAnno, bOutput):
    dict_sCancer    = {'lung': 'LUAD', 'gastric': 'STAD'}  # Distinguish between K1000E and TCGA
    sTCGA_VCF       = '%s/%s.vcf'        % (sREF_DIR, dict_sCancer[sCancer])
    sTCGA_ID_File   = '%s/%s_IDlist.txt' % (sREF_DIR, dict_sCancer[sCancer])


    ## Compare with TCGA and 1000GP
    list_sFile_TCGA  = get_sample_list        (sTCGA_ID_File)
    list_cVCF_tcga   = cVCF_parse_vcf_files   (len(list_sFile_TCGA), sTCGA_VCF)
    list_cAnno_tcga  = annotate_with_database (list_cVCF_tcga, list_cAnno, 'TCGA', sCancer)

    print(red('TCGA Match :\t%d' % len(list_cAnno_tcga), 1))

    ## Get MAF from TCGA and 1000GP
    #list_sOutput   =  calculate_ranksum (sWorkDir, sTag, sCancer, list_cAnno_tcga, dict_sCancer, 'EUR', bOutput)
    #list_sOutput   =  calculate_ranksum (sWorkDir, sTag, sCancer, list_cAnno_tcga, dict_sCancer, 'AFR', bOutput)
    list_sOutput   =  calculate_ranksum (sWorkDir, sTag, sCancer, list_cAnno_tcga, dict_sCancer, 'EAS', bOutput)

    dict_cData = {cVCF.nPos: cVCF for cVCF in list_cVCF_tcga}

    sOutFile   = '%s/vcf_info_%s.txt' % (sWorkDir, sCancer)
    OutFile    = open(sOutFile, 'w')
    for cAnno in list_sOutput:
        cVCF = dict_cData[cAnno.nStartPos]

        sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (cVCF.sChrID, cVCF.nPos, cVCF.sDBSNP_ID, cVCF.sRefNuc, cVCF.sAltNuc,
                  cVCF.fQual, cVCF.sFilter, cVCF.sInfo, cVCF.sFormat, cVCF.list_sSamples,
                  cAnno.fMAF_PNorm, cAnno.fMAF_HNorm, cAnno.fMAF_tcga, cAnno.fMAF_1000)
        OutFile.write(sOut)
    #loop END: cAnno
    OutFile.close()

    return list_sOutput
#def END: tcga_1000gp_validation


def calculate_ranksum (sWorkDir, sTag, sCancer, list_cAnno, dict_sCancer, sRace, bOutput):
    sFreqDir       = '/extdata6/Kyuhong/05_K1000E_vcf_comparison/03_vcftools'
    sFreqFile_tcga = '%s/%s/%s.final.filtered.frq'    % (sFreqDir, dict_sCancer[sCancer], dict_sCancer[sCancer])
    sFreqFile_1000 = '%s/%s/%s_%s.final.filtered.frq'   % (sFreqDir, '1000GP', '1000GP', sRace)
    dict_fMAF_tcga = get_maf_from_freq_file (sFreqFile_tcga)
    dict_fMAF_1000 = get_maf_from_freq_file (sFreqFile_1000)

    print(len(dict_fMAF_tcga))
    print(len(dict_fMAF_1000))

    list_sOutput = []
    list_fMAF_p  = []
    list_fMAF_h  = []

    for cAnno in list_cAnno:
        try:
            fMAF_tcga   = dict_fMAF_tcga[cAnno.nStartPos]
            fMAF_1000   = dict_fMAF_1000[cAnno.nStartPos]
        except KeyError: continue
        cAnno.fMAF_tcga = fMAF_tcga
        cAnno.fMAF_1000 = fMAF_1000
        list_fMAF_p.append(fMAF_tcga)
        list_fMAF_h.append(fMAF_1000)
        list_sOutput.append(cAnno)
    #loop END:

    Z, fPvalue = stats.ranksums(list_fMAF_p, list_fMAF_h)
    print(red('Final Output Cnt :\t%d'    % len(list_sOutput), 1))
    print(red('Pvalue           :\t%0.2E' % fPvalue, 1))
    print(red('Median-tcga      :\t%s'    % median(list_fMAF_p), 1))
    print(red('Median-1000      :\t%s'    % median(list_fMAF_h), 1))


    if bOutput:
        OutFile = open('%s/%s_wTCGA-1000_%s_%s.txt' % (sWorkDir, sTag, sCancer, sRace), 'w')
        for cAnno in list_sOutput:
            sOut = '%s\t%s\t%s\t%s\t%s\t%s\n' \
                   % (cAnno.sNMID, cAnno.sGeneSym, cAnno.fMAF_PNorm,
                      cAnno.fMAF_HNorm, cAnno.fMAF_tcga, cAnno.fMAF_1000)
            OutFile.write(sOut)
        OutFile.close()


    return list_sOutput
#def END: calculate_ranksum


def parse_pathways (sReturn='Pathway'):

    dict_sPathways  = {}
    dict_sGeneSym   = {}
    sInFile         = '%s/CPDB_pathways.txt' % sREF_DIR
    InFile          = open(sInFile, 'r')

    for sReadLine in InFile:
        if sReadLine.startswith('pathway'): continue
        list_sColumn = sReadLine.strip('\n').split('\t')

        sPathName       = list_sColumn[0]
        sExternalID     = list_sColumn[1]
        sSource         = list_sColumn[2]
        list_sGeneSym   = [sGene.upper() for sGene in  list_sColumn[3].split(',')]

        if sPathName not in dict_sPathways:
            dict_sPathways[sPathName] = []
        dict_sPathways[sPathName] = list_sGeneSym

        for sGeneSym in list_sGeneSym:
            sKey = sGeneSym.upper()

            if sKey not in dict_sGeneSym:
                dict_sGeneSym[sKey] = []
            dict_sGeneSym[sKey].append(sPathName)
    #loop END: sReadLine
    InFile.close()

    print('Pathways', len(dict_sPathways))
    print('Genes', len(dict_sGeneSym))


    #V-S Check:
    if not dict_sPathways:
        sys.exit('Invalid Dictionary : parse_pathways : dict_sPathways size= %d' % len(dict_sPathways))

    return dict_sPathways if sReturn == 'Pathway' else dict_sGeneSym
#def END: parse_pathways


def pathway_survey (list_cAnno):

    dict_sPathways  = parse_pathways()
    dict_cAnno      = {}
    dict_sOutput    = {}
    list_sCancerSym = dict_sPathways['Pathways in cancer - Homo sapiens (human)']

    for cAnno in list_cAnno:

        if cAnno.sGeneSym not in dict_cAnno:
            dict_cAnno[cAnno.sGeneSym] = []
        dict_cAnno[cAnno.sGeneSym].append(cAnno)
    #loop END: cAnno
    print('dict_cAnno', len(dict_cAnno))

    for sPathway in dict_sPathways:

        if sPathway not in dict_sOutput:
            dict_sOutput[sPathway] = []

        for sGeneSym in dict_sPathways[sPathway]:

            if sGeneSym not in dict_cAnno: continue

            dict_sOutput[sPathway] += dict_cAnno[sGeneSym]
    #loop END: sPathway

    list_sPaths = []
    for sPathway in dict_sOutput:
        if not dict_sOutput[sPathway]: continue
        list_sPaths.append([sPathway, dict_sOutput[sPathway]])
    #loop END: sPathway

    list_sPaths  = sorted(list_sPaths, key=lambda e:len(e[1]), reverse=True)

    list_sOutput = []

    for sPathway, list_cAnno_filt in list_sPaths:
        list_fMAF_p     = []
        list_fMAF_h     = []
        list_fMAF_tcga  = []
        list_fMAF_1000  = []

        bCancer     = False

        for cAnno in list_cAnno_filt:
            list_fMAF_p.append(cAnno.fMAF_PNorm)
            list_fMAF_h.append(cAnno.fMAF_HNorm)
            list_fMAF_tcga.append(cAnno.fMAF_tcga)
            list_fMAF_1000.append(cAnno.fMAF_1000)

        # loop END: cAnno

        for cAnno in list_cAnno_filt:
            if cAnno.sGeneSym in list_sCancerSym:
                bCancer = True
                break
        #loop END: cAnno

        fMedianP     = median(list_fMAF_p)
        fMedianH     = median(list_fMAF_h)

        fMedianT     = median(list_fMAF_tcga)
        fMedian1     = median(list_fMAF_1000)
        #fMedianT     = 1
        #fMedian1     = 1


        Z, fPvalue   = stats.ranksums(list_fMAF_p, list_fMAF_h)
        Z2, fPvalue2 = stats.ranksums(list_fMAF_tcga, list_fMAF_1000)
        #fPvalue2 = 1
        list_sOutput.append([sPathway, list_cAnno_filt, fPvalue, fPvalue2,
                             fMedianP, fMedianH, fMedianT, fMedian1, bCancer])
    #loop END: sPathway, list_cAnno_filt

    print('list_sOutput', len(list_sOutput))

    #list_sOutput  = sorted(list_sOutput, key=lambda e:(e[3] - e[4], e[2], len(e[1])), reverse=True)
    #list_sOutput  = sorted(list_sOutput, key=lambda e:len(e[1]), reverse=True)
    list_sOutput  = sorted(list_sOutput, key=lambda e:(e[3]))
    #list_sOutput  = sorted(list_sOutput, key=lambda e:e[3]-e[4], reverse=True)

    nCnt = 0
    for sPathway, list_cAnno, fPval, fPval2, fMedianP, fMedianH, fMedianT, fMedian1, bCancer in list_sOutput:

        if fMedianH > fMedianP:     continue
        if not bCancer:             continue
        #if fPval2 > 0.05:           continue
        #if len(list_cAnno) > 10000: continue
        nCnt += 1

        sOut = '%s\t%s\t%0.3E\t%0.3E\t%0.4f\t%0.4f\t%0.4f\t%0.4f' % \
               (sPathway, len(list_cAnno), fPval, fPval2, fMedianP, fMedianH, fMedianT, fMedian1)
        print(sOut)
    #loop END: sPathway, list_cAnno, fPval, fMedianP, fMedianH\

    print(green('Final Pathway Count %s' % nCnt,1))
#def END: pathway_survey


def pathway_survey_ts (list_sFinal):

    dict_sPathways  = parse_pathways()
    dict_cAnno      = {}
    dict_sOutput    = {}
    list_sCancerSym = dict_sPathways['Pathways in cancer - Homo sapiens (human)']

    for sGeneSym, nVarCnt, fMAF_PN, fMAF_HN, fMAF_tcga, fMAF_1000 in list_sFinal:

        if sGeneSym not in dict_cAnno:
            dict_cAnno[sGeneSym] = []
        dict_cAnno[sGeneSym].append([sGeneSym, fMAF_PN, fMAF_HN, fMAF_tcga, fMAF_1000])

    #loop END: cAnno
    print('dict_cAnno', len(dict_cAnno))

    for sPathway in dict_sPathways:

        if sPathway not in dict_sOutput:
            dict_sOutput[sPathway] = []

        for sGeneSym in dict_sPathways[sPathway]:

            if sGeneSym not in dict_cAnno: continue

            dict_sOutput[sPathway] += dict_cAnno[sGeneSym]
    #loop END: sPathway

    list_sPaths = []
    for sPathway in dict_sOutput:
        if not dict_sOutput[sPathway]: continue
        list_sPaths.append([sPathway, dict_sOutput[sPathway]])
    #loop END: sPathway

    list_sPaths  = sorted(list_sPaths, key=lambda e:len(e[1]), reverse=True)

    list_sOutput = []

    for sPathway, list_cAnno_filt in list_sPaths:
        list_fMAF_p     = []
        list_fMAF_h     = []
        list_fMAF_tcga  = []
        list_fMAF_1000  = []

        bCancer     = False

        for sGeneSym, fMAF_PN, fMAF_HN, fMAF_tcga, fMAF_1000 in list_cAnno_filt:
            list_fMAF_p.append(fMAF_PN)
            list_fMAF_h.append(fMAF_HN)
            list_fMAF_tcga.append(fMAF_tcga)
            list_fMAF_1000.append(fMAF_1000)

        # loop END: cAnno

        for sGeneSym, fMAF_PN, fMAF_HN, fMAF_tcga, fMAF_1000 in list_cAnno_filt:
            if sGeneSym in list_sCancerSym:
                bCancer = True
                break
        #loop END: cAnno

        fMedianP     = median(list_fMAF_p)
        fMedianH     = median(list_fMAF_h)

        fMedianT     = median(list_fMAF_tcga)
        fMedian1     = median(list_fMAF_1000)
        #fMedianT     = 1
        #fMedian1     = 1

        Z, fPvalue   = stats.ranksums(list_fMAF_p, list_fMAF_h)
        Z2, fPvalue2 = stats.ranksums(list_fMAF_tcga, list_fMAF_1000)
        #fPvalue2 = 1
        list_sOutput.append([sPathway, list_cAnno_filt, fPvalue, fPvalue2,
                             fMedianP, fMedianH, fMedianT, fMedian1, bCancer])
    #loop END: sPathway, list_cAnno_filt

    print('list_sOutput', len(list_sOutput))

    #list_sOutput  = sorted(list_sOutput, key=lambda e:(e[3] - e[4], e[2], len(e[1])), reverse=True)
    #list_sOutput  = sorted(list_sOutput, key=lambda e:len(e[1]), reverse=True)
    list_sOutput  = sorted(list_sOutput, key=lambda e:(e[3]))
    #list_sOutput  = sorted(list_sOutput, key=lambda e:e[3]-e[4], reverse=True)

    nCnt = 0
    for sPathway, list_cAnno, fPval, fPval2, fMedianP, fMedianH, fMedianT, fMedian1, bCancer in list_sOutput:

        if fMedianH > fMedianP:     continue
        if not bCancer:             continue
        #if fPval2 > 0.05:           continue
        #if len(list_cAnno) > 10000: continue
        nCnt += 1

        sOut = '%s\t%s\t%0.3E\t%0.3E\t%0.4f\t%0.4f\t%0.4f\t%0.4f' % \
               (sPathway, len(list_cAnno), fPval, fPval2, fMedianP, fMedianH, fMedianT, fMedian1)
        print(sOut)
    #loop END: sPathway, list_cAnno, fPval, fMedianP, fMedianH\

    print(green('Final Pathway Count %s' % nCnt,1))
#def END: pathway_survey

## region Add Colors

def red(string, e=0):     return '\033[%s31m%s\033[m'%('' if e == 0 else '1;', string)
def green(string, e=0):   return '\033[%s32m%s\033[m'%('' if e == 0 else '1;', string)
def yellow(string, e=0):  return '\033[%s33m%s\033[m'%('' if e == 0 else '1;', string)
def blue(string, e=0):    return '\033[%s34m%s\033[m'%('' if e == 0 else '1;', string)
def magenta(string, e=0): return '\033[%s35m%s\033[m'%('' if e == 0 else '1;', string)
def cyan(string, e=0):    return '\033[%s36m%s\033[m'%('' if e == 0 else '1;', string)
def white(string, e=0):   return '\033[%s37m%s\033[m'%('' if e == 0 else '1;', string)

def get_color (cMir, sResidue):
    if sResidue == cMir.sAltNuc: sResidue = red(sResidue,1)
    elif sResidue == cMir.sRefNuc: sResidue = green(sResidue,1)
    else: sResidue = blue(sResidue,1)
    return sResidue
#def END: get_color

## endregion

## endregion


## region Classes

## region Class cVCFData

class cVCFData:
    def __init__(self):
        self.sChrID         = ''
        self.nPos           = 0
        self.sDBSNP_ID      = ''
        self.sRefNuc        = ''
        self.sAltNuc        = ''
        self.fQual          = 0.0
        self.sFilter        = ''
        self.sInfo          = ''
        self.sFormat        = ''
        self.list_sSamples  = []
    #def END: __int__


def cVCF_parse_vcf_files (nSampleCnt, sVCFFile):

    print_start_msg('Parsing VCF file', sVCFFile)

    if not os.path.isfile(sVCFFile):
        sys.exit('File Not Found %s' % sVCFFile)

    list_sOutput = []
    InFile       = open(sVCFFile, 'r')
    nCnt         = 0

    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          | 8              | 9./..
        # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
        # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph

        if sReadLine.startswith('#'): continue  # SKIP Information Headers
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cVCF                = cVCFData()
        cVCF.sChrID         = list_sColumn[0]
        cVCF.nPos           = int(list_sColumn[1])
        cVCF.sDBSNP_ID      = list_sColumn[2]
        cVCF.sRefNuc        = list_sColumn[3]
        cVCF.sAltNuc        = list_sColumn[4]
        cVCF.fQual          = float(list_sColumn[5])
        cVCF.sFilter        = list_sColumn[6]
        cVCF.sInfo          = list_sColumn[7]
        cVCF.sFormat        = list_sColumn[8]
        cVCF.list_sSamples  = list_sColumn[9:]

        #V-S Check: Sample Size
        if len(cVCF.list_sSamples) != nSampleCnt:
            sys.exit('Invalid list of Samples : list_sSample size= %d : nSampleCnt= %d'
                     % (len(cVCF.list_sSamples), nSampleCnt))

        list_sOutput.append(cVCF)
    #loop END: sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cVCF_parse_vcf_files : list_sOutput : Size = %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cVCF_parse_vcf_files
## endregion

## region Class cANNOVAR

class cANNOVAR:
    def __init__(self):
        self.sNMID      = ''
        self.sChrID     = ''
        self.nStartPos  = 0
        self.nEndPos    = 0
        self.sRefNuc    = ''
        self.sAltNuc    = ''
        self.sGeneFunc  = ''    # ex) "intergenic"
        self.sGeneSym   = ''
        self.list_sMisc = []
        self.sDBSNP_ID  = []



        self.sMAFInfo   = ''    # ex) "0.002793   937.79  31      0       937.79  6"
        self.fMAF_PNorm = 0.0   # Patient-derived normal
        self.fMAF_HNorm = 0.0   # Healthy-derived normal
    #def END: __init__

def cAnno_parse_annovar_output (sInFile):

    print_start_msg('Parsing ANNOVAR File', sInFile)

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for i, sReadLine in enumerate(InFile):

        if sReadLine.startswith('Chr'):continue

        list_sColumn      = sReadLine.strip('\n').split(',')

        '''

        if i == 0:
            list_sHeader = list_sColumn
        elif i == 1:
            for i,(a,b) in enumerate(zip(list_sHeader, list_sColumn)):
                print('%s\t%s\t%s' % (i,a,b))
        else: break

        '''
        cAnno             = cANNOVAR()
        cAnno.sChrID      = list_sColumn[0]
        cAnno.nStartPos   = int(list_sColumn[1])
        cAnno.nEndPos     = int(list_sColumn[2])
        cAnno.sRefNuc     = list_sColumn[3]
        cAnno.sAltNuc     = list_sColumn[4]
        cAnno.sGeneFunc   = list_sColumn[5].strip('"')    # ex) "intergenic"
        cAnno.sGeneSym    = list_sColumn[6].upper().strip('"')
        cAnno.list_sMisc  = [sMisc.strip('"') for sMisc in list_sColumn[7:]]

        ## Deleterious Scores, indexed from end of list as variation exists for multiple isoforms
        cAnno.fSiftScore  = cAnno.list_sMisc[-8]
        cAnno.fPolyPhen2  = cAnno.list_sMisc[-6]
        cAnno.sMA_Predict = cAnno.list_sMisc[-3]
        cAnno.fGerp       = cAnno.list_sMisc[-2]

        ## Assign DBSNP ID
        for sMiscInfo in cAnno.list_sMisc:
            if 'rs' in sMiscInfo:
                cAnno.sDBSNP_ID = sMiscInfo.strip('"')
                break
            else:
                cAnno.sDBSNP_ID = 'NA'
        #loop END: sMiscInfo

        ## Assign NMID
        sInfoCheck = [sMisc for sMisc in cAnno.list_sMisc if 'NM_' in sMisc]
        if not sInfoCheck: ## No NMID
            cAnno.sNMID = 'None'
        else:
            sGeneInfo   = [sMisc for sMisc in cAnno.list_sMisc     if 'NM_' in sMisc][0]
            sNMID       = [sInfo for sInfo in sGeneInfo.split(':') if 'NM_' in sInfo][0]
            sKey        = sNMID.strip('"')
            cAnno.sNMID = sKey
        #if END:

        ## Assign healthy-derived MAF
        # ex) "0.002793   937.79  31      0       937.79  6"
        cAnno.sMAFInfo   = list_sColumn[-1].strip('"')
        cAnno.fMAF_PNorm = float(cAnno.sMAFInfo.split('\t')[0])
        cAnno.fMAF_HNorm = float(cAnno.sMAFInfo.split('\t')[3])

        list_sOutput.append(cAnno)
    #loop END: i, sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cAnno_parse_annovar_output : list_sOutput : Size = %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cAnno_parse_annovar_output

## endregion

## region Class cCosmic
class cCOSMIC:
    def __init__(self):
        self.sGeneName   = ''
        self.sAccID      = ''
        self.nCDSLen     = 0
        self.sHGCNID     = ''   # SKIP for now
        self.sSample     = ''   # SKIP for now
        self.sSampleID   = ''   # SKIP for now
        self.sTumorID    = ''   # SKIP for now
        self.sPriSite    = ''   # primary site  ex) pancreas
        self.sSiteSub1   = ''   # SKIP for now
        self.sSiteSub2   = ''   # SKIP for now
        self.sSiteSub3   = ''   # SKIP for now
        self.sPriHist    = ''   # primary histology
        self.sHistSub1   = ''   # SKIP for now
        self.sHistSub2   = ''   # SKIP for now
        self.sHistSub3   = ''   # SKIP for now
        self.bGenomeWide = ''   # ex) y or n
        self.sMutaID     = ''   # SKIP for now
        self.sAltType    = ''   # ex) c.35G>T
        self.sAAType     = ''   # ex) p.G12V
        self.sMutaDescri = ''   # ex) Substitution - Missense
        self.sMutaZygo   = ''   # SKIP for now
        self.bLOH        = ''   # loss of heterzygosity ex) y or n
        self.sGRCh       = ''
        self.sChrom      = ''
        self.sPos        = 0
        self.sStrand     = ''
        self.bSNP        = ''   # ex) y and n
        self.sDelete     = ''   # ex) PATHOGENIC
    #def END : __init__

def cCos_parse_cosmic_consensus (sInFile):

    print_start_msg('Parsing COSMIC File', sInFile)

    dict_Check   = {}
    list_sOutput = []
    InFile       = open(sInFile, 'r')

    for i, sReadLine in enumerate(InFile):

        if sReadLine.startswith('Gene'):continue

        list_sColumn = sReadLine.strip('\n').split('\t')
        '''
        if i == 0:
            list_sHeader = list_sColumn
        elif i == 1:
            for i,(a,b) in enumerate(zip(list_sHeader, list_sColumn)):
                print('%s\t%s\t%s' % (i,a,b))
        else: break

        '''
        cCos             = cCOSMIC()
        cCos.sGeneName   = list_sColumn[0].upper()
        cCos.sAccID      = list_sColumn[1]
        cCos.nCDSLen     = int(list_sColumn[2])
        cCos.sHGCNID     = list_sColumn[3]
        cCos.sSample     = list_sColumn[4]
        cCos.sSampleID   = list_sColumn[5]
        cCos.sTumorID    = list_sColumn[6]
        cCos.sPriSite    = list_sColumn[7]
        cCos.sSiteSub1   = list_sColumn[8]
        cCos.sSiteSub2   = list_sColumn[9]
        cCos.sSiteSub3   = list_sColumn[10]
        cCos.sPriHist    = list_sColumn[11]
        cCos.sHistSub1   = list_sColumn[12]
        cCos.sHistSub2   = list_sColumn[13]
        cCos.sHistSub3   = list_sColumn[14]
        cCos.bGenomeWide = True if list_sColumn[15] == 'y' else False
        cCos.sMutaID     = list_sColumn[16]
        cCos.sAltType    = list_sColumn[17]
        cCos.sAAType     = list_sColumn[18]
        cCos.sMutaDescri = list_sColumn[19]
        cCos.sMutaZygo   = list_sColumn[20]
        cCos.bLOH        = True if list_sColumn[21] == 'y' else False
        cCos.sGRCh       = list_sColumn[22]
        cCos.sChrom      = 'chr%s' % list_sColumn[23].split(':')[0]
        list_sPosCheck   = list(set(list_sColumn[23].split(':')[1].split('-')))
        cCos.sDelete     = list_sColumn[26]
        cCos.sMutaStatus = list_sColumn[28]


        if len(list_sPosCheck) > 1:
            cCos.sPos    = list_sPosCheck[0]
        else:
            cCos.sPos    = ''.join(list_sPosCheck)

        cCos.sStrand     = list_sColumn[24]
        cCos.bSNP        = True if list_sColumn[25] == 'y' else False

        if cCos.sPriSite not in dict_Check:
            dict_Check[cCos.sPriSite] = 0
        dict_Check[cCos.sPriSite] += 1

        cCos.sDelete     = list_sColumn[26]

        if cCos.sDelete == '':
            cCos.sDelete = 'NotDefined'

        list_sOutput.append(cCos)

    #loop END: i, sReadLine

    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cCos_parse_cosmic_consensus : list_sOutput : Size = %d' % (len(list_sOutput)))

    print(red('ALL COSMIC', 1), len(list_sOutput))
    return list_sOutput
#def END: cCos_parse_cosmic_consensus
#class END: cCosmic
## endregion

## region Class cGWAS

class cGWAS:
    def __init__(self):
        self.sDateAdded   = ''  # 06-May-2015
        self.sPubmedID    = ''  # 25917933
        self.sAuthor      = ''  # Zai CC
        self.sDate        = ''  # 20-Nov-2014
        self.sJournal     = ''  # J Psychiatr Res
        self.sURL         = ''  # http://europepmc.org/abstract/MED/25917933
        self.sStudyName   = ''  # A genome-wide association study of suicide severity scores in bipolar disorder.
        self.sDisease     = ''  # Suicide in bipolar disorder
        self.sSampleInfo  = ''  # 959 European ancestry individuals
        self.sRepSample   = ''  # NA
        self.sGeneRegion  = ''  # 10p11.22
        self.sChrID       = ''  # 10
        self.nPos         = 00  # 32704340
        self.sReportGenes = ''  # C10orf68, CCDC7, ITGB1
        self.sMappedGene  = ''  # CCDC7
        self.sUpGene      = ''  # SKIPPED, Upstream gene
        self.sDownGene    = ''  # SKIPPED, Downstream gene
        self.sSNPGenes    = ''  # SKIPPED
        self.sUpDist      = ''  # SKIPPED, Upstream gene distance
        self.sDownDist    = ''  # SKIPPED, Downstream gene distance
        self.sSNP_allele  = ''  # rs7079041-A
        self.sSNP_ID      = ''  # rs7079041
    #def END: __init__

def cGWAS_parse_gwas (sInFile):

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for i, sReadLine in enumerate(InFile):

        if sReadLine.startswith('DATE'):continue

        list_sColumn = sReadLine.strip('\n').split('\t')

        # V-S Check:
        if len(list_sColumn) != nGWAS_FILE_COLUMN_SIZE:
            sys.exit('ERROR: MAID Text File Column Size= %d\tContent= %s' % (len(list_sColumn), list_sColumn))
        cGwas              = cGWAS()
        cGwas.sDateAdded   = list_sColumn[0]
        cGwas.sPubmedID    = list_sColumn[1]
        cGwas.sAuthor      = list_sColumn[2]
        cGwas.sDate        = list_sColumn[3]
        cGwas.sJournal     = list_sColumn[4]
        cGwas.sURL         = list_sColumn[5]
        cGwas.sStudyName   = list_sColumn[6]
        cGwas.sDisease     = list_sColumn[7]
        cGwas.sSampleInfo  = list_sColumn[8]
        cGwas.sRepSample   = list_sColumn[9]
        cGwas.sGeneRegion  = list_sColumn[10]
        cGwas.sChrID       = 'chr%s' % list_sColumn[11]
        cGwas.nPos         = int(list_sColumn[12]) if list_sColumn[12] else 'NA'
        cGwas.sReportGenes = list_sColumn[13].upper()
        cGwas.sMappedGene  = list_sColumn[14].upper()
        cGwas.sUpGene      = list_sColumn[15]
        cGwas.sDownGene    = list_sColumn[16]
        cGwas.sSNPGenes    = list_sColumn[17]
        cGwas.sUpDist      = list_sColumn[18]
        cGwas.sDownDist    = list_sColumn[19]
        cGwas.sSNP_allele  = list_sColumn[20].split('-')[1] if '-' in list_sColumn[20] else 'NA'
        cGwas.sSNP_ID      = list_sColumn[21]

        list_sOutput.append(cGwas)
    #loop END: i, sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cGWAS_parse_gwas : list_sOutput : Size = %d' % (len(list_sOutput)))

    print(red('ALL GWAS', 1), len(list_sOutput))

    list_sOutput = [cGWAS for cGWAS in list_sOutput if 'cancer' in cGWAS.sDisease]

    print(red('FILTERED', 1), len(list_sOutput))

    return list_sOutput
#def END: cGWAS_parse_gwas

#class END:l cGWAS


## endregion

## region Class cExAC

class cExACData:
    def __init__(self):
        self.sChrID  = ''
        self.nPos    = 0
        self.nID     = '' #DBSNP
        self.sRefNuc = ''
        self.sAltNuc = ''
        self.fQual   = 0.0
        self.sFilter = ''
        self.sInfo   = ''
    #def END: __init__

def cExAC_parse_exac_data (sInFile):

    print_start_msg('Parsing ExAc File', sInFile)

    list_sOutput = []
    InFile       = open(sInFile, 'r')
    for i, sReadLine in enumerate(InFile):
        # File Format
        # Column Number:     | 0       | 1     | 2    | 3       | 4         |
        # Column Description:| sChrID  | nPos  | nID  | sRefNuc | sAltNuc   |
        # Column Example:    | 1       | 13372 | .    | G       | C         |
        # Column Number:     | 5       | 6       | 7
        # Column Description:| fQual   | sFilter | sInfo
        # Column Example:    | 608.91  | PASS    | Additional Details "|" separated

        if sReadLine.startswith('#'):continue

        list_sColumn = sReadLine.strip('\n').split('\t')

        cExAC            = cExACData()
        cExAC.sChrID     = 'chr%s' % list_sColumn[0]
        cExAC.nPos       = int(list_sColumn[1])
        cExAC.nID        = list_sColumn[2]
        cExAC.sRefNuc    = list_sColumn[3]
        cExAC.sAltNuc    = list_sColumn[4]
        cExAC.fQual      = float(list_sColumn[5])
        cExAC.sFilter    = list_sColumn[6]
        cExAC.sInfo      = list_sColumn[7]

        list_sOutput.append(cExAC)

    #loop END: i, sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cExAC_parse_exac_data : list_sOutput : Size = %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cCos_parse_cosmic_consensus




#endregion

## region Class cClinvar

class cClinvar:
    def __init__(self):
        self.sChrID         = ''
        self.nPos           = 0
        self.sDBSNP_ID      = ''
        self.sRefNuc        = ''
        self.sAltNuc        = ''
        self.fQual          = 0.0
        self.sFilter        = ''
        self.sInfo          = ''
        self.sFormat        = ''
        self.list_sInfo     = []
        self.dict_sInfo     = {}
        self.sDelete        = ''
    #def END: __int__

def cClin_parse_clinvar_file (sClinvarFile):

    print_start_msg('Parsing Clinvar file', sClinvarFile)

    list_sOutput = []
    InFile       = open(sClinvarFile, 'r')

    dict_sDelete = {0:'uncertain', 1:'uncertain', 2:'benign', 3:'likely_benign', 4:'likely_pathogenic',
                    5:'pathogenic', 6:'other', 7:'other', 255:'other'}

    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0       | 1        | 2          | 3       | 4
        # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
        # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
        # Column Number:     | 5       | 6        | 7          |
        # Column Description:| fQual   | sFilter  | sInfo      |
        # Column Example:    | 5645.6  | PASS     | .          |
        ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance,
        # 0 - Uncertain significance,
        # 1 - not provided
        # 2 - Benign
        # 3 - Likely benign
        # 4 - Likely pathogenic
        # 5 - Pathogenic
        # 6 - drug response
        # 7 - histocompatibility,
        # 255 - other

        if sReadLine.startswith('#'): continue  #SKIP Information Headers
        list_sColumn        = sReadLine.strip('\n').split('\t')

        cClin                = cClinvar()
        cClin.sChrID         = list_sColumn[0]
        cClin.nPos           = int(list_sColumn[1])
        cClin.sDBSNP_ID      = list_sColumn[2]
        cClin.sRefNuc        = list_sColumn[3]
        cClin.sAltNuc        = list_sColumn[4]
        cClin.fQual          = float(list_sColumn[5]) if list_sColumn[5] != '.' else '.'
        cClin.sFilter        = list_sColumn[6]
        cClin.list_sInfo     = list_sColumn[7].split(';')

        for sInfo in cClin.list_sInfo:

            if len(sInfo.split('=')) < 2: continue
            sInfoName, sValue = sInfo.split('=')

            if sInfoName not in cClin.dict_sInfo:
                cClin.dict_sInfo[sInfoName] = ''
            cClin.dict_sInfo[sInfoName] = sValue
        #loop END: sInfo

        # Multiple sig values
        list_sClinSig  = cClin.dict_sInfo['CLNSIG'].split('|')

        # Most common value
        list_sClingSig = max(set(list_sClinSig), key=list_sClinSig.count).split(',')

        # Higher value
        nClingSig      = int(max(list_sClingSig))

        #V-S Check:
        if nClingSig not in dict_sDelete:
            sys.exit('Invalid nClingSig : cClin_parse_clinvar_file : nClinSig= %s' % (nClingSig))

        cClin.sDelete  = dict_sDelete[nClingSig]

        list_sOutput.append(cClin)
    #loop END: sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cClin_parse_clinvar_file : list_sOutput : Size = %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cClin_parse_clinvar_file

## endregion

## region Class KY-Output

class cKY_Output:
    def __init__(self):
        self.sNMID        = ''
        self.fMAF_ration  = 0.0
        self.fMAF_diff    = 0.0
        self.nVarCnt      = 0
        self.sVarInfo     = 0   # [<Var1 Info>],[],[]
        self.sVarDiff     = 0   # MAF difference per variant
        self.sVarRegion   = 0   # MAF difference per variant
    #def END: __init_


def cKY_parse_analysis_file (sInFile):

    list_sOutput = []
    InFile       = open(sInFile, 'r')

    for i, sReadLine in enumerate(InFile):
        # File Format
        # Column Number:     | 0         | 1          | 2         | 3           |
        # Column Description:| sNMID     | fMAF_ratio | fAbDiff   | nVariantCnt |
        # Column Example:    | NM_000309 | 2.1810959  | 0.0426700 | 1           |
        # Column Number:     | 4                          | 5               | 6
        # Column Description:| sVarinfo                   | sVarDiff        | sVarRegion*
        # Column Example:    | [0.039  0.975   .       .] | 0.07542-0.03285 | exonic
        if sReadLine.startswith('DATE'):continue

        list_sColumn = sReadLine.strip('\n').split('\t')
        cKY               = cKY_Output()
        cKY.sNMID        = list_sColumn[0]
        cKY.fMAF_ration  = float(list_sColumn[1])
        cKY.fMAF_diff    = float(list_sColumn[2])
        cKY.nVarCnt      = int(list_sColumn[3])
        cKY.sVarInfo     = list_sColumn[4]
        cKY.sVarDiff     = list_sColumn[5]
        cKY.sVarRegion   = list_sColumn[6]

        list_sOutput.append(cKY)
    #loop END: i, sReadLine
    InFile.close()

    #V-S Check:
    if not list_sOutput:
        sys.exit('Empty List : cKY_parse_analysis_file : list_sOutput : Size = %d' % (len(list_sOutput)))

    return list_sOutput
#def END: cKY_parse_analysis_file


## endregion

## endregion


## region VCF Preprocess

def vcf_filter (sData):

    # First filter after recalibration of variants
    sWorkDir    = '/extdata6/Kyuhong/05_K1000E_vcf_comparison'
    sOutDir     = '%s/00_vcf_filter'                             % sWorkDir


    if sData == 'K1000E':
        sInFile     = '%s/Kor_All_vcall_recal.preprocessed_plus.vcf' % sWorkDir
        sOutFile    = '%s/Kor_All_recal.preprocessed.filtered_plus'  % sOutDir

    else: # TCGA Data - LUAD or STAD
        sInFile     = '%s/TCGA_%s_final.vcf'        % (sWorkDir, sData)
        sOutFile    = '%s/TCGA_%s_final'            % (sOutDir,  sData)


    # Queue Options
    sQueue      = '24_730_optiplex.q'
    sLogDir     = '%s/log_%s'            % (sOutDir, sTIME_STAMP)
    bTestRun    = False

    sScript     = 'vcftools --vcf %s '   % sInFile
    sScript    += '--max-alleles 2 '
    sScript    += '--minDP 10 '
    sScript    += '--minGQ 30 '
    sScript    += '--hwe 1e-06 '
    sScript    += '--max-missing 0.3 '
    sScript    += '--out %s '            % sOutFile
    sScript    += '--recode '
    sScript    += '--remove-filtered-all '

    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.VCF-Filter.%s'
                  % (sScript, sLogDir, sQueue, sData))
    #if END:
#def END: vcf_filter


def vcftools (sSource, bFilter):

    sFilter      = 'filtered' if bFilter == 'T' else 'unfiltered'

    sWorkDir     = '/extdata6/Kyuhong/05_K1000E_vcf_comparison'
    sOutDir      = '%s/03_vcftools/%s'                      % (sWorkDir, sSource)

    sInFile_vcf  = '%s/K1056E.%s.final.out.recode.vcf'      % (sWorkDir, sFilter)
    sFileTag     = '%s/%s.final.recal.%s'                   % (sOutDir, sSource, sFilter)

    sIDFile      = '%s/%s_list_wouts.txt'                   % (sOutDir, sSource)

    # Queue Options
    sQueue       = '24_730_optiplex.q'
    sLogDir      = '%s/log_%s'                              % (sOutDir, sTIME_STAMP)
    bTestRun     = False
    sTmpScript   = copy_temp_core_script('/extdata6/Jinman/00_K1000E_DB')


    # Separate VCF by Disease
    sScript      = 'vcftools --vcf %s '            % sInFile_vcf
    sScript     += '--max-missing 0.5 '
    sScript     += '--out %s '                     % sFileTag
    sScript     += '--keep %s '                    % sIDFile
    sScript     += '--recode;'

    # Get Allele Frequency
    sScript     += 'vcftools --vcf %s.recode.vcf ' % sFileTag
    sScript     += '--freq '
    sScript     += '--out %s; '                    % sFileTag

    # Get Variants with MAF
    sScript     += '%s dist_maf '                  % sTmpScript
    sScript     += '%s.frq %s.maf.txt;'            % (sFileTag, sFileTag)

    # Filter by MAF positions
    sScript     += 'vcftools --vcf %s.recode.vcf ' % sFileTag
    sScript     += '--out %s.maf '                 % sFileTag
    sScript     += '--positions %s.maf.txt '       % sFileTag
    sScript     += '--recode; '


    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.VCFTools.%s.%s'
                  % (sScript, sLogDir, sQueue, sSource, sFilter))
    #if END:
#def END: vcftools


def vcftools_tcga (sSource, bFilter):

    sFilter      = 'filtered' if bFilter == 'T' else 'unfiltered'

    sWorkDir    = '/extdata6/Kyuhong/05_K1000E_vcf_comparison'
    sOutDir     = '%s/03_vcftools/%s'               % (sWorkDir, sSource)
    sInFile_vcf = '%s/TCGA.%s.%s.final.vcf'         % (sWorkDir, sSource, sFilter)


    sFileTag    = '%s/%s.final.%s'                  % (sOutDir, sSource, sFilter)
    sIDFile     = '%s/%s_IDlist_whites.txt'         % (sOutDir, sSource)

    # Queue Options
    sQueue      = 'optiplex.q'
    sLogDir     = '%s/log_%s'                       % (sOutDir, sTIME_STAMP)
    bTestRun    = False
    sTmpScript  = copy_temp_core_script('/extdata6/Jinman/00_K1000E_DB')

    # Separate VCF by Disease
    sScript     = 'vcftools --vcf %s '              % sInFile_vcf
    sScript     += '--out %s '                      % sFileTag
    sScript     += '--keep %s '                     % sIDFile
    sScript     += '--recode;'

    # Get Allele Frequency
    sScript     += 'vcftools --vcf %s.recode.vcf '  % sFileTag
    sScript     += '--freq '
    sScript     += '--out %s; '                     % sFileTag

    # Get Variants with MAF
    sScript     += '%s dist_maf '                   % sTmpScript
    sScript     += '%s.frq %s.maf.txt;'             % (sFileTag, sFileTag)

    # Filter by MAF positions
    sScript     += 'vcftools --vcf %s.recode.vcf '  % sFileTag
    sScript     += '--out %s.maf '                  % sFileTag
    sScript     += '--positions %s.maf.txt '        % sFileTag
    sScript     += '--recode; '

    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.VCFTools.TCGA.%s.%s'
                  % (sScript, sLogDir, sQueue, sSource, sFilter))
        # if END:
# def END: vcftools_tcga


def vcftools_1000 (sSource, bFilter):

    sFilter      = 'filtered' if bFilter == 'T' else 'unfiltered'

    sWorkDir    = '/extdata6/Kyuhong/05_K1000E_vcf_comparison'
    sOutDir     = '%s/03_vcftools/%s'               % (sWorkDir, '1000GP')
    sInFile_vcf = '%s/%s.%s.final.vcf'              % (sWorkDir, '1000GP', sFilter)


    sFileTag    = '%s/%s.final.%s'                  % (sOutDir, sSource, sFilter)
    sIDFile     = '%s/%s_IDlist.txt'                % (sOutDir, sSource)

    # Queue Options
    sQueue      = 'optiplex.q'
    sLogDir     = '%s/log_%s'                       % (sOutDir, sTIME_STAMP)
    bTestRun    = False
    sTmpScript  = copy_temp_core_script('/extdata6/Jinman/00_K1000E_DB')

    # Separate VCF by Disease
    sScript     = 'vcftools --vcf %s '              % sInFile_vcf
    sScript     += '--out %s '                      % sFileTag
    sScript     += '--max-missing 0.5 '
    sScript     += '--keep %s '                     % sIDFile
    sScript     += '--recode;'

    # Get Allele Frequency
    sScript     += 'vcftools --vcf %s.recode.vcf '  % sFileTag
    sScript     += '--freq '
    sScript     += '--out %s; '                     % sFileTag

    # Get Variants with MAF
    sScript     += '%s dist_maf '                   % sTmpScript
    sScript     += '%s.frq %s.maf.txt;'             % (sFileTag, sFileTag)

    # Filter by MAF positions
    sScript     += 'vcftools --vcf %s.recode.vcf '  % sFileTag
    sScript     += '--out %s.maf '                  % sFileTag
    sScript     += '--positions %s.maf.txt '        % sFileTag
    sScript     += '--recode; '

    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.VCFTools.TCGA.%s.%s'
                  % (sScript, sLogDir, sQueue, sSource, sFilter))
        # if END:
# def END: vcftools_tcga


def dist_maf (sInFile, sOutFile):

    print('Getting MAF Positions')

    InFile  = open(sInFile, 'r')
    OutFile = open(sOutFile, 'w')

    for sReadLine in InFile:
        ''' File Format
        CHROM   POS     N_ALLELES       N_CHR   {ALLELE:FREQ}
        chrM    59      2       632     T:1     C:0
        '''
        if sReadLine.startswith('CHROM'):
            OutFile.write(sReadLine)
            continue

        list_sColumn = sReadLine.strip('\n').split('\t')

        sChrID       = list_sColumn[0]
        nPos         = int(list_sColumn[1])
        sMAF         = float(list_sColumn[4].split(':')[1])

        if sMAF < 1.0: OutFile.write('%s\t%s\n' % (sChrID, nPos))
    #loop END: sReadLine
    InFile.close()
    OutFile.close()
    print('Getting MAF Positions...DONE')
#def END: dis_maf

## endregion


## region ANNOVAR

def avinput (sSource, bFilter):

    sFilter     = 'filtered' if bFilter == 'T' else 'unfiltered'
    sNormal     = 'normal_plus' if 'plus' in sSource else 'normal_old'

    sWorkDir    = '/extdata6/Kyuhong/05_K1000E_vcf_comparison/'

    sVCFDir_T   = '%s/03_vcftools/%s'                       % (sWorkDir, sSource)
    sVCFDir_N   = '%s/03_vcftools/%s'                       % (sWorkDir, sNormal) #Ignore if LUAD or STAD (TCGA Data)

    if sSource not in ['LUAD', 'STAD']:
        sInVCF_T    = '%s/%s.final.recal.%s.maf.recode.vcf' % (sVCFDir_T, sSource, sFilter)
        sInVCF_N    = '%s/%s.final.recal.%s.maf.recode.vcf' % (sVCFDir_N, sNormal, sFilter)

    else: # TCGA Data
        sInVCF_T = '%s/%s.final.%s.maf.recode.vcf'          % (sVCFDir_T, sSource, sFilter)
        sInVCF_N = '%s/%s.final.%s.maf.recode.vcf'          % (sVCFDir_N, sNormal, sFilter)

    sOutDir     = '%s/04_ANNOVAR/00_avinput'                % sWorkDir
    sAVInput_T  = '%s/%s_%s.avinput'                        % (sOutDir, sSource, sFilter)
    sAVInput_N  = '%s/%s_%s.avinput'                        % (sOutDir, sNormal, sFilter)
    sOutFile    = '%s/%s_%s_wNorm.avinput'                  % (sOutDir, sSource, sFilter)

    #if not os.path.isfile(sInVCF_T): sys.exit('File Not Found %s' % sInVCF_T)
    #if not os.path.isfile(sInVCF_N): sys.exit('File Not Found %s' % sInVCF_N)

    # Queue Options
    sQueue      = 'optiplex.q'
    sLogDir     = '%s/log_%s'          % (sOutDir, sTIME_STAMP)
    bTestRun    = False
    sTmpScript  = copy_temp_core_script('/extdata6/Jinman/00_K1000E_DB')

    # Format Input for ANNOVAR
    sScript      = '%s -format vcf %s ' % (sAVINPUT, sInVCF_T)
    sScript     += '--outfile %s '      % sAVInput_T
    sScript     += '-allsample '
    sScript     += '-withfreq; '

    if not os.path.isfile (sAVInput_N):
        sScript += '%s -format vcf %s ' % (sAVINPUT, sInVCF_N)
        sScript += '--outfile %s '      % sAVInput_N
        sScript += '-allsample '
        sScript += '-withfreq; '
    #if END:

    if sSource not in ['LUAD', 'STAD']:
        sScript     += '%s attach_norm %s %s %s;' % (sTmpScript, sAVInput_T, sAVInput_N, sOutFile)


    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.AVINPUT.%s.%s -hold_jid Jinman.VCFTools.*'
                  % (sScript, sLogDir, sQueue, sSource, sFilter))
    #if END:
#def END: avinpt


def attach_norm (sInFile_T, sInFile_N, sOutFile):

    dict_fNorm = {}
    InFile_N   = open(sInFile_N, 'r')

    for sReadLine in InFile_N:
        ''' File Format
        chrM    66      66      G       A       0       158.25  4
        '''
        sCol  = sReadLine.strip('\n').split('\t')
        sKey  = '%s,%s,%s,%s,%s' % (sCol[0],sCol[1],sCol[2],sCol[3],sCol[4])
        sFreq = '%s\t%s\t%s'     % (sCol[5],sCol[6],sCol[7])

        if sKey not in dict_fNorm:
            dict_fNorm[sKey] = ''
        dict_fNorm[sKey] = sFreq
    #loop END: sReadLine

    InFile_N.close()

    InFile_T = open(sInFile_T, 'r')
    OutFile  = open(sOutFile, 'w')

    for sReadLine in InFile_T:

        sCol = sReadLine.strip('\n').split('\t')
        sKey = '%s,%s,%s,%s,%s' % (sCol[0], sCol[1], sCol[2], sCol[3], sCol[4])

        try:
            sOut = '%s\t%s\n' % (sReadLine.strip('\n'), dict_fNorm[sKey])
            OutFile.write(sOut)
        except KeyError: continue
    #lop END: sReadLine

    InFile_T.close()
    OutFile.close()
#def END: attach_norm


def annovar (sSource, bFilter):

    sFilter     = 'filtered' if bFilter == 'T' else 'unfiltered'
    bAll        = False  # False = Only Deleterious i.e. Clinvar and Cosmic70 omitted

    sOutDir     = '%s/01_table_anno/%s' % (sBASE_DIR, sSource)
    os.makedirs(sOutDir, exist_ok=True)

    if sSource not in ['LUAD', 'STAD']:
        sInFile              = '%s/00_avinput/%s_%s_wNorm.avinput'  % (sBASE_DIR, sSource, sFilter)
        if bAll: sOutTag     = '%s/%s_%s_wNorm'                     % (sOutDir, sSource, sFilter)
        else:    sOutTag     = '%s/%s_%s_wNorm_score'               % (sOutDir, sSource, sFilter)

    else: # TCGA Data
        sInFile          = '%s/00_avinput/%s_%s.avinput'            % (sBASE_DIR, sSource, sFilter)
        if bAll: sOutTag = '%s/%s_%s'                               % (sOutDir, sSource, sFilter)
        else:    sOutTag = '%s/%s_%s_score'                         % (sOutDir, sSource, sFilter)

    if not os.path.isfile(sInFile): sys.exit('File Not Found %s'    % sInFile)

    # Queue Options
    sQueue      = 'optiplex.q'
    sLogDir     = '%s/log_%s' % (sOutDir, sTIME_STAMP)
    bTestRun    = False


    # Annotate with ANNOVAR
    sScript     = '%s %s '   % (sANNOVAR, sInFile)
    sScript    += '/extdata6/Daekwan/util/annovar/humandb '
    sScript    += '-buildver hg19 '
    sScript    += '-out %s ' % sOutTag

    if bAll:
        sScript += '-protocol refGene,snp137NonFlagged,ljb26_sift,ljb26_pp2hvar,ljb26_ma,gerp++gt2,clinvar_20150629,cosmic70 '
        sScript += '-operation  g,f,f,f,f,f,f,f '

    else:
        sScript += '-protocol refGene,snp137NonFlagged,ljb26_sift,ljb26_pp2hvar,ljb26_ma,gerp++gt2 '
        sScript += '-operation  g,f,f,f,f,f, '

    sScript    += '-csvout '
    sScript    += '-nastring . '
    sScript    += '-otherinfo; '
    sScript    += 'cp %s.*_multianno.csv %s/01_table_anno' % (sOutTag, sBASE_DIR)


    if bTestRun:
        print(sScript)
    else:
        os.makedirs(sLogDir, exist_ok=True)
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.ANNOVAR.%s.%s.%s'
                  % (sScript, sLogDir, sQueue, sSource, 'All' if bAll else 'Delete', sFilter))
    #if END:
#def END: annovar
## endregion


## region Cosmic consensus data

def cosmic_anal(sCancer):

    bValidate        = False  # Check pathogenic or deleterious variants from database with ours

    ## MAF filter conditions
    fMAF_min         = 0.01  # Minimum patient-derived normal MAF
    fMAF_ratio       = 1.05   # PN MAF > HN MAF * Ratio
    nMinVarforTS     = 5     # Minimum number of variants in gene for TS candidacy

    sWorkDir         = '%s/05_COSMIC_census_gene'                  % sBASE_DIR
    sCosmicFile      = '%s/cosmic.v74.txt'                         % sREF_DIR
    sAnnoFile        = '%s/anno_output_%s.csv'                     % (sREF_DIR, sCancer)

    ## Sample Lists
    sSampleFile_PN   = '%s/%s_plus_IDlist.txt'                     % (sREF_DIR, sCancer)
    sSampleFile_HN   = '%s/normal_plus_IDlist.txt'                 % sREF_DIR

    ## VCFs from Lung or Gastric with the Healthy Norm VCF
    sVCFFile_PN      = '%s/%s.vcf'                                 % (sREF_DIR, sCancer)
    sVCFFile_HN      = '%s/normal.vcf'                             % sREF_DIR

    # Read sample files
    list_sSamples_PN = get_sample_list (sSampleFile_PN)
    list_sSamples_HN = get_sample_list (sSampleFile_HN)

    print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))

    # Load COSMIC data
    list_cCosmic     = cCos_parse_cosmic_consensus   (sCosmicFile)

    # Load ANNOVAR data
    list_cAnno       = cAnno_parse_annovar_output    (sAnnoFile)

    # Find variants annotated in COSMIC
    list_cAnno       = annotate_with_database        (list_cCosmic, list_cAnno, 'COSMIC', sCancer)

    # Parse VCF file
    list_cVCF_PN     = cVCF_parse_vcf_files         (len(list_sSamples_PN), sVCFFile_PN)
    list_cVCF_HN     = cVCF_parse_vcf_files         (len(list_sSamples_HN), sVCFFile_HN)

    print(magenta('*********************************TUMOR SUPPRESS**************************************', 1))
    list_sFinal, list_cAnno_wNMID = genotype_analysis ('TS', list_cAnno, list_sSamples_PN, list_sSamples_HN, list_cVCF_PN,
                                                      list_cVCF_HN, fMAF_min, fMAF_ratio, nMinVarforTS, bValidate)
    # Outfile for IPA
    #output_for_ipa (sWorkDir, list_sFinal, 'Cosmic', sCancer, 'TS', fMAF_min, fMAF_ratio)

    #if bValidate: evaluate_db      (list_sFinal, list_cAnno_wNMID, list_cCosmic, 'COSMIC')


    #print(magenta('**********************************ONCOGENE*************************************', 1))
    #list_cAnno_onc   = check_oncogene                (list_cAnno, fMAF_min, fMAF_ratio)
    # Convert list structure for genotype reassignment
    #list_sFinal      = genotype_analysis             (list_cAnno_onc, nSampleCnt_PN, nSampleCnt_HN, list_cVCF_PN, list_cVCF_HN, fMAF_min, fMAF_ratio)
#def END: cosmic_anal


def get_coverage (sCancer, sLabel, list_db, list_input):

    dict_sKey = {}
    for sKey in list_db:
        if sLabel == 'Pathway': continue

        sGeneSym, sPos = sKey.split(',')

        if sGeneSym not in dict_sKey:
            dict_sKey[sGeneSym] = ''
        dict_sKey[sGeneSym] = sPos
    #loop END: sKey

    list_sKey       = list(dict_sKey.keys()) if sLabel != 'Pathway' else list_db

    list_intersect  = list(set(list_sKey) & set(list_input))

    if list_intersect:
        #OutFile = open('/extdata6/Kyuhong/05_K1000E_vcf_comparison/04_ANNOVAR/08_coverage/%s_%s_tcga.txt', 'w')
        print('Intersect', len(list_intersect))

        if sLabel == 'Pathway':
            return (len(list_intersect) / len(list_input)), list_intersect
        else:
            return (len(list_intersect) / len(list_input)), [[sGeneSym, dict_sKey[sGeneSym]] for sGeneSym in list_intersect]
    else:
        return 0, []
#def END: get_coverage


def genotype_analysis (sAnalysis, list_cAnno, list_sSamples_PN, list_sSamples_HN, list_cVCF_PN, list_cVCF_HN, fMAF_min, fMAF_ratio, nMinTS, bValid):

    # Convert list structure for genotype reassignment
    list_sGenePos, list_cAnno = convert_data_structure        (list_cAnno)

    # Reassign genotypes for gene based on the genotypes of the variants within gene
    list_sGT_PN, dict_sIDs_P  = reassign_genotype             ('Patient', list_sSamples_PN, list_cVCF_PN, list_sGenePos)
    list_sGT_HN, dict_sIDs_H  = reassign_genotype             ('Healthy', list_sSamples_HN, list_cVCF_HN, list_sGenePos)

    # Recalculate MAF on gene level
    list_fMAF_PN              = recalculate_maf               (list_sGT_PN)
    list_fMAF_HN              = recalculate_maf               (list_sGT_HN)

    # Filter Genes by MAF and printout
    list_sFinal               = filter_by_maf_and_output      (sAnalysis, list_sGenePos, dict_sIDs_P, dict_sIDs_H, list_fMAF_PN, list_fMAF_HN, fMAF_min, fMAF_ratio, nMinTS, bValid)

    return list_sFinal, list_cAnno
#def END: genotype_analysis


def convert_data_structure (list_cAnno):

    dict_cAnno = {}  # Group by NMID

    for cAnno in list_cAnno:
        sKey = cAnno.sNMID

        if cAnno.sGeneFunc not in ['UTR3','UTR5','exonic','intronic']: continue
        if sKey == 'None': continue

        if sKey not in dict_cAnno:
            dict_cAnno[sKey] = []
        dict_cAnno[sKey].append(cAnno)
    #loop END: cAnno

    list_sOutput = []


    for sNMID in dict_cAnno:

        list_sChrPos = ['%s_%s' % (cAnno.sChrID,cAnno.nStartPos) for cAnno in dict_cAnno[sNMID]]
        list_sOutput.append([sNMID, len(dict_cAnno[sNMID]), list_sChrPos])
    #loop END: sNMID

    return list_sOutput, list_cAnno
#def END: convert_data_structure


def check_oncogene (list_cAnno, fMAF_min, fMAF_ratio):

    nMinCluster = 1
    nWindowSize = 1

    print(len(list_cAnno))

    dict_nPos = {}

    for cAnno in list_cAnno:

        if cAnno.fMAF_PNorm < fMAF_min:                      continue
        if cAnno.fMAF_PNorm < fMAF_ratio * cAnno.fMAF_HNorm: continue

        if cAnno.sGeneFunc not in ['UTR3','UTR5','exonic','intronic']: continue

        sKey = cAnno.nStartPos

        if sKey not in dict_nPos:
            dict_nPos[sKey] = ''
        dict_nPos[sKey] = cAnno
    #loop END: cAnno

    list_sOutput = []

    list_nPos    = list(dict_nPos.keys())
    list_sGroups = cluster(list_nPos, nWindowSize)
    list_sGroups = [list_nCluster for list_nCluster in list_sGroups if len(list_nCluster) > nMinCluster]

    if not list_sGroups:
        print(green('No Clusters', 1))
        return []
    else:
        print(green('Clusters', 1), len(list_sGroups),mean([len(list_nCluster) for list_nCluster in list_sGroups]))

        for list_nCluster in list_sGroups:
            for nPos in list_nCluster:list_sOutput.append(dict_nPos[nPos])
        #loop END: list_nCluster

        return list_sOutput
    #if END:
#def END: check_oncogene


def cluster (data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups
#def END: cluster


def evaluate_db (list_sFinalGenes, list_cAnno, list_cDatabase, sDatabase):

    if sDatabase == 'COSMIC':

        dict_cData       = {int(cCos.sPos) : cCos  for cCos in list_cDatabase}
        dict_cAnno       = {cAnno.sNMID    : cAnno for cAnno in list_cAnno if cAnno.sNMID != '' }


        for sNMID, nVarCnt_PN, fMAF_PN, fMAF_HN in list_sFinalGenes:

            sKey = dict_cAnno[sNMID].nStartPos
            print(dict_cData[sKey].sPos, dict_cData[sKey].sGeneName)

            print(sNMID, nVarCnt_PN, fMAF_PN, fMAF_HN, dict_cData[sKey].sFATHMM)
# def END: evaluate_db


## endregion


## region GWAS data

def gwas_anal(sCancer):

    bValidate        = True  # Check pathogenic or deleterious variants from database with ours


    ## MAF filter conditions
    fMAF_min         = 0.01  # Minimum patient-derived normal MAF
    fMAF_ratio       = 1.5   # PN MAF > HN MAF * Ratio
    nMinVarforTS     = 5     # Minimum number of variants in gene for TS candidacy

    sWorkDir         = '%s/06_GWAS_catalog'                        % sBASE_DIR
    sGwasFile        = '%s/gwas.v1.txt'                            % sREF_DIR

    sAnnoFile        = '%s/anno_output_%s.csv'                     % (sREF_DIR, sCancer)

    ## Sample Lists
    sSampleFile_PN   = '%s/%s_IDlist.txt'                          % (sREF_DIR, sCancer)
    sSampleFile_HN   = '%s/normal_IDlist.txt'                      % sREF_DIR

    ## VCFs from Lung or Gastric with the Healthy Norm VCF
    sVCFFile_PN      = '%s/%s.vcf'                                 % (sREF_DIR, sCancer)
    sVCFFile_HN      = '%s/normal.vcf'                             % sREF_DIR


    # Read sample files
    list_sSamples_PN = get_sample_list (sSampleFile_PN)
    list_sSamples_HN = get_sample_list (sSampleFile_HN)

    list_cGWAS       = cGWAS_parse_gwas             (sGwasFile)

    list_cAnno       = cAnno_parse_annovar_output   (sAnnoFile)

    # Find variants annotated in COSMIC
    list_cAnno       = annotate_with_database       (list_cGWAS, list_cAnno, 'GWAS', sCancer)

    # Parse VCF file
    list_cVCF_PN     = cVCF_parse_vcf_files         (len(list_sSamples_PN), sVCFFile_PN)
    list_cVCF_HN     = cVCF_parse_vcf_files         (len(list_sSamples_HN), sVCFFile_HN)

    print(magenta('*********************************TUMOR SUPPRESS**************************************', 1))
    list_sFinal      = genotype_analysis             ('TS', list_cAnno, list_sSamples_PN, list_sSamples_HN, list_cVCF_PN,
                                                      list_cVCF_HN, fMAF_min, fMAF_ratio, nMinVarforTS)
    # Outfile for IPA
    output_for_ipa (sWorkDir, list_sFinal, 'ExAC', sCancer, 'TS', fMAF_min, fMAF_ratio)
#def END: gwas_anal

## endregion


## region ExAC data

def exac_anal (sCancer):

    bValidate        = False  # Check pathogenic or deleterious variants from database with ours

    ## MAF filter conditions
    #fMAF_min         = 0.01  # Minimum patient-derived normal MAF
    #fMAF_ratio       = 1.25   # PN MAF > HN MAF * Ratio
    nMinVarforTS     = 5     # Minimum number of variants in gene for TS candidacy
    sWorkDir         = '%s/07_ExAC'                                % sBASE_DIR
    sExACFile        = '%s/exac.vcf'                               % sREF_DIR
    sAnnoFile        = '%s/anno_output_%s.csv'                     % (sREF_DIR, sCancer)


    ## Sample Lists
    sSampleFile_PN   = '%s/%s_IDlist.txt'                          % (sREF_DIR, sCancer)
    sSampleFile_HN   = '%s/normal_IDlist.txt'                      % sREF_DIR

    ## VCFs from Lung or Gastric with the Healthy Norm VCF
    sVCFFile_PN      = '%s/%s.vcf'                                 % (sREF_DIR, sCancer)
    sVCFFile_HN      = '%s/normal.vcf'                             % sREF_DIR


    # Read sample files
    list_sSamples_PN = get_sample_list (sSampleFile_PN)
    list_sSamples_HN = get_sample_list (sSampleFile_HN)

    #print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))


    # Temp pickle file for faster loading after inital parsing of ExAC file (~20gb)
    sTempFile        = '%s/temp/%s.ExAC.nonpsych.passonly.EAS.temp.data' % (sBASE_DIR, sCancer)
    if os.path.isfile(sTempFile):

        print('Found temp file')
        list_cAnno = pickle.load(open(sTempFile, 'rb'))

    else:

        # Load ExAC data
        list_cExAC       = cExAC_parse_exac_data        (sExACFile)

        # Load ANNOVAR data
        list_cAnno       = cAnno_parse_annovar_output   (sAnnoFile)

        # Find variants annotated in COSMIC
        list_cAnno       = annotate_with_database       (list_cExAC, list_cAnno, 'ExAC', sCancer)

        print(len(list_cAnno))

        TempFile        = open(sTempFile, 'wb')
        pickle.dump(list_cAnno, TempFile)
        TempFile.close()
    #if END:

    # Parse VCF file
    list_cVCF_PN     = cVCF_parse_vcf_files         (len(list_sSamples_PN), sVCFFile_PN)
    list_cVCF_HN     = cVCF_parse_vcf_files         (len(list_sSamples_HN), sVCFFile_HN)
    for fMAF_min, fMAF_ratio in [[0.05,1.25],[0.01,1.25],[0.05,1.5],[0.01,1.5]]:
        print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))

        print(magenta('*********************************TUMOR SUPPRESS******************************', 1))
        list_sFinal_TS, list_cAnno_TS   = genotype_analysis             ('TS', list_cAnno, list_sSamples_PN, list_sSamples_HN, list_cVCF_PN,
                                                          list_cVCF_HN, fMAF_min, fMAF_ratio, nMinVarforTS, bValidate)
        # Outfile for IPA
        output_for_ipa (sWorkDir, list_sFinal_TS, 'ExAC', sCancer, 'TS', fMAF_min, fMAF_ratio)

        print(green('*********************************ONCOGENE**************************************', 1))
        list_cAnno_ON                   = check_oncogene                (list_cAnno, fMAF_min, fMAF_ratio)
        list_sFinal_ON, list_cAnno_ON   = genotype_analysis             ('ON', list_cAnno_ON, list_sSamples_PN, list_sSamples_HN, list_cVCF_PN,
                                                          list_cVCF_HN, fMAF_min, fMAF_ratio, 1, bValidate)

        # Outfile for IPA
        output_for_ipa (sWorkDir, list_sFinal_ON, 'ExAC', sCancer, 'ON', fMAF_min, fMAF_ratio)
    #loop END: fMAF_min, fMAF_ratio
#def END:exac_anal

## endregion


## region Reproduction Validation

def reproduce ():

    sWorkDir         = '%s/02_MAF_comparison/03_oncogene_candidates' % sBASE_DIR
    sInFile_before   = '%s/02_lung_2x_0.05.analysis.txt'             % sWorkDir
    sInFile_after    = '%s/02_lung_wo_outs_2x_0.05.analysis.txt'     % sWorkDir

    list_cBefore     = cKY_parse_analysis_file (sInFile_before)
    list_cAfter      = cKY_parse_analysis_file (sInFile_after)

    dict_cBefore     = {cKY.sNMID: cKY for cKY in list_cBefore}
    dict_cAfter      = {cKY.sNMID: cKY for cKY in list_cAfter}

    list_sKey_b      = list(dict_cBefore.keys())
    list_sKey_a      = list(dict_cAfter.keys())

    list_sIntersect  = list(set(list_sKey_b) & set(list_sKey_a))
    list_sBeforeOnly = [sNMID for sNMID in list_sKey_b if sNMID not in list_sIntersect]
    list_sAfterOnly  = [sNMID for sNMID in list_sKey_a if sNMID not in list_sIntersect]

    print(len(list_sIntersect))
    print(len(list_sBeforeOnly))
    print(len(list_sAfterOnly))
#def END: reproduce

## endregion


## region Score-based (KH code clean-up)

def score (sCancer):

    sWorkDir             = '%s/03_Score_based_comparison'   % sBASE_DIR

    ## MAF filter conditions
    fMAF_min             = 0.01  # Minimum patient-derived normal MAF
    fMAF_ratio           = 1.5   # PN MAF > HN MAF * Ratio
    print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))

    ## Deleterious prediction conditions
    fSift_cutoff         = 0.05          # Less    = Deleterious
    fPolyPhen2_cutoff    = 0.5           # Greater = Deleterious
    list_sMutA_predict   = ['H', 'M']    # H or M  = Deleterious
    fGerp_cutoff         = 5             # Greater = Deleterious

    '''
    ## Temp Check Overlap with Variant-level analysis
    sFile1      = '%s/output/output_for_ipa_Score_Test_%s_%s_%s.txt' % (sWorkDir, sCancer, fMAF_min, fMAF_ratio)
    sFile2      = '%s/SCORE_%s_ipa.txt'                         % (sWorkDir, sCancer)
    list_sFile1 = []
    list_sFile2 = []

    for sReadLine in open(sFile1):
        list_sCol = sReadLine.strip('\n').split('\t')
        list_sFile1.append('%s' % (list_sCol[0]))

    for sReadLine in open(sFile2):
        list_sCol = sReadLine.strip('\n').split('\t')
        list_sFile2.append('%s' % (list_sCol[0]))


    print('list_sFile1', len(list(set(list_sFile1))))
    print('list_sFile2', len(list(set(list_sFile2))))

    list_sOverlap = list(set(list_sFile1) & set(list_sFile2))

    print('list_sOverlap', len(list_sOverlap))

    print([e for e in list_sFile1 if e not in list_sFile2])

    sys.exit()
    '''

    ## Annovar File
    sAnnoFile            = '%s/anno_output_%s_filtered.csv' % (sREF_DIR, sCancer)
    list_cAnno           = cAnno_parse_annovar_output(sAnnoFile)

    ## COSMIC File
    sCosmicFile          = '%s/cosmic.v74.txt'           % sREF_DIR
    list_cCosmic         = cCos_parse_cosmic_consensus   (sCosmicFile)
    '''

    list_cAnno           = tcga_1000gp_validation(sCancer, sWorkDir, list_cAnno_out, False)

    #OutFile = open('%s/list_cAnno_%s.data' % (sWorkDir, sCancer), 'wb')
    #pickle.dump(list_cAnno, OutFile)

    #InFile          = open('%s/list_cAnno_%s.data' % (sWorkDir, sCancer), 'rb')

    #list_cAnno      = pickle.load(InFile)
    #list_cAnno_out  = []
    #for cAnno in list_cAnno:

    #    fMAF_PN     = cAnno.fMAF_PNorm
    #    fMAF_HN     = cAnno.fMAF_HNorm

        # Apply MAF filter
    #    if fMAF_PN < fMAF_min:               continue
    #    if fMAF_PN < (fMAF_ratio * fMAF_HN): continue
    #    list_cAnno_out.append(cAnno)
    #loop END: cAnno


    #print(red('Total    :\t%d' % len(list_cAnno_out), 1))
    #pathway_survey(list_cAnno_out)
    '''

    list_cAnno_out  = []
    for cAnno in list_cAnno:

        if cAnno.sGeneFunc != 'exonic': continue

        fMAF_PN     = cAnno.fMAF_PNorm
        fMAF_HN     = cAnno.fMAF_HNorm

        ## Apply MAF filter
        if fMAF_PN < fMAF_min:               continue
        if fMAF_PN < (fMAF_ratio * fMAF_HN): continue

        ## Check Scores
        fSiftScore  = float(cAnno.fSiftScore) if cAnno.fSiftScore  != '.' else 1
        fPolyPhen2  = float(cAnno.fPolyPhen2) if cAnno.fPolyPhen2  != '.' else 0
        sMA_Predict = cAnno.sMA_Predict       if cAnno.sMA_Predict != '.' else 'A'
        fGerp       = float(cAnno.fGerp)      if cAnno.fGerp       != '.' else 0

        nExon       = 0  # Deleterious counter for exon region
        nNonExon    = 0  # Deleterious counter for non-exon region

        if fSiftScore < fSift_cutoff:         nExon += 1;
        if fPolyPhen2 > fPolyPhen2_cutoff:    nExon += 1;
        if sMA_Predict in list_sMutA_predict: nExon += 1;
        if fGerp > fGerp_cutoff:              nExon += 1; nNonExon += 1

        ## For exonic region, at least 2 out of 4 deleterious prediction
        if cAnno.sGeneFunc == 'exonic' and nExon >= 2:
            list_cAnno_out.append(cAnno)

        ## For non-exonic region, just GERP
        if cAnno.sGeneFunc != 'exonic' and nNonExon > 0:
            list_cAnno_out.append(cAnno)
    #loop END: cAnno

    print(red('Filtered :\t%d' % len(list_cAnno_out), 1))

    list_cAnno_tcga = tcga_1000gp_validation    ('Score', sCancer, sWorkDir, list_cAnno_out, True)
    #list_cAnno_cos  = annotate_with_database    (list_cCosmic, list_cAnno_out, 'COSMIC', sCancer)

    #pathway_survey(list_cAnno_tcga)
    #output_for_ipa (sWorkDir, list_cAnno_out, sCancer, 'Score', fMAF_min, fMAF_ratio)

#def END: score

## endregion


## region Gene-level (KH code clean-up)

def ts(sCancer):

    bValid           = False  # Check pathogenic or deleterious variants from database with ours

    sWorkDir         = '%s/02_MAF_comparison/02_tumor_suppressor_candidates' % sBASE_DIR
    dict_sCancer     = {'lung':'LUAD', 'gastric':'STAD'}  # Distinguish between K1000E and TCGA

    sTCGA_VCF        = '%s/%s.vcf'                      % (sREF_DIR, dict_sCancer[sCancer])
    sTCGA_ID_File    = '%s/%s_IDlist.txt'               % (sREF_DIR, dict_sCancer[sCancer])

    ## MAF filter conditions
    fMAF_min         = 0.01  # Minimum patient-derived normal MAF
    fMAF_ratio       = 1.5  # PN MAF > HN MAF * Ratio
    nMinVarforTS     = 5     # Minimum number of variants in gene for TS candidacy  # Read sample files

    ## Sample Lists
    list_sSamples_PN = get_sample_list ('%s/%s_plus_IDlist.txt' % (sREF_DIR, sCancer))
    list_sSamples_HN = get_sample_list ('%s/%s_plus_IDlist.txt' % (sREF_DIR, 'normal'))
    print('Patient', len(list_sSamples_PN))
    print('Healthy', len(list_sSamples_HN))

    ## Annovar File
    sAnnoFile        = '%s/anno_output_%s.csv' % (sREF_DIR, sCancer)
    list_cAnno       = cAnno_parse_annovar_output(sAnnoFile)
    list_cAnno_out   = [cAnno for cAnno in list_cAnno if cAnno.sGeneFunc == 'exonic' and cAnno.fMAF_PNorm]
    print(green('Just Exonic', 1))
    print(green('Total    :\t%d' % len(list_cAnno), 1))
    print(green('Filtered :\t%d' % len(list_cAnno_out), 1))

    # Read gene list
    list_sGenePos    = get_gene_pos (list_cAnno_out)
    #sNMIDFile        = '%s/output/output_%s.txt' % (sWorkDir, sCancer)
    #list_sNMIDCheck  = [sReadLine.strip('\n').split('\t') for sReadLine in open(sNMIDFile)]
    #dict_sNMIDCheck  = {sNMID: [fMAF_pn, fMAF_hn] for sNMID, nCnt, fMAF_pn, fMAF_hn in list_sNMIDCheck}
    #list_sGenePos    = [e for e in list_sGenePos if e[0] in dict_sNMIDCheck]

    print(len(list_sGenePos))
    #list_sGenePos    = parse_gene_poslist ('%s/%s_plus.gene.txt' % (sWorkDir, sCancer))

    ## VCFs from Lung or Gastric with the Healthy Norm VCF
    sVCFFile_PN      = '%s/%s.vcf' % (sREF_DIR, sCancer)
    sVCFFile_HN      = '%s/%s.vcf' % (sREF_DIR, 'normal')

    # Parse VCF file
    list_cVCF_PN     = cVCF_parse_vcf_files (len(list_sSamples_PN), sVCFFile_PN)
    list_cVCF_HN     = cVCF_parse_vcf_files (len(list_sSamples_HN), sVCFFile_HN)

    # Reassign genotypes for gene based on the genotypes of the variants within gene
    list_sGT_PN      = reassign_genotype ('Patient', list_sSamples_PN, list_cVCF_PN, list_sGenePos)
    list_sGT_HN      = reassign_genotype ('Healthy', list_sSamples_HN, list_cVCF_HN, list_sGenePos)

    # Recalculate MAF on gene level
    list_fMAF_PN     = recalculate_maf (list_sGT_PN)
    list_fMAF_HN     = recalculate_maf (list_sGT_HN)


    # Filter Genes by MAF
    list_sFinal      = filter_by_maf_and_output ('TS', list_fMAF_PN,
                                                 list_fMAF_HN, fMAF_min, fMAF_ratio, nMinVarforTS, bValid)


    print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))
    print(yellow('Before Filtering', 1), len(list_sGenePos))
    print(yellow('Final Gene Cnt',   1), len(list_sFinal))


    ## Get MAF from TCGA and 1000GP
    sFreqDir       = '/extdata6/Kyuhong/05_K1000E_vcf_comparison/03_vcftools'
    sFreqFile_tcga = '%s/%s/%s.final.unfiltered.frq'     % (sFreqDir, dict_sCancer[sCancer], dict_sCancer[sCancer])
    sFreqFile_1000 = '%s/%s/%s_EUR.final.filtered.frq'   % (sFreqDir, '1000GP', '1000GP')
    dict_fMAF_tcga = get_maf_from_freq_file (sFreqFile_tcga)
    dict_fMAF_1000 = get_maf_from_freq_file (sFreqFile_1000)
    print(len(dict_fMAF_tcga))
    print(len(dict_fMAF_1000))

    list_sOutput   = []
    list_fMAF_p    = []
    list_fMAF_h    = []
    for sGeneSym, nVarCnt_PN, fMAF_PN, fMAF_HN, list_sPosInfo in list_sFinal:

        list_fMAF_tcga = []
        list_fMAF_1000 = []
        for sPosInfo in list_sPosInfo:

            sChrID, sPos = sPosInfo.split('_')

            try:
                fMAF_tcga   = dict_fMAF_tcga[int(sPos)]
                fMAF_1000   = dict_fMAF_1000[int(sPos)]
            except KeyError: continue
            list_fMAF_tcga.append(fMAF_tcga)
            list_fMAF_1000.append(fMAF_1000)
            list_fMAF_p.append(fMAF_tcga)
            list_fMAF_h.append(fMAF_1000)

        #loop END: sPosInfo

        #fMAF_PN, fMAF_HN = dict_sNMIDCheck[sNMID]

        #if len(list_fMAF_tcga) > 4:
        list_sOutput.append([sGeneSym, len(list_fMAF_tcga), fMAF_PN, fMAF_HN, median(list_fMAF_tcga), median(list_fMAF_1000)] )
    #loop END:

    OutFile = open('%s/TS_wTCGA-1000_%s_updated.txt' % (sWorkDir, sCancer), 'w')
    for sNMID, nVarCnt,fMAF_PN, fMAF_HN, fMAF_tcga, fMAF_1000 in list_sOutput:
        sOut = '%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (sNMID, nVarCnt,fMAF_PN, fMAF_HN, fMAF_tcga, fMAF_1000)
        OutFile.write(sOut)
    OutFile.close()

    print(len(list_fMAF_p), len(list_fMAF_h))
    Z, fPvalue = stats.ranksums(list_fMAF_p, list_fMAF_h)
    print(red('Pvalue           :\t%0.2E'% fPvalue, 1))

    pathway_survey_ts(list_sOutput)

    # Outfile for IPA
    #output_for_ipa_ts (sWorkDir, list_sFinal, sCancer, 'TS', fMAF_min, fMAF_ratio)
#def END: recalc_maf



def get_sample_list (sSampleFile):

    if not os.path.isfile(sSampleFile):
        sys.exit('FileNotFound: %s' % sSampleFile)

    list_sSamples   = []
    InFile          = open(sSampleFile, 'r')

    for sReadLine in InFile:
        #File Format: Single Line with X amount of columns
        list_sSamples.append(sReadLine.strip('\n'))
    #loop END: sReadLine
    InFile.close()

    if not list_sSamples:
        sys.exit('Empty List : get_sample_list : list_sSample Size= %d' % len(list_sSamples))

    return list_sSamples
#def END: get_sample_list


def parse_gene_poslist (sGeneFile):

    print_start_msg('Parsing Gene List')

    if not os.path.isfile(sGeneFile):
        sys.exit('FileNotFound: %s' % sGeneFile)

    list_sGenePos   = []
    InFile          = open(sGeneFile, 'r')
    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0         | 1       | 2
        # Column Description:| sNMID     | nVarCnt | list_sPosInfo
        # Column Example:    | NM_000059 | 6       | chr13_32906558,chr13_32910842,chr13_32911042...

        if sReadLine.startswith('#'): continue

        list_sColumn  = [sColumn for sColumn in sReadLine.strip('\n').split(' ') if sColumn] # No empty elements

        sNMID         = list_sColumn[0]
        nVarCnt       = int(list_sColumn[1])
        list_sPosInfo = list_sColumn[2].split(',')

        list_sGenePos.append([sNMID, nVarCnt, list_sPosInfo])

    #loop END: sReadLine
    InFile.close()

    #V-S Check: Empty List
    if not list_sGenePos:
        sys.exit('Empty List : get_sample_list : list_sSample Size= %d' % len(list_sGenePos))

    return list_sGenePos
#def END: parse_gene_poslist


def get_gene_pos (list_cAnno):

    dict_sGeneSym = {}

    for cAnno in list_cAnno:

        if cAnno.sNMID == 'None': continue

        sKey = cAnno.sGeneSym
        sPos = '%s_%s' % (cAnno.sChrID, cAnno.nStartPos)
        if sKey not in dict_sGeneSym:
            dict_sGeneSym[sKey] = []
        dict_sGeneSym[sKey].append(sPos)
    #loop END: cAnno

    #V-S Check:
    if not dict_sGeneSym:
        sys.exit('Empty Dictionary : get_gene_pos : dict_sOutput size=%d' % (len(dict_sGeneSym)))

    list_sOutput = []
    for sGeneSym in dict_sGeneSym:
        list_sOutput.append([sGeneSym, len(dict_sGeneSym[sGeneSym]), dict_sGeneSym[sGeneSym]])

    return list_sOutput
#def END: get_gene_pos


def reassign_genotype (sSampleSource, list_sSampleIDs, list_cVCF, list_sGenePos):

    print_start_msg('Assigning gene level genotypes %s-derived Normal' % sSampleSource)

    dict_cVCF         = {'%s:%s' % (cVCF.sChrID, cVCF.nPos): cVCF for cVCF in list_cVCF}
    dict_sGTScore     = {'./.':-1, '0/0':0, '0/1':1, '1/0':1, '1/1':2} # Based on VCFtools criteria
    list_sOutput      = []
    nTotVarCnt        = 0
    nSampleCnt        = len(list_sSampleIDs)

    for sGeneSym, nVarCnt, list_sPosinfo in list_sGenePos:

        nTotVarCnt += len(list_sPosinfo)

        list_sGT_gene    = [[] for i in range(nSampleCnt)]
        list_sUpdatedPos = []
        for sPosInfo in list_sPosinfo:
            # sPosInfo format: chr13_32906558
            sChrID, sPos   = sPosInfo.split('_')
            sKey           = '%s:%s' % (sChrID, sPos)

            if sKey not in dict_cVCF: continue
            list_sSamples  = dict_cVCF[sKey].list_sSamples

            assert len(list_sGT_gene) == nSampleCnt == len(list_sSamples)

            for i, sSampleInfo in enumerate(list_sSamples):
                # sSampleInfo format: GT:AD:DP:GQ:PL
                ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
                ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
                ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
                ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
                sGT, sAD, sDP, sGQ, sPL = sSampleInfo.split(':')

                if sGT not in dict_sGTScore: continue

                list_sGT_gene[i].append(dict_sGTScore[sGT])

            #loop END: sSampleInfo

            #print(blue('PatientCnt',1), sNMID, sPosInfo, nPatCnt, len(list_sIDs))
            list_sUpdatedPos.append(sPosInfo)
        #loop END: sPosInfo

        # Update variant count with annotated variants
        nUpdatedVarCnt = [len(list_sGT_pos) for list_sGT_pos in list_sGT_gene][0]

        #print(nUpdatedVarCnt, len(list_sUpdatedPos))
        #assert nUpdatedVarCnt == len(list_sUpdatedPos)

        # Cases where there is only 1 variant and not in VCF
        if nUpdatedVarCnt == 0: continue
        try:
            list_sGT_gene  = [max(list_sGT_pos) for list_sGT_pos in list_sGT_gene if list_sGT_pos]
        except ValueError:
            print(list_sGT_gene)
            sys.exit()

        list_sOutput.append([sGeneSym, nUpdatedVarCnt, list_sGT_gene, list_sUpdatedPos])
    #loop END: sNMID, nVarCnt, list_sPosinfo

    #V-S Check: List size
    if not list_sOutput:
        sys.exit('Invalid lists : reassign_genotype : list_sOutput size= %d' % len(list_sOutput))

    print(red('TotalVarCnt', 1), nTotVarCnt)

    return list_sOutput
#def END: reassign_genotype


def recalculate_maf (list_sGenePos):

    print_start_msg('Recalculating MAF based on genotypes')

    list_sOutput = []
    for sGeneSym, nUpdatedVarCnt, list_sGT_gene, list_sPosInfo in list_sGenePos:

        fTotal  = sum([        2 for sGT_score in list_sGT_gene if sGT_score != -1])
        fAllele = sum([sGT_score for sGT_score in list_sGT_gene if sGT_score > 0])
        fMAF    = fAllele/fTotal
        list_sOutput.append([sGeneSym, nUpdatedVarCnt, fMAF, list_sPosInfo])
    #loop END: sNMID, nUpdatedVarCnt, list_sGT_gene

    #V-S Check: Empty List
    if not list_sOutput:
        sys.exit('Empty List : recalculate_maf : list_sGenePos Size= %d' % len(list_sGenePos))

    return list_sOutput
#def END: recalculate_maf


def filter_by_maf_and_output (sAnalysis, list_fMAF_PN, list_fMAF_HN, fMAF_min, fMAF_ratio, nMinTS, bValid):

    list_sFinalGenes = []

    dict_fMAF_HN     = {sGeneSym : [nVarCnt, fMAF] for sGeneSym, nVarCnt, fMAF, list_sPos in list_fMAF_HN}

    for sGeneSym, nVarCnt_PN, fMAF_PN, list_sPosInfo in list_fMAF_PN:

        if sGeneSym not in dict_fMAF_HN: continue  # No NMID
        if sAnalysis == 'TS':
            if nVarCnt_PN <= nMinTS:  continue  # Not enough variants

        nVarCnt_HN, fMAF_HN = dict_fMAF_HN[sGeneSym]


        if bValid:  # Check Normal-derived MAF as Greater than Patient-derived MAF
            if fMAF_HN < fMAF_min:               continue
            if fMAF_HN < (fMAF_ratio * fMAF_PN): continue

        else:
            if fMAF_PN < fMAF_min:               continue
            if fMAF_PN < (fMAF_ratio * fMAF_HN): continue
        #if END:

        list_sFinalGenes.append([sGeneSym, nVarCnt_PN, fMAF_PN, fMAF_HN, list_sPosInfo])

    #loop END: sNMID, nVarCnt, list_sPosInfo


    return list_sFinalGenes
#def END: filter_by_maf_and_output


def filter_by_maf_and_output2 (sAnalysis, list_sGenePos, list_fMAF_PN, list_fMAF_HN, fMAF_min, fMAF_ratio, nMinTS, bValid):

    list_sFinalGenes = []

    dict_fMAF_PN     = {sNMID : [nVarCnt, fMAF, list_sPos] for sNMID, nVarCnt, fMAF, list_sPos in list_fMAF_PN}
    dict_fMAF_HN     = {sNMID : [nVarCnt, fMAF, list_sPos] for sNMID, nVarCnt, fMAF, list_sPos in list_fMAF_HN}

    for sNMID, nVarCnt, list_sPosInfo in list_sGenePos:

        if sNMID not in dict_fMAF_PN: continue  # Not found
        if sNMID not in dict_fMAF_HN: continue  # Not found

        nVarCnt_PN, fMAF_PN, list_sPos = dict_fMAF_PN[sNMID]
        nVarCnt_HN, fMAF_HN, list_sPos = dict_fMAF_HN[sNMID]

        if sAnalysis == 'TS':
            if nVarCnt_PN < nMinTS:      continue  # Not enough variants

        if bValid:  # Check Normal-derived MAF as Greater than Patient-derived MAF
            if fMAF_HN < fMAF_min:               continue
            if fMAF_HN < (fMAF_ratio * fMAF_PN): continue

        else:
            if fMAF_PN < fMAF_min:               continue
            if fMAF_PN < (fMAF_ratio * fMAF_HN): continue
        #if END:


        list_sFinalGenes.append([sNMID, nVarCnt_PN, fMAF_PN, fMAF_HN, list_sPos])

    #loop END: sNMID, nVarCnt, list_sPosInfo

    return list_sFinalGenes
#def END: filter_by_maf_and_output2


def output_for_ipa (sWorkDir, list_cAnno, sCancer, sAnalysis, fMAF_min, fMAF_ratio):
    print_start_msg('Output for IPA %s %s %s' % (sCancer, sAnalysis, fMAF_ratio))

    sOutDir = '%s/output' % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)

    OutFile = open('%s/output_for_ipa_%s_%s_%s_%s.txt' % (sOutDir, sAnalysis, sCancer, fMAF_min, fMAF_ratio), 'w')

    for cAnno in list_cAnno:
        sOut = '%s\t%0.5f\t%0.5f\n' % (cAnno.sNMID, cAnno.fMAF_PNorm, cAnno.fMAF_HNorm)
        OutFile.write(sOut)
        # print(sOut[:-1])
    # loop END: sNMID, nVarCnt, fMAF_PN, fMAF_HN
    OutFile.close()
# def END: output_for_ipa

def output_for_ipa_ts (sWorkDir, list_sOutput, sCancer, sAnalysis, fMAF_min, fMAF_ratio):

    print_start_msg('Output for IPA %s %s %s' % (sCancer, sAnalysis, fMAF_ratio))

    sOutDir = '%s/output' % sWorkDir
    os.makedirs(sOutDir, exist_ok=True)

    OutFile = open('%s/output_for_ipa_%s_%s_%s_%s.txt' % (sOutDir, sAnalysis, sCancer, fMAF_min, fMAF_ratio), 'w')

    for sNMID, nVarCnt, fMAF_PN, fMAF_HN in list_sOutput:

        sOut = '%s\t%s\t%0.5f\t%0.5f\n' % (sNMID, nVarCnt, fMAF_PN, fMAF_HN)
        OutFile.write(sOut)
        #print(sOut[:-1])
    #loop END: sNMID, nVarCnt, fMAF_PN, fMAF_HN
    OutFile.close()
#def END: output_for_ipa_ts

## endregion


## region Variant-level (KH code clean-up)

def onco (sCancer):

    sWorkDir             = '%s/02_MAF_comparison/03_oncogene_candidates' % sBASE_DIR
    dict_sCancer         = {'lung':'LUAD', 'gastric':'STAD'}  # Distinguish between K1000E and TCGA
    sTCGA_VCF            = '%s/%s.vcf'                      % (sREF_DIR, dict_sCancer[sCancer])
    sTCGA_ID_File        = '%s/%s_IDlist.txt'               % (sREF_DIR, dict_sCancer[sCancer])


    ## MAF filter conditions
    fMAF_min             = 0.05  # Minimum patient-derived normal MAF
    fMAF_ratio           = 1.75   # PN MAF > HN MAF * Ratio

    ## Annovar File
    sAnnoFile            = '%s/anno_output_%s.csv' % (sREF_DIR, sCancer)
    list_cAnno           = cAnno_parse_annovar_output(sAnnoFile)

    list_cAnno_out  = []
    dict_nPos       = {}

    for cAnno in list_cAnno:

        #if cAnno.sNMID == 'None': continue
        if cAnno.sGeneFunc != 'exonic': continue
        fMAF_PN     = cAnno.fMAF_PNorm
        fMAF_HN     = cAnno.fMAF_HNorm

        ## Apply MAF filter
        if fMAF_PN < fMAF_min:               continue
        if fMAF_PN < (fMAF_ratio * fMAF_HN): continue

        if cAnno.nStartPos not in dict_nPos:
            dict_nPos[cAnno.nStartPos] = []
        dict_nPos[cAnno.nStartPos].append(cAnno)
        list_cAnno_out.append(cAnno)
    #loop END: cAnno

    print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))
    print(red('Total cAnno :\t%d' % len(list_cAnno), 1))
    print(red('Positions   :\t%d' % len(dict_nPos), 1))

    #list_cAnno_tcga = tcga_1000gp_validation('Onco', sCancer, sWorkDir, list_cAnno_out, True)
    #pathway_survey(list_cAnno_tcga)
    output_for_ipa(sWorkDir, list_cAnno_out, sCancer, 'Onco', fMAF_min, fMAF_ratio)

    # Variant Window   - Scrapped
    '''
    list_cAnno_out = []
    nMinCluster    = 1
    nWindowSize    = 5
    list_nPos      = list(dict_nPos.keys())
    list_sGroups   = cluster(list_nPos, nWindowSize)
    list_sGroups   = [list_nCluster for list_nCluster in list_sGroups if len(list_nCluster) > nMinCluster]

    print(green('Clusters= %d Win= %d' % (len(list_sGroups), nWindowSize), 1))

    list_sGroups = sorted(list_sGroups, key=lambda e:len(e), reverse=True)
    for list_nCluster in list_sGroups:

        for nPos in list_nCluster:
            list_cAnno_out += dict_nPos[nPos]

    #loop END: list_nCluster
    print(red('Filtered :\t%d' % len(list_cAnno_out), 1))
    '''

#def END: onco

## endregion


## region Calculate COSMIC, GWAS, and TCGA Coverage

def calc_coverage (sCancer):

    ## MAF filter conditions
    fMAF_min         = 0.01  # Minimum patient-derived normal MAF
    fMAF_ratio       = 1.5   # PN MAF > HN MAF * Ratio
    print(magenta('MAF_PN > %s and %s * MAF_HN' % (fMAF_min, fMAF_ratio), 1))

    sWorkDir         = '%s/08_coverage'                            % sBASE_DIR
    sCosmicFile      = '%s/cosmic.v74.txt'                         % sREF_DIR
    sGwasFile        = '%s/gwas.v1.txt'                            % sREF_DIR
    sAnnoFile        = '%s/anno_output_%s.csv'                     % (sREF_DIR, sCancer)

    dict_sCancer     = {'lung':'LUAD', 'gastric':'STAD'}  # Distinguish between K1000E and TCGA
    sTCGA_VCF        = '%s/%s.vcf'                            % (sREF_DIR, dict_sCancer[sCancer])
    sTCGA_ID_File    = '%s/%s_IDlist.txt'                          % (sREF_DIR, dict_sCancer[sCancer])

    list_sFile_TCGA  = get_sample_list (sTCGA_ID_File)

    list_score       = [line.strip('\n').split('\t')[3].upper() for line in open('%s/score_%s.txt' % (sWorkDir, sCancer))]
    list_ts          = [line.strip('\n').split('\t')[3].upper() for line in open('%s/ts_%s.txt'     % (sWorkDir, sCancer))]
    list_on          = [line.strip('\n').split('\t')[3].upper() for line in open('%s/on_%s.txt'     % (sWorkDir, sCancer))]
    #list_epith       = [line.strip('\n').split('\t')[0].upper() for line in open('%s/EpitSignaling_%s2.txt'  % (sWorkDir, sCancer))]
    #list_leukemia    = [line.strip('\n').split('\t')[0].upper() for line in open('%s/leukemia_%s.txt'  % (sWorkDir, sCancer))]
    #list_wnt         = [line.strip('\n').split('\t')[0].upper() for line in open('%s/wntpath_%s.txt'  % (sWorkDir, sCancer))]

    #print(list_epith)
    #print(list_on)
    #coverage_score, list_inter = get_coverage (sCancer, 'Pathway', list_epith, list_on)
    #print(list_inter)

    # Load ANNOVAR data
    list_cAnno       = cAnno_parse_annovar_output    (sAnnoFile)
    print('list_cAnno', len(list_cAnno))

    '''
    dict_sGeneSyms   = {}
    for cAnno in list_cAnno:

        fMAF_PN = cAnno.fMAF_PNorm
        fMAF_HN = cAnno.fMAF_HNorm

        if fMAF_PN < fMAF_min:               continue
        if fMAF_PN < (fMAF_ratio * fMAF_HN): continue

        if cAnno.sGeneSym == 'NONE': continue
        if cAnno.sGeneSym not in dict_sGeneSyms:
            dict_sGeneSyms[cAnno.sGeneSym] = 0
        dict_sGeneSyms[cAnno.sGeneSym] += 1

    print('Genes', len(dict_sGeneSyms))
    OutFile = open('%s/%s_fullgenelist.txt' % (sWorkDir, sCancer), 'w')
    for sGeneSym in dict_sGeneSyms:
        OutFile.write('%s\n' % sGeneSym)
    OutFile.close()

    sys.exit()
    '''

    list_sFiltered   = []
    for cAnno in list_cAnno:

        fMAF_PN = cAnno.fMAF_PNorm
        fMAF_HN = cAnno.fMAF_HNorm

        if fMAF_PN < fMAF_min:               continue
        if fMAF_PN < (fMAF_ratio * fMAF_HN): continue

        #if cAnno.sGeneSym in ['ACTN2', 'MYH2', 'MYO7A', 'NOTCH1', 'NOTCH2',
        #                      'PVRL3', 'TCF3', 'TUBA8', 'TUBB8', 'TUBB4B']:
        #if cAnno.sGeneSym in ['TUBA8', 'MYO7A']:
        #if cAnno.sGeneSym in list_inter:
            #print(cAnno.nStartPos, cAnno.sGeneSym, cAnno.sRefNuc, cAnno.sAltNuc, cAnno.sGeneFunc, cAnno.fMAF_PNorm, cAnno.fMAF_HNorm, cAnno.list_sMisc)

        list_sFiltered.append(cAnno)
    #loop END:

    if not list_sFiltered: sys.exit()

    # Load COSMIC, GWAS, and TCGA data
    list_cCosmic     = cCos_parse_cosmic_consensus   (sCosmicFile)
    list_cGWAS       = cGWAS_parse_gwas              (sGwasFile)
    list_cVCF_tcga   = cVCF_parse_vcf_files          (len(list_sFile_TCGA), sTCGA_VCF)

    # Find variants annotated in database
    list_cAnno_cosm  = annotate_with_database        (list_cCosmic, list_sFiltered, 'COSMIC', sCancer)
    list_cAnno_gwas  = annotate_with_database        (list_cGWAS, list_sFiltered, 'GWAS', sCancer)
    list_cAnno_tcga  = annotate_with_database        (list_cVCF_tcga, list_sFiltered, 'TCGA', sCancer)

    print('Anno Cosm', len(list_cAnno_cosm))
    print('Anno Gwas', len(list_cAnno_gwas))
    print('Anno TCGA', len(list_cAnno_tcga))

    sys.exit()

    # Convert list structure for genotype reassignment
    list_sGenes_cos  = ['%s,%s' % (cAnno.sGeneSym, cAnno.nStartPos) for cAnno in list_cAnno_cosm]
    list_sGenes_gwas = ['%s,%s' % (cAnno.sGeneSym, cAnno.nStartPos) for cAnno in list_cAnno_gwas]
    list_sGenes_tcga = ['%s,%s' % (cAnno.sGeneSym, cAnno.nStartPos) for cAnno in list_cAnno_tcga]

    list_all_three   = list(set(list_sGenes_cos) & set(list_sGenes_tcga))

    print('Cosmic')
    coverage_score, l1 = get_coverage (sCancer, 'Score', list_sGenes_cos, list_score)
    coverage_ts, l2    = get_coverage (sCancer, 'TS', list_sGenes_cos, list_ts)
    coverage_on, l3    = get_coverage (sCancer, 'ON', list_sGenes_cos, list_on)
    print('Score', coverage_score, l1)
    print('TS', coverage_ts, l2)
    print('ON', coverage_on, l3)

    print('Gwas')
    coverage_score, l1  = get_coverage(sCancer, 'Score', list_sGenes_gwas, list_score)
    coverage_ts, l2     = get_coverage(sCancer, 'TS', list_sGenes_gwas, list_ts)
    coverage_on, l3     = get_coverage(sCancer, 'ON', list_sGenes_gwas, list_on)
    print('Score', coverage_score, l1)
    print('TS', coverage_ts, l2)
    print('ON', coverage_on, l3)

    sys.exit()


    print('TCGA')
    coverage_score, l1  = get_coverage(sCancer, 'Score', list_sGenes_tcga, list_score)
    coverage_ts, l2     = get_coverage(sCancer, 'TS',    list_sGenes_tcga, list_ts)
    coverage_on, l3     = get_coverage(sCancer, 'ON',    list_sGenes_tcga, list_on)
    print('Score', coverage_score, l1)
    print('TS', coverage_ts, l2)
    print('ON', coverage_on, l3)

    print('Cosmic + TCGA')
    coverage_score, l1  = get_coverage(sCancer, 'Score', list_all_three, list_score)
    coverage_ts, l2     = get_coverage(sCancer, 'TS',    list_all_three, list_ts)
    coverage_on, l3     = get_coverage(sCancer, 'ON',    list_all_three, list_on)
    print('Score', coverage_score, l1)
    print('TS', coverage_ts, l2)
    print('ON', coverage_on, l3)


    print(l1)
    print(l2)
    print(l3)
#def END: calc_coverage
## endregion


## region Validate Databases

def valid_db (sDatabase):

    sWorkDir         = '%s/09_validate_dbs'                        % sBASE_DIR
    sCosmicFile      = '%s/cosmic.v74.txt'                         % sREF_DIR
    sClinvarFile     = '%s/clinvar.vcf'                            % sREF_DIR

    ## Sample Lists
    sSampleFile_HN   = '%s/normal_IDlist.txt'                      % sREF_DIR

    ## VCFs from Lung or Gastric with the Healthy Norm VCF
    sVCFFile_HN      = '%s/normal.vcf'                             % sREF_DIR

    # Read sample files
    list_sSamples_HN = get_sample_list (sSampleFile_HN)

    list_cData       = []
    if sDatabase == 'cosmic':
        # Load COSMIC data
        list_cData   = cCos_parse_cosmic_consensus   (sCosmicFile)

    elif sDatabase == 'clinvar':
        # Load ClinVar data
        list_cData   = cClin_parse_clinvar_file      (sClinvarFile)
    #if END:

    # Parse VCF file
    list_cVCF        = cVCF_parse_vcf_files         (len(list_sSamples_HN), sVCFFile_HN)

    print('Normal', len(list_cVCF))
    print(sDatabase, len(list_cData))

    list_sFinal      = get_intersection             (sDatabase, list_cVCF, list_cData)

    output_for_boxplot                              (sWorkDir, sDatabase, list_sFinal)
    output_for_histogram                            (sWorkDir, sDatabase, list_sFinal)

    print_done_msg('Validate Database')
#def END: valid_db


def get_intersection (sDatabase, list_cVCF, list_cData):

    print_start_msg('Getting Intersection %s' % sDatabase)

    if sDatabase == 'cosmic':

        # By Position
        dict_cData       = {int(cCos.sPos) : cCos  for cCos in list_cData if '-' not in cCos.sPos}
        list_sKeys_data  = list(dict_cData.keys())

    elif sDatabase == 'clinvar':

        # By Position
        dict_cData       = {cClin.nPos: cClin for cClin in list_cData}
        list_sKeys_data  = list(dict_cData.keys())

    #if END:
    dict_cVCF        = {cVCF.nPos:    cVCF for cVCF  in list_cVCF}
    list_sKeys_vcf   = list(dict_cVCF.keys())

    print('%s-PosCount' % sDatabase, len(list_sKeys_data))
    print('%s-PosCount' % 'Normal',  len(list_sKeys_vcf))

    list_sIntersection = list(set(list_sKeys_data) & set(list_sKeys_vcf))
    print('IntersectionCount', len(list_sIntersection))

    list_sOutput     = []
    for sKey in list_sIntersection:

        cData          = dict_cData[sKey]
        cVCF           = dict_cVCF[sKey]
        fMAF, nPatCnt  = calculate_maf(cVCF.list_sSamples)

        list_sOutput.append([sKey, cData.sDelete, fMAF])
    #loop END: sKey

    #V-S Check: Empty List
    if not list_sOutput:
        sys.exit('Empty List : get_intersection : list_sOutput Size= %d' % len(list_sOutput))

    return list_sOutput
#def END: get_intersection


def calculate_maf (list_sSamples):

    dict_sGTScore = {'./.':-1, '0/0':0, '0/1':1, '1/0':1, '1/1':2} # Based on VCFtools criteria
    list_sGT_gene = []

    for sSampleInfo in list_sSamples:
        # sSampleInfo format: GT:AD:DP:GQ:PL
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
        sGT, sAD, sDP, sGQ, sPL = sSampleInfo.split(':')

        list_sGT_gene.append(dict_sGTScore[sGT])
    #loop END: sSampleInfo

    fTotal  = sum([        2 for sGT_score in list_sGT_gene if sGT_score != -1])
    fAllele = sum([sGT_score for sGT_score in list_sGT_gene if sGT_score > 0])
    nPatCnt = sum([        1 for sGT_score in list_sGT_gene if sGT_score > 0])



    return fAllele/fTotal, nPatCnt
#def END: calculate_maf


def output_for_boxplot (sWorkDir, sDatabase, list_sFinal):

    sOutDir ='%s/output' % (sWorkDir)
    os.makedirs(sOutDir, exist_ok=True)

    dict_sOutput = {}
    for nPos, sPathogen, fMAF in list_sFinal:
        if sPathogen not in dict_sOutput:
            dict_sOutput[sPathogen] = []
        dict_sOutput[sPathogen].append(fMAF)
    #loop END:

    for sPathogen in dict_sOutput:
        print(sPathogen, len(dict_sOutput[sPathogen]))
    #loop END: sPathogen

    OutFile = open('%s/output_for_boxplot_%s.txt' % (sOutDir, sDatabase), 'w')

    for sKey in dict_sOutput:
        for fMAF in dict_sOutput[sKey]:

            sOut = '%s\t%0.5f\n' % (sKey, fMAF)
            #print(sOut[:-1])
            OutFile.write(sOut)
        #loop END: fMAF
    #loop END: sKey
    OutFile.close()

    print_done_msg('Output for Boxplot', sOutDir)
#def END: output_for_boxplot


def output_for_histogram (sWorkDir, sDatabase, list_sFinal):

    nBins = 100

    dict_sOutput = {}
    for nPos, sPathogen, fMAF in list_sFinal:
        if sPathogen not in dict_sOutput:
            dict_sOutput[sPathogen] = []
        dict_sOutput[sPathogen].append(fMAF)
    #loop END:

    sOutDir ='%s/output' % (sWorkDir)
    os.makedirs(sOutDir, exist_ok=True)

    for sPathogen in dict_sOutput:

        list_fMAF = dict_sOutput[sPathogen]

        ## MAF Histogram ##
        dict_fMAF = {i:[] for i in range(nBins)}

        fBinRange = max(list_fMAF) / nBins

        for nGroup, lCnt in itertools.groupby(list_fMAF, key=lambda n: n//fBinRange):
            dict_fMAF[nGroup] = dict_fMAF[nGroup] + list(lCnt)
        #loop END: nGroup, lCnt

        nSumCheck = 0
        OutFile     = open('%s/output_for_histogram_%s_%s.txt' % (sOutDir, sDatabase, sPathogen), 'w')
        sHeader     = 'Bin-%s\t%s\t%s\t%s\n' % (fBinRange, 'Count', 'Mean', 'SEM')
        OutFile.write(sHeader)

        for nGroup in dict_fMAF:

            nSumCheck += len(dict_fMAF[nGroup])

            if sum(dict_fMAF[nGroup]) != 0:
                sOut       = '%s\t%s\t%0.5f\t%0.5f\n' % \
                             (nGroup, len(dict_fMAF[nGroup]), mean(dict_fMAF[nGroup]),
                              stderror(dict_fMAF[nGroup]) if stderror(dict_fMAF[nGroup])  != 'NULL' else 0)
            else:
                sOut       = '%s\t0\t0\t0\n' % nGroup
            OutFile.write(sOut)
        #loop END: nGroup
        OutFile.close()
    #loop END: sPathogen
    print_done_msg('Output for Histogram %s Bins' % nBins, sOutDir)
#def END: output_for_histogram

## endregion


def main():
    pass
#def END: main()



if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__