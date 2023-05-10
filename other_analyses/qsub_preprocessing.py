#!/extdata6/Jinman/opt/python3/bin/python3


import os, sys

sBASE_DIR       = '/extdata6/Jinman/00_K1000E_DB'
sFILEID_DIR     = '%s/00_raw_fastq' % sBASE_DIR
sREF_DIR        = '/extdata6/Jinman/reference_genome/hg19'

## Chr Segment Size
nCHR_WINDOW     = 10000000 #1000000 # chromosome running size

## Constants
fAlpha          = '0.0'
lN              = 'None'  # constant for K1000E
bVCF            = True
bPLP            = True
list_sChrIDs    = ['chrM'] + ['chr'+str(x) for x in range(1,23)] + ['chrX','chrY']

def main():

    ## Queue Options
    sQueue          = 'optiplex.q'
    bTestRun        = False 
    sLabName        = 'stomach_normal'

    list_sFileIDs   = [sReadLine.strip('\n')for sReadLine in open('%s/%s.txt' % (sFILEID_DIR, sLabName), 'r')]

    ## Step 1: Run preprocessing script
    qsub_preprocessing(sLabName, list_sFileIDs, sQueue, bTestRun)

    ## Step 2: Merge window bam files into one
    #qsub_merge_bam    (sLabName, list_sFileIDs, sQueue, bTestRun)

    ## Step 3: Sort and index
    #qsub_sort_index   (sLabName, list_sFileIDs, sQueue, bTestRun)
#def END: main


def qsub_preprocessing (sLabName, list_sFileIDs, sQueue, bTestRun):

    sJobName        = 'Jinman.SV_Preproc'
    sProgram        = 'Only_Preprocessing.py'
    sInDir          = '%s/06_reindex_samtools/out'       % sBASE_DIR
    sOutDir         = '%s/07_preprocess/out_preproc'     % sBASE_DIR
    sLogDir         = '%s/07_preprocess/log'             % sBASE_DIR
    sTmpDir         = '%s/07_preprocess/tmp'             % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)
    for sFileID in list_sFileIDs[50:]:
        
        if os.path.isfile('%s/07_preprocess/out/%s.bam.bai' % (sBASE_DIR, sFileID)): continue
        sInFile       = '%s/%s.sort2.bam' % (sInDir, sFileID)

        nJobNo        = 0

        for sChrID in list_sChrIDs:
            lViewRange = lSetViewRange(sChrID)

            for iViewRangeIdx, sViewRange in enumerate(lViewRange):

                nJobNo  += 1
                sOutFile = '%s/%s.%s.%s.bam'     % (sOutDir, sFileID, sChrID, sViewRange)
                sVCFFile = '%s/%s.%s.%s.vcf '    % (sOutDir, sFileID, sChrID, iViewRangeIdx)

                sScript  = '/home/lab/bin/python2.7 %s %s %s %s %s %s %s %s %s %s %s' \
                           % (sProgram, sChrID, sInFile, lN, sOutDir, sViewRange, sVCFFile, bVCF, bPLP, fAlpha, sOutFile)
                sScript  = 'echo "%s" | qsub -cwd -q %s -j y -o %s -N %s.%s.%s.%s' \
                           % (sScript, sQueue, sLogDir, sJobName, sLabName, sChrID, sFileID)

                if bTestRun: print(sScript)
                else:
                    os.system(sScript)
            #loop END: iViewRangeIdx, sViewRange
        #loop END: sChrID
    #loop END: sFileID
#def END: qsub_preprocessing


def qsub_merge_bam (sLabName, list_sFileIDs, sQueue, bTestRun):

    sJobName        = 'Jinman.SV_MergeBam'
    sInDir          = '%s/07_preprocess/out_preproc'    % sBASE_DIR
    sOutDir         = '%s/07_preprocess/out_merge'      % sBASE_DIR
    sLogDir         = '%s/07_preprocess/log'            % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sFileID in list_sFileIDs:

        sInFile  = '%s/%s.*.bam'  % (sInDir,  sFileID)
        sOutFile = '%s/%s.bam'    % (sOutDir, sFileID)

        sScript  = '/home/lab/bin/samtools-0.1.19/samtools merge %s %s' % (sOutFile, sInFile)
        sScript  = 'echo "%s" | qsub -cwd -q %s -j y -o %s -N %s.%s.%s' \
                           % (sScript, sQueue, sLogDir, sJobName, sLabName, sFileID)

        if bTestRun: print(sScript)
        else:
            os.system(sScript)
    #loop END: sFileID
#def END: qsub_merge_bam


def qsub_sort_index (sLabName, list_sFileIDs, sQueue, bTestRun):

    sJobName        = 'Jinman.SV_SortIndex'
    sInDir          = '%s/07_preprocess/out_merge'    % sBASE_DIR
    sOutDir         = '%s/07_preprocess/out_sort'     % sBASE_DIR
    sLogDir         = '%s/07_preprocess/log'          % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    for sFileID in list_sFileIDs:

        sInFile  = '%s/%s.bam'  % (sInDir,  sFileID)
        sOutFile = '%s/%s'      % (sOutDir, sFileID)

        sScript  = '/home/lab/bin/samtools-0.1.19/samtools sort -@ 1 %s %s;' % (sInFile, sOutFile)
        sScript += '/home/lab/bin/samtools-0.1.19/samtools index %s.bam;'     % sOutFile
        sScript  = 'echo "%s" | qsub -cwd -q %s -j y -o %s -N %s.%s.%s'\
                    % (sScript, sQueue, sLogDir, sJobName, sLabName, sFileID)

        if bTestRun: print(sScript)
        else:
            os.system(sScript)
    #loop END: sFileID
#def END: qsub_merge_bam

## region Util Functions

def get_chrom_sizes ():
    sRefFile      = '%s/hg19.fa.fai' % sREF_DIR
    dict_sChrSize = {}

    #V-S Check: File Path
    if not os.path.isfile(sRefFile):
        sys.exit('File Not Found: %s' % sRefFile)

    for sReadLine in open(sRefFile, 'r'):
        list_sCol = sReadLine.split('\t')
        sChrID    = list_sCol[0]
        nChrSize  = int(list_sCol[1])

        if sChrID not in dict_sChrSize:
            dict_sChrSize[sChrID] = 0
        dict_sChrSize[sChrID] = nChrSize
    #loop END: sReadLine
    return dict_sChrSize
#def END: get_chrom_sizes


def lSetViewRange(sChrID):
    nChrSize     =  get_chrom_sizes()[sChrID]
    list_nWindow = []
    nRange       = int(nChrSize/nCHR_WINDOW) + 1

    for i in range(nRange):

        x = i * nCHR_WINDOW + 1
        y = (i+1) * nCHR_WINDOW

        if y > nChrSize: y = nChrSize

        list_nWindow.append('%s-%s' % (x, y))
    #loop END: i
    return list_nWindow
#end lSetViewRange

## endregion

main()
