#!/extdata6/Jinman/opt/python3/bin/python3

import os, sys, time

sTIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_'))


## region Globals

sBASE_DIR       = '/extdata6/Jinman/00_K1000E_DB'
sBAM_DIR        = '/extdata6/Kyuhong/04_merged_bam/02_sorted/out'
sFILEID_DIR     = '%s/00_fileids'   % sBASE_DIR
sKIT_DIR        = '%s/ref'          % sBASE_DIR
sBWA            = '/home/lab/bin/bwa-0.7.5a/bwa'
sSAMTOOLS       = '/home/lab/bin/samtools-0.1.19/samtools'
sSAMTOOLS13     = '/home/lab/bin/samtools-1.3.1/samtools'

sJAVA           = '/home/lab/bin/jre1.6.0_45/bin/java'
sJAVA2          = '/home/lab/bin/jre1.8.0_91/bin/java'

sGATK           = '/home/lab/bin/GenomeAnalysisTK-2.4-7/GenomeAnalysisTK.jar'
sGATK2          = '/home/lab/bin/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar'

sMARKDUPS       = '/home/lab/bin/picard-tools-1.93/MarkDuplicates.jar'
sADDOREPLACEHG  = '/home/lab/bin/picard-tools-1.93/AddOrReplaceReadGroups.jar'
sMERGESAM       = '/home/lab/bin/picard-tools-1.93/MergeSamFiles.jar'
sPICARD         = '/home/lab/bin/picard-tools-2.3.0/picard.jar'

sVCFSORT        = '/home/lab/bin/vcftools_0.1.12b/bin/vcf-sort'

sSICKLE         = '/home/lab/bin/sickle-master/sickle'


# HG19
sREF            = '/extdata6/Jinman/reference_genome/hg19/hg19.fa'
sREFINDEX       = '/extdata6/Jinman/reference_genome/hg19/hg19.fa.fai'
sREFSEQDICT     = '/extdata6/Jinman/reference_genome/hg19/hg19.dict'

# Known Sites DBs
sKNOWNSITES     = '/extdata6/Jinman/reference_genome/SGI_collaboration'
sKNOWN_OMNI     = '%s/1000G_omni2.5.hg19.sites.vcf'                         % sKNOWNSITES
sKNOWN_HAPMAP   = '%s/hapmap_3.3.hg19.sites.vcf'                            % sKNOWNSITES
sKNOWN_SNPS     = '%s/dbsnp_137.hg19.vcf'                                   % sKNOWNSITES
sKNOWN_INDELS1  = '%s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf ' % sKNOWNSITES
sKNOWN_INDELS2  = '%s/1000G_phase1.indels.hg19.sites.vcf'                   % sKNOWNSITES

list_sCHRIDs    = ['chrM'] + ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']

nMAX_QSUB_CHAR  = 131071


## endregion

## region Util Functions

def get_file_ids (sAnalysis):

    sInFile         = '%s/%s.txt' % (sFILEID_DIR, sAnalysis)
    InFile          = open(sInFile, 'r')
    list_sFileIDs   = []

    for sReadLine in InFile:
        list_sFileIDs.append (sReadLine.strip('\n'))
    #loop END:

    return list_sFileIDs
#def END: get_file_ids


def make_log_dir (sJobName):

    sLogDir = '%s/log/%s/%s' % (sBASE_DIR, sJobName, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)
    return sLogDir
#def END: make_log_dir


def get_chrom_sizes():

    dict_sChrSize = {}

    for sReadLine in open(sREFINDEX, 'r'):
        list_sCol   = sReadLine.split('\t')
        sChrID      = list_sCol[0].replace('chr', '')
        nChrSize    = int(list_sCol[1])

        if sChrID not in dict_sChrSize:
            dict_sChrSize[sChrID] = 0
        dict_sChrSize[sChrID] = nChrSize
    # loop END: sReadLine

    return dict_sChrSize
#def END: get_chrom_sizes


def set_chr_range(sChrID, nWinSize):

    nChrSize     =  get_chrom_sizes()[sChrID[3:]]
    list_nWindow = []
    nRange       = int(nChrSize/nWinSize) + 1

    for i in range(nRange):

        x = i * nWinSize + 1
        y = (i+1) * nWinSize

        if y > nChrSize: y = nChrSize

        list_nWindow.append('%s-%s' % (int(x), int(y)))
    #loop END: i
    return list_nWindow
#def END: set_chr_range


def get_fullbamlist (sAnalysis):

    if sAnalysis != 'K1056':
        sBamDir         = '%s/07_preprocess/%s/final' % (sBASE_DIR, sAnalysis)
        list_sFileIDs   = get_file_ids (sAnalysis)
        list_sOutput    = []
        for sFile in list_sFileIDs:
            sBamFile = '%s/%s.bam' % (sBamDir, sFile)

            if not os.path.isfile(sBamFile): continue

            list_sOutput.append(sBamFile)
        #loop END: sFile

        return list_sOutput

    else:
        return [sReadLine.strip('\n') for sReadLine in open('%s/%s_bam_list.txt' % (sFILEID_DIR, sAnalysis))]
#def END: get_fullbamlist


def get_capturekit_list (sAnalysis):

    sInFile      = '%s/%s_capturekit.txt' % (sFILEID_DIR, sAnalysis)
    InFile       = open(sInFile, 'r')
    list_sFileID = [sReadLine.strip('\n').split('\t') for sReadLine in InFile]
    InFile.close()
    return dict(list_sFileID)
#def END: get_capturekit_list

## endregion





def mpileup_analysis (sQueue, bTestRun):

    list_sAnalysis  = ['Park', 'Jang', 'Lee', 'Choi']
    #sAnalysis       = 'Lee'

    sOutDir         = '%s/12_for_ewha/pileup' % sBASE_DIR

    for sAnalysis in list_sAnalysis:

        sJobName        = 'Jinman.PileupAnal.%s'     % sAnalysis
        sLogDir         = make_log_dir        (sJobName)

        if sAnalysis == 'Lee':   #ERCSB used two different capture kits
            dict_sCapture = get_capturekit_list(sAnalysis)
        else:
            dict_sCapture = {'Baek':'Agilent_SureSelect_50Mb',
                             'Park':'Agilent_SureSelect_V5',
                             'Jang':'Agilent_SureSelect_50Mb',
                             'Choi':'Roche_V2'}

        sInDir          = '%s/%s'            % (sBAM_DIR, sAnalysis)
        sTempDir        = '%s/%s'            % (sOutDir, sAnalysis)
        os.makedirs(sTempDir, exist_ok=True)

        list_sFileID    = get_file_ids(sAnalysis)
        list_sCompile   = []
        for sFileID in list_sFileID:
            sInBam      = '%s/%s.bam'        % (sInDir, sFileID)
            sCaptureKit = dict_sCapture[sFileID] if sAnalysis == 'Lee' else dict_sCapture[sAnalysis]
            sKitBedFile = '%s/%s.bed'        % (sKIT_DIR, sCaptureKit)
            sOutFile    = '%s/%s.pileup.txt' % (sTempDir, sFileID)
            sOutFile2   = '%s/%s.calc.txt'   % (sTempDir, sFileID) #file with total and average depth
            sLogFile    = '%s/%s.log'        % (sLogDir, sFileID)

            list_sCompile.append(sOutFile2)
            #VS-Check
            if not os.path.isfile(sInBam):
                sys.exit('%s **Not Found**' % sInBam)
            if not os.path.isfile(sKitBedFile):
                sys.exit('%s **Not Found**' % sKitBedFile)

            sScript     = 'samtools mpileup -l %s %s > %s;' % (sKitBedFile, sInBam, sOutFile)
            sScript    += 'echo \\\"%s\\\" '                % sFileID
            sScript    += '\` awk \'{ sum += \$4 } END { if (NR > 0) print sum, sum / NR }\' %s\` > %s' % (sOutFile, sOutFile2)

            sScript     = 'echo "date;%s;date" | qsub -V -q %s -N %s.%s -cwd -j y -o %s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile)

            if bTestRun: print(sScript)
            else:        os.system(sScript)
        #loop END: sFileID

        sCompileFile = '%s/%s.compiled.txt'     % (sOutDir, sAnalysis)
        sJobName     = 'Jinman.PileupCompile.%s' % sAnalysis
        sLogDir      = make_log_dir (sJobName)
        sLogFile     = '%s/%s.compile.log'      % (sLogDir, sAnalysis)
        sScript      = 'cat %s > %s'        % (' '.join(list_sCompile), sCompileFile)
        sScript      = 'echo "date;%s;date" | qsub -V -q %s -N %s -cwd -j y -o %s -hold_jid Jinman.DepthAnal.%s.*' \
                        % (sScript, sQueue, sJobName, sLogFile, sAnalysis)

        if bTestRun: print(len(sScript))
        else:
            os.system(sScript)
            if len(sScript) > nMAX_QSUB_CHAR:
                sys.exit('Script Length Must Be <=131071 | Current Length= %d' % (len(sScript)))
    #loop END: sAnalysis
#def END: depth_analysis


def depth_analysis (sQueue, bTestRun):

    list_sAnalysis  = ['Baek', 'Park', 'Jang', 'Lee', 'Choi']
    list_sAnalysis  = ['Baek']

    sOutDir         = '%s/12_for_ewha/depth' % sBASE_DIR

    for sAnalysis in list_sAnalysis:

        sJobName        = 'Jinman.DepthAnal.%s' % sAnalysis
        sLogDir         = make_log_dir (sJobName)

        if sAnalysis == 'Lee':   #ERCSB used two different capture kits
            dict_sCapture = get_capturekit_list(sAnalysis)
        else:
            dict_sCapture   = {'Baek':'Agilent_SureSelect_50Mb',
                               'Park':'Agilent_SureSelect_V5',
                               'Jang':'Agilent_SureSelect_50Mb',
                               'Choi':'Roche_V2'}

        sInDir          = '%s/%s'            % (sBAM_DIR, sAnalysis)
        sTempDir        = '%s/%s'            % (sOutDir, sAnalysis)
        os.makedirs(sTempDir, exist_ok=True)

        list_sFileID    = get_file_ids(sAnalysis)
        list_sCompile   = []
        list_sTarFiles  = []
        for sFileID in list_sFileID:
            sInBam      = '%s/%s.bam'        % (sInDir, sFileID)
            sCaptureKit = dict_sCapture[sFileID] if sAnalysis == 'Lee' else dict_sCapture[sAnalysis]
            sKitBedFile = '%s/%s.bed'        % (sKIT_DIR, sCaptureKit)
            sOutFile    = '%s/%s.depth.txt'  % (sTempDir, sFileID)
            sOutFile2   = '%s/%s.calc.txt'   % (sTempDir, sFileID) #file with total and average depth
            sOutFile3   = '%s/%s.chr19.txt'  % (sTempDir, sFileID) #just chr19, for sangmoon
            sLogFile    = '%s/%s.log'        % (sLogDir, sFileID)

            if os.path.isfile(sOutFile3):continue

            list_sCompile.append(sOutFile2)
            list_sTarFiles.append(sOutFile3)
            #VS-Check
            if not os.path.isfile(sInBam):
                sys.exit('%s **Not Found**' % sInBam)
            if not os.path.isfile(sKitBedFile):
                sys.exit('%s **Not Found**' % sKitBedFile)

            #Original
            sScript     = '%s depth -a -b %s %s > %s;' % (sSAMTOOLS13, sKitBedFile, sInBam, sOutFile)
            sScript    += 'echo \\\"%s\\\" '           % sFileID
            sScript    += '\` awk \'{ sum += \$3 } END { if (NR > 0) print sum, sum / NR }\' %s\` > %s;' % (sOutFile, sOutFile2)

            ## Temp for Sangmoon, chr19 only ##
            sScript    += 'egrep %s %s > %s' % ('chr19', sOutFile, sOutFile3)
            ###################################

            sScript = 'echo "date;%s;date" | qsub -V -q %s -N %s.%s -cwd -j y -o %s' \
                      % (sScript, sQueue, sJobName, sFileID, sLogFile)

            if bTestRun: print(sScript)
            else:        os.system(sScript)
        #loop END: sFileID
        sJobName     = 'Jinman.DepthCompile.%s' % sAnalysis
        sLogDir      = make_log_dir (sJobName)
        sLogFile     = '%s/%s.compile.log'      % (sLogDir, sAnalysis)

        #Original
        sCompileFile = '%s/%s.compiled.txt'     % (sOutDir, sAnalysis)
        sScript      = 'cat %s > %s;'        % (' '.join(list_sCompile), sCompileFile)
        #sScript      = 'echo "date;%s;date" | qsub -V -q %s -N %s -cwd -j y -o %s -hold_jid Jinman.DepthAnal.%s.*' \
        #                % (sScript, sQueue, sJobName, sLogFile, sAnalysis)

        ## Temporary for Sangmoon, chr19 only ##
        sTarFile     = '%s/%s.chr19.tar.gz' % (sOutDir, sAnalysis)
        sScript      += 'tar -zcf %s %s'     % (sTarFile, ' '.join(list_sTarFiles))
        sScript      = 'echo "date;%s;date" | qsub -V -q %s -N %s -cwd -j y -o %s -hold_jid Jinman.DepthAnal.%s.*' \
                  % (sScript, sQueue, sJobName, sLogFile, sAnalysis)
        ########################################
        if bTestRun: print(sScript)
        else:
            os.system(sScript)
            if len(sScript) > nMAX_QSUB_CHAR:
                sys.exit('Script Length Must Be <=131071 | Current Length= %d' % (len(sScript)))
    #loop END: sAnalysis
#def END: depth_analysis


def sickle (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.Sickle.%s'        % sAnalysis
    sLogDir     = make_log_dir (sJobName)

    sInDir      = '%s/00_raw_fastq/%s'      % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/00_raw_fastq/%s'      % (sBASE_DIR, sAnalysis)

    sScore      = 'sanger'  # quality scoring scheme (must check!)

    for sFileID in list_sFileIDs:

        sInFastq1  = '%s/%s.1.fastq'        % (sInDir, sFileID)
        sInFastq2  = '%s/%s.2.fastq'        % (sInDir, sFileID)
        sOutFastq1 = '%s/%s.1.qt.fastq'     % (sOutDir, sFileID)
        sOutFastq2 = '%s/%s.2.qt.fastq'     % (sOutDir, sFileID)
        sOutFastqS = '%s/%s.2.s.fastq'      % (sOutDir, sFileID)
        sLogFile   = '%s/%s.qt.log'         % (sLogDir, sFileID)

        sScript  = '%s pe -t %s '           % (sSICKLE, sScore)
        sScript += '-f %s -r %s '           % (sInFastq1, sInFastq2)
        sScript += '-o %s -p %s '           % (sOutFastq1, sOutFastq2)
        sScript += '-s %s '                 % (sOutFastqS)

        sScript  = 'echo "date;%s;date" | qsub -V -q %s -N %s.%s -cwd -j y -o %s' \
                    % (sScript, sQueue, sJobName, sFileID, sLogFile)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop END: sFileID
#def END: sickle


def bwa_align (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.BWA.Align.%s'     % sAnalysis
    sPrevJob    = 'Jinman.Sickle.%s'        % sAnalysis
    sLogDir     = make_log_dir (sJobName)

    sInDir      = '%s/00_raw_fastq/%s'      % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/01_bwa_aln/%s'        % (sBASE_DIR, sAnalysis)

    os.makedirs(sOutDir, exist_ok=True)

    nThreads    = 4

    for sFileID in list_sFileIDs:
        for i in [1,2]:
            sInFastq = '%s/%s.%s.qt.fastq'  % (sInDir,  sFileID, i)
            sOutSai  = '%s/%s.%s.sai'       % (sOutDir, sFileID, i)
            sLogFile = '%s/%s.%s.log'       % (sLogDir, sFileID, i)

            sScript  = '%s aln %s '         % (sBWA, sREF)
            sScript += '%s '                % sInFastq
            sScript += '-f %s '             % sOutSai
            sScript += '-t %s '             % nThreads

            sScript  = 'echo "%s" | qsub -V -q %s -N %s.%s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, i, sLogFile, sPrevJob, sFileID)

            #if os.path.isfile(sOutSai): continue

            if bTestRun: print(sScript)
            else:        os.system(sScript)
        #loop END: i
    #loop END: sFileID
#def END: bwa_align


def bwa_sampe (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.BWA.Sampe.%s'             % sAnalysis
    sPrevJob    = 'Jinman.BWA.Align.%s'             % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sFastqDir   = '%s/00_raw_fastq/%s'              % (sBASE_DIR, sAnalysis)
    sWorkDir    = '%s/01_bwa_aln/%s'                % (sBASE_DIR, sAnalysis)

    for sFileID in list_sFileIDs:
        sInSai_1    = '%s/%s.1.sai'                 % (sWorkDir,  sFileID)
        sInSai_2    = '%s/%s.2.sai'                 % (sWorkDir,  sFileID)
        sInFastq_1  = '%s/%s.1.qt.fastq'            % (sFastqDir, sFileID)
        sInFastq_2  = '%s/%s.2.qt.fastq'            % (sFastqDir, sFileID)
        sOutSam     = '%s/%s.sam'                   % (sWorkDir,  sFileID)
        sLogFile    = '%s/%s.log'                   % (sLogDir,   sFileID)

        sScript     = '%s sampe %s '                        % (sBWA, sREF)
        sScript    += '-r \'@RG\\tID:%s\\tSM:%s\\tPL:%s\' ' % (sFileID, sFileID, 'illumina')
        sScript    += '%s '                         % sInSai_1
        sScript    += '%s '                         % sInSai_2
        sScript    += '%s '                         % sInFastq_1
        sScript    += '%s '                         % sInFastq_2
        sScript    += '> %s '                       % sOutSam

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s.*' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
        91
    # loop END: sFileID
#def END: bwa_sampe


def samtools_sort_index (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.Samtools.SortIndex.%s'    % sAnalysis
    sPrevJob    = 'Jinman.BWA.Sampe.%s'             % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/01_bwa_aln/%s'                % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/02_samtools/%s'               % (sBASE_DIR, sAnalysis)

    os.makedirs(sOutDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInSam      = '%s/%s.sam'                   % (sInDir,  sFileID)
        sOutBam     = '%s/%s.bam'                   % (sOutDir, sFileID)
        sOutSorted  = '%s/%s.sorted'                % (sOutDir, sFileID)
        sLogFile    = '%s/%s.log'                   % (sLogDir, sFileID)

        sScript     = '%s view -bt %s '             % (sSAMTOOLS, sREFINDEX)
        sScript    += '%s '                         % sInSam
        sScript    += '-o %s ;'                     % sOutBam

        sScript    += '%s sort -T %s %s '           % (sSAMTOOLS, sLogDir, sOutBam)
        sScript    += '-o %s ;'                     % sOutSorted

        sScript    += '%s index '                   % sSAMTOOLS
        sScript    += '%s.bam ;'                    % sOutSorted

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop END: sFileID
#def END: samtools_sort_index


def headgroups (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.Picard.AddorReplaceHeads.%s'              % sAnalysis
    sPrevJob    = 'Jinman.Samtools.SortIndex.%s'                    % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/02_samtools/%s'                               % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/03_markdups/%s'                               % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/03_markdups/tmp'                              % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.sorted.bam'                            % (sInDir,  sFileID)
        sOutBam     = '%s/%s.ar.bam'                                % (sOutDir, sFileID)
        sLogFile    = '%s/%s.ar.log'                                % (sLogDir, sFileID)

        sScript     = 'java -Xmx8g -Djava.io.tmpdir=%s -jar %s '    % (sTmpDir, sADDOREPLACEHG)
        sScript    += 'I=%s '                                       % sInBam
        sScript    += 'O=%s '                                       % sOutBam
        sScript    += 'RGLB=TCGA '      # TCGA Mixture samples
        sScript    += 'RGPL=illumina '  # don't change
        sScript    += 'RGPU=TCGA '      # TCGA
        sScript    += 'RGSM=%s '                                    % sFileID
        sScript    += 'VALIDATION_STRINGENCY=LENIENT '
        sScript    += 'SORT_ORDER=coordinate '

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: markdups


def mergesam (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.Picard.MergeSam.%s'                       % sAnalysis
    sPrevJob    = 'Jinman.Picard.AddorReplaceHeads.%s'              % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/03_markdups/%s'                               % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/03_markdups/%s'                               % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/03_markdups/tmp'                              % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.ar.bam'                                % (sInDir,  sFileID)
        sOutBam     = '%s/%s.ms.bam'                                % (sOutDir, sFileID)
        sLogFile    = '%s/%s.merge.log'                             % (sLogDir, sFileID)

        sScript     = 'java -Xmx8g -Djava.io.tmpdir=%s -jar %s '    % (sTmpDir, sMERGESAM)
        sScript    += 'I=%s '                                       % sInBam
        sScript    += 'O=%s '                                       % sOutBam
        sScript    += 'VALIDATION_STRINGENCY=LENIENT '
        sScript    += 'ASSUME_SORTED=true '

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: markdups


def markdups (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.Picard.Markdups.%s'                       % sAnalysis
    sPrevJob    = 'Jinman.Picard.MergeSam.%s'                       % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/03_markdups/%s'                               % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/03_markdups/%s'                               % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/03_markdups/tmp'                              % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.ms.bam'                                % (sInDir,  sFileID)
        sOutBam     = '%s/%s.markdup.bam'                           % (sOutDir, sFileID)
        sMetFile    = '%s/%s.markdup.dup'                           % (sOutDir, sFileID)
        sLogFile    = '%s/%s.md.log'                                % (sLogDir, sFileID)

        sScript     = 'java -Xmx8g -Djava.io.tmpdir=%s -jar %s '    % (sTmpDir, sMARKDUPS)
        sScript    += 'I=%s '                                       % sInBam
        sScript    += 'O=%s '                                       % sOutBam
        sScript    += 'M=%s '                                       % sMetFile
        sScript    += 'VALIDATION_STRINGENCY=LENIENT '
        sScript    += 'ASSUME_SORTED=true '
        sScript    += 'REMOVE_DUPLICATES=true ;'

        sScript    += '%s index %s ;'                               % (sSAMTOOLS, sOutBam)

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: markdups


def gatk_rtc (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.GATK.RealignerTargetCreator.%s'       % sAnalysis
    sPrevJob    = 'Jinman.Picard.Markdups.%s'                   % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/03_markdups/%s'                           % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/04_realn_GATK/%s'                         % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/04_realn_GATK/tmp'                        % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.markdup.bam'                       % (sInDir,  sFileID)
        sOutTmp     = '%s/%s.realn.intervals'                   % (sOutDir, sFileID)
        sLogFile    = '%s/%s.rtc.log'                           % (sLogDir, sFileID)

        sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '  % (sJAVA, sTmpDir, sGATK)
        sScript    += '-T RealignerTargetCreator '
        sScript    += '-R %s '                                  % sREF
        sScript    += '-I %s '                                  % sInBam
        sScript    += '-o %s '                                  % sOutTmp
        sScript    += '-log %s '                                % sLogFile
        sScript    += '-known %s '                              % sKNOWN_INDELS1
        sScript    += '-known %s '                              % sKNOWN_INDELS2

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: gatk_rtc


def gatk_ir (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.GATK.IndelRealigner.%s'               % sAnalysis
    sPrevJob    = 'Jinman.GATK.RealignerTargetCreator.%s'       % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/03_markdups/%s'                           % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/04_realn_GATK/%s'                         % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/04_realn_GATK/tmp'                        % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.markdup.bam'                       % (sInDir,  sFileID)
        sInTmp      = '%s/%s.realn.intervals'                   % (sOutDir, sFileID)
        sOutBam     = '%s/%s.realn.bam'                         % (sOutDir, sFileID)
        sLogFile    = '%s/%s.ir.log'                            % (sLogDir, sFileID)

        sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '  % (sJAVA, sTmpDir, sGATK)
        sScript    += '-T IndelRealigner '
        sScript    += '-R %s '                                  % sREF
        sScript    += '-I %s '                                  % sInBam
        sScript    += '-o %s '                                  % sOutBam
        sScript    += '-log %s '                                % sLogFile
        sScript    += '-targetIntervals %s '                    % sInTmp

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: gatk_ir


def gatk_br (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.GATK.BaseRecalibrator.%s'             % sAnalysis
    sPrevJob    = 'Jinman.GATK.IndelRealigner.%s'               % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/04_realn_GATK/%s'                         % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/05_recal_GATK/%s'                         % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/05_recal_GATK/tmp'                        % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.realn.bam'                         % (sInDir,  sFileID)
        sOutGrp     = '%s/%s.recal.grp'                         % (sOutDir, sFileID)
        sLogFile    = '%s/%s.br.log'                            % (sLogDir, sFileID)

        sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '  % (sJAVA, sTmpDir, sGATK)
        sScript    += '-T BaseRecalibrator '
        sScript    += '-R %s '                                  % sREF
        sScript    += '-I %s '                                  % sInBam
        sScript    += '-o %s '                                  % sOutGrp
        sScript    += '-knownSites %s '                         % sKNOWN_SNPS
        sScript    += '-knownSites %s '                         % sKNOWN_INDELS1
        sScript    += '-knownSites %s '                         % sKNOWN_INDELS2

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def EMD: gatk_br


def gatk_pr (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.GATK.PrintReads.%s'                   % sAnalysis
    sPrevJob    = 'Jinman.GATK.BaseRecalibrator.%s'             % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/04_realn_GATK/%s'                         % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/05_recal_GATK/%s'                         % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/05_recal_GATK/tmp'                        % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.realn.bam'                         % (sInDir,  sFileID)
        sInGrp      = '%s/%s.recal.grp'                         % (sOutDir, sFileID)
        sOutGrp     = '%s/%s.recal2.grp'                        % (sOutDir, sFileID)
        sOutBam     = '%s/%s.recal.bam'                         % (sOutDir, sFileID)
        sLogFile    = '%s/%s.pr.log'                            % (sLogDir, sFileID)

        sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '  % (sJAVA, sTmpDir, sGATK)
        sScript    += '-T BaseRecalibrator '
        sScript    += '-R %s '                                  % sREF
        sScript    += '-I %s '                                  % sInBam
        sScript    += '-BQSR %s '                               % sInGrp
        sScript    += '-o %s '                                  % sOutGrp
        sScript    += '-knownSites %s '                         % sKNOWN_SNPS
        sScript    += '-knownSites %s '                         % sKNOWN_INDELS1
        sScript    += '-knownSites %s ;'                        % sKNOWN_INDELS2

        sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '  % (sJAVA, sTmpDir, sGATK)
        sScript    += '-T PrintReads '
        sScript    += '-R %s '                                  % sREF
        sScript    += '-I %s '                                  % sInBam
        sScript    += '-BQSR %s '                               % sInGrp
        sScript    += '-o %s ;'                                 % sOutBam

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: gatk_pr


def samtools_resort_reindex (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    sJobName    = 'Jinman.Samtools.ReSortIndex.%s'  % sAnalysis
    sPrevJob    = 'Jinman.GATK.PrintReads.%s'       % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/05_recal_GATK/%s'             % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/06_reindex_samtools/%s'       % (sBASE_DIR, sAnalysis)

    os.makedirs(sOutDir, exist_ok=True)

    for sFileID in list_sFileIDs:
        sInBam      = '%s/%s.recal.bam'             % (sInDir,  sFileID)
        sOutBam     = '%s/%s.resort'                % (sOutDir, sFileID)
        sLogFile    = '%s/%s.log'                   % (sLogDir, sFileID)

        sScript     = '%s sort %s '                 % (sSAMTOOLS, sInBam)
        sScript    += '%s ;'                        % sOutBam

        sScript    += '%s index %s.bam '            % (sSAMTOOLS, sOutBam)

        sScript     = 'echo "%s" | qsub -V -q %s -N %s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                        % (sScript, sQueue, sJobName, sFileID, sLogFile, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: samtools_resort_reindex


def sovar_preprocessing (sAnalysis, list_sFileIDs, sQueue, bTestRun):

    ## Chr Segment Size
    nWinSize     = 25000000  # 1000000 # chromosome running size

    ## SoVar Constants
    fAlpha          = '0.0'

    sProgram        = '%s/07_preprocess/Only_Preprocessing.py' % sBASE_DIR
    sTmpDir         = '%s/07_preprocess/%s/tmp'                % (sBASE_DIR, sAnalysis)
    os.makedirs(sTmpDir, exist_ok=True)

    for sFileID in list_sFileIDs[1:]:

        ## Part 1: Preprocess - Scatter
        sJobName        = 'Jinman.SoVar.Preprocessing.%s'   % sAnalysis
        sPrevJob        = 'Jinman.Samtools.ReSortIndex.%s'  % sAnalysis
        sLogDir         = make_log_dir(sJobName)

        sInDir          = '%s/06_reindex_samtools/%s'       % (sBASE_DIR, sAnalysis)
        sOutDir         = '%s/07_preprocess/%s/scatter'     % (sBASE_DIR, sAnalysis)
        os.makedirs(sOutDir, exist_ok=True)

        sInBam          = '%s/%s.resort.bam'                % (sInDir, sFileID)
        lN              = 'None'  # Not neccessary for K1000E
        bVCF, bPLP      = True, True
        nJobNo          = 0

        list_sOutFiles  = []

        for sChrID in list_sCHRIDs:

            list_nWin = set_chr_range(sChrID, nWinSize)

            for i, sViewRange in enumerate(list_nWin):

                nJobNo += 1

                sOutFile = '%s/%s.%s.%s.bam'  % (sOutDir, sFileID, sChrID, sViewRange)
                sLogFile = '%s/%s.%s.%s.log'  % (sLogDir, sFileID, sChrID, sViewRange)

                sVCFFile = '%s/%s.%s.%s.vcf ' % (sOutDir, sFileID, sChrID, i)

                sScript  = '/home/lab/bin/python2.7 %s %s %s %s %s %s %s %s %s %s %s' \
                          % (sProgram, sChrID, sInBam, lN, sOutDir, sViewRange, sVCFFile, bVCF, bPLP, fAlpha, sOutFile)

                sScript  = 'echo "%s" | qsub -V -q %s -N %s.%s.%s.%s -cwd -j y -o %s -hold_jid %s.%s' \
                          % (sScript, sQueue, sJobName, sFileID, sChrID, nJobNo, sLogFile, sPrevJob, sFileID)

                list_sOutFiles.append(sOutFile)

                if bTestRun: print(sScript)
                else:        os.system(sScript)
            # loop END: i, sViewRange
        # loop END: sChrID

        ## Part 2: Preprocess - Gather

        sJobName    = 'Jinman.SoVar.Merge.%s'           % sAnalysis
        sPrevJob    = 'Jinman.SoVar.Preprocessing.%s'   % sAnalysis
        sLogDir     = make_log_dir(sJobName)

        sOutDir2    = '%s/07_preprocess/%s/gather'      % (sBASE_DIR, sAnalysis)
        sOutDir3    = '%s/07_preprocess/%s/final'       % (sBASE_DIR, sAnalysis)

        os.makedirs(sOutDir2, exist_ok=True)
        os.makedirs(sOutDir3, exist_ok=True)

        sInFile     = ' '.join(list_sOutFiles)
        sOutFile2   = '%s/%s.bam'                       % (sOutDir2, sFileID)
        sOutFile3   = '%s/%s'                           % (sOutDir3, sFileID)

        sScript     = '%s merge %s %s;'                 % (sSAMTOOLS, sOutFile2, sInFile)
        sScript    += '%s sort %s %s;'                  % (sSAMTOOLS, sOutFile2, sOutFile3)
        sScript    += '%s index %s.bam;'                % (sSAMTOOLS, sOutFile3)

        sScript     = 'echo "%s" | qsub -cwd -q %s -j y -o %s -N %s.%s -hold_jid %s.%s.*' \
                        % (sScript, sQueue, sLogDir, sJobName, sFileID, sPrevJob, sFileID)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    # loop END: sFileID
#def END: sovar_preprocessing


def gatk_unifiedgenotyper (sAnalysis, list_sBam, sQueue, bTestRun):

    nWinSize    = 5000000
    sOutDir     = '%s/08_variantcall/%s'             % (sBASE_DIR, sAnalysis)
    sInterDir   = '%s/intervals2'                    % sOutDir  # output files
    sChromDir   = '%s/chrom2'                        % sOutDir  # combined vcf file per chr
    sChromDir2  = '%s/chrom_sort2'                   % sOutDir  # sorted vcf files per chr
    sFinalDir   = '%s/final'                         % sOutDir  # combined output chr files (FINAL)
    sTmpDir     = '%s/08_variantcall/tmp2'           % sBASE_DIR

    os.makedirs(sOutDir,    exist_ok=True)
    os.makedirs(sInterDir,  exist_ok=True)
    os.makedirs(sChromDir,  exist_ok=True)
    os.makedirs(sChromDir2, exist_ok=True)
    os.makedirs(sFinalDir,  exist_ok=True)
    os.makedirs(sTmpDir,    exist_ok=True)

    for sChrID in list_sCHRIDs:

        list_nWindow    = set_chr_range(sChrID, nWinSize)
        list_sOutFile   = []

        for i, sViewRange in enumerate(list_nWindow):

            sJobName    = 'Jinman.GATK.UG.CALL.%s'                   % sAnalysis
            sLogDir     = make_log_dir(sJobName)

            sInterval   = '%s.%s' % (i, sViewRange)
            sLogFile    = '%s/%s.preproc.%s.%s.raw.log'              % (sLogDir, sAnalysis, sChrID, sInterval)

            sOutFile    = '%s/%s.preproc.%s.%s.raw.vcf'              % (sInterDir, sAnalysis, sChrID, sInterval)
            list_sOutFile.append(sOutFile)

            sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '   % (sJAVA, sTmpDir, sGATK)
            sScript    += '-T UnifiedGenotyper '

            for sFileID in list_sBam:
                sScript += '-I %s.bam ' % sFileID

            sScript += '-R %s '                                      % sREF
            sScript += '--dbsnp %s '                                 % sKNOWN_SNPS
            sScript += '-glm BOTH '
            sScript += '-stand_call_conf 30.0 '
            sScript += '-stand_emit_conf 10.0 '
            sScript += '-o %s '                                      % sOutFile
            sScript += '-L %s:%s '                                   % (sChrID, sViewRange)

            sScript = 'echo "%s" | qsub -V -q %s -N %s.%s.%s -cwd -j y -o %s' \
                        % (sScript, sQueue, sJobName, sChrID, sViewRange, sLogFile)

            if bTestRun: print(sScript)
            #else:        os.system(sScript)
        #loop END:i, sViewRange

        sJobName    = 'Jinman.GATK.CatVariants.%s'     % sAnalysis
        sPrevJob    = 'Jinman.GATK.UG.CALL.%s'         % sAnalysis
        sLogDir     = make_log_dir(sJobName)

        sLogFile    = '%s/%s.preproc.%s.raw.log'       % (sLogDir,    sAnalysis, sChrID)
        sOutFile2   = '%s/%s.preproc.%s.raw.vcf'       % (sChromDir,  sAnalysis, sChrID)
        sOutFile3   = '%s/%s.preproc.%s.raw.vcf'       % (sChromDir2, sAnalysis, sChrID)

        sScript     = 'java -cp %s org.broadinstitute.sting.tools.CatVariants ' % sGATK2
        sScript    += '-R %s '                                                  % sREF

        for sOutFile in list_sOutFile:
            sScript += '-V %s '  % sOutFile

        sScript    += '-out %s ' % sOutFile2
        sScript    += '-assumeSorted ;'

        ## SortVcf - VcfTools
        sScript    += 'date;cat %s | %s -t %s> %s ;date;'  % (sOutFile2, sVCFSORT, sTmpDir, sOutFile3)

        sScript     = 'echo "%s" | qsub -V -q %s -hold_jid %s.%s.* -N %s.%s -cwd -j y -o %s' \
                        % (sScript, sQueue, sPrevJob, sChrID, sJobName, sChrID, sLogFile)

        if bTestRun: print(sScript)
        #else:        os.system(sScript)
    #loop END: sChrID

    sJobName     = 'Jinman.GATK.CombineVariants.%s'              % sAnalysis
    sPrevJob     = 'Jinman.GATK.CatVariants.%s'                  % sAnalysis
    sLogDir      = make_log_dir(sJobName)

    sLogFile     = '%s/%s.preproc.raw.log'                       % (sLogDir,   sAnalysis)
    sOutFile3    = '%s/%s.preproc.combined.vcf'                  % (sFinalDir, sAnalysis)

    sScript      = 'java -Xmx24g -Djava.io.tmpdir=%s -jar %s '   % (sTmpDir, sGATK)
    sScript     += '-T CombineVariants '

    for sChrID in list_sCHRIDs:
        sInVCF   = '%s/%s.preproc.%s.raw.vcf'                    % (sChromDir2, sAnalysis, sChrID)
        sScript += '--variant:%s %s '                            % (sChrID, sInVCF)

    sScript     += '-R %s '                                      % sREF
    sScript     += '-o %s '                                      % sOutFile3
    sScript     += '-genotypeMergeOptions PRIORITIZE '
    sScript     += '-priority %s ' % (','.join(list_sCHRIDs))

    sScript      = 'echo "%s" | qsub -V -q %s -hold_jid %s.* -N %s -cwd -j y -o %s' \
                    % (sScript, sQueue, sPrevJob, sJobName, sLogFile)

    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: gatk_unifiedgenotyper


def gatk_haplotypecaller (sAnalysis, list_sBam, sQueue, bTestRun):

    nWinSize    = 5000000
    sOutDir     = '%s/08_variantcall/%s_HC'         % (sBASE_DIR, sAnalysis)
    sInterDir   = '%s/intervals'                    % sOutDir  # output files
    sChromDir   = '%s/chrom'                        % sOutDir  # combined vcf file per chr
    sChromDir2  = '%s/chrom_sort'                   % sOutDir  # sorted vcf files per chr
    sFinalDir   = '%s/final'                        % sOutDir  # combined output chr files (FINAL)
    sTmpDir     = '%s/08_variantcall/tmp'           % sBASE_DIR

    os.makedirs(sOutDir,    exist_ok=True)
    os.makedirs(sInterDir,  exist_ok=True)
    os.makedirs(sChromDir,  exist_ok=True)
    os.makedirs(sChromDir2, exist_ok=True)
    os.makedirs(sFinalDir,  exist_ok=True)
    os.makedirs(sTmpDir,    exist_ok=True)

    for sChrID in list_sCHRIDs:

        list_nWindow    = set_chr_range(sChrID, nWinSize)
        list_sOutFile   = []

        for i, sViewRange in enumerate(list_nWindow):

            sJobName    = 'Jinman.GATK.HC.CALL.%s'                   % sAnalysis
            sLogDir     = make_log_dir(sJobName)

            sInterval   = '%s.%s' % (i, sViewRange)
            sLogFile    = '%s/%s.preproc.%s.%s.raw.log'              % (sLogDir, sAnalysis, sChrID, sInterval)

            sOutFile    = '%s/%s.preproc.%s.%s.raw.vcf'              % (sInterDir, sAnalysis, sChrID, sInterval)
            list_sOutFile.append(sOutFile)

            if os.path.isfile(sOutFile): continue

            sScript     = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s '   % (sJAVA, sTmpDir, sGATK)
            sScript    += '-T HaplotypeCaller '

            for sFileID in list_sBam:
                sScript += '-I %s.bam ' % sFileID

            sScript += '-R %s '                                      % sREF
            sScript += '--dbsnp %s '                                 % sKNOWN_SNPS
            sScript += '-stand_call_conf 30.0 '
            sScript += '-stand_emit_conf 10.0 '
            sScript += '-out_mode EMIT_ALL_SITES ' # Testing
            sScript += '-o %s '                                      % sOutFile
            sScript += '-L %s:%s '                                   % (sChrID, sViewRange)

            sScript = 'echo "%s" | qsub -V -q %s -N %s.%s.%s -cwd -j y -o %s' \
                        % (sScript, sQueue, sJobName, sChrID, sViewRange, sLogFile)

            if bTestRun: print(sScript); sys.exit()
            else:        os.system(sScript)
        #loop END:i, sViewRange

        sJobName    = 'Jinman.GATK.CatVariants.HC.%s'  % sAnalysis
        sPrevJob    = 'Jinman.GATK.HC.CALL.%s'         % sAnalysis
        sLogDir     = make_log_dir(sJobName)

        sLogFile    = '%s/%s.preproc.%s.raw.log'       % (sLogDir,    sAnalysis, sChrID)
        sOutFile2   = '%s/%s.preproc.%s.raw.vcf'       % (sChromDir,  sAnalysis, sChrID)
        sOutFile3   = '%s/%s.preproc.%s.raw.vcf'       % (sChromDir2, sAnalysis, sChrID)

        sScript     = 'java -cp %s org.broadinstitute.sting.tools.CatVariants ' % sGATK2
        sScript    += '-R %s '                                                  % sREF

        for sOutFile in list_sOutFile:
            sScript += '-V %s '  % sOutFile

        sScript    += '-out %s ' % sOutFile2
        sScript    += '-assumeSorted ;'

        ## SortVcf - VcfTools
        sScript    += 'date;cat %s | %s -t %s> %s ;date;'  % (sOutFile2, sVCFSORT, sTmpDir, sOutFile3)

        sScript     = 'echo "%s" | qsub -V -q %s -hold_jid %s.%s.* -N %s.%s -cwd -j y -o %s' \
                        % (sScript, sQueue, sPrevJob, sChrID, sJobName, sChrID, sLogFile)

        if bTestRun: print(sScript)
        else:        os.system(sScript)
    #loop END: sChrID

    sJobName     = 'Jinman.GATK.CombineVariants.HC.%s'           % sAnalysis
    sPrevJob     = 'Jinman.GATK.CatVariants.HC.%s'               % sAnalysis
    sLogDir      = make_log_dir(sJobName)

    sLogFile     = '%s/%s.preproc.raw.log'                       % (sLogDir,   sAnalysis)
    sOutFile3    = '%s/%s.preproc.combined.vcf'                  % (sFinalDir, sAnalysis)

    sScript      = 'java -Xmx24g -Djava.io.tmpdir=%s -jar %s '   % (sTmpDir, sGATK)
    sScript     += '-T CombineVariants '

    for sChrID in list_sCHRIDs:
        sInVCF   = '%s/%s.preproc.%s.raw.vcf'                    % (sChromDir2, sAnalysis, sChrID)
        sScript += '--variant:%s %s '                            % (sChrID, sInVCF)

    sScript     += '-R %s '                                      % sREF
    sScript     += '-o %s '                                      % sOutFile3
    sScript     += '-genotypeMergeOptions PRIORITIZE '
    sScript     += '-priority %s ' % (','.join(list_sCHRIDs))

    sScript      = 'echo "%s" | qsub -V -q %s -hold_jid %s.* -N %s -cwd -j y -o %s' \
                    % (sScript, sQueue, sPrevJob, sJobName, sLogFile)

    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: gatk_haplotypecaller



def gatk_vr (sAnalysis, sQueue, bTestRun):

    sJobName    = 'Jinman.GATK.VariantRecalibrator.%s'      % sAnalysis
    sPrevJob    = 'Jinman.GATK.CombineVariants.%s'          % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sInDir      = '%s/08_variantcall/%s/final'              % (sBASE_DIR, sAnalysis)
    sOutDir     = '%s/09_recal_variant/%s'                  % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/09_recal_variant/tmp'                 % sBASE_DIR

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sTmpDir, exist_ok=True)

    sInVCF      = '%s/%s.preproc.combined.vcf'              % (sInDir,  sAnalysis)
    sOutTable   = '%s/%s.preproc.combined.recal'            % (sOutDir, sAnalysis)
    sOutTranch  = '%s/%s.preproc.combined.tranches'         % (sOutDir, sAnalysis)
    sLogFile    = '%s/%s.recal.log'                         % (sLogDir, sAnalysis)

    sScript     = '%s -Xmx24g -Djava.io.tmpdir=%s -jar %s ' % (sJAVA, sTmpDir, sGATK)
    sScript    += '-T VariantRecalibrator '
    sScript    += '-input %s '                              % sInVCF
    sScript    += '-R %s '                                  % sREF
    sScript    += '-resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 %s '     % sKNOWN_SNPS
    sScript    += '-resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 %s '      % sKNOWN_OMNI
    sScript    += '-resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 %s '    % sKNOWN_HAPMAP
    sScript    += '-resource:Mills,known=false,training=true,truth=true,prior=12.0 %s '         % sKNOWN_INDELS1
    sScript    += '-resource:1000G,VCF,known=false,training=true,truth=false,prior=10.0 %s '    % sKNOWN_INDELS2
    sScript    += '-an QD '
    sScript    += '-an MQRankSum '
    sScript    += '-an ReadPosRankSum '
    sScript    += '-an FS '
    sScript    += '-an DP '
    sScript    += '-mode BOTH '
    sScript    += '--maxGaussians 4 '
    sScript    += '--percentBadVariants 0.05 '
    sScript    += '-recalFile %s '                          % sOutTable
    sScript    += '-tranchesFile %s '                       % sOutTranch

    sScript     = 'echo "%s" | qsub -V -q %s -N %s -cwd -j y -o %s -hold_jid %s' \
                    % (sScript, sQueue, sJobName, sLogFile, sPrevJob)

    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: gatk_vr


def gatk_ar (sAnalysis, sQueue, bTestRun):

    sJobName    = 'Jinman.GATK.ApplyRecalibration.%s'       % sAnalysis
    sPrevJob    = 'Jinman.GATK.VariantRecalibrator.%s'      % sAnalysis
    sLogDir     = make_log_dir(sJobName)

    sVCFDir     = '%s/08_variantcall/%s/final'              % (sBASE_DIR, sAnalysis)
    sInDir      = '%s/09_recal_variant/%s'                  % (sBASE_DIR, sAnalysis)
    sFinalDir   = '%s/09_recal_variant/%s/final'            % (sBASE_DIR, sAnalysis)
    sTmpDir     = '%s/09_recal_variant/tmp'                 % sBASE_DIR

    os.makedirs(sVCFDir,   exist_ok=True)
    os.makedirs(sFinalDir, exist_ok=True)
    os.makedirs(sTmpDir,   exist_ok=True)

    sInVCF      = '%s/%s.preproc.combined.vcf'              % (sVCFDir,   sAnalysis)
    sInTable    = '%s/%s.preproc.combined.recal'            % (sInDir,    sAnalysis)
    sInTranch   = '%s/%s.preproc.combined.tranches'         % (sInDir,    sAnalysis)
    sOutVCF     = '%s/%s.preproc.final.recal.vcf'           % (sFinalDir, sAnalysis)
    sLogFile    = '%s/%s.recal.log'                         % (sLogDir,   sAnalysis)

    sScript     = '%s -Xmx24g -Djava.io.tmpdir=%s -jar %s ' % (sJAVA, sTmpDir, sGATK)
    sScript    += '-T ApplyRecalibration '
    sScript    += '-input %s '                              % sInVCF
    sScript    += '-R %s '                                  % sREF
    sScript    += '--ts_filter_level 99.0 '
    sScript    += '-mode BOTH '
    sScript    += '-recalFile %s '                          % sInTable
    sScript    += '-tranchesFile %s '                       % sInTranch
    sScript    += '-o %s; '                                 % sOutVCF
    sScript    += '-tar zcf %s.tar.gz %s %s.idx; '          % (sOutVCF, sOutVCF, sOutVCF)
    sScript     = 'echo "%s" | qsub -V -q %s -N %s -cwd -j y -o %s -hold_jid %s' \
                    % (sScript, sQueue, sJobName, sLogFile, sPrevJob)

    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: gatk_ar


def main():

    #list_sAnalysis  = ['Baek', 'SGI', 'GSK_Kor', 'ERCSB', 'ERCSB_Gastric', 'ERSCB_Lung']
    sAnalysis       = 'Baek'
    list_sFileIDs   = get_file_ids (sAnalysis)

    #Queue Options
    sQueue   = 'optiplex.q'
    bTestRun = False

    ## For Dr. Seo (ERCSB) ##
    #mpileup_analysis (sQueue, bTestRun)
    #depth_analysis (sQueue, bTestRun)
    #########################

    ## Pre-step: Sickle
    #sickle               (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 1: BWA - Alignment
    #bwa_align               (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 2: BWA - Make SAMPE
    bwa_sampe               (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 3: Samtools - Sort and Index
    #samtools_sort_index     (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 4: PicardTools - Mark Duplicates
    headgroups              (sAnalysis, list_sFileIDs, sQueue, bTestRun)
    #mergesam                (sAnalysis, list_sFileIDs, sQueue, bTestRun)
    #markdups                (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 5: GATK - RealignerTargetCreator
    #gatk_rtc                (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 6: GATK - IndelRealigner
    #gatk_ir                 (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 7: GATK - BaseRecalibrator
    #gatk_br                 (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 8: GATK - PrintReads
    #gatk_pr                 (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 8: Samtools - Sort and Index
    #samtools_resort_reindex (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    ## Step 9: SoVar - Preprocessing
    #sovar_preprocessing     (sAnalysis, list_sFileIDs, sQueue, bTestRun)

    list_sAnalysis  = ['STAD', 'LUAD']
    sAnalysis       = 'K1056'
    #sAnalysis       = 'LUAD_blood_normal'
    #sAnalysis       = 'STAD_blood_normal'

    ## Step 10: GATK - UnifiedGenotyper
    list_sBams      = get_fullbamlist (sAnalysis)
    gatk_unifiedgenotyper   (sAnalysis, list_sBams, sQueue, bTestRun)
    gatk_haplotypecaller    (sAnalysis, list_sBams, sQueue, bTestRun)

    ## Step 11: GATK - VariantRecalibrator
    #gatk_vr                 (sAnalysis, sQueue, bTestRun)

    ## Step 12: GATK - ApplyRecalibration
    #gatk_ar                 (sAnalysis, sQueue, bTestRun)
#def END: main

if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__