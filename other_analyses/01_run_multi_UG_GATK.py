#!/extdata6/Jinman/opt/python3/bin/python3

import os, sys, time

sTIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_') )

nWINSIZE    = 1000000
nBUFFER     = 1000
sBASE_DIR   = '/extdata6/Kyuhong/04_merged_bam/'
sJAVA       = '/home/lab/bin/jre1.6.0_45/bin/java'
sJAVA2      = '/home/lab/bin/jre1.8.0_91/bin/java'
sGATK       = '/home/lab/bin/GenomeAnalysisTK-2.4-7/GenomeAnalysisTK.jar'
sGATK2      = '/home/lab/bin/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar'
sPICARD     = '/home/lab/bin/picard-tools-2.3.0/picard.jar'
sVCFSORT    = '/home/lab/bin/vcftools_0.1.12b/bin/vcf-sort'
sIGVTOOLS   = '/home/lab/bin/IGVTools/igvtools'

sREF        = '/extdata6/Jinman/reference_genome/hg19/hg19.fa'
sREFINDEX   = '/extdata6/Jinman/reference_genome/hg19/hg19.fa.fai'
sREFSEQDICT = '/extdata6/Jinman/reference_genome/hg19/hg19.dict'

def get_chrom_sizes ():

    dict_sChrSize = {}

    for sReadLine in open(sREFINDEX, 'r'):
        list_sCol = sReadLine.split('\t')
        sChrID    = list_sCol[0].replace('chr', '')
        nChrSize  = int(list_sCol[1])

        if sChrID not in dict_sChrSize:
            dict_sChrSize[sChrID] = 0
        dict_sChrSize[sChrID] = nChrSize
    #loop END: sReadLine

    return dict_sChrSize

def lSetViewRange(sChrID):
    nChrSize     = get_chrom_sizes()[sChrID[3:]]
    list_nWindow = []
    nRange       = int(nChrSize/nWINSIZE) + 1

    for i in range(nRange):
        x = (i) * nWINSIZE + 1
        y = (i+1) * nWINSIZE

        if y > nChrSize: y = nChrSize

        if i == 0 or i == nRange-1:
            list_nWindow.append('%s-%s' % (int(x), int(y)))
        else:
            list_nWindow.append('%s-%s' % (int(x)-nBUFFER, int(y)+nBUFFER))

        #list_nWindow.append('%s-%s' % (int(x), int(y)))
    #loop END: i
    return list_nWindow

def main():

    sData     = 'K1062'
    sIDFile   = '%s/%s_bam_list.txt'                              % (sBASE_DIR, sData)
    sTmpDir   = '%s/03_Unified_Genotyper/tmp'                     % sBASE_DIR
    sOutDir   = '%s/03_Unified_Genotyper/out_%s'                  % (sBASE_DIR, sData) # output files
    sOutChr   = '%s/chrom2'                                       % (sOutDir) # combined output file per chr
    sOutChr2  = '%s/chrom_sort2'                                  % (sOutDir) # combined output file per chr
    sOutUnion = '%s/03_Unified_Genotyper/out_final_%s_re'            % (sBASE_DIR, sData) # combined output chr files (FINAL)
    sLogDir   = '%s/03_Unified_Genotyper/log_%s'                  % (sBASE_DIR, sTIME_STAMP)  # log file directory

    sKnownDir = '/extdata6/Jinman/reference_genome/SGI_collaboration'
    sDBSNP    = '%s/dbsnp_137.hg19.vcf'                                  % sKnownDir
    sINDEL1   = '%s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf' % sKnownDir
    sINDEL2   = '%s/1000G_phase1.indels.hg19.sites.vcf'                  % sKnownDir

    sQueue    = 'optiplex.q'
    bTestRun  = False


    os.makedirs(sOutDir,   exist_ok=True)
    os.makedirs(sOutChr,   exist_ok=True)
    os.makedirs(sOutChr2,  exist_ok=True)
    os.makedirs(sOutUnion, exist_ok=True)
    os.makedirs(sLogDir,   exist_ok=True)
    os.makedirs(sTmpDir,   exist_ok=True)

    list_sChrIDs =['chrM'] + ['chr%s' % str(x) for x in range(1, 23)] + ['chrX', 'chrY']

    for sChrID in list_sChrIDs:

        list_nWindow  = lSetViewRange(sChrID)
        list_sOutFile = []

        for i, sViewRange in enumerate(list_nWindow):  # ['1-1000000', '1000001-2000000', ...]

            sInterval = '%s.%s' % (i, sViewRange)
            sLogFile  = '%s/K1000E.preprocessed_read.%s.%s.raw.log' % (sLogDir, sChrID, sInterval)
            sOutFile  = '%s/K1000E.preprocessed_read.%s.%s.raw.vcf' % (sOutDir, sChrID, sInterval)
            list_sOutFile.append(sOutFile)

            sCmd      = '%s -Xmx8g -Djava.io.tmpdir=%s -jar %s ' % (sJAVA, sTmpDir, sGATK)
            sCmd     += '-T UnifiedGenotyper '

            for sFileID in [sLine.strip('\n') for sLine in open(sIDFile)]:
                sCmd += '-I %s.bam ' % sFileID

            sCmd     += '-R %s '      % sREF
            sCmd     += '--dbsnp %s ' % sDBSNP
            sCmd     += '-glm BOTH '
            sCmd     += '-stand_call_conf 30.0 '
            sCmd     += '-stand_emit_conf 10.0 '
            sCmd     += '-o %s '      % sOutFile
            sCmd     += '-L %s:%s '   % (sChrID, sViewRange)

            sCmd      = 'echo "%s" | qsub -V -q %s -N Jinman.CALL.%s.%s.%s -cwd -j y -o %s' \
                        % (sCmd, sQueue, sData, sChrID, sViewRange, sLogFile)

            #if bTestRun: print(sCmd)
            #else:        os.system(sCmd)
        #loop END:i, sViewRange

        sLogFile  = '%s/K1000E.preprocessed_read.%s.raw.merge.log' % (sLogDir, sChrID)
        sOutFile2 = '%s/K1000E.preprocessed_read.%s.raw.vcf'       % (sOutChr, sChrID)
        sOutFile3 = '%s/K1000E.preprocessed_read.%s.raw.vcf'       % (sOutChr2, sChrID)

        sCmd = 'java -cp %s org.broadinstitute.sting.tools.CatVariants ' % sGATK2
        sCmd += '-R %s '        % sREF

        for sOutFile in list_sOutFile:
            sCmd += '-V %s '    % sOutFile

        sCmd += '-out %s '      % sOutFile2
        sCmd += '-assumeSorted ;'

        ## TEST SortVcf - Picard
        #sCmd += '%s -Xmx80g -XX:-UseGCOverheadLimit '  % sJAVA2
        #sCmd += '-jar %s SortVcf '                     % sPICARD
        #sCmd += 'I=%s '                                % sOutFile2
        #sCmd += 'O=%s '                                % sOutFile3
        #sCmd += 'SEQUENCE_DICTIONARY=%s ;'             % sREFSEQDICT

        ## SortVcf - VcfTools
        sCmd = 'date;cat %s | %s -t %s> %s ;date;'      % (sOutFile2, sVCFSORT, sTmpDir, sOutFile3)


        sCmd = 'echo "%s" | qsub -V -q %s -hold_jid Jinman.CALL.%s.%s.* -N Jinman.MERGE.RE.%s.%s -cwd -j y -o %s' \
               % (sCmd, sQueue, sData, sChrID, sData, sChrID, sLogFile)

        #if bTestRun: print(sCmd)
        #else:        os.system(sCmd)
    #loop END: sChrID

    sLogFile  = '%s/K1000E.preprocessed_read.union.log'    % sLogDir
    sOutFile3 = '%s/K1000E.preprocessed_read.combined.vcf' % sOutUnion

    sCmd = 'java -Xmx24g -Djava.io.tmpdir=%s -jar %s ' % (sTmpDir, sGATK)
    sCmd += '-T CombineVariants '

    for sChrID in list_sChrIDs:
        sInVCF = '%s/K1000E.preprocessed_read.%s.raw.vcf' % (sOutChr2, sChrID)
        sCmd  += '--variant:%s %s ' % (sChrID, sInVCF)

    sCmd += '-R %s ' % sREF
    sCmd += '-o %s ' % sOutFile3
    sCmd += '-genotypeMergeOptions PRIORITIZE '
    sCmd += '-priority %s ' % (','.join(list_sChrIDs))


    sCmd = 'echo "%s" | qsub -V -q %s -hold_jid Jinman.MERGE.RE.%s.* -N Jinman.UNION.%s.RE -cwd -j y -o %s' \
           % (sCmd, sQueue, sData, sData,sLogFile)

    if bTestRun: print(sCmd)
    else:        os.system(sCmd)
#def END: main()


main()