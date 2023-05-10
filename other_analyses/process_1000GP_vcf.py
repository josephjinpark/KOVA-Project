#!/extdata6/Jinman/opt/python3/bin/python3



import os, sys, time

bTestRun    = False

sTIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_') )

sAnalysis      = '1000GP'
# Tools
sJAVA       = '/home/lab/bin/jre1.8.0_91/bin/java'
sCRAMTOOLS  = '/home/lab/bin/cramtools/cramtools-3.0.jar'   # for cram to bam conversion
sSAMTOOLS   = '/home/lab/bin/samtools-1.2/samtools'          # for sort and index bam file
sPICARD     = '/home/lab/bin/picard-tools-2.3.0/picard.jar'     # for bam to fastq conversion (1000GP aligned to hg38)
sGATK       = '/home/lab/bin/GenomeAnalysisTK-2.4-7/GenomeAnalysisTK.jar'

# Ref
sREF19      = '/extdata6/Jinman/reference_genome/hg19/hg19.fa'
sREF38      = '/extdata6/Jinman/reference_genome/hg38/hg38_wAlt.fa'
sREF_100GP  = '/extdata6/Jinman/reference_genome/hg19/hs37d5.fa'

# Directories
sBase_Dir   = '/extdata6/Jinman/00_K1000E_DB/11_vcftools/1000GP_vcf'
sVCF_Dir    = '%s/VCF'                              % sBase_Dir
sID_File    = '%s/%s.txt'                           % (sBase_Dir, sAnalysis)
sTmp_Dir    = '%s/tmp'                              % sBase_Dir

# for Qsub
sSymbol     = 'combine_1000GP' # qsub name symbol
sQueue      = 'optiplex.q'     # SGE queue name

os.makedirs(sTmp_Dir, exist_ok=True)


def get_file_ids (sWorkDir, sInFile):
    if not os.path.isfile(sInFile):
        list_AllFiles = os.listdir(sWorkDir)
        return [sFile for sFile in list_AllFiles if sFile.endswith('.gz')]
    else:
        return [sReadLine.strip('\n') for sReadLine in open(sInFile)]
#def END: get_file_ids



def main ():
    lists_sChrIDs   = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']

    #list_sFiles     = get_file_ids(sVCF_Dir, sID_File)

    sJobName        = 'Jinman.GATK.CombineVariants.%s'  % sAnalysis
    sLog_Dir        = '%s/log_%s'                       % (sBase_Dir, sTIME_STAMP)
    os.makedirs(sLog_Dir, exist_ok=True)

    sLogFile        = '%s/%s.combine.log'               % (sLog_Dir,  sAnalysis)
    sOutFile        = '%s/%s.combined.vcf'              % (sBase_Dir, sAnalysis)



    sScript = 'java -Xmx24g -Djava.io.tmpdir=%s -jar %s ' % (sTmp_Dir, sGATK)
    sScript += '-T CombineVariants '

    for sChrID in lists_sChrIDs:
        sInVCF = '%s/%s.%s.vcf.gz'      % (sVCF_Dir, sAnalysis, sChrID)
        sScript += '--variant:%s %s '   % (sChrID, sInVCF)

    sScript += '-R %s ' % sREF_100GP
    sScript += '-o %s ' % sOutFile
    sScript += '-genotypeMergeOptions PRIORITIZE '
    sScript += '-priority %s ' % (','.join(lists_sChrIDs))

    sScript = 'echo "%s" | qsub -V -q %s -N %s -cwd -j y -o %s' \
              % (sScript, sQueue, sJobName, sLogFile)

    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: main

main()

