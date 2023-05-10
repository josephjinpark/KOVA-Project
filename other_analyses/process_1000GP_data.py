#!/extdata6/Jinman/opt/python3/bin/python3



import os, sys, time

bTestRun    = False

sTIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_') )

sGroup      = 'CDX'
# Tools
sJAVA       = '/home/lab/bin/jre1.8.0_91/bin/java'
sCRAMTOOLS  = '/home/lab/bin/cramtools/cramtools-3.0.jar'   # for cram to bam conversion
sSAMTOOLS   = '/home/lab/bin/samtools-1.2/samtools'          # for sort and index bam file
sPICARD     = '/home/lab/bin/picard-tools-2.3.0/picard.jar'     # for bam to fastq conversion (1000GP aligned to hg38)

# Ref
sREF19      = '/extdata6/Jinman/reference_genome/hg19/hg19.fa'
sREF38      = '/extdata6/Jinman/reference_genome/hg38/hg38_wAlt.fa'

# Directories
sBase_Dir   = '/extdata4/Jinman/project/K1000E'
sID_File    = '%s/00_fileids/%s.txt'                % (sBase_Dir, sGroup)
sBAM_Dir    = '%s/01_bamfiles/1000GP_bam/%s'        % (sBase_Dir, sGroup)
sFASTQ_Dir  = '%s/01_bamfiles/1000GP_bam/%s/fastq'  % (sBase_Dir, sGroup)
sTmpDir     = '%s/01_bamfiles/tmp'                  % sBase_Dir
sLog_Dir    = '%s/log_%s'                           % (sBase_Dir, sTIME_STAMP)

# for Qsub
sSymbol     = 'process_1000GP' # qsub name symbol
sQueue      = 'optiplex.q'       # SGE queue name

os.makedirs(sBAM_Dir,   exist_ok=True)
os.makedirs(sFASTQ_Dir, exist_ok=True)
os.makedirs(sTmpDir,   exist_ok=True)
os.makedirs(sLog_Dir,   exist_ok=True)


def get_file_ids (sBAM_Dir, sInFile):

    if not os.path.isfile(sInFile):
        list_AllFiles = os.listdir(sBAM_Dir)
        return [sFile[:-5] for sFile in list_AllFiles if sFile.endswith('.cram')]

    else:
        return [sReadLine.strip('\n') for sReadLine in open(sInFile)]
#def END: get_file_ids



def main ():

    list_sFiles = get_file_ids (sBAM_Dir, sGroup)

    for sFile in list_sFiles[:1]:
        sIn_cram    = '%s.cram'        % sFile
        sOut_bam    = '%s.bam'         % sFile
        sOut_sorted = '%s.sorted'      % sFile
        sFastq_1    = '%s.1.fastq'     % sFile
        sFastq_2    = '%s.2.fastq'     % sFile
        sLogFile    = '%s/%s.log'      % (sLog_Dir, sFile)

        sScript1    = '%s view '                     % sSAMTOOLS
        #sScript1   += '--input-fmt-option ignore_md5=1 '
        sScript1   += '-T %s '                       % sREF38
        sScript1   += '-b -o %s/%s '                 % (sBAM_Dir, sOut_bam)
        sScript1   += ' %s/%s '                      % (sBAM_Dir, sIn_cram)

        #sScript1   += '--skip-md5-check '
        #sScript1   += '--ignore-md5-mismatch;date '

        sScript2    = '%s sort %s/%s '               % (sSAMTOOLS, sBAM_Dir, sOut_bam) # In
        sScript2   += '%s/%s '                       % (sBAM_Dir, sOut_sorted)         # Out

        sScript3    = '%s index %s/%s.bam '          % (sSAMTOOLS, sBAM_Dir, sOut_sorted)

        sScript4    = '%s -jar %s SamToFastq '       % (sJAVA, sPICARD)
        sScript4   += 'I=%s/%s.bam '                 % (sBAM_Dir, sOut_sorted)
        sScript4   += 'FASTQ=%s/%s '                 % (sFASTQ_Dir, sFastq_1)
        sScript4   += 'SECOND_END_FASTQ=%s/%s '      % (sFASTQ_Dir, sFastq_2)
        sScript4   += 'INCLUDE_NON_PF_READS=true '

        sScript     = '%s;%s;%s;%s;' % (sScript1, sScript2, sScript3, sScript4)
        sScript     = 'echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.%s.%s.%s'\
                      % (sScript, sLogFile, sQueue, sSymbol, sGroup, sFile.split('.')[0])

        if bTestRun == True: print(sScript)
        else:                os.system(sScript)
    #loop END: sFile
#def END: main

main()

