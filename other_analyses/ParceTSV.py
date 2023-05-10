#!/extdata6/Jinman/opt/python3/bin/python3

# For the Purpose of comparing TCGA data received from SGI and summary file downloaded from TCGA site
# From SGI :            LUAD (394)            and STAD (332) files
# TCGA as of 05/18/16 : LUAD (462 NB, 228 NT) and STAD (467 NB, 95 NT)  NB: Blood NT: Tissue


from collections import Counter
import sys, os, subprocess

sBASE_DIR   = '/extdata6/Jinman/00_K1000E_DB'
sTSV_DIR    = '%s/cancer_tsvs'  % sBASE_DIR
sFILEID_DIR = '%s/00_raw_fastq' % sBASE_DIR

## region Class Summary
class Summary:
    def __init__(self):
        self.sStudy              = ''
        self.sBarcode            = ''
        self.sDisease            = ''
        self.sDisease_name       = ''
        self.sSample_type        = ''
        self.sSample_type_name   = ''
        self.sAnalyte_typ        = ''
        self.sLibrary_type       = ''
        self.sCenter             = ''
        self.sCenter_name        = ''
        self.sPlatform           = ''
        self.sPlatform_name      = ''
        self.sAssembly           = ''
        self.sFilename           = ''
        self.sFiles_size         = ''
        self.sChecksum           = ''
        self.sAnalysis_id        = ''
        self.sAliquot_id         = ''
        self.sPatient_id         = ''
        self.sSample_id          = ''
        self.sTss_id             = ''
        self.sSample_accession   = ''
        self.sPublished          = ''
        self.sUploaded           = ''
        self.sModified           = ''
        self.sState              = ''
        self.sFlag               = ''
    #def END: __init__


def Summary_parse_tsv(sInFile):

    InFile             = open(sInFile, 'r')
    list_cSum          = []
    for sReadLine in InFile:

        cSum = Summary()

        if sReadLine.startswith('study'): continue

        list_sColumn = sReadLine.strip('\n').split('\t')

        cSum.sStudy             = list_sColumn[0]
        cSum.sBarcode           = list_sColumn[1]
        cSum.sDisease           = list_sColumn[2]
        cSum.sDisease_name      = list_sColumn[3]
        cSum.sSample_type       = list_sColumn[4]
        cSum.sSample_type_name  = list_sColumn[5]
        cSum.sAnalyte_typ       = list_sColumn[6]
        cSum.sLibrary_type      = list_sColumn[7]
        cSum.sCenter            = list_sColumn[8]
        cSum.sCenter_name       = list_sColumn[9]
        cSum.sPlatform          = list_sColumn[10]
        cSum.sPlatform_name     = list_sColumn[11]
        cSum.sAssembly          = list_sColumn[12]
        cSum.sFilename          = list_sColumn[13]
        cSum.sFiles_size        = list_sColumn[14]
        cSum.sChecksum          = list_sColumn[15]
        cSum.sAnalysis_id       = list_sColumn[16]
        cSum.sAliquot_id        = list_sColumn[17]
        cSum.sPatient_id        = list_sColumn[18]
        cSum.sSample_id         = list_sColumn[19]
        cSum.sTss_id            = list_sColumn[20]
        cSum.sSample_accession  = list_sColumn[21]
        cSum.sPublished         = list_sColumn[22]
        cSum.sUploaded          = list_sColumn[23]
        cSum.sModified          = list_sColumn[24]
        cSum.sState             = list_sColumn[25]

        list_cSum.append(cSum)
    #loop END: sReadLine
    InFile.close()

    #V-S Check
    if not list_cSum:
        sys.exit('Invalid List : Summary_parse_tsv : list_cSum size= %d' % len(list_cSum))

    return list_cSum

    #def END: parse_tsv
#class END: Summary
## endregion Class Summary

def load_fileid_txt (sInFile): return [sReadLine.strip('\n') for sReadLine in open(sInFile)]

def checksum (sCancerType, list_sFiles, list_cSum):

    # MD5 Checksum of random 5 files
    sBam_Dir        = '%s/TCGA/%s_normal' % (sFILEID_DIR, sCancerType)
    dict_sCheckSum  = {cSum.sFilename[:-4] : cSum.sChecksum for cSum in list_cSum}

    for sCommonFile in list_sFiles[:5]:
        sBamFile    = '%s/%s.bam' % (sBam_Dir, sCommonFile)

        sSTout      = subprocess.Popen('md5sum %s' % sBamFile, stdout=subprocess.PIPE, shell=True).stdout

        sCheckSum   = str(list(sSTout)[0], 'UTF-8').strip('\n').split(' ')[0]

        print(sCheckSum, dict_sCheckSum[sCommonFile])
    #loop END: sCommonFile
#def END: checksum



def main():

    list_sCancers   = ['LUAD', 'STAD']

    for sCancerType in list_sCancers:

        # Determine how many NBs and NTs from SGI by comparing it with summary (.tsv) file from TCGA
        sIDFile          = '%s/%s_normal.txt' % (sFILEID_DIR, sCancerType)

        list_SGI_files   = load_fileid_txt (sIDFile)


        dict_cSum        = {sNormType: [] for sNormType in ['NB', 'NT']}
        for sNormType in dict_cSum:

            sTSVFile             = '%s/%s_%s.tsv' % (sTSV_DIR, sCancerType, sNormType)
            dict_cSum[sNormType] = Summary_parse_tsv(sTSVFile)
        #loop END: sNormType

        list_sCommonNB  = list(set(list_SGI_files) & set([cSum.sFilename[:-4] for cSum in dict_cSum['NB']]))
        list_sCommonNT  = list(set(list_SGI_files) & set([cSum.sFilename[:-4] for cSum in dict_cSum['NT']]))
        list_sUpdateNB  = list(set([cSum.sFilename[:-4] for cSum in dict_cSum['NB']]) - set(list_SGI_files))

        dict_sBarCodeNB = {cSum.sFilename[:-4]: cSum.sAnalysis_id for cSum in dict_cSum['NB']}

        OutFile = open('%s/%s_forDL.txt' % (sBASE_DIR, sCancerType), 'w')
        for sFileName in list_sUpdateNB:
            sOut = '%s\n' % (dict_sBarCodeNB[sFileName])
            OutFile.write(sOut)
        #loop END: sFile
        OutFile.close()

        print(sCancerType)
        print('SGI Files', len(list_SGI_files))
        print('NB Files',  len(list_sCommonNB))
        print('NT Files',  len(list_sCommonNT))

        #checksum(sCancerType, list_sCommonNB, dict_cSum['NB'])
        #checksum(sCancerType, list_sCommonNT, dict_cSum['NT'])

    #loop END: sCancerType


main()
