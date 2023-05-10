#!/usr/bin/python
import os,sys
home = os.getcwd()
fn_input = sys.argv[1]
dic_type={}
dic_nm = {}
dic_nm_maf = {}
dic_nm_nn_maf={}
dic_nm_pn_maf={}
dic_nm_sc={}
dic_nm_mut={}
fn_out = file(sys.argv[2]+'.csv','wrt')
fn_out_anl=file(sys.argv[2]+'.analysis.txt','wrt')

nSift = 0
nPP2 = 0
nMA = 0
nGerp = 0
nCheck = 0

for i in file(fn_input):
    i = i.strip()
    line  = i.split('",')
    li_quat = i.split('"')
    li_com=i.split(',')
    if i.startswith('Ch'): #assign type such as exonic, intronic, etc
        s_type = None
    else:
        s_type=li_quat[1]
        s_geneSym=li_quat[3]

    if i.startswith('Chr,'): #print first row
        fn_out.write('%s\n'%i)
    else:

        try: #make s_type dictionary
            dic_type[s_type]+=1
        except:
            dic_type[s_type]=1

#chr1,866437,866437,C,T,"exonic","SAMD11",".","synonymous SNV","SAMD11:NM_152486:exon4:c.C273T:p.D91D","rs139076934",.,.,.,.,.,.,.,"0.002024     299.09  16"
#Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,snp137NonFlagged,SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationAssessor_score,MutationAssessor_pred,gerp++gt2,Otherinfo
            
        i_st_nm=i.find('NM_')# pick up the NM id
        if i_st_nm > 1:
            s_nm=i[i_st_nm:].split(':')[0]
        else:

            s_nm ='None'

        s_nm_sc='[%s\t%s\t%s\t%s]'%(li_com[-8],li_com[-6],li_com[-3],li_com[-2]) # compile deleterious score, SIFT=li_com[-8], POLYPHEN-2=li_com[-6], MutationAssessor=li_com[-3], GERP++=li_com[-2]
        nCheck += 1

        # CHECKING FOR NM by string length
        if len(s_nm) > 2: # check existing NMid, 2 is heuristic number of NM-id length cutoff
            if not li_com[-8]=='.': #attach sift score
                #print li_com
                f_sift=float(li_com[-8]) #sift score
            else:
                f_sift = 1 # non-deleterious



            if not li_com[-6]=='.': #attach polyphen2 score
                f_pp2 = float(li_com[-6])
            else:
                f_pp2=0 # non-deleterious



            if not li_com[-2]=='.': #attach gerp score
                f_gerp = float(li_com[-2].replace('"',''))
            else:
                f_gerp = 0 # non-deleterious



            if not li_com[-3]=='.': #attatch mutation assessor
                s_ma = li_com[-3].replace('"','')
            else:
                s_ma = 'A' # It means non deleterious




            temp_li = li_com[-1].replace('"','').split() # read MAF information
            f_maf = float(temp_li[0])

            if len(temp_li) > 3: # 3 is heuristic cutoff. if temp_li leght is larger than 3, it has MAF information
                f_norm_maf = float(li_com[-1].replace('"','').split()[3])
            else:
                continue

            temp_criteria   = 0 #the number of satisfying exonic deleterious score
            temp_criteria_2 = 0 #the number of satisfying non-exonic deleterious score

            ##options
            if f_maf > 2* f_norm_maf and f_maf > 0.01: # MAF cut-off
            #if f_maf > 0:
            #if f_maf:
                if f_sift < 0.05: #deleteriou SIFT cutoff
                    temp_criteria += 1
                    nSift +=1
                if f_pp2 > 0.5: #deleteriou PolyPhen2 cutoff
                    temp_criteria += 1
                    nPP2 +=1
                if s_ma =='H' or s_ma =='M': #deleterious MA cutoff
                    temp_criteria +=1
                    nMA += 1
                if f_gerp > 5: #deleterious GERP++ cutoff
                    temp_criteria_2 += 1
                    temp_criteria   += 1
                    nGerp +=1

            if s_type=='exonic': # exonic region use 4 deleterious score
                #MAF information attach step (similar with oncogene candidates detector)
                if temp_criteria > 1: 
                    fn_out.write('%s\n'%i)                  
                    try:
                        dic_type[s_type]+=1
                    except:
                        dic_type[s_type]=1
                    try:
                        dic_nm[s_nm]+=1
                    except:
                        dic_nm[s_nm]=1
                    try:
                        dic_nm_maf[s_nm]+=",%s-%s"%(f_maf,f_norm_maf)
                    except:
                        dic_nm_maf[s_nm]='%s-%s'%(f_maf,f_norm_maf)

                    try:
                        dic_nm_nn_maf[s_nm]+= f_norm_maf
                    except:
                        dic_nm_nn_maf[s_nm]= f_norm_maf

                    try:
                        dic_nm_pn_maf[s_nm]+= f_maf
                    except:
                        dic_nm_pn_maf[s_nm]= f_maf
                    try:
                        dic_nm_sc[s_nm]+=',%s-%s'%(s_nm_sc,temp_criteria)
                    except:                       
                        dic_nm_sc[s_nm]='%s-%s'%(s_nm_sc,temp_criteria)
                    
                    try:
                        dic_nm_mut[s_nm]+=',%s'%(s_type)
                    except:                       
                        dic_nm_mut[s_nm]='%s'%(s_type)

            else: #non-exonic region use only GERP++ score
                #MAF information attach step (similar with oncogene candidates detector)

                if temp_criteria_2 > 0:
                    fn_out.write('%s\n'%i)                  
                    try:
                        dic_type[s_type]+=1
                    except:
                        dic_type[s_type]=1
                    try:
                        dic_nm[s_nm]+=1
                    except:
                        dic_nm[s_nm]=1
                    try:
                        dic_nm_maf[s_nm]+=",%s-%s"%(f_maf,f_norm_maf)
                    except:
                        dic_nm_maf[s_nm]='%s-%s'%(f_maf,f_norm_maf)
                    try:
                        dic_nm_nn_maf[s_nm]+= f_norm_maf
                    except:
                        dic_nm_nn_maf[s_nm]= f_norm_maf
                    try:
                        dic_nm_pn_maf[s_nm]+= f_maf
                    except:
                        dic_nm_pn_maf[s_nm]= f_maf
                    try:
                        dic_nm_sc[s_nm]+=',%s-%s'%(s_nm_sc,temp_criteria)
                    except:                       
                        dic_nm_sc[s_nm]='%s-%s'%(s_nm_sc,temp_criteria)
                    try:
                        dic_nm_mut[s_nm]+=',%s'%(s_type)
                    except:                       
                        dic_nm_mut[s_nm]='%s'%(s_type)


print nSift
print nPP2
print nMA
print nGerp
print nCheck

fn_out_anl.write('###dic_nm\n#nm_id\tratio_MAF\tMAF_PN_Sum\tMAF_NN_sum\tpos\teach_score\t\t\teach_freq\teach_mutation_type\n')
li_nm=sorted(dic_nm.keys())
for j in xrange(len(li_nm)):
#    if dic_nm_pn_maf[li_nm[j]] > dic_nm_nn_maf[li_nm[j]]: 
#        if dic_nm_pn_maf[li_nm[j]] > 0.1 and dic_nm[li_nm[j]] >= 4:
#        if dic_nm_pn_maf[li_nm[j]] > 0.05 and dic_nm[li_nm[j]] >= 4:
    try:
        f_ratio_maf = dic_nm_pn_maf[li_nm[j]]/dic_nm_nn_maf[li_nm[j]]
    except:
        f_ratio_maf = dic_nm_pn_maf[li_nm[j]]/0.0001
    if not li_nm[j]=='None':
        fn_out_anl.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(li_nm[j],f_ratio_maf,dic_nm_pn_maf[li_nm[j]],dic_nm_nn_maf[li_nm[j]],dic_nm[li_nm[j]],dic_nm_sc[li_nm[j]],dic_nm_maf[li_nm[j]],dic_nm_mut[li_nm[j]]))
