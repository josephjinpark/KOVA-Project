#!/usr/bin/python
import os,sys
home = os.getcwd()
fn_input = sys.argv[1]
dic_geneSymbol={}
dic_type={}
#dic_mut_type={}
dic_nm = {}
dic_nm_maf = {}
dic_nm_nn_maf={}
dic_nm_pn_maf={}
dic_nm_sc={}
dic_nm_mut={}
dic_crd={}
fn_out = file(sys.argv[2]+'.csv','wrt')
fn_out_anl=file(sys.argv[2]+'.analysis.txt','wrt')

nCheck1 = 0
nCheck2 = 0
nCheck3 = 0
dict_sType = {}

for i in file(fn_input):
    i = i.strip()
    line  = i.split('",')
    li_quat = i.split('"')
    li_com=i.split(',')
    s_crd='%s_%s'%(li_com[0],li_com[1])
    if i.startswith('Ch'): #assign type such as exonic, intronic, etc
        s_type = None
    else:
        s_type=li_quat[1]
        s_geneSym=li_quat[3]

    if i.startswith('Chr,'): #print first row
        fn_out.write('%s\n'%i)
    else:
        try:
            dic_geneSymbol[s_geneSym]+=1 #make genesymbol dictionary
        except:
            dic_geneSymbol[s_geneSym]=1
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
            s_nm='None'

        s_nm_sc='[%s\t%s\t%s\t%s]'%(li_com[-8],li_com[-6],li_com[-3],li_com[-2])
        
        nCheck1 += 1
        if s_type not in dict_sType:
            dict_sType[s_type] = 0
        dict_sType[s_type] += 1

        if s_nm == 'None': continue
        if s_type=='exonic':
            
            nCheck2 += 1

            f_norm_maf = float(li_com[-1].replace('"','').split()[3])
            f_maf=float(li_com[-1].replace('"','').split()[0])
            if f_maf:
                nCheck3 += 1

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
                    dic_nm_sc[s_nm]+=',%s'%s_nm_sc
                except:
                    dic_nm_sc[s_nm]='%s'%s_nm_sc

                try:
                    dic_nm_pn_maf[s_nm]+= f_maf
                except:
                    dic_nm_pn_maf[s_nm]= f_maf
                try:
                    dic_nm_mut[s_nm]+=',%s'%(s_type)
                except:
                    dic_nm_mut[s_nm]='%s'%(s_type)
                try:
                    dic_crd[s_nm]+=',%s'%(s_crd)
                except:
                    dic_crd[s_nm]='%s'%(s_crd)



#li_mt= sorted(dic_type.keys())
#for i in xrange(len(li_mt)):
#    fn_out_anl.write('%s\t%s\n'%(li_mt[i],dic_type[li_mt[i]]))

print nCheck1
print nCheck2
print nCheck3
print len(dic_crd)





fn_out_anl.write('###dic_nm\n#nm_id\tratio_MAF\tdelta_MAF\t\t#pos\teach_coordinates\t\t\teach_freq\teach_mutation_type\n')
li_nm=sorted(dic_nm.keys())
for j in xrange(len(li_nm)):# j is just id

    if li_nm[j] == 'NM_001080396':
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
              (li_nm[j],s_ratio_maf,s_delta_maf,dic_nm[li_nm[j]],dic_crd[li_nm[j]],dic_nm_maf[li_nm[j]],dic_nm_mut[li_nm[j]])


    li_maf = dic_nm_maf[li_nm[j]].split(',')
    s_len_maf = len(li_maf)
    for k in xrange(s_len_maf): #the number on variants in the gene
        f_maf_group, f_maf_norm = li_maf[k].split('-')
        f_delta_maf =0
        f_g_maf =0
        f_n_maf =0

        if float(f_maf_norm) ==0:
            f_maf_norm = 0.0000001
        f_g_maf+= (float(f_maf_group))
        f_n_maf+= float(f_maf_norm)
    if f_n_maf==0:
        f_n_maf =0.00000001
    s_delta_maf = '%2.7f'%(f_g_maf-f_n_maf)
    f_ratio_maf =(f_g_maf/f_n_maf)
    s_ratio_maf = '%2.7f'%f_ratio_maf

    if int(dic_nm[li_nm[j]]) > 5 and li_nm[j]!='None':
        fn_out_anl.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(li_nm[j],s_ratio_maf,s_delta_maf,dic_nm[li_nm[j]],dic_crd[li_nm[j]],dic_nm_maf[li_nm[j]],dic_nm_mut[li_nm[j]]))
