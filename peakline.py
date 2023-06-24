import pandas as pd
import numpy as np
import sys
import math
import os

#raw_gwas file
df = pd.read_csv(sys.argv[1],sep = '\t')

df['p_wald'] = df['p_wald'].astype('float')

filter_df = df[df['p_wald'] < 1e-5]

gene_lst = filter_df.iloc[:,1].values.tolist()

p_value = filter_df.iloc[:,2].values.tolist()

#gene list for extract line infos
col_info = pd.read_csv('all_extract_info.csv',sep = ',',usecols =[0,1,2,3,4])

extract_num = [i-1 for i in gene_lst]

select_rows = col_info.iloc[extract_num]

add_lst = np.column_stack((p_value,np.zeros(len(p_value))))

merge_lst = np.hstack((add_lst,select_rows))

seri = pd.DataFrame(merge_lst)

seri.drop_duplicates(subset = [2],inplace = True)

seri.drop(seri.columns[1],axis = 1, inplace = True)

seri.iloc[:,2] = seri.iloc[:,2].astype(int)

#column names  p_value len_fusion_gene breakpoint1 breakpoint2 fusion_seq
seri.rename(columns = {0:'p_value',2:'fusion_id',3:'len_fusion_gene',4:'breakpoint1',5:'breakpoint2',6:'fusion_seq'},inplace = True)

#create fasta format file
def extract_orf_frame(seri):
    fasta = seri.set_index('fusion_id')['fusion_seq'].to_dict()

    out_path1 = sys.argv[2]

    if not os.path.exists(out_path1):
        os.makedirs(out_path1)

    for key,value in fasta.items():
        each_orf_filename = key + '.txt'
        filepath = os.path.join(out_path1,each_orf_filename)
        with open(filepath,'w') as f:
            f.write(">"+key+'\n'+value)

#storage gene id
#tmp_gene_id = seri.iloc[:,1]
#tmp_gene_id.to_csv(sys.argv[2],sep = '\t',index = False,header=False)

#orffinder input

#intermediate result
res1 = seri.iloc[:,:5]
res1.to_csv(sys.argv[2],sep = '\t',index = False,header=True)
