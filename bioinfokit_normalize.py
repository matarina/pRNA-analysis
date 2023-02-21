import os
os. getcwd() 
os.chdir('H:/')
from bioinfokit.analys import norm, get_data
import pandas as pd
featurecount = pd.read_csv('H:/prna_transcript_analysis/epd_promoter/chromesome_distribution/SRR8131644mrna.count',sep='\t',skiprows = 1)
featurecount.head(2)

##get RPKM
df = featurecount
df = df.drop(['Chr','Start','End','Strand'], axis=1)
df = df.set_index('Geneid')

df.head(2)
nm = norm()
nm.rpkm(df=df, gl='Length')
# get RPKM normalized dataframe
rpkm_df = nm.rpkm_norm
rpkm_df.head(2)
rpkm_df.to_csv('H:/single_cell_sepcific_tpm.csv',index=True)


##get RPM
# now, normalize raw counts using CPM method 
df = featurecount
df = df.set_index('Geneid')

df.head(2)

nm = norm()
nm.cpm(df=df)
# get CPM normalized dataframe
cpm_df = nm.cpm_norm
cpm_df.head(2)
cpm_df.to_csv('H:/hepg2_2_prna.count.csv',index=True)




###get TPMdf = featurecount
df = featurecount
df = df.drop(['Chr','Start','End','Strand'], axis=1)
df = df.set_index('Geneid')

length = [1000]*11147
df['length'] = length

df.head(2)
df.tail(2)

nm = norm()
nm.tpm(df=df, gl='Length')
# get TPM normalized dataframe
tpm_df = nm.tpm_norm
tpm_df.head(2)
tpm_df.to_csv('H:/single_cell_sepcific_tpm.csv',index=True)
