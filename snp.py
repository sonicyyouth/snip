import glob
import numpy as np
import os
import pandas as pd
fdir = r'/users/liuqun/Dropbox/Rworkspace/projects'
os.getcwd()

os.chdir(fdir+r'/snp')
rppadir = fdir+"/tcgafile/Level_3"
rppalist =[]
namemap = pd.read_table(''.join(glob.glob(rppadir+"/FILE*.txt")),index_col = 0,header = 0)
for i in glob.iglob(rppadir+"/mdanderson*.txt"):
    rpa = pd.read_table(i,index_col=0,header=1)
    rpa.columns = [namemap.ix[i[len(rppadir)+1:]][0]]
    rppalist.append(rpa)

len(rppalist[1].index)
len(rppalist[3].index)
x02 = rppalist[1].index.difference(rppalist[3].index)
x03 = rppalist[3].index.difference(rppalist[1].index)
x04 = x02 | x03
x041 = x04.delete([2,17,18,21,24,25,28,31])
x042 = dict(zip(x041[0::2],x041[1::2]))
for i in range(len(rppalist)): rppalist[i] = rppalist[i].rename(index =x042)
rppa = pd.DataFrame()
rppa = pd.concat(rppalist,axis=1)
rppana = rppa.dropna(axis=0)
rppa.to_pickle('lggrppa.pkl')
rppa.to_csv('lggrppa.csv')

# rppalen = []

# for i in range(len(rppalist)):
#     rlen = len(rppalist[i])
#     print(rlen)
#     print(i)
#     rppalen.append(rlen)


# x01 = rppalist[1].index.values
# x02 = rppalist[0].index.values
# x00 = np.setdiff1d(x01,x02)
# len(x01)
dates = pd.date_range('1/1/2000', periods= 8)
df = pd.DataFrame(np.random.randn(8,4),index = dates,columns= list('ABCD'))
df1 = df.copy()
df2 = pd.DataFrame({'a' : ['one', 'one', 'two', 'three', 'two', 'one', 'six'],'b' : ['x', 'y', 'y', 'x', 'y', 'x', 'x'],'c' : np.random.randn(7)})
df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar','foo', 'bar', 'foo', 'foo'],'B' : ['one', 'one', 'two', 'three','two', 'two', 'one', 'three'],'C' : np.random.randn(8), 'D' : np.random.randn(8)})
grouped = df.groupby('A')

index = pd.date_range('10/1/1999', periods=1100)
ts = pd.Series(np.random.normal(0.5, 2, 1100), index)
ts = pd.rolling_mean(ts, 100, 100).dropna()
key = lambda x: x.year
zscore = lambda x: (x - x.mean()) / x.std()
transformed = ts.groupby(key).transform(zscore)
compare = pd.DataFrame({'Original': ts, 'Transformed': transformed})
compare.plot()
countries = np.array(['US', 'UK', 'GR', 'JP'])
key = countries[np.random.randint(0, 4, 1000)]
grouped = data_df.groupby(key)
# display large dataframes in an html iframe
def ldf_display(df, lines=500):
    txt = ("<iframe " +
           "srcdoc='" + df.head(lines).to_html() + "' " +
           "width=1000 height=500>" +
           "</iframe>")
    return IPython.display.HTML(txt)

lggclin = pd.read_table('../tcgafile/LGG.clin.merged.txt',index_col = 0,header = 0)
lggclin1 = pd.read_table('../tcgafile/GBMLGG.merged_only_clinical_clin_format.txt',index_col = 0,header = 0)
lggclin1.shape
x = lggclin1.index.tolist()
x1 = [item for i, item in enumerate(x) if re.search('-[0-9]+\.',item)]
lggclin1 = lggclin1.drop(x1).T
lggclin1 = lggclin1.replace({",":";"})

lggclin1.to_csv('gbmclin1.csv')

import re
sdel = re.compile('-1[0-9]+\.')
for i in range(len(x)): = sdel.findall(x.tolist())
gbmclin = pd.read_csv('gbmclin1.csv',index_col = 0,header = 0)
gbmclin.shape
x = gbmclin.columns.values
x[42]  = 'vitual_status'
gbmclin.columns = x
gbmclin[x[3]].notnull().sum()
gbmclin[x[5]].notnull().sum()
gbmclin['time'] = gbmclin.iloc[:,5]
gbmclin.loc[gbmclin.iloc[:,3].notnull(),'time'] = gbmclin.loc[gbmclin.iloc[:,3].notnull(),x[3]]
gbmclin.iloc[:,42].value_counts()
gbmclin = gbmclin.replace({",":";"})
gbmclin[x[42]] = gbmclin[x[42]].map({'alive':0,'dead':1})
gbmclin['histology'] = gbmclin['histology'].replace({'untreated primary (de novo) gbm':'primary GBM','glioblastoma multiforme (gbm)':'primary GBM'})
gbmclin.to_csv('gbmclin1.csv')
lggrna = pd.read_csv('../tcgafile/GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt',sep = '\t',index_col = 0)
lggrna = lggrna.drop(lggrna.index[0])
lggrna.columns.str[13:15].value_counts()
lggrna.columns = lggrna.columns.str[:15]
x = lggrna.columns.str[:12]
lggrna1 = lggrna.drop(lggrna.columns[lggrna.columns.str[13:15] == "02"],axis = 1)
lggrna1.to_csv('primyglmrna.csv')

rembanno = pd.read_table('A-AFFY-44.adf.txt')
rembanno.to_csv('rembanno.csv')
