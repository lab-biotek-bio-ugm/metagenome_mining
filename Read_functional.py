#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from pandas import DataFrame
import pandas as pd

import matplotlib.pyplot as plt
plt.close('all')

import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')


# # Visualize Functional Analysis - Summary

# In[2]:


df_slim = pd.DataFrame(columns=['GO_ID', 'Category', 'Pathway'])
file = os.listdir()

for a in range(len(file)):
    if 'slim.csv' in file[a]:
        sample_ID = file[a].replace('_FASTQ_GO_slim.csv', '')
        x = pd.read_csv(file[a], header=None)
        x.columns = ['GO_ID', 'Category', 'Pathway', sample_ID]
        x[sample_ID] = (x[sample_ID] / x[sample_ID].sum()) * 100
        #print(x.head())
        df_slim = pd.merge(x, df_slim, on=['GO_ID','Category','Pathway'], how='left')

df_slim.head()


# In[58]:


df_slim_sum = df_slim.drop(columns=['GO_ID','Pathway'])
df_slim_sum = df_slim_sum.set_index('Category')
df_slim_sum.describe()


# In[59]:


df_slim_sum['total'] = df_slim_sum.sum(axis=1)
df_slim_sum['total'] = df_slim_sum['total']/2
df_slim_sum = df_slim_sum.sort_values('total', ascending=False)
df_slim_sum.head()


# In[63]:


df_heatmap = df_slim_sum.drop(columns=['total'])
sns.heatmap(df_heatmap[:20])#, annot=True, linewidths=.25)


# # Detailed GO

# In[64]:


df_GO = pd.DataFrame(columns=['GO_ID', 'Category', 'Pathway'])
file = os.listdir()

for a in range(len(file)):
    if 'FASTQ_GO.csv' in file[a]:
        sample_ID = file[a].replace('_FASTQ_GO.csv', '')
        x = pd.read_csv(file[a], header=None)
        x.columns = ['GO_ID', 'Category', 'Pathway', sample_ID]
        x[sample_ID] = (x[sample_ID] / x[sample_ID].sum()) * 100
        #print(x.head())
        df_GO = pd.merge(x, df_GO, on=['GO_ID','Category','Pathway'], how='left')

df_GO.head()


# In[65]:


df_GO_sum = df_GO.drop(columns=['GO_ID','Pathway'])
df_GO_sum = df_GO_sum.set_index('Category')
df_GO_sum.describe()


# In[66]:


df_GO_sum['total'] = df_GO_sum.sum(axis=1)
df_GO_sum['total'] = df_GO_sum['total']/2
df_GO_sum = df_GO_sum.sort_values('total', ascending=False)
df_GO_sum.head()


# In[70]:


df_heatmap = df_GO_sum.drop(columns=['total'])
sns.heatmap(df_heatmap[:30])#, annot=True, linewidths=.25)


# # Interpro

# In[123]:


df = pd.read_csv('ERR476421_FASTQ_I5.tsv', header=None, sep='\n')
df = df[0].str.split('\t', expand=True)
df.columns = ['Protein_accession',
              'Sequence_MD5_digest',
              'Sequence_length',
              'Analysis',
              'Signature_accession',
              'Signature_description',
              'Start_location',
              'Stop_location',
              'Score',
              'Status',
              'Date',
              'Interpro_accession',
              'Interpro_description',
              'GO_annotations']
df.head()


# In[ ]:




