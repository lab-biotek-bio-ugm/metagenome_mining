#!/usr/bin/env python
# coding: utf-8

# # Functional Analysis from EBI-MGNify
# 
# _Matin Nuhamunada_<sup>1*</sup>
# 
# <sup>1</sup>Department of Tropical Biology, Universitas Gadjah Mada;   
# Jl. Teknika Selatan, Sekip Utara, Bulaksumur, Yogyakarta, Indonesia, 55281;   
# 
# *Correspondence: [matin_nuhamunada@ugm.ac.id](mailto:matin_nuhamunada@mail.ugm.ac.id)  
# 
# ## Deskripsi Data
# Ada tiga file yang dapat diakses dari hasil analisis fungsional EBI-MGnify: (https://emg-docs.readthedocs.io/en/latest/portal.html#description-of-functional-annotation-files-available-to-download)
# 
# 1. *InterPro matches file*: it is a tab-delimited file containing 15 columns.
# 2. *Complete GO annotation file*: it is a comma-separated file containing 4 columns. The first column lists the GO terms (labelled GO:XXXXXXX) having been associated with the predicted CDSs. The second gives the GO term description while the third indicates which category the GO term belong to. There is 3 category: ‘biological process’ (higher biological process such as ‘rRNA modification’) , ‘molecular function’ (individual catalytic activity such as ‘mannosyltransferase activity’) and ‘cellular component’ (cellular localisation of the activity such as ‘mitochondrion’). The last column give the number of predicted CDSs having been annotated with the GO terms for the run.
# 3. *GO slim annotation file*: this file is derived from the ‘Complete GO annotation file’ and has the same format. The GO slim set is a cut-down version of the GO terms containing a subset of the terms in the whole GO. They give a broad overview of the ontology content without the details of the specific fine grained terms. Go slim terms are used for visualisation on the website. To illustrate how the GO slim terms relates to the GO terms, the different metal binding GO terms present in the ‘Complete GO annotation’ file are summarized as one generic metal binding term in the ‘GO slim annotation’ file. The last column give the number of predicted CDSs having been annotated with the GO slim terms for the run.
# 
# ## Challenge
# 1. Web Scraping menggunakan RESTful API dari EBI MGnify untuk mencari studi dengan kata kunci "human", "skin", dan kategori  sampel "metagenome" (hasil shotgun sequencing)
# 2. Melakukan (a) filtering studi dan pemilihan sampel untuk studi komparatif, dan (b) mengambil metadata dari sampel terpilih
# 3. Web Scraping menggunakan RESTful API dari EBI MGnify untuk mengunduh data hasil analisis fungsional (ketiga file di atas)
# 4. Melakukan data cleaning dan eksplorasi data (visualisasi)
# 5. Data mining dengan menggunakan analisis statistik: dengan dimensional reduction (PCoA, UMAP, T2SNE), dan teknik multivariat lainnya
# 

# ## Load Library

# In[7]:


import os
from pandas import DataFrame
import pandas as pd

import matplotlib.pyplot as plt
plt.close('all')

import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')


# ## Challenge 1

# In[8]:


#study data
study_data = pd.read_csv("data_study.csv")
study_data


# In[14]:


sample_data = pd.read_csv("data_sample.csv")
sample_data.sort_values('Description', ascending=False)
description = sample_data.Description.unique()
print(description)


# In[17]:


filtered_data = sample_data.loc[sample_data['Description'] == 'Human Skin Metagenome']
filtered_data = filtered_data.reset_index(drop=True)
filtered_data


# ## Challenge 2

# In[18]:


from jsonapi_client import Session, Filter
import pycurl
import html

API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/'


# In[19]:


# Updated
def get_study(search, filename):
    if not os.path.isfile(filename):
        with open(filename, 'wb') as f:
            c = pycurl.Curl()
            c.setopt(c.URL, 'https://www.ebi.ac.uk/metagenomics/api/v1/studies?lineage=root&ordering=-last_update&search='+search+'&format=csv')
            c.setopt(c.WRITEDATA, f)
            c.perform()
            c.close()
    return filename

def get_metadata(metadata, key):
    for m in metadata:
        if m['key'].lower() == key.lower():
            value = m['value']
            unit = html.unescape(m['unit']) if m['unit'] else ""
            return "{value} {unit}".format(value=value, unit=unit)
    return None

def get_analysis_result(run, extension):
    API_BASE_RUN = 'https://www.ebi.ac.uk/metagenomics/api/latest/runs'
    with Session(API_BASE_RUN) as s:
        study = s.get(run,'analysis').resource
        for i in study.downloads:
            if extension in i.file_format['name']:
                link = i.url
    return link

#perlu diupdate utk rhizosfer
def random_sampling(dataframe, amount):
    df_random = DataFrame(columns=('Sample_ID','Run_ID','Release_version','Sex','Body_site', 'Description'))
    df_random.index.name = 'No'
    a = 0
    while a < amount:
        i = np.random.choice(dataframe.index.values, 1)
        container = df_random.loc[:, 'Sample_ID']
        if not container.isin([dataframe.loc[i[0], 'Sample_ID']]).any():
            df_random.loc[i[0]] = [dataframe.loc[i[0], 'Sample_ID'],                                    dataframe.loc[i[0], 'Run_ID'],                                    dataframe.loc[i[0], 'Release_version'],                                    dataframe.loc[i[0], 'Sex'],                                    dataframe.loc[i[0], 'Body_site'],                                    dataframe.loc[i[0], 'Description']                                  ]
            a = a + 1
    return df_random


# In[ ]:





# ## Challenge 4-5 
# ### 5.1 Visualize Functional Analysis - Summary

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


# ### 5.2 Detailed GO

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


# ### 5.3 Interpro

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




