#!/usr/bin/env python
# coding: utf-8

# In[243]:


import os
from pandas import DataFrame
import pandas as pd


# # Data Cleaning

# In[134]:


# Load samples
sample_data = pd.read_csv("data_sample.csv")

# Rapikan label
label = []
for i in sample_data.Description:
    y = " ".join(i.split())
    label.append(y)

# Replace label
for x in range(len(sample_data)):
    sample_data.Description[x] = label[x]
sample_data.head()


# # Persiapan Scraping

# In[135]:


# Ambil unique value dari study id untuk filter
study_id = sample_data['MGnify ID'].unique()
study_id


# In[136]:


#filter data untuk diunduh per study
filtered_data = sample_data.loc[sample_data['MGnify ID'] == study_id[4]] 
#ganti nomor ini dari study id, untuk kasus ini ane ambil yg paling sedikit dulu

filtered_data = filtered_data.reset_index(drop=True)
filtered_data


# # Data Scraping

# In[137]:


from jsonapi_client import Session, Filter
import html
import pycurl


# In[152]:


# create get function for sample
def ebi_sample(sample):
    API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/samples'
    with Session(API_BASE) as s:
        params = {'page_size': 100}
        f = Filter(urlencode(params))
        sample = s.get(sample).resource
    return sample


# In[148]:


filtered_data.Sample[0]


# ### mari kita lihat data apa saja yg bisa diperoleh dari sample

# In[141]:


sampled = ebi_sample('ERS805755')

print(sampled.accession)
print(sampled.analysis_completed)
print(sampled.biosample)
print(sampled.collection_date)
print(sampled.environment_biome)
print(sampled.environment_feature)
print(sampled.environment_material)
print(sampled.geo_loc_name)
print(sampled.host_tax_id)
print(sampled.last_update)
print(sampled.latitude)
print(sampled.longitude)
print(sampled.runs)
print(sampled.sample_alias)
print(sampled.sample_desc)
print(sampled.sample_metadata)
print(sampled.species)


# ## kita tertarik untuk mengambil data: 1. accession, 2. runs, dan 3. sample_metadata

# In[154]:


#create containers
scraping = []
for i in filtered_data.Sample:
    sample = ebi_sample(i)
    a = sample.accession
    r = sample.runs
    m = sample.sample_metadata
    data = [a, r, m]
    scraping.append(data)
scraping


# In[155]:


# create dataframe from list
df_scrape = pd.DataFrame(scraping, columns = ["Sample", "Runs", "Metadata"]) 
df_scrape 


# # merge and write

# In[156]:


#merge based on similar value
df_scrapingresult = pd.merge(filtered_data, df_scrape, on="Sample")
df_scrapingresult


# In[157]:


#write to csv jangan lupa study yg mana
df_scrapingresult.to_csv('scraping_result_study_MGYS00000518.csv', index=False)


# # Download Analysis result

# In[236]:


import os
def download_GO_Slim(run_id):
    outname = run_id+'_FASTQ_GO_slim.csv'
    outdir = 'result_GO_slim'
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fullname = os.path.join(outdir, outname)    

    df = DataFrame(columns=('category', 'description', 'annotation counts'))
    df.index.name = 'GO term'
    
    with Session(API_BASE) as s:
        run = s.get('runs', run_id).resource
        for a in run.analyses:
            for ann in a.go_slim:
                df.loc[ann.accession] = [
                    ann.lineage, ann.description, ann.count
                ]
    return df.to_csv(fullname)


# In[242]:


for i in df_scrapingresult.Runs:
    download_GO_Slim(i[0].id)


# In[ ]:




