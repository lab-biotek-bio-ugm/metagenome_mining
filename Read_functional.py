#!/usr/bin/env python
# coding: utf-8

# # Functional Analysis from EBI-MGNify
# 
# _Matin Nuhamunada_<sup>1*</sup>, _Ahmad Ardi_<sup>1*</sup>
# 
# <sup>1</sup>Department of Tropical Biology, Universitas Gadjah Mada;   
# Jl. Teknika Selatan, Sekip Utara, Bulaksumur, Yogyakarta, Indonesia, 55281;   
# 
# *Correspondence: [matin_nuhamunada@ugm.ac.id](mailto:matin_nuhamunada@mail.ugm.ac.id)  
# 
# ## Daftar Isi
# 1. Deskripsi Data
# 2. Challenge
# 3. Load Library
# 4. Web Scraping & Data Cleaning
# 5. Exploratory Data Analysis (EDA)
# 6. Data Mining
# 
# ## A. Deskripsi Data
# Ada tiga file yang dapat diakses dari hasil analisis fungsional EBI-MGnify: (https://emg-docs.readthedocs.io/en/latest/portal.html#description-of-functional-annotation-files-available-to-download)
# 
# 1. *InterPro matches file*: it is a tab-delimited file containing 15 columns.
# 2. *Complete GO annotation file*: it is a comma-separated file containing 4 columns. The first column lists the GO terms (labelled GO:XXXXXXX) having been associated with the predicted CDSs. The second gives the GO term description while the third indicates which category the GO term belong to. There is 3 category: ‘biological process’ (higher biological process such as ‘rRNA modification’) , ‘molecular function’ (individual catalytic activity such as ‘mannosyltransferase activity’) and ‘cellular component’ (cellular localisation of the activity such as ‘mitochondrion’). The last column give the number of predicted CDSs having been annotated with the GO terms for the run.
# 3. *GO slim annotation file*: this file is derived from the ‘Complete GO annotation file’ and has the same format. The GO slim set is a cut-down version of the GO terms containing a subset of the terms in the whole GO. They give a broad overview of the ontology content without the details of the specific fine grained terms. Go slim terms are used for visualisation on the website. To illustrate how the GO slim terms relates to the GO terms, the different metal binding GO terms present in the ‘Complete GO annotation’ file are summarized as one generic metal binding term in the ‘GO slim annotation’ file. The last column give the number of predicted CDSs having been annotated with the GO slim terms for the run.
# 
# ## B. Challenge
# 1. Web Scraping menggunakan RESTful API dari EBI MGnify untuk mencari studi dengan kata kunci "human", "skin", dan kategori  sampel "metagenome" (hasil shotgun sequencing)
# 2. Melakukan (a) filtering studi dan pemilihan sampel untuk studi komparatif, dan (b) mengambil metadata dari sampel terpilih
# 3. Web Scraping menggunakan RESTful API dari EBI MGnify untuk mengunduh data hasil analisis fungsional (ketiga file di atas)
# 4. Melakukan data cleaning dan eksplorasi data (visualisasi)
# 5. Data mining dengan menggunakan analisis statistik: dengan dimensional reduction (PCoA, UMAP, T2SNE), dan teknik multivariat lainnya
# 

# ## C. Load Library

# In[1]:


import os
from pandas import DataFrame
import pandas as pd
import pycurl

from scipy import stats

import matplotlib.pyplot as plt
plt.close('all')

import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


from jsonapi_client import Session, Filter
import html

API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/'


# ## D. Webscraping & Data Cleaning
# Web Scraping menggunakan RESTful API dari EBI MGnify untuk mencari studi dengan kata kunci "human", "skin", dan kategori sampel "metagenome" (hasil shotgun sequencing)

# ### D.1 Study
# Pada part ini kita akan mencoba melihat ada berapa banyak studi terkait dengan Human Skin di EBI MGnify

# In[3]:


#obtaining study data - TO DO!
def get_biome(lineage, exp_type):
    API_BASE_BIOME = 'https://www.ebi.ac.uk/metagenomics/api/latest/biomes'
    with Session(API_BASE_BIOM) as s:
        study = s.get(run,'analysis').resource
        for i in study.downloads:
            if extension in i.file_format['name']:
                link = i.url
    return link


# In[15]:


API_BASE_BIOME = 'https://www.ebi.ac.uk/metagenomics/api/latest/biomes'
with Session(API_BASE_BIOME) as s:
    biome = s.get('root:Host-associated:Human:Skin').resource


# In[21]:


biome.samples_count


# In[22]:


biome.studies


# In[23]:


biome.samples


# In[4]:


#load study data
study_data = pd.read_csv("data_study.csv")
study_data


# In[5]:


print(study_data.describe())
print('\n' + 'Total sampel: ' + str(study_data.Samples.sum()) + '\n')
print('Centers : ' + str(len(study_data['Centre name'].unique())))
x = 0
for i in study_data['Centre name'].unique():
    x = x + 1
    print(str(x) + ' ' + i)


# ### D.2 Kesimpulan Studi
# Diperoleh 15 studi di EBI yang dilakukan dari 10 institusi, dengan total sampel terkait human skin microbiome sebanyak 4053

# ## D.3 Samples
# Pada part ini kita akan melihat sampel apa saja yang tersedia untuk human skin metagenome

# In[6]:


# Load samples
sample_data = pd.read_csv("data_sample.csv")


# In[7]:


# Rapikan label
label = []
for i in sample_data.Description:
    y = " ".join(i.split())
    label.append(y)


# In[8]:


# Replace label
for x in range(len(sample_data)):
    sample_data.Description[x] = label[x]
sample_data.head()


# In[9]:


x = 0
for i in sample_data.Description.unique():
    x = x + 1
    print(str(x) + ' ' + i)


# In[10]:


sample_data.groupby(['MGnify ID']).describe()


# In[11]:


sample_data['MGnify ID'].count()


# In[12]:


study_id = sample_data['MGnify ID'].unique()
study_id_clean = []
for x in study_id:
    y = x.split(',')
    for z in y:
        study_id_clean.append(z)
study_id_clean = list(dict.fromkeys(study_id_clean))


# In[13]:


x = 0
for i in study_id_clean:
    x = x + 1
    print(str(x) + ' ' + i)


# In[14]:


df = study_data.loc[study_data['MGnify ID'].isin(study_id_clean)]
df.reset_index(drop=True)


# ### Kesimpulan
# Hanya ada 9 Studi terkait data metagenome

# ### Ambil Metadata

# In[15]:


study_id


# In[16]:


filtered_data = sample_data.loc[sample_data['MGnify ID'] == study_id[0]] #ganti nomor ini dari study id
filtered_data = filtered_data.reset_index(drop=True)
filtered_data


# ### Getting Sample Metadata

# In[17]:


def ebi_sample(sample):
    API_BASE = 'https://www.ebi.ac.uk/metagenomics/api/latest/samples'
    with Session(API_BASE) as s:
        sample = s.get(sample).resource
    return sample


# In[18]:


sampled = ebi_sample('SRS731606')


# In[19]:


print(sampled.accession)
print(sampled.analysis_completed)
#print(sampled.biome)
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
#print(sampled.studies)


# In[20]:


sampled = ebi_sample('SRS731606')
for x in sampled.runs:
    print(x.id)
    for y in x.analyses:
        print(y.id)
        for z in y.downloads:
            if 'FASTQ_GO_slim.csv' in z.url:
                print(z.url)


# In[21]:


filename = y.id + '_FASTQ_GO_slim.csv'
link = z.url

if not os.path.isfile(filename):
    with open(filename, 'wb') as f:
        c = pycurl.Curl()
        c.setopt(c.URL, link)
        c.setopt(c.WRITEDATA, f)
        c.perform()
        c.close()


# In[21]:


with open('SRR1631558_MERGED_FASTQ_GO_slim.csv', 'wb') as f:
    c = pycurl.Curl()
    c.setopt(c.URL, 'https://www.ebi.ac.uk/metagenomics/api/v1/analyses/MGYA00023807/file/SRR1631558_MERGED_FASTQ_GO_slim.csv')
    c.setopt(c.WRITEDATA, f)
    c.perform()
    c.close()


# In[23]:


t = y.downloads[7]
t.url


# ## Challenge 2

# In[24]:


def get_analysis_result(run, extension):
    API_BASE_RUN = 'https://www.ebi.ac.uk/metagenomics/api/latest/runs'
    with Session(API_BASE_RUN) as s:
        study = s.get(run,'analysis').resource
        for i in study.downloads:
            if extension in i.file_format['name']:
                link = i.url
    return link


# In[ ]:


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


API_BASE_SAMPLE = 'https://www.ebi.ac.uk/metagenomics/api/v1/samples'
API_BASE_RUNS = 'https://www.ebi.ac.uk/metagenomics/api/v1/runs'
with Session(API_BASE_SAMPLE) as s:
    sample = s.get('SRS731606','runs').resources
    for a in sample:
        print(a)


# In[ ]:


a.accession


# In[ ]:


API_BASE_RUNS = 'https://www.ebi.ac.uk/metagenomics/api/v1/runs'
with Session(API_BASE_RUNS) as r:
    runs = r.get('SRR1631560','analyses').resources
    for x in runs:
        print(x)


# In[ ]:


y = x.downloads
z = x.go_slim


# ## Challenge 4-5 
# ### 5.1 Visualize Functional Analysis - Summary

# In[22]:


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


# In[25]:


df_slim.to_csv('df_slim.csv')


# In[23]:


df_slim_sum = df_slim.drop(columns=['GO_ID','Pathway'])
df_slim_sum = df_slim_sum.set_index('Category')
df_slim_sum.head()


# In[27]:


df_slim_sum['total'] = df_slim_sum.sum(axis=1)
df_slim_sum['total'] = df_slim_sum['total']/3
df_slim_sum = df_slim_sum.sort_values('total', ascending=True)
df_plot = df_slim_sum.drop(columns=['total'])
plt.figure()
df_plot.plot.barh(figsize=(8, 35))


# In[28]:


sns.heatmap(df_plot[:100])


# In[29]:


df_plot_z = df_plot.apply(stats.zscore)
plt.figure()
df_plot_z.plot.barh(figsize=(8, 35))


# In[30]:


sns.heatmap(df_plot_z[:100])


# ### 5.2 Detailed GO

# In[80]:


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


# In[81]:


df_GO_sum = df_GO.drop(columns=['GO_ID','Pathway'])
df_GO_sum = df_GO_sum.set_index('Category')
df_GO_sum.describe()


# In[82]:


df_GO_sum['total'] = df_GO_sum.sum(axis=1)
df_GO_sum['total'] = df_GO_sum['total']/2
df_GO_sum = df_GO_sum.sort_values('total', ascending=False)
df_GO_sum.head()


# In[83]:


df_heatmap = df_GO_sum.drop(columns=['total'])
sns.heatmap(df_heatmap[:30])#, annot=True, linewidths=.25)


# ### 5.3 Interpro

# In[ ]:


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


df.Sequence_length.describe()


# In[ ]:


df.Sequence_length.unique()


# In[ ]:


x = df.Interpro_accession.unique()
x


# In[ ]:


df.Interpro_accession.describe()


# In[ ]:


df.Interpro_accession


# In[ ]:




