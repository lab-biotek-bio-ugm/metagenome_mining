#!/usr/bin/env python
# coding: utf-8

# # Load Library

# In[2]:


import os
from pandas import DataFrame
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn import preprocessing

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from itertools import repeat


import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')


# # Data Preparation
# ## Input
# List all summarized Gene Ontology summary from interproscan result. The data is FASTQ_GO_slim files from studies.
# 
# data should be extracted to folder "data"

# In[5]:


file = []
for name in os.listdir("data"):
    if name.endswith("_FASTQ_GO_slim.csv"):
        file.append(name)
file


# ## Merge GO_Slim data
# Merge all data to single dataframe and extract target = study name \
# Output: (1) Dataframe ("df_raw") and (2) target ("df_raw_target") 

# In[3]:


# merge and extract target
for num, i in enumerate(file):
    if num == 0:
        df_raw = pd.read_csv(i)
        len_target = len(list(df_raw.columns)) - 3
        study=(list(repeat(i.replace('_FASTQ_GO_slim.csv',''), len_target)))
    else:
        df = pd.read_csv(i)
        len_target = len(list(df.columns)) - 3
        study = study + list(repeat(i.replace('_FASTQ_GO_slim.csv',''), len_target))
        df_raw = pd.merge(df_raw, df, on=['GO term','category','description'], how='left')
df_raw = df_raw.fillna(0)

# the target list for samples
sample = list(df_raw.drop(columns=["GO term", "category", "description"]).T.index)
study = [study, sample]
df_raw_target = pd.DataFrame(study, index=['Study', 'Sample']).T


# In[4]:


# Raw Data containing reads of sequences found in metagenome related to GO pathways
df_raw.head(2)


# In[5]:


df_raw_target.head(2)


# ## Write Output

# In[6]:


# Export raw data to csv
if not os.path.exists('output'):
    os.mkdir('my_folder')
df_raw.to_csv("./output/df_raw.csv")
df_raw_target.to_csv("./output/df_raw_target.csv")


# # Exploratory Data Analysis
# ## Sort by Sum of Feature

# In[7]:


#plt.figure(figsize = (20,10))
#sns.heatmap(df_raw.drop(columns=['GO term','category']).set_index('description'))


# In[8]:


x = df_raw.set_index(['GO term', 'category', 'description'])
x.index.get_level_values(0)
df_raw.drop(columns='category')


# In[9]:


# Sum feature
#df_raw_sorted = df_raw.drop(columns=['category'])
df_raw_sorted = df_raw.set_index(['GO term', 'category', 'description'])
df_raw_sorted['total'] = df_raw_sorted.sum(axis=1)
df_raw_sorted['total'] = df_raw_sorted['total']/(len(df_raw.index))

# Sort feature by sum of numbers
df_raw_sorted = df_raw_sorted.sort_values('total', ascending=True)
df_raw_sorted = df_raw_sorted.drop(columns=['total'])
#df_raw_sorted.head(2)

plt.figure(figsize = (10,18))
sns.set(font_scale=0.7) 
#sns.set_context("paper")
sns.heatmap(df_raw_sorted, 
            yticklabels=df_raw_sorted.index.get_level_values(2),
           cbar_kws={"orientation": "horizontal"})


# # Preprocessing

# In[10]:


# normalizing by MinMax Scaler
min_max_scaler = preprocessing.MinMaxScaler()
scaled_array = min_max_scaler.fit_transform(df_raw_sorted)
df_norm_minmax = pd.DataFrame(scaled_array, columns=df_raw_target.Sample, index=df_raw_sorted.index)
plt.figure(figsize = (10,18))
sns.set(font_scale=0.7) 
sns.heatmap(df_norm_minmax,
           yticklabels=df_raw_sorted.index.get_level_values(2),
           cbar_kws={"orientation": "horizontal"})


# In[20]:


from scipy.stats import zscore
df_norm_zscore = df_raw_sorted.apply(zscore)

plt.figure(figsize = (10,18))
sns.set(font_scale=0.7) 
sns.heatmap(df_norm_zscore,
           yticklabels=df_norm_zscore.index.get_level_values(1),
           cbar_kws={"orientation": "horizontal"})


# ## Prepare Feature for Machine Learning

# In[12]:


data = df_norm_minmax.T
data.head(2)


# In[17]:


eda = data.describe()
dropfeature = []
for i in range(116):
    x = eda.iloc[:,i]
    if x['25%'] == x['50%'] == 0:
    #if x['25%'] == 0:
        #print(x.name, x['25%'], x['50%'], x['75%'])
        dropfeature.append(str(x.name))
len(dropfeature)
data = data.drop(index = dropfeature)


# In[ ]:


data = data.groupby(level=0).mean()


# In[ ]:


plt.figure(figsize = (8,16))
sns.heatmap(data.T)


# In[ ]:


#target
#y = {'target': data.index.values}
y = {'target': study}

#feature
x = data.values

# Standardizing the features
x = StandardScaler().fit_transform(x)


# In[ ]:


pca = PCA(n_components=3)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2', 'principal component 3'])
principalDf


# In[ ]:


finalDf = pd.concat([principalDf, pd.DataFrame(y)], axis=1)
finalDf


# In[ ]:


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = finalDf.target.unique()
colors = ['b', 'w', 'r', 'c', 'm', 'y', 'k', 'g']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50
               , alpha=0.5)
ax.legend(targets)
ax.grid()
#bg_color = 'grey'
#ax.patch.set_facecolor(bg_color)

plt.show()


# In[ ]:


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('Principal Component 1', fontsize = 8)
ax.set_ylabel('Principal Component 2', fontsize = 8)
ax.set_zlabel('Principal Component 3', fontsize = 8)
ax.set_title('3 component PCA', fontsize = 20)
targets = finalDf.target.unique()
colors = ['b', 'w', 'r', 'c', 'm', 'y', 'k', 'g']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 50
               , alpha=0.3)
ax.legend(targets)
ax.grid()
#bg_color = 'grey'
#ax.patch.set_facecolor(bg_color)

plt.show()


# In[ ]:




