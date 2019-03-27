
# coding: utf-8

# In[14]:


import pandas as pd
import numpy as np
import sklearn
import scipy
import gzip


# In[15]:


get_ipython().magic(u'pylab inline')


# In[16]:


Sm = pd.read_csv("../../Pipeline_compound/similar3d_Wvender_polyols_SMILES.csv")


# In[17]:


Sm.info()


# In[18]:


Sm['CID'] = Sm['CID'].apply(str)


# In[19]:


import rdkit


# In[20]:


print('The rdkit version is {}.'.format(rdkit.__version__))


# In[21]:


from rdkit import Chem
mols = []
outfile = open('similar3D_Wvendor_polyols_cleaned.dat','w')
count = 0
#Sm['CID'] = Sm['CID'].apply(str)
for i in range(0,len(Sm.CanonicalSMILES)):    
    Canonm = Chem.MolFromSmiles(Sm.CanonicalSMILES[i])
    mols.append((Canonm, Sm.CID[i]))
    if Canonm is not None:
        outfile.write("%s\n" % ("\t".join(Sm.loc[i])))
    else:
        count += 1
outfile.close()
num_cids = len(Sm.CanonicalSMILES)
num_cids -= count # decrease the number of actives accordingly
print "number of cids with RDKit-invalid SMILES =", count


# In[22]:


Sm_canon = pd.read_table("similar3D_Wvendor_polyols_cleaned.dat", names = ["CID", "CanonicalSMILES", "IsomericSMILES"])


# In[23]:


Sm_canon.head()


# In[24]:


mols = []
count = 0
for i in range(0,len(Sm_canon.CanonicalSMILES)):    
    Canonm = Chem.MolFromSmiles(Sm_canon.CanonicalSMILES[i])
    mols.append((Canonm, Sm_canon.CID[i]))


# In[27]:


from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
# calculate fingerprints
fps = []
for Canonm,idx in mols:
    fps.append(Chem.RDKFingerprint(Canonm, maxPath=5))
# generate distance matrix
dist_matrix = []
num_fps = len(fps)
for i in range(1, num_fps):
    similarities = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
    dist_matrix.extend([1-x for x in similarities])
    #print i
# cluster
clusters = Butina.ClusterData(dist_matrix, num_fps, 0.5, isDistData=True) # distance cutoff = 0.5
print "number of clusters =", len(clusters)
num_clust_g5 = len([c for c in clusters if len(c) > 100])
print "number of clusters with more than 100 compounds =", num_clust_g5


# In[ ]:


# plot the size of the first 150 clusters
fig = plt.figure(1, figsize=(18, 5))
plt1 = plt.subplot(111)
plt.axis([0, 151, 0, len(clusters[0])+1])
plt.xlabel('Cluster index')
plt.ylabel('Number of molecules')
plt1.bar(range(1, 151), [len(c) for c in clusters[:150]], lw=0)
plt.show()


# In[ ]:


# take the cluster centre of each cluster
final_mols = [mols[c[0]] for c in clusters]

# sort the molecules within a cluster based on their similarity 
# to the cluster centre and sort the clusters based on their size
clusters2 = []
for c in clusters:
    if len(c) < 2: continue
    fps_clust = [Chem.RDKFingerprint(mols[i][0], maxPath=5) for i in c]
    simils = DataStructs.BulkTanimotoSimilarity(fps_clust[0], 
                                                fps_clust[1:])
    simils = [(s,i) for s,i in zip(simils, c[1:])]
    simils.sort(reverse=True)
    clusters2.append((len(simils), [i for s,i in simils]))
clusters2.sort(reverse=True)

# take 5 molecules (or 50%) of each cluster starting with the 
# largest one
idx = 0
diff = 100 - len(final_mols)
while diff > 0:
    c = clusters2[idx][1]
    if clusters2[idx][0] > 5:
        num_cmps = 5
    else:
        num_cmps = int(0.5*len(c))+1
    if num_cmps > diff: num_cmps = diff
    final_mols += [mols[i] for i in c[:num_cmps]]
    idx += 1
    diff = 100 - len(final_mols)
print "number of selected molecules =", len(final_mols)


# In[ ]:


num_clust = len(clusters)
mols_cluster0 = [final_mols[0]] + final_mols[num_clust:num_clust+5]
mcs = MCS.FindMCS([m[0] for m in mols_cluster0]) # get scaffold
scaffold = Chem.MolFromSmarts(mcs.smarts)
AllChem.Compute2DCoords(scaffold) # generate 2D coordinates for scaffold
for m,idx in mols_cluster0: # align molecules to scaffold
    AllChem.GenerateDepictionMatching2DStructure(m, scaffold)
Draw.MolsToGridImage([m[0] for m in mols_cluster0], 
                     legends=[m[1] for m in mols_cluster0], 
                     molsPerRow=4)


# In[ ]:


# write them to file
outfile = open('results/compounds_for_confirmatory_assay.txt', 'w')
outfile.write("#SAMPLE\tSMILES\n") # header
outfile.writelines("%s\t%s\r\n" % 
    (idx,Chem.MolToSmiles(m, isomericSmiles=True)) for m,i in final_mols)
outfile.close()

