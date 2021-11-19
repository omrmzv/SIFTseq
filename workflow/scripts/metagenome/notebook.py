#!/usr/bin/env python
# coding: utf-8

# In[178]:


from pyfaidx import Fasta
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import sys

# In[14]:

genome_fasta = sys.argv[1]
unfiltered_grammy = sys.argv[2]
unfiltered_tblat1_file = sys.argv[3]
filt_grammy = sys.argv[4]
filt_tblat1_file = sys.argv[5]
gi_to_tax = sys.argv[6]
outdir = sys.argv[7]
outfile = sys.argv[8]
genome = Fasta(genome_fasta)


# In[75]:


unfilt = pd.read_csv(unfiltered_grammy, sep = '\t')
unfilt = unfilt[["species", "AdjustedBlast", "RelCoverage"]]
unfilt = unfilt.dropna()
unfilt = unfilt.groupby('species').sum()
unfilt.reset_index(inplace=True)
unfilt.columns = ["species", "unfiltAB", "unfiltRC"]

filt = pd.read_csv(filt_grammy, sep = '\t')
filt = filt[["species", "AdjustedBlast", "RelCoverage"]]
filt = filt.dropna()
filt = filt.groupby('species').sum()
filt.reset_index(inplace=True)
filt.columns = ["species", "filtAB", "filtRC"]

df = pd.merge(left=unfilt, right = filt, on = 'species').astype({'species':'int'}).astype({'species':'str'})

fig = plt.figure()
ax = plt.gca()
ax.plot([min(df['unfiltAB']), max(df['unfiltAB'])], [min(df['unfiltAB']), max(df['unfiltAB'])])
ax.scatter(df['unfiltAB'], df['filtAB'])
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("Unfiltered AdjustedBlast")
ax.set_ylabel("Filtered AdjustedBlast")
plt.savefig(f'{outdir}/figures/unfilt_vs_filt_AdjBlast.png', dpi = 300, )
plt.close()

# In[419]:

filt_tblat1 = pd.read_csv(filt_tblat1_file, sep = '\t', header = None, names = ['read_id', 'ggi', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'strand'])
tmp = filt_tblat1['ggi'].str.split("|", expand = True)
filt_tblat1['gi'] = tmp[1]
filt_tblat1['refseq'] = tmp[3]
unfilt_tblat1 = pd.read_csv(unfiltered_tblat1_file, sep = '\t', header = None, names = ['read_id', 'ggi', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'strand'])
tmp = unfilt_tblat1['ggi'].str.split("|", expand = True)
unfilt_tblat1['gi'] = tmp[1]
unfilt_tblat1['refseq'] =tmp[3]

tmp = None
gi_to_tax = pd.read_csv(gi_to_tax, sep = '\t', header=None, names=['gi', 'Taxid'], dtype = {'gi': 'str', 'Taxid':'int'})

def get_species_hits(species, grammy, tblat1, gi_to_tax):
    grammy = pd.read_csv(grammy, sep = '\t')

    grammy = grammy[["Taxid", "species"]]

    grammy = grammy[grammy["species"]==int(species)]

    df = pd.merge(left=grammy, right = gi_to_tax, on = 'Taxid')

    C_content = []

    gis = np.unique(df["gi"].values).tolist()
    tblat1_sub = tblat1[tblat1['gi'].isin(gis)]
    #print(tblat1_sub)
    tblat_rows = tblat1_sub.index.tolist()
    #print(tblat_rows)
    #tblat_rows = []
    #j = 0

    for _, line in tblat1_sub.iterrows():
        refseq = line['refseq']
        sstart = int(line['sstart'])
        ssend = int(line['send'])
        strand = line['strand']

        if strand == 'Plus/Plus':
            start = sstart
            end = ssend
        else:
            start = ssend
            end = sstart
        seq = genome[refseq][start:end].seq
        C_content.append(float(seq.count('C'))/len(seq))
        #tblat_rows.append(j)
        #j+=1
    return(C_content, tblat_rows)

# In[420]:

def adjust_hits(species, filt_grammy, filt_tblat1, unfilt_grammy, unfilt_tblat1, gi_to_tax):

    C_content_filt, tblat_rows = get_species_hits(species, filt_grammy, filt_tblat1, gi_to_tax)
    C_content_unfilt, _ = get_species_hits(species, unfilt_grammy, unfilt_tblat1, gi_to_tax)
    a= stats.ks_2samp(C_content_unfilt, C_content_filt)
    print(a)
    print(len(C_content_filt))
    print(len(C_content_unfilt))
    if (a[0] < 0.05) or (a[1] > 0.01): #Dstat and pv
        C_content_filt2 = C_content_filt.copy()
        bad_rows = []
    else:
        to_remove = []
        to_rem = []

        maxks = a
        i = 0
        bad_rows = []
        for value, hit_row in zip(C_content_filt, tblat_rows):
            to_rem = to_remove.copy()
            to_rem.append(i)
            if (maxks[0]<0.01) and (maxks[1]>0.01):
                break
            dutie = C_content_filt.copy()
            to_rem = sorted(to_rem, reverse=True)
            for idx in to_rem:
                dutie.pop(idx)

            if len(dutie) > 0:

                newks = stats.ks_2samp(C_content_unfilt, dutie)

                if newks[0] < maxks[0]:
                    maxks = newks
                    to_remove.append(i)
                    bad_rows.append(hit_row)
                i+=1

        C_content_filt2 = C_content_filt.copy()
        to_remove = sorted(to_remove, reverse=True)
        for idx in to_remove:
            C_content_filt2.pop(idx)

    dens=False
    fig, axs = plt.subplots(3)
    axs[0].hist(C_content_unfilt, density=dens, bins=100, label = f"unfiltered {str(len(C_content_unfilt))} hits", alpha =0.75)
    axs[1].hist(C_content_filt,density = dens, bins=100, label = f"filtered {str(len(C_content_filt))} hits", alpha = 0.75)
    axs[2].hist(C_content_filt2, density=dens, bins=100, label = f"refiltered {str(len(C_content_filt2))} hits", alpha = 0.75)

    for ax in axs:
        ax.set(xlabel = "cytosine fraction of reference at mapped locations", ylabel = "count")
        ax.set_xlim(0,1)
        ax.legend()

    plt.savefig(f'{outdir}/figures/{str(species)}_filter.png', dpi = 300)
    plt.close()
    if (float(len(C_content_filt2))/len(C_content_filt) <0.95):
            print(f"Removed more than 5% of reads for {str(species)}")
    print(len(C_content_filt2))
    return(bad_rows)

rows_to_remove = []
for species in df['species'].to_list():
    print(species)
    a = adjust_hits(species, filt_grammy, filt_tblat1, unfiltered_grammy, unfilt_tblat1, gi_to_tax)
    rows_to_remove += a

refiltered = filt_tblat1.drop(rows_to_remove)
refiltered.to_csv(outfile, sep = '\t', header=False, index=False)
