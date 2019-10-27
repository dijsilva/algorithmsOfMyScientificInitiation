import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from diego_tools import blast_out
import os
import numpy as np

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/psi_blast/first_systems_select_after_psi_blast.txt", "r")
tot_species = []
for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
tot_a = []
tot_i = []
for sist in ref:
    path = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(sist)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    chains = (chain1, chain2)
    path_ncbi_1 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt'.format(sist,sist,chains[0])
    path_tax_1 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt'.format(sist, sist, chains[0])
    path_ncbi_2 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt'.format(sist, sist, chains[1])
    path_tax_2 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt'.format(sist, sist, chains[1])


    blast_A = blast_out(path_ncbi_1, path_tax_1)
    tuple_of_species_A = blast_A[6]
    blast_I = blast_out(path_ncbi_2, path_tax_2)
    tuple_of_species_I = blast_I[6]

    cont_A = {}
    cont_I = {}

    for lists in tuple_of_species_A:
        if lists[1] not in cont_A:
            cont_A[lists[1]] = 1
        else:
            cont_A[lists[1]] = cont_A[lists[1]] + 1

    for listsI in tuple_of_species_I:
        if listsI[1] not in cont_I:
            cont_I[listsI[1]] = 1
        else:
            cont_I[listsI[1]] = cont_I[listsI[1]] + 1
    species_A = []
    species_I = []

    for lists in tuple_of_species_A:
        if lists[1] == 1:
            species_A.append(lists[0])
    for listsI in tuple_of_species_I:
        if listsI[1] == 1:
            species_I.append(listsI[0])
    cont_species = 0
    for sp1 in species_A:
        for sp2 in species_I:
            if sp1 == sp2:
                cont_species += 1
    tot_species.append(cont_species)
    tot_a.append(len(species_A))
    tot_i.append(len(species_I))

_ref = np.arange(len(ref))
plt.figure(figsize=(38.0, 8.0))
plt.rcParams.update({'font.size': 18})
ax = plt.subplot(111)
w = 0.3
ax.bar(_ref - w, tot_a,width=w,color='k',align='center', label='k')
ax.bar(_ref, tot_i,width=w,color='r',align='center', label='r')
ax.bar(_ref + w, tot_species,width=w,color='gray',align='center', label='gray')
plt.xticks(_ref, ref)
plt.xlabel('Sistemas')
plt.ylabel('Número de espécies')
plt.legend(["Proteína A", "Proteína B", "Em comum nas duas proteínas"])
plt.savefig('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/plots_lbtc/psi_plots_all_systems_alignment.png', dpi=300) 
plt.show()
