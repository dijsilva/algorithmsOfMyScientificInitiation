import matplotlib.pyplot as plt
from matplotlib import gridspec
from diego_tools import blast_out
import os

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/domains_approved.txt", "r")
n_species = {}
for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
for sist in ref[105:120]:
    path = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(sist)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    chains = (chain1, chain2)
    path_ncbi_1 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}.txt'.format(sist,sist,chains[0])
    path_tax_1 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}.txt'.format(sist, sist, chains[0])
    path_ncbi_2 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}.txt'.format(sist, sist, chains[1])
    path_tax_2 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}.txt'.format(sist, sist, chains[1])


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
    cont_sp = 0
    for sp1 in species_A:
        for sp2 in species_I:
            if sp1 == sp2:
                cont_sp += 1

    n_species[str(sist)] = cont_sp

fig = plt.figure(figsize=(14.0, 8.0))
fig.suptitle("Quantidade de sequências", fontsize=16)
gs = gridspec.GridSpec(len(n_species), 1)
for n in range(len(n_species)):
    ax1 = fig.add_subplot(gs[:2])
    ax1.set_title("Distribuição das espécies no blast")
    ax1.plot(cont_A.keys(), cont_A.values(), 'b+', label='Chain {}'.format(chains[0]))
    ax1.set_ylabel('Número de espécies')

plt.show()
