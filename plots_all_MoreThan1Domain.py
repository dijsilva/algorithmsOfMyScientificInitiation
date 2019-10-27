import matplotlib.pyplot as plt
from matplotlib import gridspec
from diego_tools import blast_out_Normal
import os
import glob

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/moreThan1Domain_estructureApproved_without4y61.txt", "r")

for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
for sist in ref:
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files".format(sist), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)
    path_ncbi_1 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt'.format(sist,sist,chains[0])
    path_tax_1 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt'.format(sist, sist, chains[0])
    path_ncbi_2 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt'.format(sist, sist, chains[1])
    path_tax_2 = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt'.format(sist, sist, chains[1])


    blast_A = blast_out_Normal(path_ncbi_1, path_tax_1)
    tuple_of_species_A = blast_A[6]
    blast_I = blast_out_Normal(path_ncbi_2, path_tax_2)
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
    con = 0
    for sp1 in species_A:
        for sp2 in species_I:
            if sp1 == sp2:
                con += 1

    fig = plt.figure(figsize=(12.0, 9.0))
    plt.rcParams.update({'font.size': 14})
    #fig.suptitle("{}".format(sist), fontsize=24)
    gs = gridspec.GridSpec(2,6)

    ax1 = fig.add_subplot(gs[0,:4])
    #ax1.set_title("Distribuição das espécies no blast")
    ax1.plot(cont_A.keys(), cont_A.values(), 'k+', label='Proteína A') #{}'.format(chains[0]))
    ax1.set_ylabel('Número de espécies')

    ax2 = fig.add_subplot(gs[1,:4])
    ax2.plot(cont_I.keys(), cont_I.values(), 'r+', label='Proteína B') #{}'.format(chains[1]))
    ax2.set_xlabel('Quantidade de vezes que as espécies se repetem')
    ax2.set_ylabel('Número de espécies')
    ax1.legend(loc=1)
    ax2.legend(loc=1)

    ax3 = fig.add_subplot(gs[:,5])
    ax3.set_ylabel('Número de espécies')
    y = [len(species_A), len(species_I), con]
    x = (['A', 'B', 'A ∩ B'])#.format(chains[0], chains[1])]
    #x = ['{}'.format(chains[0]), '{}'.format(chains[1]), 'A ∩ B')#.format(chains[0], chains[1])]
    barlist = ax3.bar(x, y)
    barlist[0].set_color('k')
    barlist[1].set_color('r')
    barlist[2].set_color('gray')
    plt.savefig('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/plots_lbtc/MoreThan1Domain/{}_MoreThan1Domain.svg'.format(sist), dpi=300)
    #plt.show()
