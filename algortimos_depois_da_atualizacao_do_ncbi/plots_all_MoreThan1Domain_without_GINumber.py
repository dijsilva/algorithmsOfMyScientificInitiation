import matplotlib.pyplot as plt
from matplotlib import gridspec
import sys 
sys.path.insert(0, '/home/dsilva/Dropbox/lbtc/pibic_2018/scripts')
from diego_tools import blast_out_Normal

from diego_tools_withoutGINumber import blast_out as blast_out_WithoutGiNumber
import os
import glob

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/systems_after_blast_MoreThan1Domain_withoutGINumber.txt", "r")

for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
#mudança na quantidade de sequencias porque o computador draco nao conseguiu fazer o download de todas as sequencias
ref = ref[0:244]
for sist in ref:
    #coletar as cadeias de cada sistema
    system = sist[0:4]
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files".format(system), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)
    path_ncbi_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt'.format(system,system,chains[0])
    path_tax_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt'.format(system, system, chains[0])
    path_ncbi_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt'.format(system, system, chains[1])
    path_tax_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt'.format(system, system, chains[1])

    if (len(sist) > 4):
        if ((sist[-2] == '_') and (sist[-1] == chain1)):
            blast_A = blast_out_WithoutGiNumber(path_ncbi_A, path_tax_A)
            tuple_of_species_A = blast_A[6]
            blast_B = blast_out_Normal(path_ncbi_B, path_tax_B) 
            tuple_of_species_B = blast_B[6]
        elif ((sist[-2] == '_') and sist[-1] == chain2):
            blast_A = blast_out_Normal(path_ncbi_A, path_tax_A)
            tuple_of_species_A = blast_A[6]
            blast_B = blast_out_WithoutGiNumber(path_ncbi_B, path_tax_B) 
            tuple_of_species_B = blast_B[6]
        else:
            print('Algum problema na formatação do sistema {}'.format(system))
    if (len(sist) == 4):
        blast_A = blast_out_WithoutGiNumber(path_ncbi_A, path_tax_A)
        tuple_of_species_A = blast_A[6]
        blast_B = blast_out_WithoutGiNumber(path_ncbi_B, path_tax_B) 
        tuple_of_species_B = blast_B[6]

    cont_A = {}
    cont_B = {}

    for lists in tuple_of_species_A:
        if lists[1] not in cont_A:
            cont_A[lists[1]] = 1
        else:
            cont_A[lists[1]] = cont_A[lists[1]] + 1

    for listsB in tuple_of_species_B:
        if listsB[1] not in cont_B:
            cont_B[listsB[1]] = 1
        else:
            cont_B[listsB[1]] = cont_B[listsB[1]] + 1
    species_A = []
    species_B = []

    for lists in tuple_of_species_A:
        if lists[1] == 1:
            species_A.append(lists[0])
    for listsB in tuple_of_species_B:
        if listsB[1] == 1:
            species_B.append(listsB[0])
    con = 0
    for sp1 in species_A:
        for sp2 in species_B:
            if sp1 == sp2:
                con += 1

    fig = plt.figure(figsize=(12.0, 9.0))
    plt.rcParams.update({'font.size': 14})
    #fig.suptitle("{}".format(system), fontsize=24)
    gs = gridspec.GridSpec(2,6)

    ax1 = fig.add_subplot(gs[0,:4])
    #ax1.set_title("Distribuição das espécies no blast")
    ax1.plot(cont_A.keys(), cont_A.values(), 'k+', label='Proteína A') #{}'.format(chains[0]))
    ax1.set_ylabel('Número de espécies')

    ax2 = fig.add_subplot(gs[1,:4])
    ax2.plot(cont_B.keys(), cont_B.values(), 'r+', label='Proteína B') #{}'.format(chains[1]))
    ax2.set_xlabel('Quantidade de vezes que as espécies se repetem')
    ax2.set_ylabel('Número de espécies')
    ax1.legend(loc=1)
    ax2.legend(loc=1)

    ax3 = fig.add_subplot(gs[:,5])
    ax3.set_ylabel('Número de espécies')
    y = [len(species_A), len(species_B), con]
    x = (['A', 'B', 'A ∩ B'])#.format(chains[0], chains[1])]
    #x = ['{}'.format(chains[0]), '{}'.format(chains[1]), 'A ∩ B')#.format(chains[0], chains[1])]
    barlist = ax3.bar(x, y)
    barlist[0].set_color('k')
    barlist[1].set_color('r')
    barlist[2].set_color('gray')
    plt.savefig('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/plots_lbtc/MoreThan1Domain/{}_MoreThan1Domain.svg'.format(system), dpi=300)
    #plt.show()
