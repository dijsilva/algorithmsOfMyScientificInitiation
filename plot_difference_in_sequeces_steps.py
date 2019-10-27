import matplotlib.pyplot as plt 
import os
from matplotlib.pyplot import figure

systems = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/blast/systems_filtred_by_lenght0.3_evalue.txt", "r")

for s in sis:
    s.split()
    s = s[0:4]
    systems.append(s)
for system in systems:
    path = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(system)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    chains = [chain1, chain2]


    firstSequences = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}.txt".format(system, system, chains[0]), "r")

    contFirst = 0
    for line in firstSequences:
        contFirst += 1

    
    secondSequences = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_l_filtred.fas".format(system, system, chains[0]), "r")

    #valor de -1 para compensar a sequencia original do pdb que foi colocada na primeira posição de todas os MSAs
    contSecond = (-1)
    for line in secondSequences:
        contSecond += 1
    contSecond = contSecond / 2


    thirdSequences = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_le_filtred.fas".format(system, system, chains[0]), "r")

    #valor de -1 para compensar a sequencia original do pdb que foi colocada na primeira posição de todas os MSAs
    thirdSecond = (-1)
    for line in thirdSequences:
        thirdSecond += 1
    thirdSecond = thirdSecond / 2
    

    secondSequences03 = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e_filtred.fas".format(system, system, chains[0]), "r")

    #valor de -1 para compensar a sequencia original do pdb que foi colocada na primeira posição de todas os MSAs
    contSecond03 = (-1)
    for line in secondSequences03:
        contSecond03 += 1
    contSecond03 = contSecond03 / 2


    thirdSequences03 = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_el0.3_filtred.fas".format(system, system, chains[0]), "r")

    #valor de -1 para compensar a sequencia original do pdb que foi colocada na primeira posição de todas os MSAs
    contThird03 = (-1)
    for line in thirdSequences03:
        contThird03 += 1
    contThird03 = contThird03 / 2

    x = ['Todas as sequências', 'Filtro pelo tamanho', 'Filtro do tamanho + e-value', 'Filtro e-value', 'Filtro e-value + tamanho0.3']
    y = [contFirst, contSecond, thirdSecond, contSecond03, contThird03]
    
    figure(num=None, figsize=(14, 9))
    plt.bar(x[0],y[0], color='black')
    plt.bar(x[1],y[1], color='blue')
    plt.bar(x[2],y[2], color='blue')
    plt.bar(x[3],y[3], color='red')
    plt.bar(x[4],y[4], color='red')
    plt.xlabel('Passos de filtragem')
    plt.ylabel('Número de sequências')
    plt.savefig('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/plots_lbtc/sequencias_perdidas/{}_usado_nos_primeiros_alinhamentos.svg'.format(system))
    plt.close()