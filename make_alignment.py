import os.path
import os

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/blast/systems_filtred_by_lenght0.3_evalue.txt", "r")
for system in sis:
    system.split()
    system = system[0:4]
    ref.append(system)
for system in ref:
    path = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(system)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    chains = (chain1, chain2)
    for chain in chains:
        os.system("cd /home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ && mafft --globalpair --maxiterate 16 --inputorder sequences_{}_{}_el0.3_filtred.fas > sequences_{}_{}_el0.3_filtred.fasta".format(system, system, chain, system, chain))
