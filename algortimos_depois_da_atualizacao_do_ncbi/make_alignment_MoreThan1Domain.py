import glob 
import os

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/system_selected_afer_blast_MoreThan1Domain.txt", "r")
for system in sis:
    system.split()
    system = system[0:4]
    ref.append(system)
for system in ref:
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files".format(system), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)
    for chain in chains:
        os.system("cd /home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ && clustalw -infile=sequences_{}_{}_MoreThan1Domain.fas  -type=protein -OUTORDER=INPUT -OUTFILE=sequences_{}_{}_MoreThan1Domain.fasta".format(system, system, chain, system, chain))
