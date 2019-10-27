from Bio import SeqIO
import os
import glob

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/moreThan1Domain_estructureApproved_without4y61.txt", "r")

#output
out_A = open('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/system_selected_afer_blast_filtred_MoreThan1Domain.txt', 'w')

for system in sis:
    system.split()
    system = system[0:4]
    ref.append(system)
for system in ref:
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/".format(system), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)
    pdb_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}'.format(system, system, chains[0])
    pdb_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}'.format(system, system, chains[1])
    sequences = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}_MoreThan1Domain.txt'.format(system, system, chains[0])
    id_ = open(sequences, "r")

    lenght = {}
    for record in SeqIO.parse(pdb_A, 'fasta'):
        lenght[chains[0]] = len(record.seq)

    for record in SeqIO.parse(pdb_B, 'fasta'):
        lenght[chains[1]] = len(record.seq)

    cont = 0
    for line in id_:
        cont += 1

    if (cont < lenght[chains[0]]) or (cont < lenght[chains[1]]):
        pass
    else:
        out_A.write('{}\n'.format(system))
        print(lenght, cont, system)
out_A.close()
