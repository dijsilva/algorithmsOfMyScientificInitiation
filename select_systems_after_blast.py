from Bio import SeqIO
import os

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/domains_approved.txt", "r")
out_A = open('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/blast/verificando.txt', 'w')

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
    pdb_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}'.format(system, system, chains[0])
    pdb_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}'.format(system, system, chains[1])
    sequences = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}.txt'.format(system, system, chains[0])
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
