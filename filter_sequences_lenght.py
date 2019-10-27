from Bio import SeqIO
import numpy as np
import os

min_ = 0.7
max_ = 1.3

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/blast/systems_first_alignment.txt", "r")

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
    average = {}
    for chain in chains:
        file_ = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e_filtred.fas'.format(system, system, chain), 'r')
        tot = 0
        cont = 0
        for record in SeqIO.parse(file_, "fasta"):
            tot += len(record.seq)
            cont += 1
        media = int(tot / cont)
        average[chain] = media

    filea = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e_filtred.fas'.format(system, system, chains[0]), 'r')
    fileb = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e_filtred.fas'.format(system, system, chains[1]), 'r')
    idsa = []
    idsb = []
    seqsa = []
    seqsb = [] 
    newida = []
    newidb = []
    newseqa = []
    newseqb= []
    for record in SeqIO.parse(filea, "fasta"):
        idsa.append(record.id)
        seqsa.append(record.seq)

    for record in SeqIO.parse(fileb, "fasta"):
        idsb.append(record.id)
        seqsb.append(record.seq)

    for pos, i in enumerate(seqsa):
        if len(i) > int((average[chains[0]]*min_)) and len(i) < int((average[chains[0]]*max_)):
            newida.append(idsa[pos])
            newseqa.append(seqsa[pos])
            newidb.append(idsb[pos])
            newseqb.append(seqsb[pos])
    new2ida = []
    new2idb = []
    new2seqa = []
    new2seqb = []
    for pos, j in enumerate(newseqb):
        if len(j) > int((average[chains[1]]*min_)) and len(j) < int((average[chains[1]]*max_)):
            new2ida.append(newida[pos])
            new2seqa.append(newseqa[pos])
            new2idb.append(newidb[pos])
            new2seqb.append(newseqb[pos])
    out_a = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_el0.3_filtred.fas'.format(system, system, chains[0]), 'w')
    out_b = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_el0.3_filtred.fas'.format(system, system, chains[1]), 'w')
    if idsa[0] not in new2ida:
        out_a.write('>{}\n{}\n'.format(idsa[0], seqsa[0]))
    if idsb[0] not in new2idb:
        out_b.write('>{}\n{}\n'.format(idsb[0], seqsb[0]))
    for pos, id_a in enumerate(new2ida):
        out_a.write('>{}\n{}\n'.format(id_a, str(new2seqa[pos])))
    out_a.close()
    for pos, id_b in enumerate(new2idb):
        out_b.write('>{}\n{}\n'.format(id_b, str(new2seqb[pos])))
    out_b.close()
