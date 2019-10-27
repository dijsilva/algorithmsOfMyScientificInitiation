from Bio import SeqIO
import os


cutoff = 0.00001

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/domains_approved.txt", "r")

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
    print('{}'.format(system))
    ncbi_a = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}.txt'.format(system, system, chains[0]), 'r')
    ncbi_b = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}.txt'.format(system, system, chains[1]), 'r')
    ref_a = {}
    ref_b = {}
    for line in ncbi_a:
        val_a = (float(line.split(' ')[1]))
        id_a = str(line.split('|')[3])
        ref_a[id_a] = val_a
    for line in ncbi_b:
        val_b = (float(line.split(' ')[1]))
        id_b = str(line.split('|')[3])
        ref_b[id_b] = val_b

    sequences_a = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}.fas'.format(system, system, chains[0]), 'r')
    sequences_b = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}.fas'.format(system, system, chains[1]), 'r')

    id_a = []
    seq_a = []
    id_b = []
    seq_b = []
    for record in SeqIO.parse(sequences_a, "fasta"):
        id_a.append(record.id)
        seq_a.append(record.seq)

    for record in SeqIO.parse(sequences_b, "fasta"):
        id_b.append(record.id)
        seq_b.append(record.seq)

    newid_a = []
    newseq_a = []
    newid_b = []
    newseq_b = []
    for pos, i in enumerate(id_a):
        if pos > 0:
            if 'pir' not in i.split('(')[0]:
                if ref_a[i.split('(')[0]] < cutoff:
                    newid_a.append(id_a[pos])
                    newseq_a.append(seq_a[pos])
                    newid_b.append(id_b[pos])
                    newseq_b.append(seq_b[pos])
    new2id_a = []
    new2seq_a = []
    new2id_b = []
    new2seq_b = []
    for pos, j in enumerate(newid_b):
        if 'pir' not in j.split('(')[0]:
            if ref_b[j.split('(')[0]] < cutoff:
                new2id_a.append(newid_a[pos])
                new2seq_a.append(newseq_a[pos])
                new2id_b.append(newid_b[pos])
                new2seq_b.append(newseq_b[pos])

    out_a = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e03_filtred.fas'.format(system, system, chains[0]), 'w')
    out_b = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e03_filtred.fas'.format(system, system, chains[1]), 'w')
    out_a.write(">{}\n{}\n".format(id_a[0], seq_a[0]))
    out_b.write(">{}\n{}\n".format(id_b[0], seq_b[0]))
    for pos, id_ in enumerate(new2id_a):
        out_a.write(">{}\n{}\n".format(id_, new2seq_a[pos]))
    out_a.close()
    for pos, id_b in enumerate(new2id_b):
        out_b.write(">{}\n{}\n".format(id_b, new2seq_b[pos]))
    out_b.close()
