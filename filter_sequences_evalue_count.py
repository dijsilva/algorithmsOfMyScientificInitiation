from Bio import SeqIO
import os


cutoff = 0.00001

ref = []

lenComplex, qtdSeq, Complex = [], [], []

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
    in_a = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e_filtred.fas'.format(system, system, chains[0])
    in_b = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_e_filtred.fas'.format(system, system, chains[1])
    seqA, seqB = [], []
    for record in SeqIO.parse(in_a, 'fasta'):
        seqA.append(record)
    qtdSeq.append(len(seqA))
    for record in SeqIO.parse(in_b, 'fasta'):
        seqB.append(record)

    lenComplex.append(len(seqA[0].seq) + len(seqB[0].seq))
    Complex.append(system)

    print('{} {} {} {} {}'.format(system, len(seqA), len(seqA[0].seq), len(seqB[0].seq), (len(seqA[0].seq) + len(seqB[0].seq))))
