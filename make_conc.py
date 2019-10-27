from Bio import SeqIO
import os

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/systems_filtred_by_lenght_evalue.txt", "r")

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

    sequence_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_le_filtred.fasta'.format(system, system, chains[0])
    sequence_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_le_filtred.fasta'.format(system, system, chains[1])

    id_A = []
    id_B = []
    seq_A = []
    seq_B = []

    for record in SeqIO.parse(sequence_A, 'fasta'):
        seq_A.append(record.seq)
        id_A.append(record.id)

    for record in SeqIO.parse(sequence_B, 'fasta'):
        seq_B.append(record.seq)
        id_B.append(record.id)

    conc_seq = []
    conc_id = []
    for pos, i in enumerate(id_A):
        conc_id.append(id_A[pos]+'__'+id_B[pos])
        conc_seq.append(seq_A[pos]+seq_B[pos])

    out_a  = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}{}_le_filtred.fasta'.format(system, system, chains[0], chains[1]), 'w')
    for pos, i in enumerate(conc_id):
        out_a.write('>{}\n{}\n'.format(conc_id[pos], conc_seq[pos]))
    out_a.close()
