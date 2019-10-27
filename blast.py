from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time 
import numpy as np
import os

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/domains_approved2.txt", "r")

for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
for sist in ref:
    path = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(sist)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    chains = (chain1, chain2)
    for chain in chains:
        try:
            for record in SeqIO.parse("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}".format(sist, sist, chain), "fasta"):
                inicio = time.time()
                print('Fazendo blasp para o arquivo fas_{}_ch_{}'.format(sist, chain))
                result = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=100000, service="psi")
                fim = time.time()
                print('{:.0f} min de duração'.format((fim - inicio)/ 60))
                print('\n')
                records = NCBIXML.read(result)
            pre_tax = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt".format(sist, sist, chain), "w")
            ncbi_info = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt".format(sist, sist, chain), "w")
            for aln in records.alignments:
                ncbi_info.write('{} {} {} {}\n'.format(aln.hit_id, aln.hsps[0].expect, aln.hsps[0].score, aln.hsps[0].identities))
                pre_tax.write('{}\n'.format(aln.hit_def))
        except:
            try:
                for record in SeqIO.parse("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}".format(sist, sist, chain), "fasta"):
                    inicio = time.time()
                    print('Fazendo blasp para o arquivo fas_{}_ch_{}'.format(sist, chain))
                    result = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=100000, service="psi")
                    fim = time.time()
                    print('{:.0f} min de duração'.format((fim - inicio)/ 60))
                    print('\n')
                    records = NCBIXML.read(result)
                pre_tax = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt".format(sist, sist, chain), "w")
                ncbi_info = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt".format(sist, sist, chain), "w")
                for aln in records.alignments:
                    ncbi_info.write('{} {} {} {}\n'.format(aln.hit_id, aln.hsps[0].expect, aln.hsps[0].score, aln.hsps[0].identities))
                    pre_tax.write('{}\n'.format(aln.hit_def))
            except:
                try:
                    for record in SeqIO.parse("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}".format(sist, sist, chain), "fasta"):
                        inicio = time.time()
                        print('Fazendo blasp para o arquivo fas_{}_ch_{}'.format(sist, chain))
                        result = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=100000, service="psi")
                        fim = time.time() 
                        print('{:.0f} min de duração'.format((fim - inicio)/ 60))
                        print('\n')
                        records = NCBIXML.read(result)
                    pre_tax = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt".format(sist, sist, chain), "w")
                    ncbi_info = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt".format(sist, sist, chain), "w")
                    for aln in records.alignments:
                        ncbi_info.write('{} {} {} {}\n'.format(aln.hit_id, aln.hsps[0].expect, aln.hsps[0].score, aln.hsps[0].identities))
                        pre_tax.write('{}\n'.format(aln.hit_def))
                except:
                    break
        pre_tax.close()
        ncbi_info.close()
