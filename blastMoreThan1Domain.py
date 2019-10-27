from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time 
import numpy as np
import os
import glob

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/RepeatMoreThan1Domain_estructureApproved.txt", "r")

for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
for sist in ref:
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files".format(sist), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)
    for chain in chains:
        try:
            for record in SeqIO.parse("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}".format(sist, sist, chain), "fasta"):
                inicio = time.time()
                print('Fazendo blasp para o arquivo fas_{}_ch_{}'.format(sist, chain))
                result = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=100000)
                fim = time.time()
                print('{:.0f} min de duração'.format((fim - inicio)/ 60))
                print('\n')
                records = NCBIXML.read(result)
            pre_tax = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain), "w")
            ncbi_info = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain), "w")
            for aln in records.alignments:
                ncbi_info.write('{} {} {} {}\n'.format(aln.hit_id, aln.hsps[0].expect, aln.hsps[0].score, aln.hsps[0].identities))
                pre_tax.write('{}\n'.format(aln.hit_def))
        except:
            try:
                for record in SeqIO.parse("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}".format(sist, sist, chain), "fasta"):
                    inicio = time.time()
                    print('Fazendo blasp para o arquivo fas_{}_ch_{}'.format(sist, chain))
                    result = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=100000)
                    fim = time.time()
                    print('{:.0f} min de duração'.format((fim - inicio)/ 60))
                    print('\n')
                    records = NCBIXML.read(result)
                pre_tax = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain), "w")
                ncbi_info = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain), "w")
                for aln in records.alignments:
                    ncbi_info.write('{} {} {} {}\n'.format(aln.hit_id, aln.hsps[0].expect, aln.hsps[0].score, aln.hsps[0].identities))
                    pre_tax.write('{}\n'.format(aln.hit_def))
            except:
                try:
                    for record in SeqIO.parse("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}".format(sist, sist, chain), "fasta"):
                        inicio = time.time()
                        print('Fazendo blasp para o arquivo fas_{}_ch_{}'.format(sist, chain))
                        result = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=100000)
                        fim = time.time() 
                        print('{:.0f} min de duração'.format((fim - inicio)/ 60))
                        print('\n')
                        records = NCBIXML.read(result)
                    pre_tax = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain), "w")
                    ncbi_info = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain), "w")
                    for aln in records.alignments:
                        ncbi_info.write('{} {} {} {}\n'.format(aln.hit_id, aln.hsps[0].expect, aln.hsps[0].score, aln.hsps[0].identities))
                        pre_tax.write('{}\n'.format(aln.hit_def))
                except:
                    break
        pre_tax.close()
        ncbi_info.close()
