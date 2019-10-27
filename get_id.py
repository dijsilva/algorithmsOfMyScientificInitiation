import os
from diego_tools import blast_out


ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/domains_approved.txt", "r")

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
    path_ncbi_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt'.format(sist,sist,chains[0])
    path_tax_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt'.format(sist, sist, chains[0])
    path_ncbi_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_ncbi_info_{}_ch_{}.txt'.format(sist, sist, chains[1])
    path_tax_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_pre_tax_{}_ch_{}.txt'.format(sist, sist, chains[1])
    print('Coletando os ids do sistema {}'.format(sist))
    blast_A = blast_out(path_ncbi_A, path_tax_A)
    tuple_of_species_A = blast_A[6]
    blast_B = blast_out(path_ncbi_B, path_tax_B)
    tuple_of_species_B = blast_B[6]

    species_A = []
    species_B = []

    for lists in tuple_of_species_A:
        if lists[1] == 1:
            species_A.append(lists[0])
    for listsB in tuple_of_species_B:
        if listsB[1] == 1:
            species_B.append(listsB[0])
    conc = []
    for sp1 in species_A:
        for sp2 in species_B:
            if sp1 == sp2:
                conc.append(sp1)

    sequences_A = []
    sequences_B = []
    for species in conc:
        for i_lists, lists in enumerate(blast_A[3]):
            for j_sp, sp in enumerate(lists):
                if species == sp:
                    pos = (blast_A[2][i_lists])
                    seq_A = (blast_A[0][pos].split(' ')[0].split('|')[1])
                    sequences_A.append(seq_A)
        for i_lists, lists in enumerate(blast_B[3]):
            for j_sp, sp in enumerate(lists):
                if species == sp:
                    pos = (blast_B[2][i_lists])
                    seq_B = (blast_B[0][pos].split(' ')[0].split('|')[1])
                    sequences_B.append(seq_B)
    out_A = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_id_seq_and_species_{}_{}.txt'.format(sist, sist, chains[0]),'w')
    out_B = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/psi_blast_out/psi_id_seq_and_species_{}_{}.txt'.format(sist, sist, chains[1]),'w')
    a_temp = [] 
    for pos, seqs_A in enumerate(sequences_A):
        if seqs_A in a_temp:
            a_temp.append(seqs_A)
            out_A.write('{} {}|{}\n'.format(seqs_A, a_temp.count(seqs_A), conc[pos]))
        else:
            a_temp.append(seqs_A)
            out_A.write('{} 1|{}\n'.format(seqs_A, conc[pos]))
    out_A.close()
    b_temp = []
    for pos, seqs_B in enumerate(sequences_B):
        if seqs_B in b_temp:
            b_temp.append(seqs_B)
            out_B.write('{} {}|{}\n'.format(seqs_B, b_temp.count(seqs_B), conc[pos]))
        else:
            a_temp.append(seqs_B)
            out_B.write('{} 1|{}\n'.format(seqs_B, conc[pos]))
    out_B.close()
