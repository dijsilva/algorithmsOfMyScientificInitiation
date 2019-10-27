import os
import sys

sys.path.insert(0, '/home/dsilva/Dropbox/lbtc/pibic_2018/scripts')
from diego_tools import blast_out_Normal

from diego_tools_withoutGINumber import blast_out as blast_out_WithoutGiNumber
import glob


ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/systems_after_blast_MoreThan1Domain_withoutGINumber.txt", "r")

for s in sis:
    s = s.split("\n")[0]
    ref.append(s)
for sist in ref:
    #coletar as cadeias de cada sistema
    system = sist[0:4]
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files".format(system), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)
    path_ncbi_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt'.format(system,system,chains[0])
    path_tax_A = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt'.format(system, system, chains[0])
    path_ncbi_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt'.format(system, system, chains[1])
    path_tax_B = '/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/pre_tax_{}_ch_{}_MoreThan1Domain.txt'.format(system, system, chains[1])
    print('Coletando os ids do sistema {}'.format(system))

    #armazenar as tuplas que contem a informação do nome da espécie a qtd de vezes que ela aparece
    if (len(sist) > 4):
        if ((sist[-2] == '_') and (sist[-1] == chain1)):
            blast_A = blast_out_WithoutGiNumber(path_ncbi_A, path_tax_A)
            tuple_of_species_A = blast_A[6]
            blast_B = blast_out_Normal(path_ncbi_B, path_tax_B) 
            tuple_of_species_B = blast_B[6]
        elif ((sist[-2] == '_') and sist[-1] == chain2):
            blast_A = blast_out_Normal(path_ncbi_A, path_tax_A)
            tuple_of_species_A = blast_A[6]
            blast_B = blast_out_WithoutGiNumber(path_ncbi_B, path_tax_B) 
            tuple_of_species_B = blast_B[6]
        else:
            print('Algum problema na formatação do sistema {}'.format(system))
    if (len(sist) == 4):
        blast_A = blast_out_WithoutGiNumber(path_ncbi_A, path_tax_A)
        tuple_of_species_A = blast_A[6]
        blast_B = blast_out_WithoutGiNumber(path_ncbi_B, path_tax_B) 
        tuple_of_species_B = blast_B[6]

    species_A = []
    species_B = []

    #saber quais espécies aparecem apenas uma vez no sistema
    for lists in tuple_of_species_A:
        if lists[1] == 1:
            species_A.append(lists[0])
    for listsB in tuple_of_species_B:
        if listsB[1] == 1:
            species_B.append(listsB[0])

    #conc é uma lista com as espécies que aparecem apenas uma vez e são comuns nos dois sistemas
    conc = []
    for sp1 in species_A:
        for sp2 in species_B:
            if sp1 == sp2:
                conc.append(sp1)
    sequences_A = []
    sequences_B = []
    for species in conc:
        #blast_A[3] é a lista com os organismos correspondentes as sequencias selecionadas
        for i_lists, lists in enumerate(blast_A[3]):
            for j_sp, sp in enumerate(lists):
                if species == sp:
                    #blast_A[2] é a lista com os códigos do ncbi
                    pos = (blast_A[2][i_lists])
                    #blast_A[0] é a lista com códigos do ncbi retirados do output do blast... seq_A então é o código do ncbi
                    seq_A = (blast_A[0][pos].split(' ')[0].split('|')[1])
                    sequences_A.append(seq_A)
        #mesmo procedimento feito para A é então feito para B
        for i_lists, lists in enumerate(blast_B[3]):
            for j_sp, sp in enumerate(lists):
                if species == sp:
                    pos = (blast_B[2][i_lists])
                    seq_B = (blast_B[0][pos].split(' ')[0].split('|')[1])
                    sequences_B.append(seq_B)
    out_A = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}_MoreThan1Domain.txt'.format(system, system, chains[0]),'w')
    out_B = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}_MoreThan1Domain.txt'.format(system, system, chains[1]),'w')

    #a lista a_temp é criada para controlar quais sequencias já foram inseridas no output
    a_temp = [] 
    for pos, seqs_A in enumerate(sequences_A):
        if seqs_A in a_temp:
            """o método .count() é usado para fornecer a informação de quantas vezes um código 
            aparece na lista já que um mesmo código pode parecer mais de uma vez, mas referente 
            a espécies diferentes"""
            a_temp.append(seqs_A)
            out_A.write('{} {}|{}\n'.format(seqs_A, a_temp.count(seqs_A), conc[pos]))
        else:
            a_temp.append(seqs_A)
            out_A.write('{} 1|{}\n'.format(seqs_A, conc[pos]))
    out_A.close()

    #a lista a_temp é criada para controlar quais sequencias já foram inseridas no output
    b_temp = []
    for pos, seqs_B in enumerate(sequences_B):
        if seqs_B in b_temp:
            """o método .count() é usado para fornecer a informação de quantas vezes um código 
            aparece na lista já que um mesmo código pode parecer mais de uma vez, mas referente 
            a espécies diferentes"""
            b_temp.append(seqs_B)
            out_B.write('{} {}|{}\n'.format(seqs_B, b_temp.count(seqs_B), conc[pos]))
        else:
            a_temp.append(seqs_B)
            out_B.write('{} 1|{}\n'.format(seqs_B, conc[pos]))
    out_B.close()
