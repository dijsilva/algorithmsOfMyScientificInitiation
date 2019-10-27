from Bio import Entrez
import os
import glob

ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/moreThan1Domain_estructureApproved_without4y61.txt", "r")

for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
for sist in ref:
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files".format(sist), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)

    ValueOfVerificA = os.path.isfile('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_MoreThan1Domain.fas'.format(sist, sist, chains[0]))
    ValueOfVerificB = os.path.isfile('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_MoreThan1Domain.fas'.format(sist, sist, chains[1]))
    
    verific = None

    if((ValueOfVerificA == True) and (ValueOfVerificB == True)):
        lenOfFileA = int(os.path.getsize('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_MoreThan1Domain.fas'.format(sist, sist, chains[0])))
        lenOfFileB = int(os.path.getsize('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_MoreThan1Domain.fas'.format(sist, sist, chains[1])))

        if((lenOfFileA > 0) and (lenOfFileB > 0)):
            verific = True
        else:
            verific = False
    else:
        verific = False

    if (verific == False):
        print('Sistema {} com arquivos ainda não gerados. Trabalhando neles agora.'.format(sist))
        file_A = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}_MoreThan1Domain.txt'.format(sist, sist, chains[0]), 'r')
        file_B = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/id_seq_and_species_{}_{}_MoreThan1Domain.txt'.format(sist, sist, chains[1]), 'r')
        
        print('Baixando as sequências do sistema {}'.format(sist))
        print('Baixando as sequências da cadeia {}'.format(chains[0]))
        ids_a = []
        for line in file_A:
            ids_a.append(line.split(' ')[0])

        out_sequences_A = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_MoreThan1Domain.fas'.format(sist, sist, chains[0]), 'w')

        pdb_A = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}'.format(sist, sist, chains[0]), 'r')
        pdb_A = pdb_A.read()
        sequence_pdb_A = pdb_A.split('\n')
        title_pdb_A = sequence_pdb_A[0:1][0]
        title_pdb_A = title_pdb_A.split('|')[0]
        seq_pdb_A = sequence_pdb_A[1:]
        seq_pdb_A = ''.join(seq_pdb_A)
        out_sequences_A.write('{}\n{}\n'.format(title_pdb_A, seq_pdb_A))
        
        a_temp = []   
        for id_A in ids_a:
            Entrez.email = "diegosilva.unb@gmail.com"
            handle = Entrez.efetch(db="protein", id=id_A, rettype="fasta", retmode="text")
            sequences = handle.read()
            sequences = sequences.split('\n')
            title = sequences[0:1][0]
            title = title.split(' ')[0]
            seq = sequences[1:]
            seq = ''.join(seq)
            if id_A in a_temp:
                a_temp.append(id_A)
                out_sequences_A.write('{}({})\n{}\n'.format(title, a_temp.count(id_A), seq))
            else:
                a_temp.append(id_A)
                out_sequences_A.write('{}(1)\n{}\n'.format(title, seq))
        out_sequences_A.close()
        
        print('Baixando as sequências da cadeia {}'.format(chains[1]))
        ids_b = []
        for line in file_B:
            ids_b.append(line.split('|')[0])
        out_sequences_B = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_MoreThan1Domain.fas'.format(sist, sist, chains[1]), 'w')

        pdb_B = open('/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/fas_{}_ch_{}'.format(sist, sist, chains[1]), 'r')
        pdb_B = pdb_B.read()
        sequence_pdb_B = pdb_B.split('\n')
        title_pdb_B = sequence_pdb_B[0:1][0]
        title_pdb_B = title_pdb_B.split('|')[0]
        seq_pdb_B = sequence_pdb_B[1:]
        seq_pdb_B = ''.join(seq_pdb_B)
        out_sequences_B.write('{}\n{}\n'.format(title_pdb_B, seq_pdb_B))
        
        b_temp = []
        for id_B in ids_b:
            Entrez.email = "diegosilva.unb@gmail.com"
            handle = Entrez.efetch(db="protein", id=id_B, rettype="fasta", retmode="text")
            sequences = handle.read()
            sequences = sequences.split('\n')
            title = sequences[0:1][0]
            title = title.split(' ')[0]
            seq = sequences[1:]
            seq = ''.join(seq)
            if id_B in b_temp:
                b_temp.append(id_B)
                out_sequences_B.write('{}({})\n{}\n'.format(title, b_temp.count(id_B), seq))
            else:
                b_temp.append(id_B)
                out_sequences_B.write('{}(1)\n{}\n'.format(title, seq))
        out_sequences_B.close()
    else:
        print('Tudo certo com o sistema {}, partindo para o próximo.\n'.format(sist))
        pass
