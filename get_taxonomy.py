from urllib import request
import pickle
import os

#carrega os códigos PDB dos sistemas
pdbCodeFile = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/domains_approved.txt")
codeSystems = []
for line in pdbCodeFile:
    line = line.split()
    line = line[0]
    codeSystems.append(line)


#carrega o dicionário com os códigos do UniProt dos sistemas
uniprotCodeFile = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/codigos_uniprot_e_pdb.txt", "rb")
codesUniProt = pickle.load(uniprotCodeFile)
uniprotCodeFile.close()

#dicionario para armazenar os códigos do PDB_UniProt e a Taxonomia
codeUniprotAndTaxonomy = {}

#baixa as informações do banco de dados do UniProt usando a biblioteca url
for system in codesUniProt:
    codes = codesUniProt[system]
    for code in codes:
        #lista com as informações de taxonomia que serão baixadas
        taxonomyList = []
        #baixa as informações de taxonomia e cria um arquivo de log com essas informações
        data = request.urlretrieve("https://www.uniprot.org/uniprot/{}.txt".format(code), "log.txt")
        log = open('log.txt', 'r')
        for line in log:
            if line.startswith("OC"):
                line = (line.split()[1:])
                for item in line:
                    taxonomyList.append(item)
        #apaga o arquivo de log onde foram salvos as informações de taxonomia
        os.remove("log.txt")
        #grava a taxonomia no dicionario
        codeUniprotAndTaxonomy['{}_{}'.format(system, code)] = taxonomyList