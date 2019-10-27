import matplotlib.pyplot as plt 
import matplotlib
from matplotlib.pyplot import figure 
import pickle

#le o dicionario com os dados de taxonomia

file = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/codeUniprotAndTaxonomy.txt", "rb")
codeUniprotAndTaxonomy = pickle.load(file)

#dicionario com as informações de taxonomia e o codigo PDB
taxonomy = {}

for item in codeUniprotAndTaxonomy:
    system = item[0:4]
    uniprotCode = item[5:]
    if item[0:4] not in taxonomy:
        #escrito desta forma para selecionar apenas o primeiro termo sem a presença do ponto e vírgula
        taxonomy[item[0:4]] = codeUniprotAndTaxonomy[item][1][0:-1]
    else:
        if taxonomy[item[0:4]] == codeUniprotAndTaxonomy[item][1][0:-1]:
            pass
        else:
            taxonomy['{}'.format(item[0:4])] = '{}&{}'.format(taxonomy[item[0:4]], codeUniprotAndTaxonomy[item][1][0:-1])

#dicionario para contar a quantidade de cada um dos tipos de termos
taxonomyQuantitie = {}

for item in taxonomy:
    if taxonomy[item] not in taxonomyQuantitie:
        taxonomyQuantitie[taxonomy[item]] = 1
    else:
        taxonomyQuantitie[taxonomy[item]] += 1

figure(num=None, figsize=(22, 16))
matplotlib.rcParams['xtick.labelsize'] = 6
plt.bar(taxonomyQuantitie.keys(), taxonomyQuantitie.values())
plt.xticks(rotation=120)
plt.xlabel("Taxonomia")
plt.ylabel("Quantidade")
plt.show()
