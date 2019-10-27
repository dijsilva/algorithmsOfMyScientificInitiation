import os
import glob

sis = open('/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/moreThan1Domain_estructureApproved_without4y61.txt', 'r')
ref = []


for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
controllerOfSystems = []
for sist in ref:
    fastas = glob.glob1("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/".format(sist), 'fas*')
    chain1 = fastas[0][-1]
    chain2 = fastas[1][-1]
    chains = (chain1, chain2)

    ncbiChainA = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain1), "r")
    ncbiChainB = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/ncbi_info_{}_ch_{}_MoreThan1Domain.txt".format(sist, sist, chain2), "r")
    print('Analisando {}'.format(sist))
    controllerOfSystemsTemporary = []
    for line in ncbiChainA:
        if line.startswith('gi'):
            pass
        else:
            controllerOfSystemsTemporary.append('{}_{}'.format(sist, chain1))
            break
    for line in ncbiChainB:
        if line.startswith('gi'):
            pass
        else:
            if ('{}_{}'.format(sist,chain1) in controllerOfSystemsTemporary):
                controllerOfSystems.append('{}'.format(sist))
                controllerOfSystemsTemporary = []
                break
            else:
                controllerOfSystemsTemporary.append('{}_{}'.format(sist, chain2))
    if (len(controllerOfSystemsTemporary) > 0):
        controllerOfSystems.append(controllerOfSystemsTemporary[0])
    print('{} Analisado'.format(sist))
print(controllerOfSystems)
