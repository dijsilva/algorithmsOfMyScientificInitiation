from pyvmd import *
import os
import numpy as np
import math

"""
Algoritmo para análise dos sistemas e aplicação de critérios para
escolha dos sistemas a serem estudados.
"""

"""
criar a lista de sistemas com apenas 1 domínio para cada proteína
"""

ref = []
sis = open("/home/dsilva/pibic_2018/arquivos/sist_with_1domain.txt")
out = open("/home/dsilva/pibic_2018/arquivos/domains_approved.txt", "w")
incompat = open("/home/dsilva/pibic_2018/arquivos/sist_incomp.txt", "w")

"""
Criar uma lista com os sistemas a serem estudados
"""
for s in sis:
    s.split()
    s = s[0:4]
    ref.append(s)
ref = ref[0:(len(ref)-1)]

#saber o nome dos domínios porque cada sistema tem domínios diferentes
for sistem in ref:
    path = "/home/dsilva/pibic_2018/proteins/2ch/{}/files/pfam/full".format(sistem)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    
    mol = System(name="{}".format(sistem))
    mol.load("/home/dsilva/pibic_2018/proteins/2ch/{}/files/{}.ent.pdb".format(sistem, sistem))

    lenA = mol.selectAtoms("protein and chain {} and name CA".format(chain1))
    lenA = lenA["resid"]
    len1 = []
    for l in lenA:
        if l > 0:
            len1.append(l)
    lenB = mol.selectAtoms("protein and chain {} and name CA".format(chain2))
    lenB = lenB["resid"]
    len2 = []
    for l2 in lenB:
        if l2 > 0:
            len2.append(l2)
    matriz = np.zeros((len(len1),len(len2)))
    contR = 0
    #cria uma matriz com os valores de distância entre cada aminoácido
    for id_1 in len1:
        sel1 = mol.selectAtoms("protein and chain {} and resid {} and name CA".format(chain1, id_1))
        coord1 = sel1.coords
        coord1 = coord1[0:1,:]
        contC = 0
        for id_2 in len2:
            sel2 = mol.selectAtoms("protein and chain {} and resid {} and name CA".format(chain2, id_2))
            coord2 = sel2.coords
            coord2 = coord2[0:1,:]
            a = (float((coord2[0, 0] - coord1[0, 0])**2)) + (float((coord2[0,1] - coord1[0, 1])**2)) + (float((coord2[0, 2] - coord1[0, 2])**2))
            b = math.sqrt(a)
            matriz[contR,contC] = b
            contC += 1
        contR += 1
    #aplica um teste lógico para converter os valores da matriz em valores booleanos com critério de distância estabelecido
    test = np.where(matriz < 8.0, 1, 0)
    cont_ptnA = 0
    cont_ptnB = 0
    #conta as posições que atendam o critério de distancia da proteína B
    for x in range(len(len1)):
        a = np.sum(test[x:(x+1),:])
        if a >= 1:
            cont_ptnA += 1
    #conta as posições que atendam o critério de distância da proteína B
    for y in range(len(len2)):
        g = np.sum(test[:,y:(y+1)])
        if g >= 1:
            cont_ptnB += 1
    perc_A = ((cont_ptnA/(len(len1))*100))
    perc_B = ((cont_ptnB/(len(len2))*100))
    if perc_A > 50 or perc_B > 50 or len(len1) > (len(len2)*5) or len(len2) > (len(len1)*5):
        incompat.write("{} \n".format(sistem))
    else:
        out.write("{} \n".format(sistem))
    mol.delete()

out.close()
