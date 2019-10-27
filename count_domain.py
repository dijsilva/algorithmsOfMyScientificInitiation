import os
import numpy as np

path1 = "/home/dsilva/data_pibic_2018/proteins/2ch"
list_dir = os.listdir(path1)
#out = open("systems_with_more_than_one_domain.txt", "w")
listOfSystems_1, listOfSystems_2 = [], []
for directory in list_dir:
    path2 = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(directory)
    arch = os.listdir(path2)
    dict_domain = {}
    numberOfDomains = []
    for n in range(len(arch)):
        domain = (arch[n])
        domain = domain[-1]
        if domain not in dict_domain:
            dict_domain[domain] = 1
        else:
            dict_domain[domain] = dict_domain[domain] + 1
    print(directory, dict_domain)
    if len(dict_domain) == 2:
        values = (dict_domain.values())
        for b, c in enumerate(values):
            numberOfDomains.append(c)
        if (numberOfDomains[0] == numberOfDomains[1]) and numberOfDomains[0]  == 1:
            listOfSystems_2.append(directory)
        else:
            listOfSystems_1.append(directory)
            #out.write("{} \n".format(directory))
    elif len(dict_domain) >= 0 and len(dict_domain) != 2:
        listOfSystems_1.append(directory)
        #out.write("{} \n".format(directory))
#out.close()
