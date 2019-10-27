import numpy as np
import random
x = 278
system = '5F5S'
genome = np.zeros((5, x))
lista = (list(range(x)))
for i in range(0,5):
    for pos, p in enumerate(lista):
        genome[i,pos] = p
    random.shuffle(lista)
np.save("/home/dsilva/data_pibic_2018/MI/{}_AB/files/genome.coev_and_5_rands_{}.npy".format(system, system.lower()), genome)
