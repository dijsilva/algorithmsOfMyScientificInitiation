import matplotlib.pyplot as plt
from matplotlib import gridspec
from tools import blast_out

path_ncbi_A = '/home/dsilva/pibic_2018/proteins/2ch/2few/files/blast_out/ncbi_info_2few_ch_A.txt'
path_tax_A = '/home/dsilva/pibic_2018/proteins/2ch/2few/files/blast_out/pre_tax_2few_ch_A.txt' 
path_ncbi_B = '/home/dsilva/pibic_2018/proteins/2ch/2few/files/blast_out/ncbi_info_2few_ch_B.txt'
path_tax_B = '/home/dsilva/pibic_2018/proteins/2ch/2few/files/blast_out/pre_tax_2few_ch_B.txt'

blast_A = blast_out(path_ncbi_A, path_tax_A)
tuple_of_species_A = blast_A[6]
blast_I = blast_out(path_ncbi_B, path_tax_B)
tuple_of_species_I = blast_I[6]

cont_A = {}
cont_I = {}

for lists in tuple_of_species_A:
    if lists[1] not in cont_A:
        cont_A[lists[1]] = 1
    else:
        cont_A[lists[1]] = cont_A[lists[1]] + 1

for listsI in tuple_of_species_I:
    if listsI[1] not in cont_I:
        cont_I[listsI[1]] = 1
    else:
        cont_I[listsI[1]] = cont_I[listsI[1]] + 1
species_A = []
species_I = []

for lists in tuple_of_species_A:
    if lists[1] == 1:
        species_A.append(lists[0])
for listsI in tuple_of_species_I:
    if listsI[1] == 1:
        species_I.append(listsI[0])
con = 0
for sp1 in species_A:
    for sp2 in species_I:
        if sp1 == sp2:
            con += 1

fig = plt.figure(figsize=(14.0, 8.0))
fig.suptitle("Sistema 2few", fontsize=16)
gs = gridspec.GridSpec(2,3)

ax1 = fig.add_subplot(gs[0,:2])
ax1.set_title("Distribuição das espécies no blast")
ax1.plot(cont_A.keys(), cont_A.values(), 'b+', label='Cadeia A')
ax1.set_ylabel('Número de espécies')

ax2 = fig.add_subplot(gs[1,:2])
ax2.plot(cont_I.keys(), cont_I.values(), 'r+', label='Cadeia B')
ax2.set_xlabel('N de vezes que uma espécie se repete')
ax2.set_ylabel('Número de espécies')
ax1.legend(loc=1)
ax2.legend(loc=1)

ax3 = fig.add_subplot(gs[:,2])
ax3.set_title('N de espécies que não se repetem no blast')
y = [len(species_A), len(species_I), con]
x = ['A', 'I', 'A ∩ I']
#aux3 = plt.subplot2grid(shape=(2,2), loc=(1,1), rowspan=1, colspan=2)
barlist = ax3.bar(x, y)
barlist[0].set_color('b')
barlist[1].set_color('r')
barlist[2].set_color('k')
plt.savefig('2few.png', dpi=300)
plt.show()
