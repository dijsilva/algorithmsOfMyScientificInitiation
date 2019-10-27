from Bio import AlignIO
from Bio import SeqIO
import os

####################Part of algorithm to make the subalignments##########################
ref = []
sis = open("/home/dsilva/Dropbox/lbtc/pibic_2018/arquivos/blast/systems_filtred_by_lenght_evalue.txt", "r")
for system in sis:
    system.split()
    system = system[0:4]
    ref.append(system)

for system in ref:
    path = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/pfam/full".format(system)
    arch = os.listdir(path)
    chain1 = arch[0]
    chain1 = chain1[-1]
    chain2 = arch[1]
    chain2 = chain2[-1]
    chains = (chain1, chain2)

    sequences_a = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_le_filtred.fasta".format(system, system, chains[0])
    sequences_b = "/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sequences_{}_{}_le_filtred.fasta".format(system, system, chains[1])

    alignment_sequence_a = AlignIO.read(sequences_a, "fasta")
    alignment_sequence_b = AlignIO.read(sequences_b, "fasta")

    sequence_a = alignment_sequence_a[0,:].seq
    sequence_b = alignment_sequence_b[0,:].seq

    positions_without_gap_a, positions_without_gap_b = [], []

    for pos, aa in enumerate(sequence_a):
        if aa not in '-._':
            positions_without_gap_a.append(pos)
        else:
            pass

    for pos, aa in enumerate(sequence_b):
        if aa not in '-._':
            positions_without_gap_b.append(pos)
        else:
            pass

    new_alignment_a = alignment_sequence_a[:,positions_without_gap_a[0]:positions_without_gap_a[1]]
    new_alignment_b = alignment_sequence_b[:,positions_without_gap_b[0]:positions_without_gap_b[1]]

    for position in positions_without_gap_a[1:]:
        new_alignment_a += alignment_sequence_a[:,position:position+1]
    for position in positions_without_gap_b[1:]:
        new_alignment_b += alignment_sequence_b[:,position:position+1]


    #these lines of code were for save the subalignments before select the sequences with more than 20 percent of gaps 

    #out_a = open("/home/dsilva/Desktop/sub_alignment_a.fas", "w")
    #out_b = open("/home/dsilva/Desktop/sub_alignment_b.fas", "w")
    #AlignIO.write(new_alignment_a, out_a, "fasta")
    #AlignIO.write(new_alignment_b, out_b, "fasta")
    #out_a.close()
    #out_b.close()

    #####################Part of algorithm to select the sequences with less than 20 percent of gaps##############################

    length_of_subalignment_a = len(new_alignment_a[0,:].seq)
    length_of_subalignment_b = len(new_alignment_b[0,:].seq)

    sequences_without_much_gaps_a, sequences_without_much_gaps_b = [], []

    for i in range(len(new_alignment_a)):
        counter_of_gaps = 0
        for aa in new_alignment_a[i,:].seq:
            if aa in "-._":
                counter_of_gaps += 1
            else:
                pass
        if ((counter_of_gaps / length_of_subalignment_a) < 0.20):
            sequences_without_much_gaps_a.append(i)
        else:
            pass

    for i in range(len(new_alignment_b)):
        counter_of_gaps = 0
        for aa in new_alignment_b[i,:].seq:
            if aa in "-._":
                counter_of_gaps += 1
            else:
                pass
        if ((counter_of_gaps / length_of_subalignment_b) < 0.20):
            sequences_without_much_gaps_b.append(i)
        else:
            pass


    sequences_without_much_gaps = []
    for index, position in enumerate(sequences_without_much_gaps_a):
        if position in sequences_without_much_gaps_b:
            sequences_without_much_gaps.append(position)
        else:
            pass

    idseqa, seqa = [], []

    for record in new_alignment_a:
        idseqa.append(record.id)
        seqa.append(record.seq)

    idseqb, seqb = [], []

    for recordb in new_alignment_b:
        idseqb.append(recordb.id)
        seqb.append(recordb.seq)

    newidseqa, newseqa = [], []
    newidseqb, newseqb = [], []

    for position, info in enumerate(sequences_without_much_gaps):
        newidseqa.append(idseqa[info])
        newseqa.append(seqa[info])
        newidseqb.append(idseqb[info])
        newseqb.append(seqb[info])

    out_a = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sub_{}_{}.fasta".format(system, system, chains[0]), "w")
    out_b = open("/home/dsilva/data_pibic_2018/proteins/2ch/{}/files/blast_out/sub_{}_{}.fasta".format(system, system, chains[1]), "w")

    for pos, id_a in enumerate(newidseqa):
        out_a.write('>{}\n{}\n'.format(id_a, str(newseqa[pos])))
    for pos, id_b in enumerate(newidseqb):
        out_b.write('>{}\n{}\n'.format(id_b, str(newseqb[pos])))
    out_a.close()
    out_b.close()

