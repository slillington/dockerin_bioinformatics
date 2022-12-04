import re
from Bio import SeqIO

docks = SeqIO.parse("dockerin_seqs_HMMERserver.fa","fasta")
kinase_motif_seq = open("kinmotif_dockerins.fasta",'a')
s1 = re.compile("R.S")
s2 = re.compile("K.S")
s3 = re.compile("S.R")
s4 = re.compile("S.K")

i = 0
k = 0
for d in docks:
    j = 0
    a = []
    a.append(s1.findall(str(d.seq))) 
    a.append(s2.findall(str(d.seq)))
    a.append(s3.findall(str(d.seq)))
    a.append(s4.findall(str(d.seq)))
    j = [x for x in a if x]
    if len(j) > 0:
        print(d.id + " has %i PKC recognition motifs" %len(j))
        print(d.seq)
        SeqIO.write(d, kinase_motif_seq, "fasta")
        i += 1

    k += 1

print ("%i dockerins with >= 1 PKC substrate motif out of %i" %(i,k))

