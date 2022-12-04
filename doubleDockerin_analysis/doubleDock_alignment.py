#Pairwise alignment of double dockerin sequences with Biopython

# Import libraries
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO

seqs = SeqIO.parse(r"C:\Users\steph\Box\OmalleyLabfiles\Data\bioinformatics_work\dockerin_work\doubleDockerin_analysis\double_dockerins.fasta","fasta")
print(len(seqs))