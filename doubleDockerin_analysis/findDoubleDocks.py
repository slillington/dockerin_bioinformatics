"""
This code uses the list of dockerins identified from the HMM
with their index information to find proteins with consecutive dockerin repeats.
We want to see how the linkers between dockerin domains are conserved.

Step 1 should be to read dockerinIDs_plus_subseqs.txt and 
output a dictionary with the key being the sequence ID and
the value being a list of the sequence range to pull out.

Step 2 reads the dictionary and outputs double dockerin sequences that have 
double dockerins with <= 25 residues between end of dock 1 and start of dock 2.
Outputs these sequences in FASTA format to double_dockerins.fasta.
Also, collect some statistics on linker length:
    Plot a histogram of linker length.
    Compute AA composition of linker.

"""
import re
import numpy as np
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt 
from collections import Counter

#Matplotlib settings
SMALLEST = 16
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALLEST)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#Read the file
raw = open(r'C:\Users\steph\Box\OmalleyLabfiles\Data\bioinformatics_work\dockerin_work\dockerinIDs_plus_subseqs.txt','r')
entries = raw.readlines()
#strip everything but the ID and the subsequence
keys = ['start']
values = [['placeholder']]
i = 0
for entry in entries:
    if  not (entry.find('FUNG') == -1 and entry.find('PIRSE') == -1):
        start = entry.find('|')+1
        end = entry.find(' DE')
        #start2 = entry.find('(',end)
        #end2 = entry.find(')',start2)
        parts = [x.strip() for x in re.split(r'/',entry[start:end])]
        #idN = parts[0] + ' ' + entry[start2:end2+1]
        idN = parts[0]
        #print(idN)
        seq_num = tuple([int(x) for x in parts[1].split('-')])
        #print(parts)
        if not keys[-1] == idN: #Check if sequence ID already in list
            keys.append(idN)
            v = [] #empty list to add first entry to
            v.append(seq_num)
            values.append(v)
            i += 1
        else:
            #First check if subseq is already in list
            if seq_num in values[i]:
                continue 
            else:
                values[i].append(seq_num)

subseq_dict = dict(zip(keys[1:],values[1:])) #A dictionary with Uniprot IDs and lists of subsequence ranges
#print(subseq_dict)

#Step 2: Step through list of keys and grab double dockerin sequences. Output linker stats
result_seqs = SeqIO.to_dict(SeqIO.parse(r"C:\Users\steph\Box\OmalleyLabfiles\Data\bioinformatics_work\dockerin_work\PF02013_full_length_sequences_new.fasta","fasta"))
output_file = open("double_dockerins.fasta",'a')
output_seqs = []
linkers = []
#print(result_seqs)
for key in subseq_dict:
    seq_list = subseq_dict[key]
    if len(seq_list) == 1:
        continue
    else:
        #Compare adjacent tuples to see if they are close enough
        for i in range(1,len(seq_list)):
            if min(seq_list[i]) - max(seq_list[i-1]) <= 25: #They are close enough
                start = min(seq_list[i-1])
                startL = max(seq_list[i-1])
                end = max(seq_list[i])
                endL = min(seq_list[i])
                #print(start,end)
                seq2 = str(result_seqs[key].seq)
                #print(len(seq2))
                #Create entry to add to FASTA file
                output_seq = ">" + str(key) + "_" + str(start)+"-"+str(end)+"\n" + seq2[start:end] + "\n"
                output_file.write(output_seq)

                #Add linker to separate list
                linkers.append(seq2[startL:endL])

#Plot histogram of linker lengths
fig, ax = plt.subplots(2)
lengths = np.array([len(x) for x in linkers])
ax[0].hist(lengths,bins=25,color='k')
ax[0].set_xlabel("Linker length")
ax[0].set_ylabel("Count")


#Compute average linker sequence composition
fullstr = ''.join(elem for elem in linkers)
freqs = dict(Counter(fullstr))
total = sum(freqs.values())
for k in freqs.keys():
    freqs[k] = freqs[k] / total


#Plot bar chart of composition
ax[1].bar(freqs.keys(),freqs.values(),color='k')
ax[1].set_ylabel("Frequency in linker")
#Print out compositions
print("Linker amino acid composition")

for k in sorted(freqs, key=freqs.get, reverse=True):
    print(str(k) + ':' + str(freqs[k]))


plt.show()








    

