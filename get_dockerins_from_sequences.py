
import re
from Bio import SeqIO
import pandas as pd
import os
import numpy as np
pd.set_option('display.max_colwidth',-1)
#Open the multiple sequence alignment
def get_dockerins_from_sequences(msa_file):
    msa = open("fungi_seqdb_fungi_hmm.sto",'r')
    entries = msa.readlines()
    entries = [e for e in entries if '[subseq from]' in e]

    #Now pull ID and subsequence and put into dictionary
    keys = []
    values = []

    for string in entries:
        parts = re.split(r'[/\s]',string)
        keys.append(parts[1].replace('sp','tr'))
        values.append(parts[2])
    subseq_dict =dict(zip(keys,values))
    #Parse fasta file with full length sequences
    result_seqs = list(SeqIO.parse("dockerin_containing_seqs_gutfungi_Uniprot.fasta","fasta"))

    #dockerin_seqs = []

    for s in result_seqs:
        try:
            indices = subseq_dict[s.id].split('-')
            #dockerin_seqs.append(s.seq[int(indices[0]):int(indices[1])])
            s.seq = s.seq[int(indices[0]):int(indices[1])]
            s.id = str(s.id) + "/residues_" + subseq_dict[s.id]
        except KeyError:
            print("No key " + s.id + " in multiple sequence alignment.")

    SeqIO.write(result_seqs, "dockerin_seqs.fasta", "fasta")

def parse_groups(group_analysis, idtable):
    groups = pd.read_csv(group_analysis,sep='\t')
    ids = pd.read_csv(idtable,sep='\t')
    #Make a list of the sequence names that we can add to the groups table then write to read_csv
    seq_names = []
    for node in groups['Node']:
        id = str(ids.loc[ids['Id'] == node]['Label'])
        #Isolate only part that we want
        id2 = id[id.find('tr|'):id.find('OX=')]
        #print(id)
        seq_names.append(id2)

    #print(seq_names)
    groups.insert(2,"ProteinName",seq_names, True)
    groups.to_csv('~/Bioinformatics/dockerin_work/20191009_122656/analysis/dockerin_clusters_cutoffEval15.txt',sep='\t',index=False)


def compile_subgroup_statistics(newgroupfile):
    new_groups = pd.read_csv(newgroupfile,sep='\t')
    #Take top 10 largest clusters (listed in descending order of size)
    '''
    hemicellulases = [[],[],[],[],[],[],[],[],[],[]]
    cellulases = [[],[],[],[],[],[],[],[],[],[]]
    esterases = [[],[],[],[],[],[],[],[],[],[]]
    plyases = [[],[],[],[],[],[],[],[],[],[]]
    coths = [[],[],[],[],[],[],[],[],[],[]]
    mannanases = [[],[],[],[],[],[],[],[],[],[]]
    phosphatases = [[],[],[],[],[],[],[],[],[],[]]
    cbms = [[],[],[],[],[],[],[],[],[],[]]
    serp_lect_exp = [[],[],[],[],[],[],[],[],[],[]]
    uncharacterized = [[],[],[],[],[],[],[],[],[],[]]
    others = [[],[],[],[],[],[],[],[],[],[]]
    '''
    hemicellulases = np.zeros(77)
    cellulases = np.zeros(77)
    esterases = np.zeros(77)
    plyases = np.zeros(77)
    coths = np.zeros(77)
    mannanases = np.zeros(77)
    phosphatases = np.zeros(77)
    cbms = np.zeros(77)
    serp_lect_exp = np.zeros(77)
    uncharacterized = np.zeros(77)
    others = np.zeros(77)
    neocallimastix = np.zeros(77)
    piromyces = np.zeros(77)
    anaeromyces = np.zeros(77)
    caecomyces = np.zeros(77)
    other_species = np.zeros(77)
    for i in range(1,77):
        entries = new_groups.loc[new_groups['Subgroup'] == i]['ProteinName']
        for cell in entries:
            #Add write dockerin sequence to file in FASTA format for HMM creation

            if re.search('xyl',cell,re.IGNORECASE) or re.search('galact',cell,re.IGNORECASE) or re.search('arab',cell,re.IGNORECASE):
                hemicellulases[i] += 1
            elif re.search('glu',cell,re.IGNORECASE) or re.search('cellulase',cell,re.IGNORECASE) or re.search('glycosid',cell,re.IGNORECASE) or re.search('cello',cell,re.IGNORECASE):
                cellulases[i] += 1
            elif re.search('ester',cell,re.IGNORECASE) or re.search('Alpha/beta-hydrolase',cell,re.IGNORECASE):
                esterases[i] += 1
            elif re.search('lyase',cell,re.IGNORECASE):
                plyases[i] += 1
            elif re.search('Coth',cell,re.IGNORECASE):
                coths[i] += 1
            elif re.search('mann',cell,re.IGNORECASE) or re.search('hydrolase family 26',cell,re.IGNORECASE):
                mannanases[i] += 1
            elif re.search('phosph',cell,re.IGNORECASE):
                phosphatases[i] += 1
            elif re.search('binding',cell,re.IGNORECASE):
                cbms[i] += 1
            elif re.search('serpin',cell,re.IGNORECASE) or re.search('subtilase',cell,re.IGNORECASE) or re.search('lectin',cell,re.IGNORECASE) or re.search('swollenin',cell,re.IGNORECASE) or re.search('expansin',cell,re.IGNORECASE):
                serp_lect_exp[i] += 1
            elif re.search('Uncharacterized',cell,re.IGNORECASE):
                uncharacterized[i] += 1
            else:
                others[i] += 1
        #Loop through again to count sequences from each species
            if re.search('neocallimastix',cell,re.IGNORECASE):
                neocallimastix[i] += 1
            elif re.search('piromyces',cell,re.IGNORECASE):
                piromyces[i] += 1
            elif re.search('anaeromyces',cell,re.IGNORECASE):
                anaeromyces[i] += 1
            elif re.search('caecomyces',cell,re.IGNORECASE):
                caecomyces[i] += 1
            else:
                other_species[i] += 1

        with open("20191009_122656/analysis/cutoffEval15_subgroup_stats.txt",'a+') as file:
            file.write("#################################################\nSubgroup %i:\nHemicellulases: %i\nCellulases: %i\nEsterases: %i\nPolysaccharide lyases: %i\nCotH proteins: %i\n\
Mannanases: %i\nPhosphoesterases/phosphatases: %i\nCBMs: %i\nSerpins/Lectins/Expansins/Subtilases: %i\nUncharacterized: %i\n\
Other proteins: %i\nNeocallimastix: %i\nPiromyces: %i\nAnaeromyces: %i\nCaecomyces: %i\nOther species: %i\n" %(i,hemicellulases[i] \
        ,cellulases[i],esterases[i],plyases[i],coths[i],mannanases[i],phosphatases[i],cbms[i],\
serp_lect_exp[i],uncharacterized[i],others[i],neocallimastix[i],piromyces[i],anaeromyces[i],caecomyces[i],other_species[i]))

    #Add metrics for entire sequence database
    file = open("20191009_122656/analysis/cutoffEval15_subgroup_stats.txt",'a+')
    file.write("#################################################\nTotal database:\nHemicellulases: %i\nCellulases: %i\nEsterases: %i\nPolysaccharide lyases: %i\nCotH proteins: %i\n\
Mannanases: %i\nPhosphoesterases/phosphatases: %i\nCBMs: %i\nSerpins/Lectins/Expansins/Subtilases: %i\nUncharacterized: %i\n\
Other proteins: %i\nNeocallimastix: %i\nPiromyces: %i\nAnaeromyces: %i\nCaecomyces: %i\nOther species: %i\n" %(np.sum(hemicellulases) \
,np.sum(cellulases),np.sum(esterases),np.sum(plyases),np.sum(coths),np.sum(mannanases),np.sum(phosphatases),np.sum(cbms),\
np.sum(serp_lect_exp),np.sum(uncharacterized),np.sum(others),np.sum(neocallimastix),np.sum(piromyces),np.sum(anaeromyces),np.sum(caecomyces),np.sum(other_species)))

def construct_files_for_HMM(newgroupfile,dock_fasta_path,subgroup):
    new_groups = pd.read_csv(newgroupfile,sep='\t')
    entries = new_groups.loc[new_groups['Subgroup'] == subgroup]['ProteinName']
    dock_fasta = open(dock_fasta_path,'r')
    dock_seqs = dock_fasta.read()
    dock_fasta.close()

    newfilename = "20191009_122656/analysis/subfamilies/" + "Evalcutoff15_dockerinsubgroup" + str(subgroup) + ".fasta"
    newfile = open(newfilename,'a+')
    for name in entries:
        idx1 = dock_seqs.find(name)
        idx2 = dock_seqs.find('>',idx1)
        seq_to_write = dock_seqs[idx1-1:idx2]
        newfile.write(seq_to_write)



def main():
    #msa_file = "fungi_seqdb_fungi_hmm.sto"
    #group_file = "20191009_122656/analysis/GROUPS_10_Info.txt"
    #idtable = "20191009_122656/NodeTable_cs.txt"
    #get_dockerins_from_sequences(msa_file)
    #parse_groups(group_file,idtable)
    file = "20191009_122656/analysis/dockerin_clusters_cutoffEval15.txt"
    fasta_path = "dockerin_seqs.fasta"
    #compile_subgroup_statistics(file)
    for i in range(1,11):
        construct_files_for_HMM(file,fasta_path,i)


if __name__ == "__main__":
    main()
