This project includes code to pull dockerin sequences from the proteomes of
anaerobic fungi using the dockerin Hidden Markov Model and characterize them.

get_dockerins_from_sequences.py:
Takes a multiple sequence alignment of dockerin-containing proteins and
strips out the dockerin sequences.

Uses a Sequence Similarity Network (SSN) to cluster related dockerin sequences
into subfamilies (separate program, SSNpipe)

Analyzes the distributions of catalytic domain families within dockerin subfamilies
to identify potential connections between enzyme function and fused dockerin subfamily.

kinase_seq_check.py:
Searches dockerin sequences across anaerobic fungal proteome for putative CotH kinase
sequence motifs to determine whether dockerins are potential substrates of these kinases.
Identified motifs were mapped onto the dockerin structure to predict whether dockerin
phosphorylation could serve as a post-translational binding trigger.

findDoubleDocks.py:
Strips double dockerin sequences from a FASTA file of dockerin containing proteins
and computes some statistics on linker amino acid composition + linker length
distribution.

