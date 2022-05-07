import sys
import re
from Bio import SeqIO
aa_dict = {	
				'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A',                               # Alanine
				'TGC' : 'C', 'TGT' : 'C',                                                           # Cysteine
				'GAC' : 'D', 'GAT' : 'D',                                                           # Aspartic Acid
				'GAA' : 'E', 'GAG' : 'E',                                                           # Glutamic Acid
				'TTC' : 'F', 'TTT' : 'F',                                                           # Phenylalanine
				'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G',                               # Glycine
				'CAC' : 'H', 'CAT' : 'H',                                                           # Histidine
				'ATA' : 'I', 'ATC' : 'I', 'ATT' : 'I',                                             # Isoleucine
				'AAA' : 'K', 'AAG' : 'K',                                                           # Lysine
				'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L', 'TTA' : 'L', 'TTG' : 'L',   # Leucine
				'ATG' : 'M',                                                                         # Methionine
				'AAC' : 'N', 'AAT' : 'N',                                                           # Asparagine
				'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P',                               # Proline
				'CAA' : 'Q', 'CAG' : 'Q',                                                           # Glutamine
				'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R', 'AGA' : 'R', 'AGG' : 'R',   # Arginine
				'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S', 'AGC' : 'S', 'AGT' : 'S',   # Serine
				'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T',                               # Threonine
				'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V',                               # Valine
				'TGG' : 'W',                                                                         # Tryptophan
				'TAC' : 'Y', 'TAT' : 'Y',                                                           # Tyrosine
				'TAA' : 'U', 'TAG' : 'U', 'TGA' : 'U'                                              # Stop
				}

def cds2aa(seq):
	seq = str(seq)
	length = len(seq)
	aa_seq = ""
	for i in range(0,(length - 1),3):
		start = i
		end = i + 3
		codon = seq[start:end].upper()
		if codon in  aa_dict.keys():
			aa = aa_dict[f"{codon}"]
		else: 
			aa = "N"
		aa_seq = aa_seq + aa
	return aa_seq



fa_dict = SeqIO.to_dict(SeqIO.parse(f"{sys.argv[1]}","fasta"))
for key, val in fa_dict.items():
	print(">" + key)
	print(cds2aa(val.seq))
