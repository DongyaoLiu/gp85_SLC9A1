import sys
import re
import copy
GenBank = sys.argv[1] 
mode = str(f"{sys.argv[2]}")

gb = open(f"{GenBank}","r")
aa_dict = {'GCA':'A', 'GCC' :'A', 'GCG': 'A', 'GCT': 'A',                               # Alanine
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

def nuc2aa(nuli_string, start_codon):
	aa_out = ""
	start_codon = int(start_codon) - 1
	nuli_string = "".join(nuli_string.split()).upper()
	nuli_string = str(nuli_string)[start_codon:]
	if len(nuli_string) % 3 == 0:
		for i in range(0,len(nuli_string),3):
			target_nuli = nuli_string[i:i+3]
			if target_nuli in aa_dict.keys():
				aa = aa_dict[f"{target_nuli}"]
			else:
				aa ="N"
			aa_out = f"{aa_out}{aa}" 
		return aa_out
	else:
		return ""





pep_dict={}
line_flag = 0
trans_start = False 
for line in gb.readlines():
	if re.match("LOCUS",line):
		default_length = line.split()[2]
	line_flag = line_flag + 1
	line = line.strip()
	strain = re.search(r"subgroup[:]{0,1}\W([A-Z])\W|subgroup[:]{0,1}\WALV-([A-Z])\W|subgroup.-.([A-Z])|genotype: ALV-([A-Z])", line)
	date = re.search(r"/collection_date.+([0-9]{4})",line)
	country = re.search(r"country.+(China|Japan|USA|Egypt|India|Malaysia|Nigeria|Taiwan)",line)
	gp85 = re.search(r"(note|product|gene).+(gp85|GP85|SU)\W",line)
	gp37 = re.search(r"(note|product|gene).+(gp37)|(gp37)",line)
	start_end = re.search(r"mat_peptide\W+([0-9]{1,4})..([0-9]{3,4})",line)
	pep_start = re.search(r"translation=\"([A-Z]+)",line)
	pep_end = re.search(r"([A-Z]{0,})\"",line)
	CDS = re.search(r"codon_start.+([0-9]{1,2})",line)
	Orange = re.match("ORIGIN",line)
	nuli = re.search(r"[0-9]\W(.+)",line)
	if Orange:
		nuli_begin = "start"
		nuli_str = ""
		continue
	if "nuli_begin" in globals():
		if nuli:
			nuli_str =f"{nuli_str}{nuli.group(1)}"
#			print(nuli.group(1))
	if pep_start:
		Pep = pep_start.group(1)
		trans_start = "True"
	#	print(pep_start.group(1) ,"start")
		continue
	if pep_end:
		if trans_start == "True":
			trans_start = "False"
			Pep = f"{Pep}{pep_end.group(1)}"
	#		print(pep_end.group(1),"end")
	if trans_start == "True":
		Pep = f"{Pep}{line}"
	#	print(line)
	if "first_cds" not in globals():		
		if CDS:
			Codon = re.search(r"codon_start.+([0-9]{1,2})",line).group(1)
			first_cds = "find"
	if re.search(r"CDS.+[.>]([0-9]+)",line):
		length = re.search(r"CDS.+[.>]+([0-9]+)",line).group(1)
		if "Codon" in globals():
			GP85_annotation = [Codon, length]
		else:	
			GP85_annotation	= [1, length]
	if start_end:
		Start_end = start_end.groups()
	if gp85:
		GP85 = "gp85"
		if "Start_end" in globals():
			GP85_annotation = [Start_end[0],Start_end[1]]
			del Start_end
	if gp37:
		GP37 = "gp37"
		if "Start_end" not in globals():
			GP37_annotation =["NA", "NA"]
		else:
			GP37_annotation = [Start_end[0],Start_end[1]]
			del Start_end
	if country != None:
		Country = country
	#	print(Country.group(1))
	if re.search(r"LOCUS", line):
		LOCUS = line.split()[1]
#	#	print(LOCUS)
	if date != None:
		Date = date
	#	print(date.group(1))
	#if re.match(r"VERSION\W+", line):
	#	version = re.match(r"VERSION\W+([A-Za-z0-9.]+)\W", line).group(1)	
	#	print(version)
	if strain:
		strain_info = strain	
	if re.match("//",line):
		del nuli_begin
		if "Codon" not in globals():
			Codon = 1
		if "Pep" not in globals():
			Pep = "NA"
		if "Country" not in globals():
			Country = re.match(r"(NA)", "NA")
		if "Date" not in globals():
			Date = re.match(r"(NA)", "NA")
		if "GP85" not in globals():
			GP85 = "NA"
			GP85_annotation = ["NA", "NA"]			
		if "GP37" not in globals():
			GP37 = "NA"
			GP37_annotation =["NA", "NA"]
		if "GP85_annotation" not in globals():	
			GP85_annotation = [1 ,default_length]
		if "first_cds" in globals():
			del first_cds
		else:
			Pep = ""
#			print(LOCUS, "not")
		if "strain_info" in globals(): 
			for i in [1,2,3,4]:
				if strain_info != None:
					if strain_info.group(i) != None:
						print(">" + LOCUS, strain_info.group(i), Country.group(1),Date.group(1),length,GP85,
								GP85_annotation[0],GP85_annotation[1],GP37,GP37_annotation[0],GP37_annotation[1])
						del strain_info 
						del Country
						del Date
						del GP85
						del GP37
						del GP85_annotation
						del GP37_annotation
						break
		else:
			print(">" + LOCUS, "NA", Country.group(1),Date.group(1),length,GP85,
					GP85_annotation[0],GP85_annotation[1],GP37,GP37_annotation[0], GP37_annotation[1])
			del Country 
			del Date
			del GP85
			del GP37
			del GP85_annotation
			del GP37_annotation
		if mode == "faa":
			if Pep != "":
				print(Pep)
			else:
				print(nuc2aa(nuli_str,Codon))
		del nuli_str
		del Pep	
		del default_length
