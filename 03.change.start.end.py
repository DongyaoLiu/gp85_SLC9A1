from Bio import SeqIO
import re
import sys

fa_dict = SeqIO.to_dict(SeqIO.parse(f"/Users/liudongyao/Downloads/NWAFU/ALV/{sys.argv[1]}", "fasta"))
tblastn_file = open(f"/Users/liudongyao/Downloads/NWAFU/ALV/{sys.argv[2]}","r")


rm_list = ["KC841155","KF796649"]
for line in tblastn_file.readlines():
	line_list = line.strip().split("\t")
	name = line_list[1]
	if name in rm_list:
		continue
	start = int(line_list[8]) - 1
	end = int(line_list[9])
	seq = str(fa_dict[f"{name}"].seq)[start:end]
	print(">"+f"{name}")
	print(seq)
	
