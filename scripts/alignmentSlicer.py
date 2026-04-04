from Bio import AlignIO, SeqIO
import argparse
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv, sys, os, re


def hamming_distance(seq1, seq2):
	agreeBP = sum([base1 == base2 for base1, base2 in zip(seq1, seq2) if base1 != "-" and base2 != "-"])
	disagreeBP = sum([base1 != base2 for base1, base2 in zip(seq1, seq2) if base1 != "-" and base2 != "-"])

	return (agreeBP, disagreeBP)

def window(seq, width=500, step=100):
	seqlen = len(seq)
	for i in range(0, seqlen, step):
		if i + width > seqlen:
			j = seqlen
		else:
			j = i + width
		yield (i, j)
		if j == seqlen:
			break

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("input_alignment")
	parser.add_argument("--window", nargs="?", const=1, default=1000)
	parser.add_argument("--step", nargs="?", const=1, default=100)
	parser.add_argument("--minVarSite", nargs="?", const=1, default=10)
	parser.add_argument("--bedFile", nargs="?", const=1, default="")
	parser.add_argument("--outPrefix", nargs="?", const=1, default="")
	args = parser.parse_args()

	# Load the sequence alignment.

	if args.input_alignment in ["-","stdin"]:
		infile = sys.stdin
	else:
		infile = args.input_alignment

	with open(args.bedFile) as f_bed:
		for line in f_bed:
			chrom, orig_start, orig_end = line.strip().split()

	# Load the multiple sequence alignment and sort records by name in ascending order.
	alignment = AlignIO.read(infile, "fasta")

	length_aln = alignment.get_alignment_length()

	for win in window(range(length_aln), width=int(args.window), step=int(args.step)):
		start, end = win
		win_aln = alignment[:, start:end]
		num_varSite = 0
#		print (start, end, end-start)
		for idx in range(end-start):
			set_elements = set(win_aln[:,idx])
			if set_elements == {"N"} or set_elements == {"-"}:
				continue
			else:
				if "N" in set_elements:
					set_elements.remove("N")
				elif "-" in set_elements:
					set_elements.remove("-")

				if len(set_elements) > 1:
#					print(set_elements)
					num_varSite += 1

		if num_varSite >= int(args.minVarSite):
			outfname = args.outPrefix
#			print(num_varSite)
			with open(outfname + "%s-%s.fa" % (int(orig_start) + start, int(orig_start) + end) , "w") as fout:
				SeqIO.write(win_aln, fout, "fasta")
