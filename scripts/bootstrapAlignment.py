from Bio import AlignIO, SeqIO
import argparse
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv, sys, os, re, glob, random



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("input_Path2Aln")
	parser.add_argument("--numBT", nargs="?", const=1, default=10)
	parser.add_argument("--outPrefix", nargs="?", const=1, default="")
	parser.add_argument("--seed", nargs="?", const=1, default=10000)
	args = parser.parse_args()

	if int(args.seed) == -999:
		seed = random.choice(range(1000000000)) 
	else:
		seed = int(args.seed) * 3
	
	random.seed(seed)
	sys.stderr.write(str(seed) + "\n")

	# Load the sequence alignment.

	list_sliceAlns = glob.glob(args.input_Path2Aln + "/*.fa")
	num_slices = len(list_sliceAlns)

	count = 0
	for i in range(int(args.numBT)):
		idx_bt_samples = [random.choice(range(num_slices)) for x in range(num_slices)]
		bt_sliceAlns = [list_sliceAlns[i] for i in idx_bt_samples]
		for f in bt_sliceAlns:
			alignment = AlignIO.read(f, "fasta")
			if count == 0:
				cat_align = alignment
			else:
				cat_align += alignment	
			count += 1

		fout = sys.stdout
		
		AlignIO.write(cat_align, fout, "fasta")
		fout.close()


