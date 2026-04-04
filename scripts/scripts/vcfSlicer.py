import sys, os, argparse
from pysam import VariantFile
from pysam import tabix_index

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
	parser.add_argument("input_vcf")
	parser.add_argument("--window", nargs="?", const=1, default=1000)
	parser.add_argument("--step", nargs="?", const=1, default=100)
	parser.add_argument("--minVarSite", nargs="?", const=1, default=10)
	parser.add_argument("--bedFile", nargs="?", const=1, default="")
	parser.add_argument("--outPrefix", nargs="?", const=1, default="")
	args = parser.parse_args()

	# Load the sequence vcf.

	with open(args.bedFile) as f_bed:
		for line in f_bed:
			chrom, orig_start, orig_end = line.strip().split()

	orig_start = int(orig_start)
	orig_end = int(orig_end)

	outfname = args.outPrefix
	
	# Load the multiple sequence alignment and sort records by name in ascending order.
	vcf_in = VariantFile(args.input_vcf)

	length_aln = orig_end - orig_start + 1

	for win in window(range(length_aln), width=int(args.window), step=int(args.step)):
		start, end = win
		adj_start = start + orig_start
		adj_end = end + orig_start

		fname_out = outfname + "%s-%s.vcf" % (adj_start, adj_end) 

		l_snvs = []
		for record in vcf_in.fetch(chrom, adj_start, adj_end):
			l_snvs.append(record)

		if len(l_snvs) > int(args.minVarSite):
			with VariantFile(fname_out, "w", header=vcf_in.header) as vcf_out:
				for rec in l_snvs:
					vcf_out.write(rec)
			tabix_index(fname_out, preset="vcf")


