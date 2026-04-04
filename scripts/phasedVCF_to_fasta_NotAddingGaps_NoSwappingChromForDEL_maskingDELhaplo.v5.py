
## This script is used to calculate the observed divergence between human and chimp using a sample of human genomes. The divergence is taken by averaging the divergence estimates between the chimp ref and each of the samples (haploid sequences) according to the cPickle and TPED file.
##
## Usage:
##		python new_calculate_divergence_humanchimp_seq_v3_average_div_between_chimp_and_TPED.py  _path_chimpRef  _path_humanRef   _filename_cPickle_biteMap  _option_calculation  __fname_TPEDfile  >  chimp_human_divergence.dat
##
## The input file '_filename_cPickle_biteMap' is typically based on called positions across all samples.
## This is a correct calculation for sequence divergence based on Nei 1987 (equation 5.3 and 10.20).

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import random, os, re, sys, time, gzip, vcf, pysam, copy
import numpy as np
import argparse
import pickle as cPickle

import sys, os, argparse
from pysam import VariantFile
from pysam import tabix_index


parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("--bedFile", nargs="?", const=1, default="")
parser.add_argument("--outgrp", nargs="?", const=1, default="")
parser.add_argument("--humanRef", nargs="?", const=1, default="")
parser.add_argument("--maskArray", nargs="?", const=1, default="")
parser.add_argument("--maskDELbed", nargs="?", const=1, default="")
parser.add_argument("--maskPosFile", nargs="?", const=1, default="")
parser.add_argument("--window", nargs="?", const=1, default=1000)
parser.add_argument("--step", nargs="?", const=1, default=100)
parser.add_argument("--minVarSite", nargs="?", const=1, default=10)
parser.add_argument("--outPrefix", nargs="?", const=1, default="")
args = parser.parse_args()

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


def opener(filename):
	f = open(filename,'r')
	if (f.read(2) == '\x1f\x8b'):
		f.seek(0)
		return gzip.GzipFile(fileobj=f)
	else:
		f.seek(0)
		return f


path_chimpdir = os.path.abspath(args.outgrp)
list_chimpchr = os.listdir(path_chimpdir)

path_humandir = os.path.abspath(args.humanRef)
list_humanchr = os.listdir(path_humandir)

path_bytemap = os.path.abspath(args.maskArray)

chimpfiles = {}
## get the file objs for all the chimp sequences as a dictionary. each key is a numeric chromosome id. 
for name in list_chimpchr:
	if re.search('seq(\.gz)?$', name):
		tmp = os.path.basename(name).split(".")[0].split("_")[0]
		name = os.path.join(path_chimpdir, name)
		if name.endswith('.gz'):
			f = gzip.open(name, 'rt')
		else:
			f = open(name, 'r')
		chimpfiles[tmp] = f


humanfiles = {}
## get the file objs for all the human sequences as a dictionary. each key is a numeric chromosome id. 
for name in list_humanchr:
	if re.search('seq(\.gz)?$', name):
		tmp = os.path.basename(name).split(".")[0].split("_")[0]
		name = os.path.join(path_humandir, name)
		if name.endswith('.gz'):
			f = gzip.open(name, 'rt')
		else:
			f = open(name, 'r')
		humanfiles[tmp] = f

# load bytearray by chromosome
dict_chromID_bytemap = {}
for f in os.listdir(path_bytemap):
	if re.search('cpickle.gz$', f):
		chromID = os.path.basename(f).split(".")[0].split("_")[0]
#		chromID = os.path.basename(f).split(".")[0]
		dict_chromID_bytemap[chromID] = os.path.join(path_bytemap, f)


# construct a list of tuples, in which each tuple has the structure (chromID, PhyPos)
# This list is used to check if the variable sites between human and chimp ref genomes are also variable in our samples.
# Excluding sites that are not fixed between human and chimp provide a conservative estimate for the divergence of human and chimp. 

bases = ['A','T','C','G']
pattern = re.compile(b'0+')

curr_chrom = None

if args.maskPosFile != "":
	f_maskPosFile = open(args.maskPosFile)
else:
	f_maskPosFile = None


if args.maskDELbed != "":
	f_maskDELbed = open(args.maskDELbed)
else:
	f_maskDELbed = None


dict_pos_randomAlleles = {}
dict_IDSample_DELchrom = {}
dict_idxSample_DELchrom = {}
dict_IDSample_DUPchrom = {}
dict_idxSample_DUPchrom = {}

################33

with open(args.bedFile) as f_bed:
	for line in f_bed:
		chrom, orig_start, orig_end = line.strip().split()

orig_start = int(orig_start)
orig_end = int(orig_end)

for line in humanfiles[chrom]:
	if line.startswith('#'):
		continue
	else:
		human_seq_chrom = line

for line in chimpfiles[chrom]:
	if line.startswith('#'):
		continue
	else:
		chimp_seq_chrom = line

if chrom != curr_chrom:
	curr_chrom = chrom
	try:
		curr_bytearray = cPickle.load(gzip.open(dict_chromID_bytemap[curr_chrom]))
	except KeyError:
		curr_bytearray = []

outfname = args.outPrefix

vcf_in = VariantFile(args.vcf)

length_aln = orig_end - orig_start + 1

for win in window(range(length_aln), width=int(args.window), step=int(args.step)):
	start, end = win
	adj_start = start + orig_start
	adj_end = end + orig_start

	tmp_human_seq = human_seq_chrom[adj_start + 1: adj_end + 1]
	tmp_chimp_seq = chimp_seq_chrom[adj_start + 1: adj_end + 1]

	set_tmp_human_seq = set(tmp_human_seq)
	set_tmp_chimp_seq = set(tmp_chimp_seq)


	try:
		set_tmp_human_seq.remove("-")
	except KeyError:
		pass
	try:
		set_tmp_human_seq.remove("N")
	except KeyError:
		pass
	try:
		set_tmp_chimp_seq.remove("-")
	except KeyError:
		pass
	try:
		set_tmp_chimp_seq.remove("N")
	except KeyError:
		pass

	if set_tmp_human_seq == set() or set_tmp_chimp_seq == set():
		continue

	chimp_seq = list(tmp_chimp_seq)
	fname_VCFout = outfname + "%s-%s.vcf" % (adj_start, adj_end) 
	fname_FASTAout = outfname + "%s-%s.fa" % (adj_start, adj_end) 

	l_snvs = []
	total_number_haploidSamples = 2 * len(vcf_in.header.samples)
	list_sampleID_seqs = []
	for sample in vcf_in.header.samples:
		human_seq1 = copy.deepcopy(list(tmp_human_seq))
		human_seq2 = copy.deepcopy(list(tmp_human_seq))
		list_sampleID_seqs.append([sample, human_seq1, human_seq2])

	vcf_records = vcf_in.fetch(chrom, adj_start, adj_end)
	list_outVCF = []
	countSNVs = 0 
	for record in vcf_records:
		if curr_bytearray == [] or curr_bytearray[record.pos] == 49:
			list_outVCF.append((record.pos, record))
			l_snvs.append(record)
			countSNVs += 1

	if countSNVs >= int(args.minVarSite):
		for record in l_snvs:
			idx_dash = None

			lookup_REF_ALT = [str(record.ref[0])]
			lookup_REF_ALT.extend(["-" if str(x) == "<DEL>" else "+" if str(x) == "<DUP>" else str(x) for x in record.alts])

			if '-' in lookup_REF_ALT or '+' in lookup_REF_ALT:
				avail_nucleotides = [ x for x in ["A","G","C","T"] if x not in lookup_REF_ALT and x != record.info['AA'][0]]
				# simplify the choice of nucleotide for DEL/DUP by picking the very first one in the list avail_nucleotides
				nucleotide_forDELsample = avail_nucleotides[0]
				if '-' in lookup_REF_ALT:
					idx_dash = lookup_REF_ALT.index('-')
				else:
					idx_dash = lookup_REF_ALT.index('+')
				lookup_REF_ALT[idx_dash] = nucleotide_forDELsample

			for idx, sample in enumerate(record.samples):

				alleles = list(record.samples[sample].alleles)
#				if None in alleles:
#					print (record.pos, sample, record.samples[sample].alleles)
#					exit()
				if None in alleles:
					alleles = ["N", "N"]

#				print (alleles[0], record.pos - (adj_start + 1))
#				print (alleles[1], record.pos - (adj_start + 1))
				list_sampleID_seqs[idx][1][record.pos - (adj_start + 1)] = alleles[0]
				list_sampleID_seqs[idx][2][record.pos - (adj_start + 1)] = alleles[1]
#			if record.info['AA'][0] == "<DUP>" or  record.info['AA'][0] == "<DEL>":
#				print (del_pos, record.info['AA'][0], nucleotide_forDELsample)
#				chimp_seq[record.pos - (pos_start + 1)] = nucleotide_forDELsample

		if curr_bytearray != []:
			for segment in re.finditer(pattern, curr_bytearray[adj_start + 1 : adj_end + 1]):
				seg_start = segment.start()
				seg_end = segment.end()
				chimp_seq[seg_start : seg_end] = "N" * (seg_end - seg_start)
				for idx in range(len(list_sampleID_seqs)):
					list_sampleID_seqs[idx][1][seg_start : seg_end] = "N" * (seg_end - seg_start)
					list_sampleID_seqs[idx][2][seg_start : seg_end] = "N" * (seg_end - seg_start)

		l_pos_exclude = []

		if f_maskDELbed is not None:
			# keep tracking if a position is within the BED for exclusion for later use
			for line in f_maskDELbed:
				line = line.strip().split()
				# check if the file f_maskDELbed contains the right region
				if line[0] != chromID:
					continue

				l_pos, l_record = zip(*list_outVCF)
				for pos in l_pos:
					if pos >= int(line[1]) + 1 and pos <= int(line[2]):
						l_pos_exclude.append(pos)

				try:
					maskStartPos_del = (int(line[1]) + 1) - (adj_start + 1)
					maskEndPos_del = int(line[2]) - (adj_start + 1)
					for idxSample in sorted(dict_idxSample_DELchrom):
						l_chrom_idx = dict_idxSample_DELchrom[idxSample]
						for idxChrom in l_chrom_idx:
							list_sampleID_seqs[idxSample][idxChrom + 1][maskStartPos_del: (maskEndPos_del + 1)] = "-" * ((maskEndPos_del + 1) - maskStartPos_del)
		#						print ("sample_idx%s\tchrom_idx%s" % (idxSample, idxChrom))
				except ValueError:
					continue


		if f_maskPosFile is not None:
			count_exMissingSites = 0
			for line in f_maskPosFile:
				line = line.strip().split()
				try:
					if int(line[1]) >= adj_start + 1  and int(line[1]) <= adj_end:
						maskPos = int(line[1]) - (adj_start + 1)
						chimp_seq[maskPos] = "N"
						for idx in range(len(list_sampleID_seqs)):
							list_sampleID_seqs[idx][1][maskPos] = "N" 
							list_sampleID_seqs[idx][2][maskPos] = "N"
				except ValueError:
					continue
			print ("%s missing sites were masked in the output FASTA files" % count_exMissingSites)


		out_list = []
		out_list.append(">CMP" + "\n")
		seq = ''.join(chimp_seq) + "\n"
		out_list.append(seq)
		len_seq = len(chimp_seq)

		for sample in list_sampleID_seqs:
			sampleID = sample[0]
			if sampleID in dict_IDSample_DELchrom:
				l_idx_DELchrom = dict_IDSample_DELchrom[sampleID]
				print ("l_idx_DELchrom", sampleID, l_idx_DELchrom)
			else:
				l_idx_DELchrom = []

			if sampleID in dict_IDSample_DUPchrom:
				l_idx_DUPchrom = dict_IDSample_DUPchrom[sampleID]
				print ("l_idx_DUPchrom", sampleID, l_idx_DUPchrom)
			else:
				l_idx_DUPchrom = []

			for seq in range(1,3):
				svComment = ''
				if (seq - 1) in l_idx_DELchrom:
					svComment += '-'
				if (seq - 1) in l_idx_DUPchrom:
					svComment += '+'
				out_list.append(">" + sampleID + "_" + str(seq) + svComment + "\n")
#				print (">" + sampleID + "_" + str(seq) + svComment + "\n")

				if len(sample[seq]) != len_seq:
					print(sampleID, len_seq, len(sample[seq]))
					exit("WTF!")

				seq1 = ''.join(sample[seq]) + "\n"
				out_list.append(seq1)

		with open(fname_FASTAout, "w") as fout:
			for entry in out_list:
				fout.write(entry)


		with VariantFile(fname_VCFout, "w", header=vcf_in.header) as vcf_out:
			for rec in l_snvs:
				vcf_out.write(rec)
		tabix_index(fname_VCFout, preset="vcf", force=True)


