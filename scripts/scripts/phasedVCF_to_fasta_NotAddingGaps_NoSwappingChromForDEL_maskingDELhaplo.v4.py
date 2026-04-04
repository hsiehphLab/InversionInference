
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

parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("bedfile")
parser.add_argument("fout")
parser.add_argument("--outgrp", nargs="?", const=1, default="/net/eichler/vol26/home/hsiehph/genomicData/humanHg19_chimpClint/dir_clint/")
parser.add_argument("--humanRef", nargs="?", const=1, default="/net/eichler/vol26/home/hsiehph/genomicData/humanHg19_chimpClint/dir_human/")
parser.add_argument("--maskArray", nargs="?", const=1, default="/net/eichler/vol26/home/hsiehph/genomicData/bytearray_HapMapGenMapGRC37_exHGDPmyMask_maskTNFRSF10CDdup/")
parser.add_argument("--maskDELbed", nargs="?", const=1, default="")
parser.add_argument("--maskPosFile", nargs="?", const=1, default="")
args = parser.parse_args()


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
	if re.search('seq$', name):
		tmp = os.path.basename(name).split(".")[0].split("_")[0]
#		tmp = os.path.basename(name).split(".")[0]
		name = os.path.join(path_chimpdir, name)
		f = open(name, 'r')
		chimpfiles[tmp] = f


humanfiles = {}
## get the file objs for all the human sequences as a dictionary. each key is a numeric chromosome id. 
for name in list_humanchr:
	if re.search('seq$', name):
		tmp = os.path.basename(name).split(".")[0].split("_")[0]
#		tmp = os.path.basename(name).split(".")[0]
		name = os.path.join(path_humandir, name)
		f = open(name, 'r')
		humanfiles[tmp] = f

# load bytearray by chromosome
dict_chromID_bytemap = {}
for f in os.listdir(path_bytemap):
	if re.search('cpickle.gz$', f):
		chromID = os.path.basename(f).split(".")[0].split("_")[0]
#		chromID = os.path.basename(f).split(".")[0]
		dict_chromID_bytemap[chromID] = os.path.join(path_bytemap, f)



# dict_chroms = cPickle.load(open(sys.argv[3],'rb'))

# construct a list of tuples, in which each tuple has the structure (chromID, PhyPos)
# This list is used to check if the variable sites between human and chimp ref genomes are also variable in our samples.
# Excluding sites that are not fixed between human and chimp provide a conservative estimate for the divergence of human and chimp. 

bases = ['A','T','C','G']
f_loci = open(args.bedfile)
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



for locus in f_loci:
	vcf_reader = vcf.Reader(filename=args.vcf)
#	print(locus.strip())
	locus = locus.strip().split()
	chromID, start, end = locus[:3]

	chrom = chromID

	pos_start = int(start)
	pos_end = int(end)
#	pos_start, pos_end = [int(x) for x in pos.split('-')]

	for line in humanfiles[chrom]:
		if line.startswith('#'):
			continue
		else:
			tmp_human_seq = line[pos_start + 1: pos_end + 1]
	humanfiles[chrom].seek(0)

	for line in chimpfiles[chrom]:
		if line.startswith('#'):
			continue
		else:
			tmp_chimp_seq = line[pos_start + 1 : pos_end + 1]
	chimpfiles[chrom].seek(0)

	chimp_seq = list(tmp_chimp_seq)

	if chrom != curr_chrom:
		curr_chrom = chrom
		try:
			curr_bytearray = cPickle.load(gzip.open(dict_chromID_bytemap[curr_chrom]))
		except KeyError:
			curr_bytearray = []

	vcf_records = vcf_reader.fetch(chrom, pos_start , pos_end)
	total_number_haploidSamples = 2 * len(vcf_records.samples)
	list_sampleID_seqs = []
	for sample in vcf_records.samples:
		human_seq1 = copy.deepcopy(list(tmp_human_seq))
		human_seq2 = copy.deepcopy(list(tmp_human_seq))
		list_sampleID_seqs.append([sample, human_seq1, human_seq2])

	dict_pos_randomAlleles = {}
	dict_IDSample_DELchrom = {}
	dict_idxSample_DELchrom = {}
	dict_IDSample_DUPchrom = {}
	dict_idxSample_DUPchrom = {}
#	list_idx_SwappingChrom = []

	list_outVCF = []
	for record in vcf_records:
		idx_dash = None
		if curr_bytearray == [] or curr_bytearray[record.POS] == 49:
			list_outVCF.append((record.POS, record))

		lookup_REF_ALT = [str(record.REF[0])]
		lookup_REF_ALT.extend(["-" if str(x) == "<DEL>" else "+" if str(x) == "<DUP>" else str(x) for x in record.ALT])

		if '-' in lookup_REF_ALT or '+' in lookup_REF_ALT:
			avail_nucleotides = [ x for x in ["A","G","C","T"] if x not in lookup_REF_ALT and x != record.INFO['AA'][0]]
#			nucleotide_forDELsample = random.choice(avail_nucleotides)
			# simplify the choice of nucleotide for DEL/DUP by picking the very first one in the list avail_nucleotides
			nucleotide_forDELsample = avail_nucleotides[0]
			if '-' in lookup_REF_ALT:
				idx_dash = lookup_REF_ALT.index('-')
			else:
				idx_dash = lookup_REF_ALT.index('+')
			lookup_REF_ALT[idx_dash] = nucleotide_forDELsample

		for idx, sample in enumerate(record.samples):
			if sample['GT'] != '.':
				alleles = [lookup_REF_ALT[int(x)] if x != "." else "N" for x in re.split("\/|\|", sample['GT'])]
			else:
				alleles = ["N", "N"]
			try:
				if record.POS == 22981867 or record.POS == 22991200 or record.POS == 24521264 or record.POS == 866265 :
					del_pos = record.POS - (pos_start + 1)
					del_idx = [i for i,base in enumerate(alleles) if base == nucleotide_forDELsample]
					if del_idx != []:
						if "<DEL>" in [str(x) for x in record.ALT]:
							dict_IDSample_DELchrom[sample.sample] = del_idx
							dict_idxSample_DELchrom[idx] = del_idx
							print(sample.sample, idx, del_idx)
						elif "<DUP>" in [str(x) for x in record.ALT]:
							dict_IDSample_DUPchrom[sample.sample] = del_idx
							dict_idxSample_DUPchrom[idx] = del_idx
							print(sample.sample, idx, del_idx)
#						if len(del_idx) == 1 and del_idx[0] == 1:
#							list_idx_SwappingChrom.append(idx)

				list_sampleID_seqs[idx][1][record.POS - (pos_start + 1)] = alleles[0]
				list_sampleID_seqs[idx][2][record.POS - (pos_start + 1)] = alleles[1]

			except IndexError:
				print(sample)
				print(lookup_REF_ALT)
				print(list_sampleID_seqs[idx][1][record.POS - (pos_start + 1)])
				print(len(list_sampleID_seqs[idx][1]))
				print(record.POS - pos_start)
				print(alleles, record.POS, pos_start)
				exit()

		if record.INFO['AA'][0] == "<DUP>" or  record.INFO['AA'][0] == "<DEL>":
			print (del_pos, record.INFO['AA'][0], nucleotide_forDELsample)
			chimp_seq[record.POS - (pos_start + 1)] = nucleotide_forDELsample

	if curr_bytearray != []:
		for segment in re.finditer(pattern, curr_bytearray[pos_start + 1 : pos_end + 1]):
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
				maskStartPos_del = (int(line[1]) + 1) - (pos_start + 1)
				maskEndPos_del = int(line[2]) - (pos_start + 1)
				for idxSample in sorted(dict_idxSample_DELchrom):
					l_chrom_idx = dict_idxSample_DELchrom[idxSample]
					for idxChrom in l_chrom_idx:
						list_sampleID_seqs[idxSample][idxChrom + 1][maskStartPos_del: (maskEndPos_del + 1)] = "-" * ((maskEndPos_del + 1) - maskStartPos_del)
#						print ("sample_idx%s\tchrom_idx%s" % (idxSample, idxChrom))
			except ValueError:
				continue

#	print (len(chimp_seq))
#	print (len(list_sampleID_seqs[1][1]))
#	print (len(list_sampleID_seqs[1][2]))
#	exit()

	# output SNVs if not within the BED records to a VCF file.
#	vcf_writer = vcf.Writer(open('chr%s_%s_%s.exSNVs_in_BEDandMasking.vcf' % (chrom, start, end), 'w'), vcf_records)

#	for entry in list_outVCF:
#		currPos = entry[0]
#		currRecord = entry[1]
#		if currPos not in l_pos_exclude:
#			vcf_writer.write_record(currRecord)


	if f_maskPosFile is not None:
		count_exMissingSites = 0
		for line in f_maskPosFile:
			line = line.strip().split()
			try:
				if int(line[1]) >= pos_start + 1  and int(line[1]) <= pos_end:
					maskPos = int(line[1]) - (pos_start + 1)
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
#	out_list_nexus = []
#	out_list_nexus.append("Chimpanzee" + "\t" + seq)
	len_seq = len(chimp_seq)

	# Note that here I force the first chromosome (aka with idx "_1") to be the one carries OCN DEL allele
	# To do so, I switch chromosomes' indices if necessary (according to *list_idx_SwappingChrom*)
#	for idx in list_idx_SwappingChrom:
#		print(list_sampleID_seqs[idx][0], idx)
#		print(list_sampleID_seqs[idx][1][del_pos], list_sampleID_seqs[idx][2][del_pos])
#		list_sampleID_seqs[idx][1], list_sampleID_seqs[idx][2] = list_sampleID_seqs[idx][2], list_sampleID_seqs[idx][1]
#		print(list_sampleID_seqs[idx][1][del_pos], list_sampleID_seqs[idx][2][del_pos])

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
			print (">" + sampleID + "_" + str(seq) + svComment + "\n")

			if len(sample[seq]) != len_seq:
				print(sampleID, len_seq, len(sample[seq]))
				exit("WTF!")

			seq1 = ''.join(sample[seq]) + "\n"
			out_list.append(seq1)

	with open(args.fout, "w") as fout:
		for entry in out_list:
			fout.write(entry)
 
#	with open(fprefix_out + "_" + chromID + "_" + start + "_" + end + ".nexus", "wb") as fout:
#	numSeqs = len(out_list_nexus)
#	fout.write("#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=%s NCHAR=%s;\n" % (numSeqs, len_seq))
#	fout.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\nMATRIX\n")
#	for entry in out_list_nexus:
#		fout.write(entry)
#	fout.write(";\n\nEND;\n")



