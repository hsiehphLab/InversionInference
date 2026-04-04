import sys, os, re, gzip, random
import pickle as cPickle

def extract_chimp_chr(chrom, list_chimp_fileObj):
	with open(list_chimp_fileObj[chrom]) as f:
		for line in f:
			if line.startswith('#'):
				continue
			else:
				return [chrom, line]


if __name__ == "__main__":

	newINVID = sys.argv[1]
	fin_INVID = sys.argv[2]

	with open(fin_INVID) as f_invID:
		for line in f_invID:
			line = line.strip().split()
			if line[0] == newINVID:
				origINVID, chrom, start, end = line[1:5]
				break

	fin_INV = sys.argv[3]
	
	try:
		rmLowConf = sys.argv[7]
		rmLowConf = True
	except IndexError:
		rmLowConf = False

	path_bytemap_dict_chrom = os.path.abspath(sys.argv[6])
	fnames = os.listdir(path_bytemap_dict_chrom)

	dict_fnames_bytemap = {}
	for name in fnames:
		if name.startswith('chr'):
			chrID = re.search('chr([a-zA-Z0-9]+).*', name).group(1)
			dict_fnames_bytemap[chrID] = os.path.join(path_bytemap_dict_chrom, name)

	path_humandir = os.path.abspath(sys.argv[4])
	list_humanchr = os.listdir(path_humandir)
	path_chimpdir = os.path.abspath(sys.argv[5])
	list_chimpchr = os.listdir(path_chimpdir)
	humanfiles = {}
	chimpfiles = {}

	for name in list_humanchr:
		if re.search('seq$', name):
			chr = re.search('chr([a-zA-Z0-9]+).*\.seq$', name).group(1)
			name = os.path.join(path_humandir, name)
			humanfiles[chr] = name
#			f = open(name, 'r')
#			humanfiles[chr] = f
	
	for name in list_chimpchr:
		if re.search('seq$', name):
			chr = re.search('chr([a-zA-Z0-9]+).*\.seq$', name).group(1)
			name = os.path.join(path_chimpdir, name)
			chimpfiles[chr] = name
#			f = open(name, 'r')
#			chimpfiles[chr] = f


	dict_sampleID_ListCN = {}
	list_usedPos = []

	with open(fin_INV) as fin, sys.stdout as fout:
		header = fin.readline().strip().split()
		new_header = ["##fileformat=VCFv4.2", "##filedate=20170918","##INFO=<ID=AA,Number=A,Type=String,Description=\"Ancestral state based on chimpanzee seq\">",
				"##INFO=<ID=VT,Number=A,Type=String,Description=\"Variant type\">", "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated ALT Allele Frequencies\">",
				"##INFO=<ID=AR2,Number=1,Type=Float,Description=\"Allelic R-Squared: estimated squared correlation between most probable REF dose and true REF dose\">",
				"##INFO=<ID=DR2,Number=1,Type=Float,Description=\"Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">",
				"##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">",
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
				"##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + P(AA)]\">",
				"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">"]
		fout.write("\n".join(new_header))

		new_header = ["\n#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
		new_header.extend(header[7:])
		fout.write("\t".join(new_header) + "\n")

		curr_chimp_chr = None
		curr_human_chr = None
		for wholeline in fin:
			flag_nextVariant = False
			line = wholeline.strip().split()
			if line[4] != origINVID:
				continue

			if line[0].startswith("chr"):
				chrID = line[0][3:]
			else:
				chrID = line[0]
			start = int(start)
			end = int(end)
			length = end - start + 1
			pos = int(end - 0.01 * length)
#			pos = int(line[3])

			while True:
				if pos <= start or pos > end:
					sys.stderr.write(wholeline)
					flag_nextVariant = True
					break

				if curr_chimp_chr is None:
					curr_chimp_chr = extract_chimp_chr(chrID, chimpfiles)
					curr_human_chr = extract_chimp_chr(chrID, humanfiles)
					curr_bytemap_chr = cPickle.load(gzip.open(dict_fnames_bytemap[chrID]))
					try:
						chimp_base = curr_chimp_chr[1][pos]
						human_base = curr_human_chr[1][pos]
					except TypeError:
						print ("WTF_1, chr%s, %s" % (chrID, pos))
						chimp_base = "N"
						human_base = "N"
				elif curr_chimp_chr[0] != chrID:
					curr_chimp_chr = extract_chimp_chr(chrID, chimpfiles)
					curr_chimp_chr = extract_chimp_chr(chrID, chimpfiles)
					curr_bytemap_chr = cPickle.load(gzip.open(dict_fnames_bytemap[chrID]))
					try:
						chimp_base = curr_chimp_chr[1][pos]
						human_base = curr_human_chr[1][pos]
					except TypeError:
						print ("WTF_2,%s" % pos)
						chimp_base = "N"
						human_base = "N"
				else:
					chimp_base = curr_chimp_chr[1][pos]
					human_base = curr_human_chr[1][pos]


				if curr_bytemap_chr[pos] != 49 :
#					pos = int(random.uniform(start, end))
					val = curr_bytemap_chr[pos-10:pos+1]
					try:
						idx_last_availbase = (len(val) -1) - list(reversed(val)).index(49)
						pos = pos - 10 + idx_last_availbase
#						print (pos, pos-100000, idx_last_availbase)
					except ValueError:
						pos -= 10
				else:
					if human_base in ["N","-"] or chimp_base in ["N","-"] or human_base != chimp_base or pos in list_usedPos:
						pos -= 1
					else:
						list_usedPos.append(pos)
						break

			if flag_nextVariant:
				continue

			if rmLowConf:
				genos = [ "./." if x not in ["0|0","0|1","1|0","1|1"] else x for x in line[7:]]
			
			else:
				genos = [ x if x in ["0|0","0|1","1|0","1|1"] else x.replace("_lowconf","") if x.find("_lowconf") != -1 else "./." for x in line[7:]]

			CNtype = "INV"
			list_bases = ["A","C","G","T"]
			list_bases.remove(human_base)
			new_human_base = list_bases[-1]
#			list_out = ["chr"+str(chrID), pos, newINVID, human_base, "<"+CNtype+">", ".", ".", "AA="+chimp_base+";VT="+CNtype, "GT"]
			list_out = ["chr"+str(chrID), pos, newINVID, chimp_base, new_human_base, ".", ".", "AA="+chimp_base+";VT="+CNtype, "GT"]

			list_out.extend(genos)
			fout.write("\t".join([str(x) for x in list_out]) + "\n")
			
