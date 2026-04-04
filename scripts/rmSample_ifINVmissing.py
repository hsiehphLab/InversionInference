import sys, os, re

if __name__ == "__main__":

	if sys.argv[1] in ["-", "stdin"]:
		fin = sys.stdin
	else:
		exit("Error! The first argument is stdin. exit")

	keyWord = sys.argv[2]

	l_snvs = []
	indices_missing = []
	for line in fin:
		if line.startswith("#"):
			if line.startswith("#CHROM"):
				l_header = line.strip().split()
			else:
				print (line.strip())
		else:
			l_line = line.strip().split()
			if re.match(l_line[2], keyWord):
				indices_missing.extend([i for i, y in enumerate(l_line) if y in [".|.", "./."]])
#				print(l_line)
#				print (indices_missing)
			l_snvs.append(l_line)



	if indices_missing == []:
		print ("\t".join(l_header))
		for snv in l_snvs:
			print( "\t".join(snv))
				
	else:
		new_l_header = [y for i, y in enumerate(l_header) if i not in indices_missing]
		print ("\t".join(new_l_header))

		for snv in l_snvs:
			new_snv = [y for i, y in enumerate(snv) if i not in indices_missing]
#			if "INV" in snv[2]:
#				print (new_snv)
#				exit()
			print ("\t".join(new_snv))


