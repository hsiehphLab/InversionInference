import sys, os
import pandas as pd
from ete3 import Tree
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer


if __name__ == "__main__":

	fname_tree = sys.argv[1]
	tree = Phylo.read(fname_tree, "newick")
	tree.collapse("CMP")

	fname_mapping_hap_INV = sys.argv[2]

	try:
		o_prefix = sys.argv[3]
	except IndexError:
		o_prefix = os.path.basename(sys.argv[1])


	df_mapping_hap_INV = pd.read_csv(fname_mapping_hap_INV, delimiter="\t")

	list_traits = [1 if x==True else 0 for x in df_mapping_hap_INV.loc[:,"INV"]]
	list_haps = [x for x in df_mapping_hap_INV.loc[:,"orig_hapID"]]

	list_hapTrait = zip(list_haps, list_traits)

	list_trait_aln = []

	for pair in list_hapTrait:
		if pair[1] == 0:
			list_trait_aln.append(SeqRecord(Seq("A"), id=pair[0]))
		elif pair[1] == 1:
			list_trait_aln.append(SeqRecord(Seq("T"), id=pair[0]))

	trait_aln = MultipleSeqAlignment(list_trait_aln)

	scorer = ParsimonyScorer()
	minHomoplasy = (scorer.get_score(tree, trait_aln))

	print ("%s\t%s" % (o_prefix, minHomoplasy))


