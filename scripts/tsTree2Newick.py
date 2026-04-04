import sys, os
import tskit
import pandas as pd
from io import StringIO
from ete3 import Tree
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer


if __name__ == "__main__":
	
	fin = sys.argv[1]
	hap = sys.argv[2]
	inv = sys.argv[3]
	out_homoplasy1 = sys.argv[4]
	out_homoplasy2 = sys.argv[5]

	ts = tskit.load(fin)

	dict_map_node2label = {}
	list_samples = []
	with open(hap) as f:
		count = 0
		for i,line in enumerate(f):
			if line.startswith("sample"):
				continue
			idx = 1
			while idx < 3:
				tmp = line.strip().split()
				dict_map_node2label[count] = tmp[0] + "_%s" % idx
				list_samples.append(tmp[0] + "_%s" % idx)
				idx += 1
				count += 1

	df_mapping_hap_INV = pd.read_csv(inv, delimiter="\t")
	df_mapping_hap_INV = df_mapping_hap_INV.set_index("orig_hapID")
	df_mapping_hap_INV_sorted = df_mapping_hap_INV.loc[list_samples]
	list_traits = [1 if x==True else 0 for x in df_mapping_hap_INV_sorted.loc[:,"INV"]]	
	list_haps = list(df_mapping_hap_INV_sorted.index)
	list_hapTrait = zip(list_haps, list_traits)

	list_trait_aln = []

	for pair in list_hapTrait:
		if pair[1] == 0:
			list_trait_aln.append(SeqRecord(Seq("A"), id=pair[0]))
		elif pair[1] == 1:
			list_trait_aln.append(SeqRecord(Seq("T"), id=pair[0]))

	trait_aln = MultipleSeqAlignment(list_trait_aln)

	trees = ts.trees()
	count = 0

	with open(out_homoplasy1, "w") as fout1, open(out_homoplasy2, "w") as fout2:
		for t in trees:
			# output newick trees
			print(t.newick(node_labels=dict_map_node2label))

			# infer the number of mutations using the Fitch's parsimony in the tskit's implementation
			ancestral_state, mutations = t.map_mutations(alleles=[0,1], genotypes=list_traits)
			fout1.write("%s\t%s\t%.5f\n" % (count, len(mutations), t.total_branch_length))

			# infer the number of mutations using the Fitch's parsimony in the Biopython's implementation
			biopython_tree = Phylo.read(StringIO(t.newick(node_labels=dict_map_node2label)), "newick")
			scorer = ParsimonyScorer()
			minHomoplasy = (scorer.get_score(biopython_tree, trait_aln))
			fout2.write("%s\t%s\t%.5f\n" % (count, minHomoplasy, biopython_tree.total_branch_length()))

			count += 1


