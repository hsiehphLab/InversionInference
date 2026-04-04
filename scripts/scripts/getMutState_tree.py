import sys, os
import pandas as pd
from ete3 import Tree
from Bio import Phylo
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer

#sys.argv = ["","chr8-7301025-INV-5297356.treefile", "mapping_hap_INV.txt", "temp"]
#sys.argv = ["","12001024-12051024.treefile", "mapping_hap_INV.txt", "temp"]
#sys.argv = ["","8201024-8301024.treefile", "mapping_hap_INV.txt", "temp"]

fname_tree = sys.argv[1]
tree = Phylo.read(fname_tree, "newick")
tree.collapse("CMP")

fname_mapping_hap_INV = sys.argv[2]
df_mapping_hap_INV = pd.read_csv(fname_mapping_hap_INV, delimiter="\t")

list_traits = [1 if x==True else 0 for x in df_mapping_hap_INV.loc[:,"INV"]]
list_haps = [x for x in df_mapping_hap_INV.loc[:,"orig_hapID"]]
list_trait_aln = []
list_hapTrait = zip(list_haps, list_traits)
for pair in list_hapTrait:
	if pair[1] == 0:
		list_trait_aln.append(SeqRecord(Seq("A"), id=pair[0]))
	elif pair[1] == 1:
		list_trait_aln.append(SeqRecord(Seq("T"), id=pair[0]))

trait_aln = MultipleSeqAlignment(list_trait_aln)
terms = tree.get_terminals()
terms.sort(key=lambda term: term.name)
trait_aln.sort()
column_i = trait_aln[:, 0]
column_i == len(column_i) * column_i[0]
clade_states = dict(zip(terms, [{c} for c in column_i]))

score_i = 0
count = 0

for clade in tree.get_nonterminals(order="postorder"):
	clade_childs = clade.clades
	left_state = clade_states[clade_childs[0]]
	right_state = clade_states[clade_childs[1]]
	state = left_state & right_state
	if not state:
		state = left_state | right_state
		score_i = score_i + 1
	clade_states[clade] = state
	clade.name =  "Inner%s_" % count + "".join(list(state)) + ":"
	count += 1

print ("NumberMutationEvents:%s" % score_i)

Phylo.write(tree,"%s.nwk" % sys.argv[3]  ,"newick")

#Phylo.write(tree,"temp.xml","phyloxml")
#Phylo.convert("temp.nwk", "newick", "temp.nex","nexus")
#Phylo.convert("temp.nwk", "newick", "temp.xml","phyloxml")




