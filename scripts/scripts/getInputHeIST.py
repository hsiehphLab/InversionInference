import sys, os

from ete3 import Tree

if __name__ == "__main__":

	fname_contree = sys.argv[1]
	tree1 = Tree(fname_contree, format=1)
	tree2 = Tree(fname_contree, format=1)

	fname_mapping_hap_INV = sys.argv[2]
	list_hapINV = []
	with open(fname_mapping_hap_INV) as fin:
		for line in fin:
			if line.startswith("hapID"):
				continue
			tmp_line = line.strip().split()
			if tmp_line[1] == "TRUE":
				list_hapINV.append(tmp_line[2])

	count = 0
	for node in tree2.traverse():
		if not node.is_leaf():
			if node.name == "":
				continue
			else:
				count += 1
				node.name = "I%s" % count

	out_tree1 = tree1.write(format=1)
	out_tree2 = tree2.write(format=1)


	print("#NEXUS")
	print("begin trees;")
	print("tree tree_1 = " + out_tree1)
	print("tree tree_2 = " + out_tree2)
	print("end;")
	print("begin hemiplasytool;")
	for inv in list_hapINV:
		print("set derived taxon=%s" % inv)
	print("set conversion type=extend")
	print("end;")




