import numpy as np
import pandas as pd
import sys, os
from sklearn.metrics.pairwise import pairwise_distances
import itertools



if __name__ == "__main__":

	keyword = sys.argv[3]

	tped = pd.read_csv(sys.argv[1], header=None, delimiter=" ", na_values=["."])
	hapID = pd.read_csv(sys.argv[2], header=None, delimiter="\t")

	indexINV = tped.index[tped.iloc[:,1].str.contains(keyword, na=False)].tolist()[0]

	inv = pd.DataFrame(tped.iloc[indexINV,4:]).reset_index(drop=True)
	inv = inv.rename(index=hapID.iloc[:,0], columns={indexINV:"INV"})
	
	pd_haplo_matrix = tped.iloc[:,4:].transpose().reset_index(drop=True)

	# compute element-wise difference for pairs of rows in pd_haplo_matrix
	l_comb = []
	for c in itertools.combinations(range(pd_haplo_matrix.shape[0]),2): 
		try:
			l_pair_Diff = [x for x in list(abs(pd_haplo_matrix.reindex(c).diff().reset_index(drop=True).loc[1,:]))]
			l_pair_hapID = hapID.reindex(c)
			l_pair_hapID = [y[0] for y in hapID.reindex(c).values.tolist()]
			if np.nan in l_pair_hapID:
				print (l_pair_Diff)
				print (l_pair_hapID)
				print (c)
				exit("Error!")
		except KeyError:
			print (c)
			exit("Error!")
		try:
			l_out = [int(x) for x in inv.loc[l_pair_hapID,"INV"].values.tolist()]
		except ValueError:
			print ("#####")
			print( inv.loc[l_pair_hapID,"INV"].values.tolist())
			print (l_pair_Diff)
			print (l_pair_hapID)
			print (c)
			print( inv.loc["EUR_CEU_NA12878_hap1","INV"])
			exit("Error!")
			
		l_out.extend(l_pair_hapID)
		l_out.extend(l_pair_Diff)
		sys.stdout.write("\t".join([str(x) for x in l_out]) + "\n")
	

	# compute pairwise IBS
	IBS = pd.DataFrame(1 - pairwise_distances(pd_haplo_matrix, metric = "hamming", force_all_finite=False))
	IBS = IBS.rename(columns = hapID.iloc[:,0], index = hapID.iloc[:,0])

	# matrix to pairwise
	tri_IBS = IBS.where(np.triu(np.ones(IBS.shape),k=1).astype(np.bool))
	df_pw_IBS = tri_IBS.stack().reset_index()
	df_pw_IBS.columns = ["h1","h2","IBS"]

	df_pw_IBS = pd.merge(df_pw_IBS, inv, left_on="h1", right_index=True)
	df_pw_IBS = pd.merge(df_pw_IBS, inv, left_on="h2", right_index=True)
	

	df_pw_IBS.columns = ["h1","h2","IBS","h1_INVgeno","h2_INVgeno"]

	df_pw_IBS.to_csv(sys.argv[4], sep="\t", index=False)



