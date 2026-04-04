import sys, os

if __name__ == "__main__":

    dict_sampleID_INVgeno = {}

    l_rmSamples = []
    with sys.stdin as fin:
        for line in fin:
            l_line = line.strip().split()

            sampleID, hap = l_line[0].split("_")
            orientation = l_line[1]

            if sampleID not in dict_sampleID_INVgeno:
                dict_sampleID_INVgeno[sampleID] = {}

            if hap not in dict_sampleID_INVgeno[sampleID]:
                if orientation == "DIR":
                    inv = 0
                elif orientation == "INV":
                    inv = 1
                dict_sampleID_INVgeno[sampleID][hap] = inv
            else:
                l_rmSamples.append(sampleID)
                del dict_sampleID_INVgeno[sampleID]

    l_sample = ["seqnames", "start", "end", "POS", "orig_ID", "verdict", "categ"]
    seqnames =  sys.argv[1]
    start = sys.argv[2]
    end = sys.argv[3]
    pos = int( (int(start) + int(end))/2)
    orig_ID = "%s-%s-%s-%s" % (seqnames, start, "INV", int(end)-int(start))
    verdict = "pass"
    categ = "inv"

    l_geno = [ str(x) for x in [seqnames, start, end, pos, orig_ID, verdict, categ]]

    for s in dict_sampleID_INVgeno:
        l_sample.append(s)
        l_hap = []
        for h in sorted(dict_sampleID_INVgeno[s]):
            l_hap.append(dict_sampleID_INVgeno[s][h])
        l_geno.append("|".join([str(x) for x in l_hap]))

    print ("\t".join(l_sample))
    print ("\t".join(l_geno))

    with open("rmSamples.id", "w") as fout:
        for s in l_rmSamples:
            fout.write(s+"\n")
    
#    for i,v in enumerate(l_sample):
#        print("\t".join([v, l_geno[i]]))
