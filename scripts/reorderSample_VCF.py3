import gzip, sys, os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", nargs="?", const=1, default="stdin")
    parser.add_argument("--IDfile")
    args = parser.parse_args()


    list_IDs = []
	# only samples in list_IDs will be written out
    with open(args.IDfile) as f_newOrder:
        for line in f_newOrder:
            line = line.strip()
            list_IDs.append(line)

    dict_idx_ID = {}
    if args.vcf == "stdin":
        f_vcf = sys.stdin
    elif args.vcf.endswith("gz"):
        f_vcf = gzip.open(args.vcf)
    else:
        f_vcf = open(args.vcf)

    list_missingSample = []

    for line in f_vcf:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                out_line = []
                line = line.strip().split()
                for idx, entry in enumerate(line):
                    if idx < 9:
                        out_line.append(entry)
                    if idx >= 9:
                        dict_idx_ID[idx] = entry

                list_valid_samples = list(dict_idx_ID.values())
                for sample in list_IDs:
                    if sample in list_valid_samples:
                        out_line.append(sample)

                print(('\t'.join(out_line)))

            else:
                print((line.strip()))
        else:
            dict_ID_vals = {}
            line = line.strip().split()
            out_line = []
            for idx, val in enumerate(line):
                if idx < 9:
                    out_line.append(line[idx])
                else:
                    try:
                        sampleID = dict_idx_ID[idx]
                    except KeyError:
                        print((list(dict_idx_ID.keys())))
                    dict_ID_vals[sampleID] = val

            for sample in list_IDs:
                try:
                    out_line.append(dict_ID_vals[sample])
                except KeyError:
                    if sample not in list_missingSample:
                        list_missingSample.append(sample)
                        sys.stderr.write(sample + "\n")
            print(("\t".join(out_line)))



