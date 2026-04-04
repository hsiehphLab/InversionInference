import subprocess
import pandas
# import matplotlib.pyplot as plt
import numpy
import os
import sys


# mostly based on information from
# https://myersgroup.github.io/relate/


# EXTERNAL TOOLS
# relate uses 'Rscript' for plotting
# needs libraries ggplot2, gridextra
# can be obtained with conda using
# conda install -c r r
# conda install -c r r-ggplot2
# conda install -c r r-gridextra


# BINARIES
# can be downloaded at https://myersgroup.github.io/relate/
# linux binaries (x86) work on the CRI cluster
# relateBinary = "./binaries/relate_v1.1.2_x86_64_static/bin/Relate"
# relateFileFormatBinary = "./binaries/relate_v1.1.2_x86_64_static/bin/RelateFileFormats"
# relatePopSizeScript = "./binaries/relate_v1.1.2_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh"
# mac binaries for office
# relateBinary = "../../tools/officeMac/relate_v1.1.6_MacOSX/bin/Relate"
# relateFileFormatBinary = "../../tools/officeMac/relate_v1.1.6_MacOSX/bin/RelateFileFormats"
# relatePopSizeScript = "../../tools/officeMac/relate_v1.1.6_MacOSX/scripts/EstimatePopulationSize/EstimatePopulationSize.sh"
# relateCoalescenceRateBinary = "../../tools/officeMac/relate_v1.1.6_MacOSX/bin/RelateFileFormats"
#relateBinary = "/home/steinrue/labshare/projects/grantStuff/populationSplit/tools/cluster/relate_v1.1.6_x86_64_static/bin/Relate"
relateBinary = "~/bin/relate_v1.1.7_x86_64_static/bin/Relate"
relateFileFormatBinary = "~/bin/relate_v1.1.7_x86_64_static/bin/RelateFileFormats"
relatePopSizeScript = "~/bin/relate_v1.1.7_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh"
relateCoalescenceRateBinary = "~/bin/relate_v1.1.7_x86_64_static/bin/RelateFileFormats"


def prepareInput (input_vcf_name, input_haps_file, input_sample_file):

	print ("[CONVERT_VCF_TO_RELATE]")

	# convert vcf to relate input

	# get basename from vcf, because relate doesn't want extensions
	base_vcf_name = os.path.splitext(input_vcf_name)[0]

	# binaries/relate_v1.1.2_x86_64_static/bin/RelateFileFormats --mode ConvertFromVcf --haps exp_dataset1.haps --sample exp_dataset1.sample --input ../exp_dataset1
	vcfToRelateCmd = relateFileFormatBinary + \
		" --mode ConvertFromVcf " + \
		" --haps " + input_haps_file + \
		" --sample " + input_sample_file + \
		" --input " + base_vcf_name

	# show it
	print (vcfToRelateCmd)

	# run the command
	# print (os.getcwd())
	process = subprocess.Popen (vcfToRelateCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, close_fds=True)
	print (process.communicate()) 

	process.terminate()

	# relate is algergic to SNPs being listed more than once, so we have to remove duplicates
	# reading everything at once is maybe not the most efficient, but ok for now
	relateInput = pandas.read_csv (input_haps_file, header=None, sep=" ")
	# the second row has the SNP numbers
	# we drop all but the first occurence
	relateInput.drop_duplicates (subset=[2], keep='first', inplace=True)
	# and write with duplicates removed
	ofs = open (input_haps_file, "w")
	relateInput.to_csv (ofs, header=None, sep=" ", index=False)
	ofs.close ()

	print ("[DONE]")


def prepareAuxiliaryFiles (recoProb, numLoci, input_genetic_map_file, input_sample_file, input_poplabels_file, samplesInDemes):
	
	# we need a genetic map
	# so make a file for a unifom genetic map

	print ("[FAKE_GENETIC_MAP]")

	# (pos, cM/MB, cM)
	# (p, r, rdist)
	# r[i] = (rdist[i+1] - rdist[i])/(p[i+1] - p[i]) * 1e6
	ofs = open (input_genetic_map_file, "w")
	# recoProb * 100 because it is CENTI morgan
	# not sure whether the formating might cause problems for certain parameters
	ofs.write (f"0\t{recoProb*100*1e6}\t0\n")
	ofs.write (f"{int(numLoci):d}\t{recoProb*100*1e6}\t{numLoci*recoProb*100}\n")
	ofs.close ()

	print ("[DONE_FAKING_GENETIC_MAP]")

	print ("[PRODUCE_FILE_WITH_POPLABELS]")

	# first get the labels of the individuals in the vcf
	sampleFile = pandas.read_csv (input_sample_file, header=None, sep="\t")
	# first two rows are bogus, labels are in first column
	indLabels = sampleFile.iloc[2:,0]

	# now write the poplabels file
	ofs = open (input_poplabels_file, "w")
	# header
	ofs.write ("sample population group sex\n")
	# write a line for each sample
	runningCount = 0
	currPopLabel = 0
	nextPopLabelSwitch = samplesInDemes[currPopLabel]
	for ind in indLabels:
		ofs.write (f"{ind} POP_{currPopLabel} GRP_{currPopLabel} NA\n")
		runningCount += 1
		if runningCount >= nextPopLabelSwitch:
			currPopLabel += 1
			# otherwise, we are done
			if (currPopLabel < len(samplesInDemes)):
				nextPopLabelSwitch += samplesInDemes[currPopLabel]
	ofs.close ()

	print ("[DONE]")


def runRelate (input_haps_file, input_sample_file, input_genetic_map_file, input_poplabels_file, mu, effectiveSize, popSizeIterations, yearsPerGen, binParameters, treeOutputPrefix, tree_log_file, sizeOutputPrefix, size_log_file):

	print ("[RUN_RELATE_TREES]")

	# first estimate local trees using relate
	# ./binaries/relate_v1.1.2_x86_64_static/bin/Relate --mode All -m 1.25e-8 -N 20000 --haps exp_dataset1.haps --sample exp_dataset1.sample --map exp_dataset1.map --seed 4711 --output result_exp_dataset1
	# --memory can be used to control memory usage
	relateTreeCmd = relateBinary + \
		" --mode All " + \
		" -m " + f"{mu:.4e}" + \
		" -N " + f"{int(effectiveSize):d}" + \
		" --haps " + input_haps_file + \
		" --sample " + input_sample_file + \
		" --map " + input_genetic_map_file + \
		" --seed " + f"{numpy.random.randint (0,9999999)}" + \
		" --output " + treeOutputPrefix + \
		f" 2>&1 | tee -a {tree_log_file}"

	# show it
	print (relateTreeCmd)

	# run the command
	process = subprocess.Popen (relateTreeCmd, stderr=subprocess.PIPE, shell=True, close_fds=True)
	process.communicate()

	# in addition to the log-file, this produces the files treeOutputPrefix.mut/.anc which are the input for size inference

	process.terminate()

	print ("[DONE]")

	print ("[RUN_RELATE_POPSIZE]")

	# and then estimate the population size
	# ./binaries/relate_v1.1.2_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh --input result_exp_dataset1 -m 1.25e-8 --poplabels exp_dataset1.poplabels --seed 4711 --output result_size_exp_dataset1 --num_iter 2
	# --threads (could be used for parallelisation)
	relateSizeCmd = relatePopSizeScript + \
		" --input " + treeOutputPrefix + \
		" -m " + f"{mu:.4e}" + \
		" --poplabels " + input_poplabels_file + \
		" --seed " + f"{numpy.random.randint (0,9999999)}" + \
		" --num_iter " + f"{popSizeIterations:d}" + \
		" --years_per_gen " + f"{yearsPerGen:d}" + \
		" --bins " + ",".join([str(x) for x in binParameters]) + \
		" --output " + sizeOutputPrefix + \
		f" 2>&1 | tee -a {size_log_file}"

	# show it
	print (relateSizeCmd)

	# run the command
	process = subprocess.Popen (relateSizeCmd, stderr=subprocess.PIPE, shell=True, close_fds=True)
	process.communicate()
	process.wait()
	process.terminate()
	process.kill()
	print ("[DONE]")




def createHistoryFile (output_coalescent_rate_file, output_history_file, minGen, maxGen, resolution):

	print ("[CREATE_REAL_OUTPUT]")

	# read in the output	
	ifs = open(output_coalescent_rate_file)

	# skip first line
	ifs.readline()
	# second line has generation times
	generTimes = ifs.readline()
	# third line has coalescent rates
	coalRates = ifs.readline()
	# should be everything
	ifs.close ()

	# now get the numbers from the lines
	generTimes = numpy.array([float(x) for x in generTimes.split()])
	# skip the first two, cause they just indices
	coalRates = numpy.array([float(x) for x in coalRates.split()][2:])

	# convert the rates to (diploid) sizes
	# 0.5/rate according to the documentation (seems to work)
	popSizes = 0.5/coalRates

	# make a grid
	times = numpy.exp(numpy.linspace(numpy.log(minGen), numpy.log(maxGen), resolution))

	# now get the sizes
	sizes = numpy.zeros (len(times))

	# where does the grid fall in terms of the right boundaries
	# use the right boundaries
	sizeIdx = numpy.searchsorted (generTimes[1:], times)

	# get the corresponding sizes
	sizesOnGrid = popSizes[sizeIdx]

	# and write it to file
	ofs = open (output_history_file, "w")
	ofs.write ("gener,size\n")
	for i in range(len(times)):
		ofs.write (f"{times[i]},{sizesOnGrid[i]}\n")
	ofs.close()

	print ("[DONE]")


def main():

	# <script_name> seed vcfFile numLoci outDir
	assert (len(sys.argv) == 5)

	# seeding
	numpy.random.seed (int(sys.argv[1]))

	# vcf-input-file
	input_vcf_file = sys.argv[2]

	# num loci
	numLoci = int(sys.argv[3])

	# the output dir
	outDir = sys.argv[4]


	# popgen parameters
	mu = 1.25e-8
	recoProb = 1.25e-8

	# relate parameters
	# how many years per generation
	# default: 28
	# If we use 28 here, then the output is in actual years
	# If we use 1 here, then the output seems to be in actual number of generations
	yearsPerGen = 1

	# effective size (kinda like an arbitrary scaling factor, haploid)
	effectiveSize = 2 * 100

	# num iterations when estimating the popsize
	# default: 10
	popSizeIterations = 10
	# popSizeIterations = 5

	# parameter for the piecewise bins
	# number are given in years (or generations, if yearsPerGen = 1)
	# format is [lower, upper, stepsize]
	# changepoints for piecewise history: [0, 10^(lower), 10^(lower+stepsize), 10^(lower+2*stepsize), ..., 10^(upper), infty)
	# i have no idea what the default is
	binParameters = [2,7,0.25]


	# just the name of the file, without extensions or pathname
	baseFilename = os.path.splitext (os.path.basename (input_vcf_file))[0]

	# for some stupid reason we do have to get the absolute path to the vcf
	# OR MAYBE NOT, BUT CURRENTLY IT ONLY WORKS IF IT IS EXECUTED FROM OUTDIR

	# also switch to outDir, so maye relate stuff lands there
	# os.chdir (outDir)
	# relate input file
	analysisDir = "./"
	# I think we need to make it outside
	# os.makedirs (analysisDir)
	input_haps_file = os.path.join (analysisDir, f"{baseFilename}.haps")
	input_sample_file = os.path.join (analysisDir, f"{baseFilename}.sample")
	input_poplabels_file = os.path.join (analysisDir, f"{baseFilename}.poplabels")
	input_genetic_map_file = os.path.join (analysisDir, f"{baseFilename}.map")

	# relate output files (or better: prefixes)
	# treeOutputPrefix = os.path.join (analysisDir, f"result_tree_{baseFilename}")
	treeOutputPrefix = f"result_tree_{baseFilename}"
	tree_log_file = treeOutputPrefix + ".log"
	# sizeOutputPrefix = os.path.join (analysisDir, f"result_size_{baseFilename}")
	sizeOutputPrefix = f"result_size_{baseFilename}"
	size_log_file = sizeOutputPrefix + ".log"

	# this is the final output from relate that we can use to create a history
	output_coalescent_rate_file = sizeOutputPrefix + ".coal"

	# this is where we put the history
	output_history_file = sizeOutputPrefix + ".csv"


	# input parameters


	# prepare the proper input for relate from the vcf
	prepareInput (input_vcf_file, input_haps_file, input_sample_file)

	# # other files to prepare
	# prepareAuxiliaryFiles (recoProb, numLoci, input_genetic_map_file, input_sample_file, input_poplabels_file, [25, 25])
	# run smaller samples, this should work
	prepareAuxiliaryFiles (recoProb, numLoci, input_genetic_map_file, input_sample_file, input_poplabels_file, [10000, 10000])

	# do the analysis
	runRelate (input_haps_file, input_sample_file, input_genetic_map_file, input_poplabels_file, mu, effectiveSize, popSizeIterations, yearsPerGen, binParameters, treeOutputPrefix, tree_log_file, sizeOutputPrefix, size_log_file)

	# # create a history for the popsize
	# createHistoryFile (output_coalescent_rate_file, output_history_file, minGen, maxGen, resolution)


if __name__ == "__main__":
	main()
