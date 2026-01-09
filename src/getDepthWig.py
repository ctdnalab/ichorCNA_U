import sys,subprocess
import pandas as pd


def main():
	usage = """
	<Input file (bed/tsv/cram/bam>
	<Chromosome size tsv (two column)>
	<Ref genome fasta>
	<Output file path>
	<Bin size (in kb)>
	<Tabix path>
	<Samtools path>
	<Chroms to analyze (comma separated)>"
	"""
	argv = sys.argv[1:]
	try:
		inputFile = argv.pop(0)
		chromSizeFile = argv.pop(0)
		refGenome = argv.pop(0)
		outName = argv.pop(0)
		binSize = int(argv.pop(0))
		tabixPath = argv.pop(0)
		samtoolsPath = argv.pop(0)
		goodChroms = argv.pop(0).split(",")
	except:
		sys.exit(usage)

	dfChromSize = pd.read_csv(chromSizeFile,sep='\t',header=None)

	binSize = binSize * 1000
	binSets = []
	for r in dfChromSize.itertuples():
		if r._1 not in goodChroms:
			continue
		allBinStart = list(range(1,r._2,binSize))
		allBinEnd = [x + binSize -1 for x in allBinStart]
		if allBinEnd[-1] > r._2:
			allBinEnd[-1] = r._2
		bins = [ [r._1,allBinStart[x],allBinEnd[x]] for x in range(len(allBinStart)) ]
		binSets.append(bins)
	dfBins = pd.DataFrame(sum(binSets,[]),columns=["chrom","binStart","binEnd"])

	allChromBins = dfBins.groupby('chrom',sort=False)
	outF = open(outName,"w")
	for chromBins in allChromBins:
		outF.write(f"fixedStep chrom={chromBins[0]} start=1 step={binSize} span={binSize}\n")
		for i in chromBins[1].itertuples():
			region = f"{i.chrom}:{i.binStart}-{i.binEnd}"
			if inputFile.endswith(".bed.gz") or inputFile.endswith(".bedpe.gz") or inputFile.endswith(".bed.bgz") or inputFile.endswith(".bedpe.bgz"):
				cmd = f"{tabixPath} -p bed {inputFile} {region} | wc"
				readCount = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE).stdout.decode('utf-8')
				readCount = outF.write(readCount.split()[0]+"\n")
			elif inputFile.endswith(".tsv.gz") or inputFile.endswith(".tsv.bgz"):
				cmd = f"{tabixPath} -s 1 -b 2 -e 3 {inputFile} {region} | wc"
				readCount = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE).stdout.decode('utf-8')
				readCount = outF.write(readCount.split()[0]+"\n")
			elif inputFile.endswith(".cram") or inputFile.endswith(".bam"):
				cmd = [samtoolsPath, 'view', '-cT',refGenome,inputFile,region]
				readCount = subprocess.run(cmd,stdout=subprocess.PIPE).stdout.decode('utf-8')
				readCount = outF.write(readCount)
	outF.close()
	return

main()