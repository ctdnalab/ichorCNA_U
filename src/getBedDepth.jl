using QuickArgParse
using CSV,DataFrames

function main()
	req = ["InputFile","BinSize","ChromSize","OutFile"]
	reqhelp = ["Path to input BED file","Bin size in kb",
	"Tab-delimited file of chrom sizes","Output wig file path"]
	opts = ["bedtoolspath","samtoolspath","chroms"]
	optsHelp=["Path to bedtools","Path to samtools",
	"Chromosomes to analyze. Comma separated, can include ranges (ie 1:22,X). ALL for all in chrom size file"]
	optsDefault=["bedtools","samtools","ALL"]
	flags = ["p"]
	flagHelp = ["Add 'chr' prefix to all chrom names."]
	title = "BED/XAM Bin Read Counter"
	desc = "Counts total reads in fixed size bins from a bed,bam,cram etc file. 
Generates a temp bed file of the bins. Requires bedtools and (if run on bam/cram file) samtools."

	R = process_reqs(req;reqHelp=reqhelp,title=title,desc=desc,
		optTag=opts,optHelp=optsHelp,optDefault=optsDefault,
		flags=flags,flagHelp=flagHelp)
	build_usage!(R)
	A = parse_args(R)

	inputFile = A["InputFile"]
	outName = A["OutFile"]
	bedtoolsPath = A["bedtoolspath"]
	samtoolsPath = A["samtoolspath"]
	if occursin(".wig",outName) == false
		outName = outName * ".wig"
	end

	cSize = readlines(open(A["ChromSize"],"r"))
	So = String[]           #Order of chroms in the size file
	S = Dict{String,Int}()  #Size lookup dict
	for ln in cSize
		if length(ln) == 0
			continue
		end
		data = [string(x) for x in split(ln,"\t")]
		push!(So,data[1])
		S[data[1]] = parse(Int,data[2])
	end

	#Process chrom size input file and generate bins
	binSize = parse(Int,A["BinSize"]) * 1000
	fileName = replace(basename(outName),".wig"=>"")
	tempBed = fileName * ".ichorCNAtempFixedBins.bed"
	bedF = open(tempBed,"w")
	for c in So
		allS = collect(range(1,stop=S[c],step=binSize))
		allE = [x + binSize for x in allS]
		bins = [ [allS[x],allE[x]] for x in 1:length(allS) ]
		for b in bins
			write(bedF,"$c\t$(b[1])\t$(b[2])\n")
		end
	end
	close(bedF)

	#Process chromosome inputs
	if A["chroms"] == "ALL"
		allChrom = So
	else
		allChrom = String[]
		rawChrom = [string(x) for x in split(A["chroms"],",")]
		rF = r"\d+\:\d+"
		for c in rawChrom
			if isnothing(match(rF,c)) == false && length(split(c,":")) == 2
				cR = [parse(Int,x) for x in split(c,":")]
				allC = [string(x) for x in collect(range(cR[1],stop=cR[2]))]
				if A["p"] == true
					allC = ["chr" * x for x in allC]
				end
				for c in allC
					push!(allChrom,c)
				end
			else
				if A["p"] == true
					c = "chr" * c
					push!(allChrom,c)
				else
					push!(allChrom,c)
				end
			end
		end
	end

	#cmd = `$toolPath coverage -counts -sorted -g $(A["ChromSize"]) -a $tempBed -b $inputFile`
	if endswith(inputFile,".bed") || endswith(inputFile,".bedpe") || endswith(inputFile,".tsv")
		cmdFinal = `$bedtoolsPath coverage -counts -sorted -g $(A["ChromSize"]) -a $tempBed -b $inputFile`
	elseif endswith(inputFile,".cram") || endswith(inputFile,".bam") || endswith(inputFile,".sam")
		chromSet = join([x for x in allChrom]...," ")
		tempSam = fileName * ".collateTemp.samtools"
		cmdA = cmd("$samtoolsPath view -F 3852 $inputFile $chromSet")
		cmdB = `$samtoolsPath collate -T $tempSam -`
		cmdC = `$bedtoolsPath coverage -counts -sorted -g $(A["ChromSize"]) -a $tempBed -b stdin`
		cmdFinal = pipeline(cmdA,cmdB,cmdC)
	end
	colTypes = Dict(1=>String,4=>Int)
	countData = CSV.File(open(cmdFinal,"r"),delim='\t',types=colTypes,header=false) |> DataFrame
	
	outF = open(outName,"w")
	for c in allChrom
		cData = filter(x->x.Column1 == c,countData)
		write(outF,"fixedStep chrom=$c start=1 step=$binSize span=$binSize\n")
		for r in eachrow(cData)
			write(outF,"$(r[:Column4])\n")
		end
	end
	close(outF)
	run(`rm $tempBed`)

	return
end

main()
