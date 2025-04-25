# ichorCNA_U

ichorCNA_U is a modified version of ichorCNA v0.3.2 that removes hard-coded human-specific chromosome name requirments, 
enabling ichorCNA to be run on a wider range of species. It also adds a slurm-compatible wrapper for high-throughput analysis of 
many samples, and generation of a summary table with tumor fractions and stats from the highest likelihood solutions.

The original ichorCNA source code, along with ichorCNA usage information and best practices, can be found here:
https://github.com/broadinstitute/ichorCNA

## Citations
If ichorCNA_U is used in published research, please cite:  
>Favaro, Patricia F, et al. Feasibility of circulating tumor DNA analysis in dogs with naturally-occurring malignant and benign splenic lesions. (2022) Scientific Reports doi:10.1038/s41598-022-09716-6

All studies using ichorCNA_U should also cite the original ichorCNA publication:  
>Adalsteinsson, Ha, Freeman, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. (2017) Nature Communications doi:10.1038/s41467-017-00965-y


## Dependencies
R:
 - R >= 3.6.0
 - GenomicRanges
 - HMMcopy
 - optparse

Python:
 - Python3
 - YAML

Julia:
 - CSV
 - DataFrames
 - Query
 - QuickArgParse https://github.com/brmcdonald/QuickArgParse  

3rd party tools:  
Bedtools and Samtools, to calculate bin read depths  
https://bedtools.readthedocs.io/en/latest/index.html  
https://www.htslib.org/  

HMMcopy_utils, to calculate gc content for processing reference genomes  
https://github.com/shahcompbio/hmmcopy_utils  
GenMap, for reference genome mappability calculations  
https://github.com/cpockrandt/genmap  

Note that the ichorCNA package does not need to be separately installed for R, ichorCNA_U only uses the code provided in this repository.

## Usage

The tarball ichorCNA_U/src_ichorCNA/inst.tar.gz must be decompressed to use the default panel of normal reference data included with ichorCNA.

ichorCNA_U reads all inputs from params_IchorCNA.yaml; IchorCNA parameters, input and 
output paths, and options selecting steps to run are all included there. 
Sample read data can be in bam, cram or bed file formats.  

There are two major steps: calculating read depth in each window to generate wig files and running
ichorCNA on the wig files to infer tumor fraction. Slurm job scripts are generated for each file
for each step. If slurm is not available, the job scripts are treated as bash scripts and run
locally.  

To run the pipeline, set all parameters and use:
```Bash
python ichorCNA_U.py params_ichorCNA.yaml
```

For interpretation of ichorCNA outputs and selection of the best tumor fraction solution, see the documentation for ichorCNA at the link above.

## Details and Limitations
 - The primary changes to ichorCNA were removal of all references to seqInfo from the source code, which reformats chromosome names to fit
standard formats for the human genome and breaks if other chromosome names are found.
 - Bedtools is used to calculate bin read depths instead of hmmcopy_util's readCounter (used in the original ichorCNA publication), because it enables a wider range of input file types and is actively maintained.
 - Sex chromosome adjustment still assumes human chromosomes X and Y. Analysis can/should be limited to autosomal chromosomes only.
 - Since centromere location information is not often available for less well characterized genomes, windows overlapping 
them are not automatically excluded. The window mappability filter may alleviate this problem somewhat.
 - Some reference genomes still cause problems with the modified ichorCNA, for reasons that are not yet clear.  

## Reference processing using GenMap, HMMcopy_utils, and Bedtools
In order to build a new panel of normals, ichorCNA needs wig files for mappability and gc content of each window plus wig files of read depths from multiple
normal samples. A new set of files needs to be generated for each window size. IchorCNA includes 500kb and 1000kb windows for hg19 and hg38 human genome builds, 
and the following example builds files for 500kb windows. All the necessary files for a new panel of normals can be generated using the following tools:

 - hmmcopy:			https://github.com/Bioconductor/copy-number-analysis/wiki/HMMcopy
 - wigToBigWig:		http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
 - GenMap:			https://github.com/cpockrandt/genmap
 - a tab-delimited file of all chromosomes and their sizes (for wigToBigWig)
 - Aligned bam files of healthy control plasma samples

#### 1. Reference genome gc content wig file via hmmcopy_utils gcCounter
By default, gcCounter doesn't calculate gc content for windows that contain any Ns. Because ichorCNA windows are large, 
for some lower quality genomes this excludes almost all windows. The -f flag is not listed in the usage statement, 
but sets gcCounter to calculate gc content for all windows regardless of Ns. Windows with too many Ns often have poor
mappability and end up being excluded by ichorCNA from downstream analysis.
```Bash
gcCounter -f -w 500000 /path/to/refgenome.fna > refGenome.gc.500k.wig
```
#### 2. Reference genome mappability wig file via genmap and hmmcopy_utils mapCounter
```Bash
genmap index -F /path/to/refgenome.fna -I genmap.index
genmap map -K 50 -E 1 -w -T [num_cores] -I genmap.index -O genmap.k50.e1.wig
wigToBigWig genmap.k50.e1.wig chrom.sizes genmap.k50.e1.bw
mapCounter -w 500000 ./genmap.k50.e1.bw > refGenome.map.500k.wig
```
#### 3. Healthy control sample bin read depth wig files via bedtools/ichorCNA_U
This can be done by using the healthy sample bams as inputs for ichorCNA_U with "runDepthCalc" 
set to true while "runIchorCNV" and "calculateStats" are set to false. In that case all the ichorCNA-specific inputs are 
ignored and only the wig files are generated.

#### 4. Run createPanelOfNormals.R from src_ichorCNA/scripts
The inputs at this step are the reference genome files generated in steps 1 and 2 along with a text file containing the paths
to all of the normal sample read depth wig files from step 3. Chromosomes should be specified here (as an array specified in R code in quotes) 
to avoid trying to normalize any small contigs that may be in the genome assembly. Regions to be excluded (centromeres etc) 
should also be provided here if available. See src_ichorCNA/inst/exdata/GRCh37.p13_centromere_UCSC-gapTable.txt for format.
```Bash
Rscript /path/to/src_ichorCNA/scripts/createPanelOfNormals.R \  
  -f normalwig_paths.txt \  
  --gcWig refGenome.gc.500k.wig \  
  --mapWig refGenome.map.500k.wig \  
  --chrs "c(1:22)" \  
  -o refGenome.PoN \  
```
The resulting **refGenome.PoN_median.rds** file is the reference file used for the normalCtrl parameter in params_IchorCNA.yaml


