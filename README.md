# SingleCellProject
The set of scripts in this project were developed to infer recombination breakpoints from hundreds of single sperm cells, 
each sequenced at very low coverage, and the corresponding high coverage (>15x) WGS of the donor. 

The single sperm cells' VCF is used both to phase the parental genome and infer the recombination breakpoints of each single
cell. The WGS of the donor is required to know with confidence the set of heterozygous sites (hetSNPs) per chromosome to be phased by 
using the (haploid) sperm cells. 

Single sperm cells genotypes should be called in two ways:
1) bulk genotype calling (vcf file containing 1 diploid call)

   Individual's variant calling from single cells as a bulk is done with the multiple replicates option in GATK. This produces a set of heterozygous SNPs (hetSNPs) which will then be used to subset the single cell .vcf. In this step one can apply filters like minimum number of cells (or minor allele frequency) where the allele was observed.  

2) genotype calling per single cell (vcf file containing X single cells)



The inference of recombination breakpoints from SC involves three steps:

1) generate input files (GT output files from vcftools)
   - retrieving common sites between bulk vcf and single cell vcfs and extracting such sites from the SC vcf with: merging_SCvcf_hetSNPs.bash
   - removing sites falling into regions of low mapping quality and keepign accessible regions of the genome with: subsetting_chr_extracting.bash.
   - this last script also splits chromosomes (chromosomes: 1-7) into short and long arms.
   
2) infer the parental origin (hidden state) of every allele along the single cell genome.
   This step is done using two methods that mainly differ on the way parental chromosomes are phased. 
   
   Method 1 looks at "high quality sites" (high calling confidence and observed across multiple cells) and uses the most common pairwise combination of alleles to form haplotypes. By the end of this step there are often multiple haplotypes as there are not enough "high quality sites" to link the haplotypes to each other. In a second stepm sites with less quality, mainly sites observed in a small number of cells, are used to link the set of haplotypes built with high confidence sites. Once parental chromosomes are phased they are used to infer the chromosome of origin of each allele along the single cell genome using the HMM. Using this phasing method is it not straightforward to compute error rates needed for the HMM. Therefore, I used 0.15 as average error rate. Furthermore, this method involves comparing all pairs of sites which can lead to very long running times. See below: Running method1
   
      Method 2 compares every single cell against each other and computes the proportion of matching alleles for each pair of cells. Non-recombining cells are  show either very low (similarity < 5%) or very high (similarity > 95%) similarity with each other. Non-recombining cells are separated into two haplotypic groups based on those similarity proportions. The most common allelic combination among cells within each group represents a parental chromosome. Error rates can then be estimated by computing how many times, within a set of cells with a given haplotypic configuration, the alternative allele is present among those cells. The error rate per site is then used in the HMM. See below: Running method2
   
   Both methods result in a matrix of inferred states across hetSNPs (rows) and cells (columns) which are used to compute recombination rates. 
A visual inspection of the phasing is be done by comparing each single cell and the reconstructed parental chromosomes (before and after HMM). This is done using a sliding window approach, with windows of 10 SNPs moving one SNPs at a time and asking how many alleles within a 10 SNP window are compatible with each of the parental chromosome (SC_phasedHapHMap_<indID>_<chrTag>_QCvar25_NbCells10_NbLinks5_4_*.pdf)

   
3) infer recombination maps from the inferred states 
    
    - From the matrices inferred above one can then compute recombination rates within Xbp (default: 1Mb) windows along the chromosomes. 
    - The script recombinationRates_and_plotting.R outputs a list of recombining cells and their breakpoint intervals (each pair of chromosome coordinates indicate start and end of the breakpoint interval or more specificaly the last position to be inferred as belonging to one parental and the first position to be inferred to belong to the second parental), recombination rates along the chromosome (recombRates_SNPDensity_chr*.pdf), comparison between physical and genetic distance (genetic_vs_physicalDist_chr*.pdf), number of breakpoints per genomic window (IndID.chrID.stateSeqHMM_clean_bp_counts.txt) and number of cells per genomic window (IndID.chrID.stateSeqHMM_clean_cell_counts.txt), plots with breakpoint resolution (genomic distance between informative sites with alternative parental origin) per individual and the summary of this information across individuals (summaryBp_chr*.txt). 


Running method1
(because method2 performs better the description of method1 will come later)

The algorithm used to phase the parental chromosomes from the single sperm cells is an adaptation from Lu, Zong and Fan et al
Science 2012. Briefly, in a first step it builds up small haplotypes based on variants called with high confidence. Secondly, it
ancors variants called with low confidence to the haplotypes assembled in the first step. By doing so one is able to reconstruct
the full parental chromosomes. 
###------------

Running method2

1) run phasing and inference of parental allelic origin:
   
   Rscript ./phasing_recomb_method2_august.R indID chrID
   example: Rscript ./phasing_recomb_method2_august.R indID chr21

2) run recombination rate computation and corresponding plots

   Rscript ./recombinationRates_and_plotting.R chrNb
   example: Rscript ./recombinationRates_and_plotting.R 21

Note: to run both scripts one need to have installed the packages: parallel, filesstrings and RColorBrewer. At the moment some function called by 
recombinationRates_and_plotting.R are contained in the breakpoint_functions.R script whose location needs to be provided. 

####----------------------




