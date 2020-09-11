# SingleCellProject
The set of scripts in this project were developed to infer recombination breakpoints from hundreds of single sperm cells, 
each sequenced at very low coverage, and the corresponding high coverage (>15x) WGS of the donor. 

The single sperm cells' VCF is used both to phase the parental genome and infer the recombination breakpoints of each single
cell. The WGS of the donor is required to know with confidence the set of heterozygous sites (hetSNPs) per chromosome to be phased by 
using the (haploid) sperm cells. 


The inference of recombination breakpoints from SC involves three steps:

1) generate input files (GT output files from vcftools)
   - retrieving common sites between bulk vcf and single cell vcfs and extracting such sites from the SC vcf with: merging_SCvcf_hetSNPs.bash
   - removing sites falling into regions of low mapping quality and keepign accessible regions of the genome with: subsetting_chr_extracting.bash.
   - this last script also splits chromosomes (chromosomes: 1-7) into short and long arms.
   
2) infer the parental origin (hidden state) of every allele along the single cell genome.
   This step is done using two methods that mainly differ on the way parental chromosomes are phased. 
   
   Method 1 looks at "high quality sites" (high calling confidence and observed across multiple cells) and uses the most common pairwise combination of alleles to form haplotypes. There are often multiple haplotypes as there are not enough "high quality sites" to link the haplotypes to each other. In a second step sites with less quality, mainly sites observed in a small number of cells, are used to link these initial haplotypes. Once parental chromosomes are phased they are used to infer the chromosome of origin of each allele along the single cell genome using the HMM. Using this method is it not straightforward to compute error rates needed for the HMM. Therefore, I used 0.15 as average error rate. Furthermore, this method involves comparing all pairs of sites which can lead to very long running times. See below: Running method1
   
   Method 2 compares every single cell against each other and for each pair of cells the proportion of sites that are similar. Non-recombining cells are those for which most of the comparisons show either very low  (similarity < 5%) or very high (similarity > 95%) similarity. These cells are then separated them into two groups. The most common allelic combination among cells within each group represents a parental chromosome. Error rates can then be computed by computing how many times a non-parental allele is seen among cells within each group of cells. The error rate per site is then used in the HMM. See below: Running method2
   
   Both methods result in a matrix of inferred states across hetSNPs (rows) and cells (columns) which are used to compute recombination rates. 
   
3) infer recombination maps from the inferred states 
    
    - From the matrices inferred above one can then compute recombination rates within Xbp (default: 1Mb) windows along the chromosomes. 
    - The script recombinationRates_and_plotting.R outputs a list of recombining cells and their breakpoint intervals (each pair of chromosome coordinates indicate start and end of the breakpoint interval or more specificaly the last position to be inferred as belonging to parental One and the first position to be inferred to belong to parental Two), recombination rates along the chromosome (recombRates_SNPDensity_chr*.pdf), comparison between physical and genetic distance (genetic_vs_physicalDist_chr*.pdf), number of breakpoints per genomic window (IndID.chrID.stateSeqHMM_clean_bp_counts.txt) and number of cells per genomic window (IndID.chrID.stateSeqHMM_clean_cell_counts.txt), plots with breakpoint resolution (genomic distance between informative sites with alternative parental origin) per individual and the summary of this information across individuals (summaryBp_chr*.txt). 


Running method1
(because method2 performs better the description of method1 will come later)

The algorithm used to phase the parental chromosomes from the single sperm cells is an adaptation from Lu, Zong and Fan et al
Science 2012. Briefly, in a first step it builds up small haplotypes based on variants called with high confidence. Secondly, it
ancors variants called with low confidence to the haplotypes assembled in the first step. By doing so one is able to reconstruct
the full parental chromosomes. 
###------------

Running method2

run phasing and inference of parental allelic origin:


qsub -P awadallalab -V -cwd -b y -pe smp 5 -l h_vmem=8G ./launching_sc_method2.bash <INDid> <chrNb>





2) running the phasing of parental chromosomes from the multiple SCs (longest part)

3) comparing the phased parental chromosomes with each single cell (rough estimation no HMM)

4) running an HMM model to identify the breakpoints by taking into account error rates


To run this algorithm it is required to compare the individual's variant calling from WGS with the variant calling 
from the single cells as a bulk (multiple replicates option in GATK). This produces a set of heterozygous SNPs (hetSNPs) which
will then be used to subset the single cell VCF. (futher details can be found in the README within prep_inputF). 



A visual inspection of the phasing can be done by looking at the comparison between each single cell and the reconstructed 
parental chromosomes. This is done with a sliding window approach, with windows of 10 SNPs moving one SNPs at a time and asking
how many alleles within a 10 SNP window are compatible with each of the parental chromosome. 


