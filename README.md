# SingleCellProject
The set of scripts in this folder were developed to infer recombination breakpoints from hundreds of single sperm cells, 
each sequenced at very low coverage, and the corresponding high coverage (>15x) WGS of the donor. 

The single sperm cells' VCF is used both to phase the parental genome and infer the recombination breakpoints of each single
cell. The WGS of the donor is required to know with confidence the set of heterozygous sites per chromosome to be phased by 
using the (haploid) sperm cells. 


The inference of recombination breakpoints from SC involves three steps:

1) generating input files (GT output files from vcftools)

2) running the phasing of parental chromosomes from the multiple SCs (longest part)

3) comparing the phased parental chromosomes with each single cell (rough estimation no HMM)

4) running an HMM model to identify the breakpoints by taking into account error rates

Scripts required for step ONE are in the folder : prep_inputF
Scripts required for step TWO are in the folder: phasing_vs3
Scripts required for step THREE & FOUR are in the folder: infer_recBreakPoints

To run this algorithm it is required to compare the individual's variant calling from WGS with the variant calling 
from the single cells as a bulk (multiple replicates option in GATK). This produces a set of heterozygous SNPs (hetSNPs) which
will then be used to subset the single cell VCF. (futher details can be found in the README within prep_inputF). 

The algorithm used to phase the parental chromosomes from the single sperm cells is an adaptation from Lu, Zong and Fan et al
Science 2012. Briefly, in a first step it builds up small haplotypes based on variants called with high confidence. Secondly, it
ancors variants called with low confidence to the haplotypes assembled in the first step. By doing so one is able to reconstruct
the full parental chromosomes. (further details can be found in the README within phasing_vs3). 

A visual inspection of the phasing can be done by looking at the comparison between each single cell and the reconstructed 
parental chromosomes. This is done with a sliding window approach, with windows of 10 SNPs moving one SNPs at a time and asking
how many alleles within a 10 SNP window are compatible with each of the parental chromosome. 


