This repository contains the workflows, data, and code that accompany the article 'Resolving marsupial phylogenomics by evaluating congruence across loci'
--------------------------------------------------------------------------------------------------------------------------------------------

19 September 2017

Bioinformatics
--------------

The bioinformatics directory contains data and code associated with processing of the sequencing products.

It contains the three directories, containing data and workflows as follows:

1. 	    	      		   targets

This directory contains two fasta files

marsupialTargets.fa -- these are the targets sequences that were used to design probes for sequence capture
                    -- sequence names begin with the putative Sarcophilus orthologue, followed by sample information for the corresponding transcriptome

targetExons.fa      -- the target exon sequences, from Sarcophilus only

2.		       	   	   clean

perl and shell scripts used to clean illumina sequence reads

3.   	       	       	       	   assembly

perl and shell scripts used in the assembly of cleaned illumina sequence reads
and the identification and phasing of heterozygous sites

Scripts are provided in this repository are for archival purposes. See also:
https://github.com/MVZSEQ
https://github.com/jasongbragg/exon-capture-phylo


Phylogenomics
-------------

This directory contains the data and results files associated with the phylogenomic analyses.

There are four directories as follows:

1, 2, and 3.			   loci_clean, loci_clean_12, loci_clean_3

These three directories correspond to analyses with all the alignment sites, first and second codon positions, and third codon positions, respectively.

There is one directory for the ML analysis of each locus. The .Rdata files in these directories are to be loaded in the R computing environment. Those with the loci names followed by .Rdata contain the results for PhyML analyses, while those called adeqres.Rdata contain the results for the model adequacy analyses described in the manuscript. These folders also contain the output of estimates of dN/dS from the software PAML.

There is a directory called further_analyses which contains the results of clustering and astral analyses.

loci_clean contains an additional three directories with prefix concatenated.analyses, which contain the maximum likelihood results for concatenated analyses for each of the complete data sets, with either all the data, the first and second codon positions, or the third codon position as indicated by the suffix of the folder name.

4. 	   	    	  	   Dating

The directory loci contains the loci that were used for the dating analyses. Analyses were performed with code names to facilitate preparation for MCMCtree analyses. The codes used can be found in the R file astraltrees.Rdata.

The directories with the MCMCtree suffix contain the input and output files of replicate dating analyses.