This repository contains the workflows, data, and code that accompany the article 'Resolving marsupial phylogenomics by evaluating congruence across loci'
--------------------------------------------------------------------------------------------------------------------------------------------

7 April 2017

Bioinformatics
--------------

This repository contains data and code associated with the bioinformatics in the article:

There are three directories, containing data and workflows:

1. targets -------------------------

This directory contains two fasta files

marsupialTargets.fa -- these are the targets sequences that were used to design probes for sequence capture
                    -- sequence names begin with the putative Sarcophilus orthologue, followed by sample information for the corresponding transcriptome

targetExons.fa      -- the target exon sequences, from Sarcophilus only

------------------------------------


2. clean ---------------------------

perl and shell scripts used to clean illumina sequence reads

------------------------------------


3. assembly ------------------------

perl and shell scripts used in the assembly of cleaned illumina sequence reads
and the identification and phasing of heterozygous sites

------------------------------------


Scripts are provided here primarily for archival purposes. See also:
https://github.com/MVZSEQ
https://github.com/jasongbragg/exon-capture-phylo


Phylogenomics
-------------

This repository contains the phylogenetics data and code associated with the article:

There are four directories containing data:

1. loci_clean -------------------------

This directory contains one folder for the ML analysis of each locus.

------------------------------------


2. loci_clean_12 and loci_clean_3 ---------------------------



------------------------------------


3. Dating ------------------------



------------------------------------