# map_HIV

Overview

This project focuses on mapping HIV integration sites within the human genome using BWA. It identifies insertion sites, refines them based on specific criteria, and plots them relative to chromosome lengths. This analysis helps in understanding the distribution and potential impacts of viral integration across the genome.

Objectives

Map HIV sequences to the human genome using BWA.
Output insertion sites and filter for high-quality mappings.
Visualize the distribution of insertion sites relative to chromosome lengths, with a focus on chromosomes 17 and 19.

Prerequisites

To get started, ensure you have the following:

A working installation of BWA (Burrows-Wheeler Aligner).
The FASTA file of the human genome (reference genome).
Workflow
Mapping HIV Reads:

Use BWA to map HIV sequences to the human genome. The resulting alignment file should contain all mapped reads.
Extracting High-Quality Mappings:

Filter the mapped reads to retain only those with a mapping quality (MAPQ) score ≥ 20.
Separate mappings into two files: one for forward strand and one for reverse strand.
Extracting Chromosome and Position Information:

For each mapped read, extract the chromosome and position information.
Add this information to the first or last column of the output file.
Refining Insertion Sites:

Combine insertion sites that are within a range of ±1 or ±2 base pairs into a single insertion site.
If multiple similar insertion sites are found, retain only one.
Visualization:

Generate visual plots showing the distribution of insertion sites across the genome, normalized by chromosome length.
Pay special attention to the insertion sites on chromosomes 17 and 19 and compare their distribution with other chromosomes.

