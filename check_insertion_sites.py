import pyranges as pr
import matplotlib.pyplot as plt
import sys
import pandas as pd

def extract_gene_coordinates(gtf_filename):
    gtf = pr.read_gtf(gtf_filename)
    return gtf
def extract_insertion_sites(data_filename):
    insertion_sites_df = pd.read_csv(data_filename, sep=' ', header=None)
    insertion_sites_df.columns = ['Chromosome', 'Start', 'Sequence']
    insertion_sites_df['End'] = insertion_sites_df['Start']
    insertion_sites = pr.PyRanges(insertion_sites_df)
    return insertion_sites

#count and plot insertions overlapping with a gene region and normalize based on length of chromosome
def count_insertions_in_genes(data_filename, gene_coordinates):
    insertion_sites = extract_insertion_sites(data_filename)
    genes = gene_coordinates[gene_coordinates.Feature == "gene"]
    overlap_results = insertion_sites.overlap(genes)
    # Count the number of insertions inside genes
    insertion_in_gene = overlap_results.df[['Chromosome', 'Start', 'End', 'Sequence']]
    insertion_count = insertion_in_gene.shape[0]
    print(f"Number of insertion sites in genes: {insertion_count}")
    insertion_count_by_chromosome = insertion_in_gene.groupby('Chromosome').size() 
    
    # Normalize insertion counts by chromosome length
    normalized_insertion_counts = insertion_count_by_chromosome / pd.Series(chromosome_lengths)

    # Sort chromosomes in order
    chromosomes_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    insertion_count_by_chromosome = normalized_insertion_counts.reindex(chromosomes_order, fill_value=0)

    # Plot the count of insertion sites for each chromosome
    insertion_count_by_chromosome.plot(kind='bar', color='skyblue')
    plt.xlabel('Chromosome')
    plt.ylabel('Insertion Sites Count')
    plt.title('Insertion Sites in Genes by Chromosome')
    plt.show()

#count and plot insertions overlapping with a promoter region and normalize based on length of chromosome
def count_insertions_in_promoters(data_filename, gene_coordinates):
    insertion_sites = extract_insertion_sites(data_filename)
    genes = gene_coordinates[gene_coordinates.Feature == "gene"]
    ## Define promoter regions (e.g., 1000 base pairs upstream of TSS)
    promoter_regions = genes.df.assign(Start=genes.df.Start - 1000, End=genes.df.Start)
    promoter_regions = pr.PyRanges(promoter_regions)
    # Check Overlaps with Promoters
    overlap_results_promoters = insertion_sites.overlap(promoter_regions)
    insertion_in_promoters = overlap_results_promoters.df[['Chromosome', 'Start', 'End', 'Sequence']]
    insertion_count = insertion_in_promoters.shape[0]
    print(f"Number of insertion sites in promoters: {insertion_count}")
    insertion_count_by_chromosome = insertion_in_promoters.groupby('Chromosome').size()

    # Normalize insertion counts by chromosome length
    normalized_insertion_counts = insertion_count_by_chromosome / pd.Series(chromosome_lengths)

    # Sort chromosomes in order
    chromosomes_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    insertion_count_by_chromosome = normalized_insertion_counts.reindex(chromosomes_order, fill_value=0)

    # Plot the count of insertion sites for each chromosome
    insertion_count_by_chromosome.plot(kind='bar', color='skyblue')
    plt.xlabel('Chromosome')
    plt.ylabel('Insertion Sites Count')
    plt.title('Insertion Sites in promoters by Chromosome')
    plt.show()

def count_insertions_perGeneCoverage(data_filename, gene_coordinates):
    insertion_sites = extract_insertion_sites(data_filename)
    genes = gene_coordinates[gene_coordinates.Feature == "gene"]
    #Count the number of insertions inside genes
    overlap_results = insertion_sites.overlap(genes)
    insertion_in_gene = overlap_results.df[['Chromosome', 'Start', 'End', 'Sequence']]
    insertion_count = insertion_in_gene.shape[0]
    # Normalize insertion counts by chromosome length
    # Sort chromosomes in order
    chromosomes_order = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    import pdb;pdb.set_trace()
    coverage = apply_coverage(genes, chromosomes_order)
    insertion_count_by_chromosome = insertion_in_gene.groupby('Chromosome').size()
    import pdb;pdb.set_trace()
    normalized_insertion_counts = insertion_count_by_chromosome / pd.Series(coverage)
    insertion_count_by_chromosome = normalized_insertion_counts.reindex(chromosomes_order, fill_value=0)
    # Plot the count of insertion sites for each chromosome
    insertion_count_by_chromosome.plot(kind='bar', color='skyblue')
    plt.xlabel('Chromosome')
    plt.ylabel('Insertion Sites Count')
    plt.title('Insertion Sites in Genes by Chromosome')
    plt.show()

def apply_coverage(genes, chromosomes_order):
   coverage = dict()
   for chrom in chromosomes_order: 
      cov = calculate_coverage(genes[genes.Chromosome==chrom])
      coverage[chrom] = cov
   return coverage

def calculate_coverage(clustered_genes):
   start = 0
   end = 0
   coverage = 0
   for i, j in zip(clustered_genes.Start, clustered_genes.End):
      if i > end:
         start = i
         end = j
         coverage += (end - start + 1)
      elif j > end:
         coverage += (j - end)
         end = j
      else: continue
   return coverage

# Read chromosome lengths
chromosome_lengths = {}
with open(sys.argv[3]) as lengths_file:
    for line in lengths_file:
        parts = line.strip().split()
        chromosome_lengths[parts[0]] = int(parts[2])

gene_coordinates = extract_gene_coordinates(sys.argv[1])
if sys.argv[4] == "gene":
   count_insertions_in_genes(sys.argv[2], gene_coordinates)
elif sys.argv[4] == "promoter": 
   count_insertions_in_promoters(sys.argv[2], gene_coordinates)
elif sys.argv[4] == "genecoverage":
   count_insertions_perGeneCoverage(sys.argv[2], gene_coordinates)
