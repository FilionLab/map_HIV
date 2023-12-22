import matplotlib.pyplot as plt

from collections import defaultdict
import sys

# Read chromosome lengths
chromosome_lengths = {}
with open(sys.argv[2]) as lengths_file:
    for line in lengths_file:
        parts = line.strip().split()
        chromosome_lengths[parts[0]] = int(parts[2])

# Read mapped sites
mapped_sites_count = defaultdict(int)
with open(sys.argv[1]) as sites_file:
    for line in sites_file:
        chromosome = line.strip().split()[0]
        if "random" in chromosome or "alt" in chromosome or "Un" in chromosome or "chrM" in chromosome: continue 
        mapped_sites_count[chromosome] += 1

chromosomes = sorted([item for item in list(mapped_sites_count.keys())], key=lambda x: (float('inf') if x in {"chrX", "chrY", "chrM"} else int(x[3:]), x))
ratios = [mapped_sites_count[item]/chromosome_lengths[item] for item in chromosomes]
# Create a bar plot
plt.bar(chromosomes, ratios)
plt.xlabel('Chromosome')
plt.ylabel('Insertions per Length')
plt.title('Insertions per Chromosome Length')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

