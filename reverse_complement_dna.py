import sys

def reverse_complement_dna(sequence):
    complement_dict = {'N':'N','A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement_sequence = ''.join(complement_dict[base] for base in reversed(sequence))
    return complement_sequence

file_path = sys.argv[1]

sequence_data = []
with open(file_path) as f:
    header = next(f)
    for line in f:
        line = line.strip().split()
        entry = {
            "ref": line[0],
            "pos": line[1],
            "seq": line[2]
        }
        sequence_data.append(entry)

#Reverse complement all sequences in the dictionary
for entry in sequence_data:
    entry["seq"] = reverse_complement_dna(entry["seq"])

print(header)
for entry in sequence_data:
    print(entry["ref"], entry["pos"], entry["seq"])
