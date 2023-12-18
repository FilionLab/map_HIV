import sys

def process_file(input_file, output_file):
    data = [{"chromosome": "test", "position": 0}]

    with open(input_file) as f, open(output_file, 'w') as outfile:
        next(f)
        for line in f:
            parts = line.split()
            chromosome = parts[0]
            position = int(parts[1])
            sequence = parts[2]
            
            # Check if the chromosome is already in the dictionary
            if chromosome not in [i["chromosome"] for i in data] and (chromosome != "phiX" and chromosome != "HIVmini"):
                  # Check if the position is unique or differs by +/- 1 or 2 or 3
                  if all(abs(pos - position) > 3 for pos in [i["position"] for i in data]):
                           data.append({"chromosome": chromosome, "position": int(position)})
                           outfile.write(f"{chromosome} {position} {sequence}\n")
            
            elif chromosome != "phiX" and chromosome != "HIVmini":
                           # If chromosome not in the dictionary, add it with the current position
                           outfile.write(f"{chromosome} {position} {sequence}\n")
                           data.append({"chromosome": chromosome, "position": int(position)})

input_file = sys.argv[1]
output_file = sys.argv[2]

process_file(input_file, output_file)

