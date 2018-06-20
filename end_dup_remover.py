import sys, os

if len(sys.argv) < 2:
	sys.exit("Usage: %s ends_from_delta_parser" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])

edge_file = sys.argv[1]

with open(edge_file) as e_f:
	lines = e_f.readlines()

lines = [ tuple(line.strip().split("\t")) for line in lines ]

unique_ends = set()

for line in lines:

	end1 = line[0:8]
	end2 = line[8:16]
	
	end_struct = frozenset( (end1, end2) )
	
	if end_struct not in unique_ends:

		unique_ends.add(end_struct)
		print("\t".join(line))
