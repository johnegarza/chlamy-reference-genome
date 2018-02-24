import sys, os, csv, re, fnmatch

if len(sys.argv) < 2:
	sys.exit("Usage: %s BAC_fosmid_pairs.tsv" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])

BAC_pairs = sys.argv[1]

pair_ID = re.compile(r"([VTPQ]{3})")

with open(BAC_pairs) as pairs:
	pair_data = csv.reader(pairs, delimiter="\t")

	for pair in pair_data:

		pair_type = pair_ID.match(pair[3]).group(1) #extract the VTP/PTQ part from the beginning of identifier

		if (pair_type == "VTP"): #VTP denotes a fosmid
	
			print(pair[0] + "\t" + str(pair[1]) + "\t" + str(pair[2]) + "\t" + pair[3] + "\t" + str(pair[4]) + "\t" + str(pair[5]) + "\t" + pair[6] + "\t" + pair[7])
