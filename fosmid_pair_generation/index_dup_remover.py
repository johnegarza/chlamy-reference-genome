#this script scans through the fosmid pairs and inserts an index (corresponding to a line in the block file)
#indicating which block a given fosmid end matches

import sys, os, csv, re, fnmatch

if len(sys.argv) < 2:
	sys.exit("Usage: %s indexed_pairs" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])

fosmid_pairs = sys.argv[1]


with open(fosmid_pairs) as pairs:
	pair_data = csv.reader(pairs, delimiter="\t")

	for pair in pair_data:

		line1 = int(pair[4])
		line2 = int(pair[9])

		if (line1 != line2):
			
			print(pair[0] + "\t" + str(pair[1]) + "\t" + str(pair[2]) + "\t" + pair[3] + "\t" + str(pair[4]) + "\t" + pair[5] + "\t" + str(pair[6]) + "\t" + str(pair[7]) + "\t" + pair[8] + "\t" + str(pair[9]))



