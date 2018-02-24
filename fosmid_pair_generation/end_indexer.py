#this script scans through the fosmid pairs and inserts an index (corresponding to a line in the block file)
#indicating which block a given fosmid end matches

import sys, os, csv, re, fnmatch

if len(sys.argv) < 3:
	sys.exit("Usage: %s fosmid_pairs tab_delim_results" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])
if not os.path.exists(sys.argv[2]):
	sys.exit("Error: File '%s' not found" % sys.argv[2])


fosmid_pairs = sys.argv[1]
block_file = sys.argv[2]

#pair_ID = re.compile(r"([VTPQ]{3})")

with open(fosmid_pairs) as pairs:
	pair_data = csv.reader(pairs, delimiter="\t")

	for pair in pair_data:

		end1_chr = pair[0]
		end1_start = int(pair[1])
		end1_stop = int(pair[2])
		end2_chr = pair[4]
		end2_start = int(pair[5])
		end2_stop = int(pair[6])

		with open(block_file) as blocks:
			block_data = csv.reader(blocks, delimiter="\t")
			
			end1_line = 0
			end2_line = 0

			for block in block_data:

				chrom = block[0]
				start = int(block[1])
				stop = int(block[2])
				line = int(block[8])

				if ( (end1_chr == chrom) and (end1_start >= start) and (end1_stop <= stop) ):
					end1_line = line

				if ( (end2_chr == chrom) and (end2_start >= start) and (end2_stop <= stop) ):
					end2_line = line


		if ( (end1_line != 0) and (end2_line != 0) ):
			
			print(pair[0] + "\t" + str(pair[1]) + "\t" + str(pair[2]) + "\t" + pair[3] + "\t" + str(end1_line) + "\t" + str(pair[4]) + "\t" + str(pair[5]) + "\t" + pair[6] + "\t" + pair[7] + "\t" + str(end2_line))



