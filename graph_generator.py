from collections import defaultdict
import sys, os, csv, re, fnmatch
from block_node import Node
from fosmid_edge import Edge

if len(sys.argv) < 3:
	sys.exit("Usage: %s all_assembly_alignments.tsv bac_pairs.tsv" % sys.argv[0])
if not os.path.exists(sys.argv[1]):
	sys.exit("Error: File '%s' not found" % sys.argv[1])
if not os.path.exists(sys.argv[2]):
	sys.exit("Error: File '%s' not found" % sys.argv[2])


alignment_file = sys.argv[1] #tab_delim_results.tsv
pair_file = sys.argv[2] #full_bac_matches.txt

ref_to_asm = defaultdict(list)

dup_filter = set()

with open(alignment_file) as a_f:
	alignment_data = csv.reader(a_f, delimiter="\t")

	for place, alignment in enumerate(alignment_data):

		ref_chr = str(alignment[0])
		ref_start = int(alignment[1]) #1 for all data, 2 for sanitized
		ref_stop = int(alignment[2])  #2 for all data, 3 for sanitized
		asm_scaf = str(alignment[3])  #3 for all data, 4 for santized
		asm_start = int(alignment[4])
		asm_stop = int(alignment[5])
		tracker = place + 1

		#ref_to_asm[ref_chr].append((ref_start, ref_stop, tracker, asm_scaf))
		Node(tracker, ref_chr, ref_start, ref_stop, asm_scaf, asm_start, asm_stop)

with open(pair_file) as p_f:
	pair_data = csv.reader(p_f, delimiter="\t")

	for pair in pair_data:

		left_chr = str(pair[0])
		left_start = int(pair[1])
		left_stop = int(pair[2])
		right_chr = str(pair[4])
		right_start = int(pair[5])
		right_stop = int(pair[6])

		left_tuples = ref_to_asm[left_chr]

		for l_t in left_tuples:
			
			if ( (left_start >= l_t[0]) and (left_stop <= l_t[1]) ):
				#left_side = True
				left_line = l_t[2]
				left_scaf = l_t[3]

				right_tuples = ref_to_asm[right_chr]

				for r_t in right_tuples:
				
					if ( (right_start >= r_t[0]) and (right_stop <= r_t[1]) ):

						right_line = r_t[2]
						right_scaf = r_t[3]

						if (left_scaf == right_scaf):

							if ( abs(left_line - right_line) < 4 ):

								match_string = "Line " + str(left_line) + " matches line " + str(right_line) + "\n"

								if(match_string not in dup_filter):

									dup_filter.add(match_string)
									print("Using data " + str(pair))
									print(match_string)






