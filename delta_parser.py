from __future__ import print_function #support python3 style print function, allowing me to change default terminal \n behavior
import sys, os, csv, re, argparse
from mapper_helper import RefMapper
from collections import defaultdict

alignments = defaultdict(RefMapper)

parser = argparse.ArgumentParser(description='parses a delta file into a format that can be used to map ref coordinates to query coordinates')
parser.add_argument('delta_file', metavar='*.delta', help='a delta file produced by nucmer')
#TODO design/document a format for this fosmid file
parser.add_argument('fosmids', metavar='fosmids.tsv', help = 'a fosmid paired end file')
parser.add_argument('alignment_file', metavar='tab_delim_results', help='an alignment file from nucmer')
args = parser.parse_args()

delta_file = args.delta_file
fosmid_file = args.fosmids
block_file = args.alignment_file
delta_regex = re.compile(r"^-?\d+$")

with open(delta_file) as d_f:
	lines = d_f.readlines()

	del lines[:2] #delete the first 2 lines in the lines list, which are just generic info
	#saves me writing a case in my if-elif-else block just to handle the first 2 lines of the file

	#TODO some of these can be removed
	ref_name = "you shouldn't be able to read this"
	query_name = "you shouldn't be able to read this either"
	ref_array = []
	query_array = []
	ref_start = 0 #do I really need all these initialization lines?
	ref_end = 0   #probably not but I have bad habits from other languages
	q_start = 0   #I'm used to Java scopes
	q_end = 0
	q_rev = False

	for line in lines:

		#line.strip()

		if (line[0] == ">"): #the beginning of a seq, in the format >ref_name query_name ref_length query_length
			name_info = line.split()
			ref_name = name_info[0].lstrip(">")
			query_name = name_info[1]
#			print("seq start- names")

		elif (line[0].isdigit() and line[0] == "0"): #marks end of an alignment block
			

			#to avoid writing multiple delta cases (see the else at the end of this if/elif/else block)
			#when a query is reversed, all calculations are performed as if it were forward;
			#here I map it back to reverse before dumping the data
			if(q_rev):

				for index, temp_q in enumerate(query_array):

					mapped_tuple = ( q_start-temp_q[0] , q_start-temp_q[1] )
					query_array[index] = mapped_tuple

			alignments[ref_name].update(ref_array, query_array, query_name)

		elif( delta_regex.search(line) is None): #beginning of alignment block, in the format
			
			a_info = line.split(" ")         #ref_start ref_end query_start query_end errors errors errors

			ref_start = int(a_info[0])
			ref_end = int(a_info[1])
			q_start = int(a_info[2])
			q_end = int(a_info[3])

			#for later use- check to see if the query alignment is reversed or not
			if(q_start > q_end):
				q_rev = True
			else:
				q_rev = False

			ref_array = []    #reset arrays; old info has been printed in the above elif
			query_array = []

			ref_array.append((ref_start, ref_end)) #initialize arrays with new info

			if(q_rev):
				q_end = q_start - q_end
				query_array.append((0, q_end))
			else:
				query_array.append((q_start, q_end))
#			print("seq start- alignment info >" + line.rstrip() + "<")

			ref_tracker = ref_start

			if(q_rev):
				query_tracker = 0
			else:
				query_tracker = q_start

			curr_list_index = 0

		else: #this line contains delta info
#			print("delta line")
			delta = int(line)
			abs_delta = abs(delta)

			ref_tracker = int(ref_array[curr_list_index][0])
			query_tracker = int(query_array[curr_list_index][0])

			#old_ref_tuple = ref_array[curr_list_index]
			#old_q_tuple = query_array[curr_list_index]
			del ref_array[curr_list_index]
			del query_array[curr_list_index]
			curr_list_index += 1

			if (delta == -1): #move ref up by 1, no other changes (edge case)

				ref_array.append( (ref_tracker + 1, ref_end) )
				query_array.append( (query_tracker, q_end) )
				curr_list_index -= 1
			
			elif (delta == 1): #move query up by 1, no other changes (edge case)
				
				ref_array.append( (ref_tracker, ref_end) )
				query_array.append( (query_tracker + 1, q_end) )
				curr_list_index -= 1

			elif (delta < 0): #makes a gap in the ref

				ref_array.append( (ref_tracker, ref_tracker + (abs_delta - 2)) )
				ref_array.append( (ref_tracker + abs_delta, ref_end) )
				query_array.append( (query_tracker, query_tracker + (abs_delta - 2)) )
				query_array.append( (query_tracker + (abs_delta - 1), q_end) )

			else: #makes a gap in the query

				ref_array.append( (ref_tracker, ref_tracker + (abs_delta - 2)) )
				ref_array.append( (ref_tracker + (abs_delta - 1), ref_end) )
				query_array.append( (query_tracker, query_tracker + (abs_delta - 2)) )
				query_array.append( (query_tracker + abs_delta, q_end) )

'''
for key in alignments:
	print("---------------" + str(key) + "---------------")
	print(str(alignments[key]))
'''

#for key in alignments:
#	str(alignments[key])
#os.system("pause")

real_mapper = defaultdict(RefMapper)
for key in alignments:
	#for rm in alignments[key]:
	for index, val in enumerate(alignments[key].query_names):
		new_key = val
		new_ref_coords = alignments[key].query_coords[index]
		new_query_coords = alignments[key].ref_coords[index]
		new_query_name = key
		real_mapper[new_key].update(new_ref_coords, new_query_coords, new_query_name)


'''
print("scaf 32 aligns to")
str(alignments["scaffold32_size125501_pilon"])
print("chromosome 1 aligns to")
str(real_mapper["chromosome_1"])
'''	

with open(fosmid_file) as f_f:
	lines = f_f.readlines()
	for index, line in enumerate(lines):
		line = line.split("\t")

		left_ref = str(line[0])
		left_end_start = int(line[1])
		left_end_stop = int(line[2])
		left_line = int(line[4])

		right_ref = str(line[5])
		right_end_start = int(line[6])
		right_end_stop = int(line[7])
		right_line = int(line[9])

		with open(block_file) as blocks:
			block_data = blocks.readlines()
			left_name = block_data[left_line - 1].split("\t")[3]
			right_name = block_data[right_line - 1].split("\t")[3]
			

		'''
		left_map_name = real_mapper[left_name].map(left_end_start, left_end_stop, index, 0)[3]
		right_map_name = real_mapper[right_name].map(right_end_start, right_end_stop, index, 1)[3]
		
		left_tuple = alignments[left_map_name].map(left_end_start, left_end_stop, index, 0)
		right_tuple = alignments[right_map_name].map(right_end_start, right_end_stop, index, 1)
		'''

		left_tuple = real_mapper[left_name].map(left_end_start, left_end_stop, index, 0)
		right_tuple = real_mapper[right_name].map(right_end_start, right_end_stop, index, 1)

		if( left_tuple[0] and right_tuple[0] ):

			
			#with the block at the end removed, the following block does nothing
			'''
			#TODO min/max might only be needed on lines 197 and 198- double check and remove
			#the 2 min and 2 max calls directly below if this is the case
			end1_scaf = str(left_tuple[3])
			end1_start = min( int(left_tuple[4]), int(left_tuple[5]) )
			end1_stop = max( int(left_tuple[4]), int(left_tuple[5]) )
			end2_scaf = str(right_tuple[3])
			end2_start = min( int(right_tuple[4]), int(right_tuple[5]) )
			end2_stop = max( int(right_tuple[4]), int(right_tuple[5]) )
			'''

			ans = []
			ans.append(left_name)
			ans.extend( [str(x) for x in left_tuple[1:]] )
			ans.append( str(line[3]) ) #fosmid ID, VTP*
			#ans.append( str(end1_line) ) 
			ans.append(str(left_line))
			ans.append( right_name )
			ans.extend( [str(x) for x in right_tuple[1:]] )
			ans.append( str(line[8]) ) #fosmid ID, VTP*
			#ans.append( str(end2_line) )
			ans.append(str(right_line))
			print("\t".join(ans).strip())



			#removing following section since node_list_indexed has the line data already
			#however, the following code may produce slightly more accurate results- return to this when time TODO
			'''
			with open(block_file) as blocks:
				block_data = csv.reader(blocks, delimiter="\t")
			
				end1_line = 0
				end2_line = 0

				for block in block_data:

					chrom = str(block[3])
					start = min( int(block[4]), int(block[5]) )
					stop = max( int(block[4]), int(block[5]) )
					file_line = int(block[8])

					if ( (end1_scaf == chrom) and (end1_start >= start) and (end1_stop <= stop) ):
						end1_line = file_line

					if ( (end2_scaf == chrom) and (end2_start >= start) and (end2_stop <= stop) ):
						end2_line = file_line


			if ( (end1_line != 0) and (end2_line != 0) ):
			###

				ans = []
				ans.append(left_name)
				ans.extend( [str(x) for x in left_tuple[1:]] )
				ans.append( str(line[3]) ) #fosmid ID, VTP*
				ans.append( str(end1_line) ) 
				ans.append( right_name )
				ans.extend( [str(x) for x in right_tuple[1:]] )
				ans.append( str(line[7].strip()) ) #fosmid ID, VTP*
				ans.append( str(end2_line) )
				print("\t".join(ans).strip())
			'''


