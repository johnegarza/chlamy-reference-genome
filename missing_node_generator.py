name_to_len = dict()
with open("assembly.fasta.fai") as fai:
	lines = fai.readlines()
	parsed_lines = [ x.split("\t") for x in lines ]
	for line in parsed_lines:
		name_to_len[line[0]] = int(line[1])

final_data = []

last_node = ("scaf", 0, 0, "chromosome", 0, 0, 0, 0, 0)

with open("node_list.tsv") as node_data:
	nodes = node_data.readlines()
	for node in nodes:
		split = node.split("\t")

		#detects the switch to a new assembly scaffold
		if split[0] != last_node[0]:
			stop_coord = name_to_len[last_node[0]]
			curr_stop = int(last_node[2]) + 1
			if curr_stop < stop_coord:
				temp = []
				temp.append(last_node[0])
				temp.append(str(curr_stop))
				
