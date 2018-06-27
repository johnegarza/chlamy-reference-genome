name_to_len = dict()
with open("assembly.fasta.fai") as fai:
	lines = fai.readlines()
	parsed_lines = [ x.split("\t") for x in lines ]
	for line in parsed_lines:
		name_to_len[line[0]] = int(line[1])

#for edge case- reading first line of the file
name_to_len["scaf"] = 0

final_data = []
included_scafs = set()

#assembly scaffold, start, stop, reference chromosome, start, stop, assembly length, reference length, line number
last_node = ("scaf", 0, 0, "chromosome", 0, 0, 0, 0, 0)

with open("../node_list.tsv") as node_data:
	nodes = node_data.readlines()
	for node in nodes:
		split = node.split("\t")

		#detects the switch to a new assembly scaffold
		if split[0] != last_node[0]:

			included_scafs.add(split[0])

			stop_coord = name_to_len[last_node[0]]
			curr_stop = int(last_node[2]) + 1
			if curr_stop < stop_coord:

				node_len = (stop_coord - curr_stop) + 1
				temp = [last_node[0], str(curr_stop), str(stop_coord), "no_match", "0", "0", node_len, "0", "-2"]

				new_node = "\t".join(temp)
				final_data.append(new_node)

			curr_start = int(split[1]) - 1
			if curr_start > 1:
				temp = [split[0], "1", str(curr_start), "no_match", "0", "0", str(curr_start), "0", "-2"]
				new_node = "\t".join(temp)
				final_data.append(new_node)

		last_stop = int(last_node[2]) + 1
		curr_start = int(split[1]) - 1
		elif (curr_start - last_stop > 0):
			node_len = (curr_start - last_stop) + 1
			temp = [split[0], str(last_stop), str(curr_start), "no_match", "0", "0", str(node_len), "0", "-2"]
			new_node = "\t".join(temp)
			final_data.append(new_node)

		final_data.append(node)

del name_to_len["scaf"] #remove artificially added edge case entry
for scaffold in name_to_len:
	if scaffold not in included_scafs:
		stop = int(name_to_len[scaffold])
		temp = [ scaffold, "1", str(stop), "no_match", "0", "0", str(stop), "0", "-2" ]
		new_node = "\t".join(temp)
		final_data.append(new_node)

for line in final_data:
	print(line)

