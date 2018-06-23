import pysam

samfile = pysam.AlignmentFile("../novoalign/imp3.merged.sorted.bam", "rb")

for position in samfile.pileup("scaffold1_size20154211_pilon", 0, 0, stepper="nofilter"):

	print( "base " + str(position.pos) + " with depth " + str(position.n) )

samfile.close()
