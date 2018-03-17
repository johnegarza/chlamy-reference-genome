with open("out.1delta") as file1:
	for line in file1.readlines():
		num_list = line.split()
		if(len(num_list) == 7 ):
			if(not( int(num_list[0]) < int(num_list[1]) and int(num_list[2]) < int(num_list[3]) ) ):
				
				if( int(num_list[0]) > int(num_list[1]) and int(num_list[2]) > int(num_list[3])  ):
					print("both bigger")
				elif( int(num_list[0]) > int(num_list[1]) ):
					print("first bigger")
				elif( int(num_list[2]) > int(num_list[3]) ):
					print("second bigger")
