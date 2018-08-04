import os
import random
import datetime
import csv



def mismatches(SEQ1, SEQ2, MAX = float("inf"), IGNORE_N = 0 ):
	"""Returns the number of mismatches between two strings.
	MAX sets the number of max number of mismatches that are reported.
	Lowering MAX increases performance.
	IGNORE_N = 1 will ignore mismatches with N."""
	mismatches = 0
	if IGNORE_N == 0:
		for i in range(len(SEQ1)):
			if SEQ1[i] != SEQ2[i]:
				mismatches += 1
			if mismatches >= MAX:
				break
		return mismatches
	else:
		for i in range(len(SEQ1)):
			if SEQ1[i] != 'N' and SEQ2[i] != 'N':
				if SEQ1[i] != SEQ2[i]:
					mismatches += 1
			if mismatches >= MAX:
				break
	return mismatches


def best_distance(SEQ1, LIST, MAX = float("inf"), IGNORE_N = 0, PRINT = 0 ):
	"""finds the best match for a sequence in a list of sequences.
	MAX sets the number of max number of mismatches before it moves on.
	Lowering MAX increases performance.
	IGNORE_N = 1 will ignore mismatches with N."""
	y = MAX
	for i in range(len(LIST)):
	   z = mismatches(SEQ1, LIST[i], y, IGNORE_N)
	   if z < y:
		  y = z
	   if z == 0:
		  break
	return y

def randDNA(n):
	x = ''.join(random.choice(['A','C','T','G']) for _ in range(n))
	return x




os.chdir("/Users/sasha/OneDrive - Leland Stanford Junior University/Github/Guide-to-bar-seq/")

#Make a matrix of nearest neighbors for various bc.lengths and library sizes
bc_number = [10,1000000]
bc_length = [10,15,20,25,30]

for k in bc_number:
	#start a file to write to
	x = open('nearest_neighbor' + str(k) + 'length' + str(bc_length[0]) + '_' + str(bc_length[-1]) + '.csv', 'w')
	for i in bc_length:	
		print (str(datetime.datetime.now()) + " Making random barcodes for " + str(k) + " barcodes of length " + str(i))
		bc = []
		for j in range(k):
			bc.append(''.join(random.choice(['A','C','T','G']) for _ in range(i)))	
		nn = []
		print ("Finding nearest neighbors")
		for j in range(k):
			nn.append(best_distance(bc[j],bc[:j] + bc[j+1:]))
		x.write(','.join(str(e) for e in nn))
		x.write("\n")
	x.close()




		

	

