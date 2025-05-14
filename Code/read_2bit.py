
import io
import struct
import sys
import py2bit
    
def count_pattern_matches(sequence, pattern):

    count = 0
    start_index = 0
    while True:
        index = sequence.find(pattern, start_index)
        if index == -1:
            break
        count += 1
        start_index = index + 1  
    return count

# function get all indicies of match but I didn't use use!!
def find_all_pattern_indices(sequence, pattern):
    indices = []
    start_index = 0
    while True:
        index = sequence.find(pattern, start_index)
        if index == -1:
            break
        indices.append(index)
        start_index = index + 1 # Move past the found pattern
    return indices

# ham handed way of creating all 4096 sixmers
def generate_6mers():
	bases = ['A','G','C','T']
	sixmers = []
	for b1 in bases:
		for b2 in bases:
			for b3 in bases:
				for b4 in bases:
					for b5 in bases:
						for b6 in bases:
							sixmers.append(b1 + b2 + b3 + b4 + b5 + b6)
	return sixmers

# function gets dna sequence by file + chromosome name
def get_chr_seq(file,chrom):
	with py2bit.open(file) as tb:
		sequence = tb.sequence(chrom)
		tb.close()
		return sequence

# function gets list of chromosomes from target file
def get_chrs(file):
	with py2bit.open(file) as tb:
		chr_dict = tb.chroms()
		chr_list = list(chr_dict.keys())
		tb.close()
		return chr_list


if __name__ == "__main__":
    file_path = "hs1.2bit"  # path to .2bit file
    # downloaded from https://hgdownload.soe.ucsc.edu/gbdb/hs1/
	
    # create list of all potential 6mers
    sixmers = generate_6mers()

    # get chromosome list from target .2bit file (also helps confirm all data is there)
    allchrs = get_chrs(file_path)
    
    #for each possible sixmer, loop through each chr and count matches
    for x in sixmers:
        sixmercount = 0
        for c in allchrs:
            genomic_sequence = get_chr_seq(file_path,c)
            if genomic_sequence:
                match_count = count_pattern_matches(genomic_sequence, x)
                sixmercount += match_count
        #after iterating all chrs, write out the total for the sixmer
        f = io.open("6merCounts.txt","a")
        f.write(x + "," + str(sixmercount)+ "\n")
        f.close()
