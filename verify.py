# SAMVerify in Python3
# @author: Stephen Petrides

import sys
import json
import math

## Print Reference ##

def find_in_ref(chromosome_str, position_str, seq_len_str, ref_str, ref_str_idx='.fai', rev_comp=False):
	""""""
	chr_len = len(chromosome_str)
	pos_int = int(position_str)
	seq_len = int(seq_len_str)
	if ref_str_idx == '.fai':
		ref_str_idx = ref_str + ref_str_idx
	elif ref_str_idx == False:
		rev_comp == False
		ref_str_idx = ref_str + '.fai'

	if ref_str.endswith('.fa') or ref_str.endswith('.fasta'):
		genome_file = open(ref_str, 'rt')

		with open(ref_str_idx, 'rt') as index_file:
			index_lines = index_file.readlines()
		for line in index_lines:
			line_list = line.split()
			if chromosome_str == line_list[0]:
				offset = int(line_list[2])
				whole_chr_len = int(line_list[1])
				updated_pos = pos_int + math.floor(pos_int/int(line_list[3]))

		genome_file.seek(offset+updated_pos-1, 0)
		ref_seq = genome_file.read(seq_len+2).replace('\n', '').upper()
		if len(ref_seq) != seq_len:
			ref_seq = ref_seq[:-1]
		genome_file.close()
	else:
		chromosome_str = ref_str + '/' + chromosome_str + '.fa'
		chromosome_file = open(chromosome_str, 'rt')
		chromosome_file_str = chromosome_file.read().replace('\n', '').upper()
		ref_seq = chromosome_file_str[chr_len+pos_int:chr_len+pos_int+seq_len]
		chromosome_file.close()

	if rev_comp is True:
		ref_seq = reverse_complement(ref_seq)
		return ref_seq
	else:
		return ref_seq

## Reconstruct Reference ##

def reconstruct(line):
	""""""
	sam_line = line.split()

	flag_str = str(bin(int(sam_line[1])))
	if flag_str[-3] == '1':
		return 'Read unmapped'

	reconstruct_str = sam_line[9]
	cigar_str = sam_line[5]
	for field in sam_line:
		if field[0:2] == 'MD':
			md_str = field[5:]

	num = set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
	SEQ_pos = -1
	num_str = ''

	for count in range(0, len(cigar_str)):
		if cigar_str[count] in num:
			num_str+=cigar_str[count]
		else:
			num_int = int(num_str)
			num_str = ''
			SEQ_pos += num_int
			if cigar_str[count] == 'I':
				reconstruct_str = reconstruct_str[0:SEQ_pos] + reconstruct_str[SEQ_pos+num_int:len(reconstruct_str)]

	SEQ_pos = -1
	num_str = ''

	for count in range(0, len(md_str)):
		if md_str[count] in num:
			num_str+=md_str[count]
		else:
			if num_str == '':
				break;
			num_int = int(num_str)+1
			num_str = ''
			SEQ_pos += num_int
			if md_str[count] != '^':
				reconstruct_str = reconstruct_str[0:SEQ_pos] + md_str[count].lower() + reconstruct_str[SEQ_pos+1:]
				#print(md_str[count] + '\t' + str(SEQ_pos) + '\t'  + reconstruct_str)
			else:
				count+=1
				reconstruct_str = reconstruct_str[0:SEQ_pos] + md_str[count].lower() + reconstruct_str[SEQ_pos:]
	return reconstruct_str

## Reverse Complement Sequence ##
def reverse_complement(seq):
	""""""
	seq = seq.upper()
	qes = ''
	for base in seq:
		if base is 'A':
			qes+='T'
		elif base is 'T':
			qes+='A'
		elif base is 'C':
			qes+='G'
		elif base is 'G':
			qes+='C'
		else:
			return 'Error: Incorrect bases included in sequence. Use only ATCG.'
	qes = qes[::-1]
	return qes


#### MAIN ####
def main():
	## Set Parameters ##

	argc = len(sys.argv)
	param_set = set()
	param_count = 1

	manual = ["\nCommands and Options", "\n\n\tBasic Command Structure:", "\n\n\t\t$ python3 verify.py <flag> <path>",
	"\n\t\t\t# Flags can come in any order but must be followed by corresponding information",
	"\n\tFlags", "\n\t\t-man", "\n\t\t\t# Prints Commands and Options section of MANUAL.", "\n\t\t-sam (required)",
	"\n\t\t\t# Followed by the path/to/file_name.sam", "\n\t\t-ref (required)",
	"\n\t\t\t# Either a path to a whole genome.fa or a path to directory containing chr_.fa files", "\n\t\t-ref_idx (optional)"
	"\n\t\t\t# Use if whole genome is used and .fai is in different directory or ",
	"\n\t\t\t with difference file name than [whole_genome_file.fa].fai", "\n\t\t-header (optional)",
	"\n\t\t\t# Followed by the header of target read", "\n\t\t\t# If no header is given the first read is verified"]

	while param_count <= argc-1:
		param_str = sys.argv[param_count]
		param_set.add(param_str)
		param_count +=1
		if param_str[0] is not '-':
			print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: " + param_str)
			sys.exit(0)
		else:
			if param_str == '-man':
				print(''.join(manual))
				sys.exit(0)
			elif param_str == '-sam':
				param_str = sys.argv[param_count]
				param_count +=1
				sam_str = param_str
			elif param_str == '-ref':
				param_str = sys.argv[param_count]
				param_count +=1
				ref_str = param_str
			elif param_str == '-ref_idx':
				param_str = sys.argv[param_count]
				param_count +=1
				ref_idx = param_str
			elif param_str == '-header':
				param_str = sys.argv[param_count]
				param_count +=1
				header_str = param_str

   ## Print Sequence From SAM ##
	sam_file = open(sam_str, 'rt')

	if '-header' in param_set:
		print('\nHeader:\t' + header_str)
		for line in sam_file:
			sam_line = line.split()
			if sam_line[0] == header_str:
				flag_str = str(bin(int(sam_line[1])))
				if flag_str[-3] == '1':
					print('Read is unmapped.')
					continue;
				position = sam_line[3]
				raw_seq = sam_line[9]
				print('Position:\t' + position + '\t' + sam_line[2] + '\t' + sam_line[11])
				print('Read:\t\t' + raw_seq)

				if '-ref' in param_set and '-ref_idx' in param_set:
					print('Reference:\t' + find_in_ref(sam_line[2], sam_line[3], len(sam_line[9]), ref_str, ref_idx))
					print('Reconstruction:\t' + reconstruct(line))
				elif '-ref' in param_set:
					print('Reference:\t' + find_in_ref(sam_line[2], sam_line[3], len(sam_line[9]), ref_str))
					print('Reconstruction:\t' + reconstruct(line))	
				else: 
					print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: No reference genome specified.")

	else:
		for line in sam_file:			# without a header specified, it takes the first sam record
			sam_line = line.split()
			first_word = sam_line[0]
			if first_word[0] != '@':
				flag_str = str(bin(int(sam_line[1])))
				if flag_str[-3] == '1':
					continue;
				print('Header:\t\t' + first_word)
				position = sam_line[3]
				raw_seq = sam_line[9]
				print('Position:\t' + position + '\t' + sam_line[2])
				print('Read:\t\t' + raw_seq)
				if '-ref' in param_set and '-ref_idx' in param_set:
					print('Reference:\t' + find_in_ref(sam_line[2], sam_line[3], len(sam_line[9]), ref_str, ref_idx))
					print('Reconstruction:\t' + reconstruct(line))
					break;
				elif '-ref' in param_set:
					print('Reference:\t' + find_in_ref(sam_line[2], sam_line[3], len(sam_line[9]), ref_str))
					print('Reconstruction:\t' + reconstruct(line))
					break;
				else:
					print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: No reference genome specified.")	
	sam_file.close()

if __name__ == "__main__":
    main()








