# SAMVerify in Python3
# @author: Stephen Petrides

import sys
import json

## Print Reference ##

def find_in_ref(sam_line, ref_dir_str):
	chromosome = sam_line[2]
	chromosome_str = ref_dir_str + '/' + chromosome + '.fa'
	chromosome_file = open(chromosome_str, 'rt')
	position = int(sam_line[3])
	seq_len = len(raw_seq)

	chromosome_file_str = chromosome_file.read().replace('\n', '').upper()
	ref_seq = chromosome_file_str[len(chromosome)+position:len(chromosome)+position+seq_len]

	print("Reference:\t" + ref_seq)

	else:
		print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: No reference genome specified.")

## Reconstruct Reference ##

def reconstruct(sam_line):
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
	print('Reconstruction:\t' + reconstruct_str)




#### MAIN ####
def main():
	## Set Parameters ##

	argc = len(sys.argv)
	param_set = set()
	param_count = 1

	manual = ["\nCommands and Options", "\n\n\tBasic Command Structure:", "\n\n\t\t$ python3 verify.py <flag> <path>",
	"\n\t\t\t# Flags can come in any order but must be followed by corresponding information",
	"\n\tFlags", "\n\t\t-man", "\n\t\t\t# Prints Commands and Options section of MANUAL.", "\n\t\t-sam",
	"\n\t\t\t# Followed by the path/to/file_name.sam", "\n\t\t-ref_dir",
	"\n\t\t\t# Followed by the path to the directory containing chr_.fa", "\n\t\t-header (optional)",
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
			elif param_str == '-ref_dir':
				param_str = sys.argv[param_count]
				param_count +=1
				ref_dir_str = param_str
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
				position = sam_line[3]
				raw_seq = sam_line[9]
				print('Position:\t' + position + '\t' + sam_line[2] + '\t' + sam_line[11])
				print('Read:\t\t' + raw_seq)

				if '-ref_dir' in param_set:
					find_in_ref(sam_line, ref_dir_str)			# call function to print from reference
					reconstruct(sam_line)						# call function to reconstruct from CIGAR and TAGs
				else: 
					print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: No reference genome specified.")

	else:
		for line in sam_file:			# without a header specified, it takes the first sam record
			sam_line = line.split()
			first_word = sam_line[0]
			if first_word[0] != '@':
				print('Header:\t\t' + first_word)
				position = sam_line[3]
				raw_seq = sam_line[9]
				print('Position:\t' + position + '\t' + sam_line[2])
				print('Read:\t\t' + raw_seq)
				if '-ref_dir' in param_set:
					find_in_ref(sam_line, ref_dir_str)
					reconstruct(sam_line)
					break;
				else:
					print("Parameters not entered correctly. See MANUAL for proper syntax.\nError: No reference genome specified.")	

if __name__ == "__main__":
    main()








