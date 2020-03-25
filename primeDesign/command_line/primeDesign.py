# Design pegRNAs and ngRNAs for prime editing
##### Import libraries
import os
import sys
import re
import time
import argparse
import logging
from argparse import RawTextHelpFormatter

##### Argument handeling
parser = argparse.ArgumentParser(description = '''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Software for the design of pegRNAs for flexible prime editing! Please visit ----- https://github.com/jyhsu15/pegDesigner ----- for more documentation on how to use the software.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------''', formatter_class=RawTextHelpFormatter)

# Inputs for de-novo design of pegRNAs and nicking gRNAs
parser.add_argument('-f', '--file', required = True, type = str, help = '''Input file (.txt or .csv) for pegDesigner

*** Example .TXT file *** --------------------------------------------------------------
|											|
|	target_01_substitution	ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC		|
|	target_01_insertion	ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC	|
|	target_01_deletion	ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC		|
|											|
 ---------------------------------------------------------------------------------------

*** Example .CSV file *** --------------------------------------------------------------
|											|
|	target_01_substitution,ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC		|
|	target_01_insertion,ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC		|
|	target_01_deletion,ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC			|
|											|
 ---------------------------------------------------------------------------------------

*** Formatting different DNA edits *** -------------------------------------------------
|											|
|	Substitution edit:	Format: (reference/edit)	Example:(G/A)		|
|	Insertion edit:		Format: (+insertion)		Example:(+ATCG)		|
|	Deletion edit:		Format: (-deletion)		Example:(-ATCG)		|
|											|
 ---------------------------------------------------------------------------------------

*** Combination edit example *** -------------------------------------------------------
|											|
|	Reference:			ATGCTGTGAT G TCGTGATG    A			|
|	Edit:				A--CTGTGAT C TCGTGATGatcgA			|
|	Sequence format:	A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A			|
|											|
 ---------------------------------------------------------------------------------------

''')

# Inputs for the design parameters of pegRNAs and nicking gRNAs
parser.add_argument('-pe_format', '--pe_format', type = str, default = 'NNNNNNNNNNNNNNNNN/NNN[NGG]', help = "***** Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM]. Examples: NNNNNNNNNNNNNNNNN/NNN[NGG], NNNNNNNNNNNNNNNNN/NNN[NG] *****\n\n")
parser.add_argument('-pbs', '--pbs_length_list', type = int, default = 0, nargs = '+', help = '***** List of the lengths of the primer binding sites for pegRNA extension *****\n\n')
parser.add_argument('-rt', '--rt_length_list', type = int, default = 0, nargs = '+', help = '***** List of the lengths of the reverse transcription template for pegRNA extension (default: 20 to 50nt with insertion) *****\n')
parser.add_argument('-nick_dist_min', '--nicking_distance_minimum', type = int, default = 0, nargs = '+', help = '***** Minimum nicking distance for PE3 in both directions *****\n\n')
parser.add_argument('-nick_dist_max', '--nicking_distance_maximum', type = int, default = 100, nargs = '+', help = '***** Maximum nicking distance for PE3 in both directions *****\n\n')

# Output directory
parser.add_argument('-out', '--out_dir', default = '0', type = str, help = '***** Name of output directory *****\n\n')

args = parser.parse_args()

##### Initialize arguments
file_in = args.file

pe_format = args.pe_format
pbs_length_list = args.pbs_length_list
rt_length_list = args.rt_length_list
nicking_distance_minimum = args.nicking_distance_minimum
nicking_distance_maximum = args.nicking_distance_maximum

# Default PBS and RTT lengths to design
if pbs_length_list == 0:
	pbs_length_list = list(range(11, 18))
if rt_length_list == 0:
	rt_length_list = list(range(10, 35, 5))

# Output directory date and time stamped
out_dir = args.out_dir
if out_dir == '0':
	out_dir = '%s_PrimeDesign' % str(time.strftime("%y%m%d_%H.%M.%S", time.localtime()))

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

# Initialize logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler(out_dir + '/PrimeDesign.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

##### IUPAC code map
iupac2bases_dict = {'A':'A','T':'T','C':'C','G':'G','a':'a','t':'t','c':'c','g':'g',
'R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACTG]',
'r':'[ag]','y':'[ct]','s':'[gc]','w':'[at]','k':'[gt]','m':'[ac]','b':'[cgt]','d':'[agt]','h':'[act]','v':'[acg]','n':'[actg]',
'(':'(',')':')','+':'+','-':'-','/':'/'}

def iupac2bases(iupac):

	try:
		bases = iupac2bases_dict[iupac]
	except:
		logger.error('Symbol %s is not within the IUPAC nucleotide code ...' % str(iupac))
		sys.exit(1)

	return(bases)

##### Reverse complement function
def reverse_complement(sequence):
	sequence = sequence
	new_sequence = ''
	for base in sequence:
		if base == 'A':
			new_sequence += 'T'
		elif base == 'T':
			new_sequence += 'A'
		elif base == 'C':
			new_sequence += 'G'
		elif base == 'G':
			new_sequence += 'C'
		if base == 'a':
			new_sequence += 't'
		elif base == 't':
			new_sequence += 'a'
		elif base == 'c':
			new_sequence += 'g'
		elif base == 'g':
			new_sequence += 'c'
		elif base == '[':
			new_sequence += ']'
		elif base == ']':
			new_sequence += '['
		elif base == '+':
			new_sequence += '+'
		elif base == '-':
			new_sequence += '-'
		elif base == '/':
			new_sequence += '/'
		elif base == '(':
			new_sequence += ')'
		elif base == ')':
			new_sequence += '('
	return(new_sequence[::-1])

##### Extract reference and edited sequence information
def process_sequence(input_sequence):

	# Check formatting is correct
	format_check = ''
	for i in input_sequence:
		if i == '(':
			format_check += '('
		elif i == ')':
			format_check += ')'
		elif i == '/':
			format_check += '/'
		elif i == '+':
			format_check += '+'
		elif i == '-':
			format_check += '-'

	# Check composition of input sequence
	if len(input_sequence) != sum([1 if x in ['A','T','C','G','(',')','+','-','/'] else 0 for x in input_sequence.upper()]):
		logger.error('Input sequence %s contains a character not in the following list: A,T,C,G,(,),+,-,/ ...' % str(input_sequence))
		sys.exit(1)

	# Check formatting
	if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
		if '((' not in format_check: # Checks both directions for nested parantheses
			if '()' not in format_check: # Checks for empty annotations
				if sum([1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']]) == 0:
					pass
				else:
					logger.error('Input sequence %s has more than one edit annotation per parantheses set (i.e. //,  +- , -/, etc.) ...' % str(input_sequence))
					sys.exit(1)
			else:
				logger.error('Input sequence %s has empty parantheses without an edit annotation (i.e. /,  + , -) ...' % str(input_sequence))
				sys.exit(1)
		else:
			logger.error('Input sequence %s has nested parantheses which is not allowed ...' % str(input_sequence))
			sys.exit(1)
	else:
		logger.error('Input sequence %s does not have full sets of parantheses ...' % str(input_sequence))
		sys.exit(1)

	# Create mapping between input format and reference and edit sequence
	editformat2sequence = {}
	edits = re.findall('\(.*?\)', input_sequence)
	for edit in edits:
		if '/' in edit:
			editformat2sequence[edit] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','')]
		elif '+' in edit:
			editformat2sequence[edit] = ['' , edit.split('+')[1].replace(')','')]
		elif '-' in edit:
			editformat2sequence[edit] = [edit.split('-')[1].replace(')',''), '']

	# Create mapping between edit number and reference and edit sequence
	editformat2sequence = {}
	editnumber2sequence = {}
	edit_idxs = [[m.start(), m.end()] for m in re.finditer('\(.*?\)', input_sequence)]
	edit_counter = 1
	for edit_idx in edit_idxs:
		edit = input_sequence[edit_idx[0]:edit_idx[1]]

		# Create edit format and number to sequence map
		if '/' in edit:
			editformat2sequence[edit] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower(), edit_counter]
			editnumber2sequence[edit_counter] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower()]

		elif '+' in edit:
			editformat2sequence[edit] = ['' , edit.split('+')[1].replace(')','').lower(), edit_counter]
			editnumber2sequence[edit_counter] = ['' , edit.split('+')[1].replace(')','').lower()]

		elif '-' in edit:
			editformat2sequence[edit] = [edit.split('-')[1].replace(')',''), '', edit_counter]
			editnumber2sequence[edit_counter] = [edit.split('-')[1].replace(')',''), '']

		edit_counter += 1

	edit_start = min([i.start() for i in re.finditer('\(', input_sequence)])
	edit_stop = max([i.start() for i in re.finditer('\)', input_sequence)])

	edit_span_sequence_w_ref = input_sequence[edit_start:edit_stop + 1]
	edit_span_sequence_w_edit = input_sequence[edit_start:edit_stop + 1]
	for edit in editformat2sequence:
		edit_span_sequence_w_ref = edit_span_sequence_w_ref.replace(edit, editformat2sequence[edit][0])
		edit_span_sequence_w_edit = edit_span_sequence_w_edit.replace(edit, editformat2sequence[edit][1])

	edit_start_in_ref = re.search('\(', input_sequence).start()
	edit_stop_in_ref_rev = re.search('\)', input_sequence[::-1]).start()

	edit_span_length_w_ref = len(edit_span_sequence_w_ref)
	edit_span_length_w_edit = len(edit_span_sequence_w_edit)

	reference_sequence = input_sequence
	edit_sequence = input_sequence
	editnumber_sequence = input_sequence
	for edit in editformat2sequence:
		reference_sequence = reference_sequence.replace(edit, editformat2sequence[edit][0])
		edit_sequence = edit_sequence.replace(edit, editformat2sequence[edit][1])
		editnumber_sequence = editnumber_sequence.replace(edit, str(editformat2sequence[edit][2]))

	return(editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev)

##### Dictionary for to organize different DNA targets
target_design = {}
with open(file_in, 'r') as f:
	for line1 in f:

		# Parse .txt files with space
		if file_in.endswith('.txt') or file_in.endswith('.TXT'):
			target_name, target_sequence = line1.strip('\n').split()

		# Parse .txt files with comma
		elif file_in.endswith('.csv') or file_in.endswith('.CSV'):
			target_name, target_sequence = line1.strip('\n').split(',')

		else:
			logger.error('Input file %s does not end with .txt or .csv ...' % str(file_in))
			sys.exit(1)

		target_sequence = target_sequence.upper()
		editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(target_sequence)

		# Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
		target_design[target_name] = {'target_sequence':target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}

##### Find cut index and reformat PE format parameter
if (pe_format.count('[') + pe_format.count(']')) == 2:

	if pe_format.count('/') == 1:

		# Find indices but shift when removing annotations
		cut_idx = re.search('/', pe_format).start()
		pam_start_idx = re.search('\[', pe_format).start()
		pam_end_idx = re.search('\]', pe_format).start()

		# Find pam and total PE format search length
		pam_length = pam_end_idx - pam_start_idx - 1
		pe_format_length = len(pe_format) - 3

		# Check if cut site is left of PAM
		if cut_idx < pam_start_idx:

			# Shift indices with removal of annotations
			pam_start_idx = pam_start_idx - 1
			pam_end_idx = pam_end_idx - 2
			spacer_start_idx = 0
			spacer_end_idx = pam_start_idx

		else:
			pam_end_idx = pam_end_idx - 1
			cut_idx = cut_idx - 2
			spacer_start_idx = pam_end_idx
			spacer_end_idx = len(pe_format) - 3
	
	else:
		logger.error('PE format parameter %s needs to cut site / within the spacer (i.e. NNNNNNNNNNNNNNNNN/NNN[NGG]) ...' % str(pe_format))
		sys.exit(1)

else:
	logger.error('PE format parameter %s needs to have one [PAM] present in its sequence (i.e. NNNNNNNNNNNNNNNNN/NNN[NGG]) ...' % str(pe_format))
	sys.exit(1)

# Remove annotations and convert into regex
pe_format_rm_annotation = pe_format.replace('/', '').replace('[', '').replace(']', '')
# print('---------- Prime editing spacer search parameters ----------')
# print('PE format:\t%s' % pe_format_rm_annotation)
# print('Spacer:\t\t%s' % pe_format_rm_annotation[spacer_start_idx:spacer_end_idx])
# print('PAM:\t\t%s' % pe_format_rm_annotation[pam_start_idx:pam_end_idx])

# Create PE format and PAM search sequences
pe_format_search_plus = ''
for base in pe_format_rm_annotation:
	pe_format_search_plus += iupac2bases(base)
pe_format_search_minus = reverse_complement(pe_format_search_plus)

pam_search = ''
pam_sequence = pe_format_rm_annotation[pam_start_idx:pam_end_idx]
for base in pam_sequence:
	pam_search += iupac2bases(base)

# print('PE search (+):\t%s' % pe_format_search_plus)
# print('PE search (-):\t%s' % pe_format_search_minus)
# print('\n')

##### Initialize data storage for output
pe_design = {}
logger.info('Searching for pegRNAs and nicking gRNAs for target sequences ...')
counter = 1
total_regions = len(target_design.keys())
for target_name in target_design:

	# pegRNA spacer search for (+) and (-) strands with reference sequence
	reference_sequence = target_design[target_name]['reference_sequence']
	find_guides_ref_plus = [[m.start()] for m in re.finditer('(?=%s)' % pe_format_search_plus, reference_sequence, re.IGNORECASE)]
	find_guides_ref_minus = [[m.start()] for m in re.finditer('(?=%s)' % pe_format_search_minus, reference_sequence, re.IGNORECASE)]

	# pegRNA spacer search for (+) and (-) strands with edit number sequence
	editnumber_sequence = target_design[target_name]['editnumber_sequence']
	find_guides_editnumber_plus = [[m.start()] for m in re.finditer('(?=%s)' % pam_search.replace('[', '[123456789'), editnumber_sequence, re.IGNORECASE)]
	find_guides_editnumber_minus = [[m.start()] for m in re.finditer('(?=%s)' % reverse_complement(pam_search).replace('[', '[123456789'), editnumber_sequence, re.IGNORECASE)]

	# Find pegRNA spacers targeting (+) strand
	if find_guides_ref_plus:

		for match in find_guides_ref_plus:

			# Extract matched sequences and annotate type of prime editing
			full_search = reference_sequence[match[0]:match[0] + pe_format_length]
			spacer_sequence = full_search[spacer_start_idx:spacer_end_idx]
			extension_core_sequence = full_search[:cut_idx]
			downstream_sequence_ref = full_search[cut_idx:]
			downstream_sequence_length = len(downstream_sequence_ref)
			pam_ref = full_search[pam_start_idx:pam_end_idx]

			# Check to see if the extended non target strand is conserved in the edited strand
			try:
				extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence, edit_sequence).start(), re.search(extension_core_sequence, edit_sequence).end()
				downstream_sequence_edit = edit_sequence[extension_core_end_idx:extension_core_end_idx + downstream_sequence_length]
				pam_edit = edit_sequence[extension_core_start_idx:extension_core_start_idx + pe_format_length][pam_start_idx:pam_end_idx]
				
				## Annotate pegRNA
				# Check if PAM is mutated relative to reference sequence
				if pam_ref == pam_edit.upper():
					pe_annotate = 'PAM_intact'

				else:
					# Check to see if mutation disrupts degenerate base positions within PAM
					if re.search(pam_search, pam_edit.upper()):
						pe_annotate = 'PAM_intact'

					else:
						pe_annotate = 'PAM_mutated'

				# Store pegRNA spacer
				nick_ref_idx = match[0] + cut_idx
				nick_edit_idx = extension_core_start_idx + cut_idx
				target_design[target_name]['pegRNA']['+'].append([nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

			except:
				continue

	# Find pegRNA spacers targeting (-) strand
	if find_guides_ref_minus:

		for match in find_guides_ref_minus:

			# Extract matched sequences and annotate type of prime editing
			full_search = reference_sequence[match[0]:match[0] + pe_format_length]
			spacer_sequence = full_search[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
			extension_core_sequence = full_search[pe_format_length - cut_idx:]
			downstream_sequence_ref = full_search[:pe_format_length - cut_idx]
			downstream_sequence_length = len(downstream_sequence_ref)
			pam_ref = full_search[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

			# Check to see if the extended non target strand is conserved in the edited strand
			try:
				extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence, edit_sequence).start(), re.search(extension_core_sequence, edit_sequence).end()
				downstream_sequence_edit = edit_sequence[extension_core_start_idx - downstream_sequence_length:extension_core_start_idx]
				pam_edit = edit_sequence[extension_core_end_idx - pe_format_length:extension_core_end_idx][pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]
				
				## Annotate pegRNA
				# Check if PAM is mutated relative to reference sequence
				if pam_ref == pam_edit.upper():
					pe_annotate = 'PAM_intact'

				else:
					# Check to see if mutation disrupts degenerate base positions within PAM
					if re.search(reverse_complement(pam_search), pam_edit.upper()):
						pe_annotate = 'PAM_intact'

					else:
						pe_annotate = 'PAM_mutated'

				# Store pegRNA spacer
				nick_ref_idx = match[0] + (pe_format_length - cut_idx)
				nick_edit_idx = extension_core_start_idx - downstream_sequence_length + (pe_format_length - cut_idx)
				target_design[target_name]['pegRNA']['-'].append([nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

			except:
				continue

	# Find ngRNA spacers targeting (+) strand
	if find_guides_editnumber_plus:

		for match in find_guides_editnumber_plus:

			# Extract matched sequences and annotate type of prime editing
			full_search = editnumber_sequence[:match[0] + pam_length]
			
			full_search2ref = full_search
			full_search2edit = full_search
			for edit_number in editnumber2sequence:
				full_search2ref = full_search2ref.replace(str(edit_number), editnumber2sequence[edit_number][0])
				full_search2edit = full_search2edit.replace(str(edit_number), editnumber2sequence[edit_number][1])

			if len(full_search2edit[-pe_format_length:]) == pe_format_length:

				# Identify ngRNA sequence information from edit sequence
				full_search_edit = full_search2edit[-pe_format_length:]
				spacer_sequence_edit = full_search_edit[spacer_start_idx:spacer_end_idx]
				pam_edit = full_search_edit[pam_start_idx:pam_end_idx]

				# Use reference sequence to find nick index
				full_search_ref = full_search2ref[-pe_format_length:]
				spacer_sequence_ref = full_search_ref[spacer_start_idx:spacer_end_idx]
				pam_ref = full_search_ref[pam_start_idx:pam_end_idx]

				# Annotate ngRNA
				if spacer_sequence_edit.upper()	== spacer_sequence_ref.upper():
					ng_annotate = 'PE3'
				else:
					ng_annotate = 'PE3b'

				# Store ngRNA spacer
				nick_ref_idx = re.search(full_search_ref, reference_sequence).end() - (pe_format_length - cut_idx)
				target_design[target_name]['ngRNA']['+'].append([nick_ref_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

	# Find ngRNA spacers targeting (-) strand
	if find_guides_editnumber_minus:

		for match in find_guides_editnumber_minus:

			# Extract matched sequences and annotate type of prime editing
			full_search = editnumber_sequence[match[0]:]
			
			full_search2ref = full_search
			full_search2edit = full_search
			for edit_number in editnumber2sequence:
				full_search2ref = full_search2ref.replace(str(edit_number), editnumber2sequence[edit_number][0])
				full_search2edit = full_search2edit.replace(str(edit_number), editnumber2sequence[edit_number][1])

			if len(full_search2edit[:pe_format_length]) == pe_format_length:

				# Identify ngRNA sequence information from edit sequence
				full_search_edit = full_search2edit[:pe_format_length]
				spacer_sequence_edit = full_search_edit[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
				pam_edit = full_search_edit[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

				# Use reference sequence to find nick index
				full_search_ref = full_search2ref[:pe_format_length]
				spacer_sequence_ref = full_search_ref[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
				pam_ref = full_search_ref[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

				# Annotate ngRNA
				if spacer_sequence_edit.upper()	== spacer_sequence_ref.upper():
					ng_annotate = 'PE3'
				else:
					ng_annotate = 'PE3b'

				# Store ngRNA spacer
				nick_ref_idx = re.search(full_search_ref, reference_sequence).start() + (pe_format_length - cut_idx)
				target_design[target_name]['ngRNA']['-'].append([nick_ref_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

	# Grab index information of edits to introduce to target sequence
	edit_start_in_ref = int(target_design[target_name]['edit_start_in_ref'])
	edit_stop_in_ref_rev = int(target_design[target_name]['edit_stop_in_ref_rev'])
	edit_span_length_w_ref = int(target_design[target_name]['edit_span_length'][0])
	edit_span_length_w_edit = int(target_design[target_name]['edit_span_length'][1])


	pe_design[target_name] = {}
	
	# Design pegRNAs targeting the (+) strand
	for peg_plus in target_design[target_name]['pegRNA']['+']:

		pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus
		pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))

		# See if pegRNA spacer can introduce all edits
		nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
		if nick2edit_length >= 0:

			# Loop through RTT lengths
			for rt_length in rt_length_list:

				# See if RT length can reach entire edit
				if (nick2edit_length + edit_span_length_w_edit) < rt_length:

					# Loop through PBS lengths
					for pbs_length in pbs_length_list:

						# Construct pegRNA extension to encode intended edit(s)
						pegRNA_ext = edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rt_length]

						# Check to see if pegRNA extension is within input sequence
						if len(pegRNA_ext) == (pbs_length + rt_length):

							# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
							if pegid not in pe_design[target_name]:

								# First list is for peg extension, second list is for nicking guide
								pe_design[target_name][pegid] = [[],[]]

							pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rt_length, pegRNA_ext])

			# Create ngRNAs targeting (-) strand for (+) pegRNAs
			if pegid in pe_design[target_name]:
				for ng_minus in target_design[target_name]['ngRNA']['-']:
					ng_nick_ref_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
					nick_distance = ng_nick_ref_idx - pe_nick_ref_idx
					if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
						pe_design[target_name][pegid][1].append([ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, '-', nick_distance])

	# Design pegRNAs targeting the (-) strand
	for peg_minus in target_design[target_name]['pegRNA']['-']:

		pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
		pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

		# See if pegRNA spacer can introduce all edits
		nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
		if nick2edit_length >= 0:

			# Loop through RTT lengths
			for rt_length in rt_length_list:

				# See if RT length can reach entire edit
				if (nick2edit_length + edit_span_length_w_edit) < rt_length:

					# Loop through PBS lengths
					for pbs_length in pbs_length_list:

						# Construct pegRNA extension to encode intended edit(s)
						pegRNA_ext = edit_sequence[pe_nick_edit_idx - rt_length:pe_nick_edit_idx + pbs_length]

						# Check to see if pegRNA extension is within input sequence
						if len(pegRNA_ext) == (pbs_length + rt_length):

							# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
							if pegid not in pe_design[target_name]:

								# First list is for peg extension, second list is for nicking guide
								pe_design[target_name][pegid] = [[],[]]

							pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-', pbs_length, rt_length, pegRNA_ext])

			# Create ngRNAs targeting (+) strand for (-) pegRNAs
			if pegid in pe_design[target_name]:
				for ng_plus in target_design[target_name]['ngRNA']['+']:
					ng_nick_ref_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
					nick_distance = ng_nick_ref_idx - pe_nick_ref_idx
					if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
						pe_design[target_name][pegid][1].append([ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, '+', nick_distance])

	logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter, total_regions))
	counter += 1

# print('------- Target design -------')
# for target_name in target_design:
# 	print('I:\t%s' % target_design[target_name]['target_sequence'])
# 	print('R:\t%s' % target_design[target_name]['reference_sequence'])
# 	print('E:\t%s' % target_design[target_name]['edit_sequence'])
# 	print('+')

# 	print('pegRNAs targeting (+) strand')
# 	for entry in target_design[target_name]['pegRNA']['+']:
# 		print('\t'.join(map(str, entry)))

# 	print('pegRNAs targeting (-) strand')
# 	for entry in target_design[target_name]['pegRNA']['-']:
# 		print('\t'.join(map(str, entry)))

# 	print('ngRNAs targeting (+) strand')
# 	for entry in target_design[target_name]['ngRNA']['+']:
# 		print('\t'.join(map(str, entry)))

# 	print('ngRNAs targeting (-) strand')
# 	for entry in target_design[target_name]['ngRNA']['-']:
# 		print('\t'.join(map(str, entry)))

# print('------- PE design -------')
# for target_name in pe_design:
# 	for pegid in pe_design[target_name]:
# 		print('-- pegRNA --')
# 		for entry in pe_design[target_name][pegid][0]:
# 			print('\t'.join(map(str, entry)))

# 		print('-- ngRNAs --')
# 		for entry in pe_design[target_name][pegid][1]:
# 			print('\t'.join(map(str, entry)))

# Output pegRNAs
pegRNAs_summary_f = '%s_PrimeDesign_summary.csv' % str(time.strftime("%m-%d-%y_%H.%M.%S", time.localtime()))
logger.info('Writing pegRNA and ngRNA designs into output file %s ...' % pegRNAs_summary_f)

counter = 1
with open(out_dir + '/%s' % pegRNAs_summary_f, 'w') as f:
	
	f.write(','.join(map(str, ['Target_name', 'Target_sequence', 'pegRNA_number', 'gRNA_type', 'Spacer_sequence', 'PAM_sequence', 'Extension_sequence', 'Strand', 'Annotation', 'Nick_index', 'Distance_from_PE_nick', 'PBS_length', 'RT_length','First_extension_nucleotide', 'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM', 'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM'])) + '\n')
	for target_name in pe_design:
		
		for pegid in pe_design[target_name]:

			# Write pegRNAs
			for pegRNA_entry in pe_design[target_name][pegid][0]:
				pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rt_length, pegRNA_ext = pegRNA_entry

				if pe_strand == '+':
					pegRNA_ext = reverse_complement(pegRNA_ext)
					pegRNA_ext_first_base = pegRNA_ext[0]

				elif pe_strand == '-':
					pe_spacer_sequence = reverse_complement(pe_spacer_sequence)
					pe_pam_ref = reverse_complement(pe_pam_ref)
					pegRNA_ext_first_base = pegRNA_ext[0]

				f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, pe_nick_ref_idx, '', pbs_length, rt_length, pegRNA_ext_first_base, 'caccG' + pe_spacer_sequence[1:] + 'gtttt', 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence[1:]), 'gtgc' + pegRNA_ext, 'aaaa' + reverse_complement(pegRNA_ext)])) + '\n')

			# Write ngRNAs
			for ngRNA_entry in pe_design[target_name][pegid][1]:
				ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

				if ng_strand == '-':
					ng_spacer_sequence = reverse_complement(ng_spacer_sequence_edit)
					ng_pam = reverse_complement(ng_pam_edit)

				f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence, ng_pam, '', ng_strand, ng_annotate, ng_nick_ref_idx, nick_distance, '', '', '', 'caccG' + ng_spacer_sequence[1:], 'aaac' + reverse_complement('G' + ng_spacer_sequence[1:]), '', ''])) + '\n')

			counter += 1