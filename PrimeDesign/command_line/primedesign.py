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
Software for the design of pegRNAs for flexible prime editing! Please visit ----- https://github.com/jyhsu15/PrimeDesign ----- for more documentation on how to use the software.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------''', formatter_class=RawTextHelpFormatter)

# Inputs for de-novo design of pegRNAs and nicking gRNAs
parser.add_argument('-f', '--file', required = True, type = str, help = '''Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (Required)

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
parser.add_argument('-pe_format', '--pe_format', type = str, default = 'NNNNNNNNNNNNNNNNN/NNN[NGG]', help = "***** Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM] (Default: NNNNNNNNNNNNNNNNN/NNN[NGG]). Examples: NNNNNNNNNNNNNNNNN/NNN[NGG], NNNNNNNNNNNNNNNNN/NNN[NG] *****\n\n")
parser.add_argument('-pbs', '--pbs_length_list', type = int, default = 0, nargs = '+', help = '***** List of primer binding site (PBS) lengths for the pegRNA extension (Default: 10 to 16 nt). Example: 12 13 14 15 *****\n\n')
parser.add_argument('-rt', '--rtt_length_list', type = int, default = 0, nargs = '+', help = '***** List of reverse transcription (RT) template lengths for the pegRNA extension (Default: 10 to 16 nt). Example: 10 15 20 *****\n')
parser.add_argument('-nick_dist_min', '--nicking_distance_minimum', type = int, default = 0, nargs = '+', help = '***** Minimum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 0). *****\n\n')
parser.add_argument('-nick_dist_max', '--nicking_distance_maximum', type = int, default = 100, nargs = '+', help = '***** Maximum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 100). *****\n\n')
parser.add_argument('-silent_mut', '--silent_mutation', default = 'N', choices = ['N', 'n', 'Y', 'y'], type = str, help = '***** Introduce silent mutation into PAM assuming sequence is in-frame. Currently only available with SpCas9 PE (Default: N). *****\n\n')

# Output directory
parser.add_argument('-out', '--out_dir', default = '0', type = str, help = '***** Name of output directory (Default: ./DATETIMESTAMP_PrimeDesign). *****\n\n')

args = parser.parse_args()

##### Initialize arguments
file_in = args.file

pe_format = args.pe_format
pbs_length_list = args.pbs_length_list
rtt_length_list = args.rtt_length_list
nicking_distance_minimum = args.nicking_distance_minimum
nicking_distance_maximum = args.nicking_distance_maximum
silent_mutation = args.silent_mutation

# Default PBS and RTT lengths to design
if pbs_length_list == 0:
	pbs_length_list = list(range(10, 16))
if rtt_length_list == 0:
	rtt_length_list = list(range(10, 16))

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

# GC content
def gc_content(sequence):
    sequence = sequence.upper()
    GC_count = sequence.count('G') + sequence.count('C')
    GC_content = float(GC_count)/float(len(sequence))

    return("%.2f" % GC_content)

# Reverse complement function
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

# Amino acid code
codon_dict = {
    'GGG':['Gly','G', 0.25],'GGA':['Gly','G', 0.25],'GGT':['Gly','G', 0.16],'GGC':['Gly','G', 0.34],
    'GAG':['Glu','E', 0.58],'GAA':['Glu','E', 0.42],'GAT':['Asp','D', 0.46],'GAC':['Asp','D', 0.54],
    'GTG':['Val','V', 0.47],'GTA':['Val','V', 0.11],'GTT':['Val','V', 0.18],'GTC':['Val','V', 0.24],
    'GCG':['Ala','A', 0.11],'GCA':['Ala','A', 0.23],'GCT':['Ala','A', 0.26],'GCC':['Ala','A', 0.4],
    'AGG':['Arg','R', 0.2],'AGA':['Arg','R', 0.2],'AGT':['Ser','S', 0.15],'AGC':['Ser','S', 0.24],
    'AAG':['Lys','K', 0.58],'AAA':['Lys','K', 0.42],'AAT':['Asn','N', 0.46],'AAC':['Asn','N', 0.54],
    'ATG':['Met','M', 1],'ATA':['Ile','I', 0.16],'ATT':['Ile','I', 0.36],'ATC':['Ile','I', 0.48],
    'ACG':['Thr','T', 0.12],'ACA':['Thr','T', 0.28],'ACT':['Thr','T', 0.24],'ACC':['Thr','T', 0.36],
    'TGG':['Trp','W', 1],'TGA':['End','X', 0.52],'TGT':['Cys','C', 0.45],'TGC':['Cys','C', 0.55],
    'TAG':['End','X', 0.2],'TAA':['End','X', 0.28],'TAT':['Tyr','Y', 0.43],'TAC':['Tyr','Y', 0.57],
    'TTG':['Leu','L', 0.13],'TTA':['Leu','L', 0.07],'TTT':['Phe','F', 0.45],'TTC':['Phe','F', 0.55],
    'TCG':['Ser','S', 0.06],'TCA':['Ser','S', 0.15],'TCT':['Ser','S', 0.18],'TCC':['Ser','S', 0.22],
    'CGG':['Arg','R', 0.21],'CGA':['Arg','R', 0.11],'CGT':['Arg','R', 0.08],'CGC':['Arg','R', 0.19],
    'CAG':['Gln','Q', 0.75],'CAA':['Gln','Q', 0.25],'CAT':['His','H', 0.41],'CAC':['His','H', 0.59],
    'CTG':['Leu','L', 0.41],'CTA':['Leu','L', 0.07],'CTT':['Leu','L', 0.13],'CTC':['Leu','L', 0.2],
    'CCG':['Pro','P', 0.11],'CCA':['Pro','P', 0.27],'CCT':['Pro','P', 0.28],'CCC':['Pro','P', 0.33],
}

# Create codon swap dictionaries
aa2codon = {}
for codon in codon_dict:
    if codon_dict[codon][1] not in aa2codon:
        aa2codon[codon_dict[codon][1]] = []

    aa2codon[codon_dict[codon][1]].append([codon, codon_dict[codon][2]])

codon_swap_0 = {}
codon_swap_1_1 = {}
codon_swap_1_2 = {}
codon_swap_2 = {}
for codon in codon_dict:

    codon_swap_0[codon] = []
    codon_swap_1_1[codon] = []
    codon_swap_1_2[codon] = []
    codon_swap_2[codon] = []

    for other_codon in aa2codon[codon_dict[codon][1]]:

        # Check if PAM disrupted with silent mutations
        if codon[1:] != other_codon[0][1:]:
            codon_swap_0[codon].append(other_codon)

        if codon[2:] != other_codon[0][2:]:
            codon_swap_1_1[codon].append(other_codon)

        if codon[:1] != other_codon[0][:1]:
            codon_swap_1_2[codon].append(other_codon)

        if codon[:2] != other_codon[0][:2]:
            codon_swap_2[codon].append(other_codon)

for codon in codon_dict:
    codon_swap_0[codon] = sorted(codon_swap_0[codon], key = lambda x: x[1], reverse = True)
    codon_swap_1_1[codon] = sorted(codon_swap_1_1[codon], key = lambda x: x[1], reverse = True)
    codon_swap_1_2[codon] = sorted(codon_swap_1_2[codon], key = lambda x: x[1], reverse = True)
    codon_swap_2[codon] = sorted(codon_swap_2[codon], key = lambda x: x[1], reverse = True)

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

	editnumber2sequence = target_design[target_name]['editnumber2sequence']

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

		pe_annotate_constant = pe_annotate

		# See if pegRNA spacer can introduce all edits
		nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
		if nick2edit_length >= 0:

			# Loop through RTT lengths
			for rtt_length in rtt_length_list:

				# See if RT length can reach entire edit
				nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
				if nick2lastedit_length < rtt_length:

					# Loop through PBS lengths
					for pbs_length in pbs_length_list:
						pe_pam_ref_silent_mutation = ''

						# Construct pegRNA extension to encode intended edit(s)
						# pegRNA_ext = edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length]
						# Patch for NGG PAMs - may need to build something more generalizable in the future
						if (silent_mutation.upper() == 'Y') and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):

							if pe_annotate_constant == 'PAM_intact':

								nick_aa_index = int(pe_nick_edit_idx)%3

								if nick_aa_index == 0:
									original_codon = edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 6].upper()

									if len(codon_swap_0[original_codon.upper()]) > 1:
										new_codon = codon_swap_0[original_codon][0][0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 3] + new_codon + edit_sequence[pe_nick_edit_idx + 6:pe_nick_edit_idx + rtt_length])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon
										pe_annotate = 'PAM_mutated_silent_mutation'

									else:
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

								elif nick_aa_index == 1:
									original_codon_1 = edit_sequence[pe_nick_edit_idx + 2:pe_nick_edit_idx + 5].upper()
									original_codon_2 = edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + 8].upper()

									if len(codon_swap_1_1[original_codon_1.upper()]) > 1:
										new_codon = codon_swap_1_1[original_codon_1][0][0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 2] + new_codon + edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + rtt_length])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon[1:] + original_codon_2[:1].lower()
										pe_annotate = 'PAM_mutated_silent_mutation'

									elif len(codon_swap_1_2[original_codon_2.upper()]) > 1:
										new_codon = codon_swap_1_2[original_codon_2][0][0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:pe_nick_edit_idx + rtt_length])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + original_codon_1[1:].lower() + new_codon[:1]
										pe_annotate = 'PAM_mutated_silent_mutation'

									else:
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

								elif nick_aa_index == 2:
									original_codon = edit_sequence[pe_nick_edit_idx + 4:pe_nick_edit_idx + 7].upper()

									if len(codon_swap_2[original_codon.upper()]) > 1:
										new_codon = codon_swap_2[original_codon][0][0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 4] + new_codon + edit_sequence[pe_nick_edit_idx + 7:pe_nick_edit_idx + rtt_length])
										pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 4].lower() + new_codon[:2]
										pe_annotate = 'PAM_mutated_silent_mutation'

									else:
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

							else:
								pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

						else:
							pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

						# Check to see if pegRNA extension is within input sequence
						if len(pegRNA_ext) == (pbs_length + rtt_length):

							# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
							if pegid not in pe_design[target_name]:

								# First list is for peg extension, second list is for nicking guide
								pe_design[target_name][pegid] = [[],[]]

							if pe_pam_ref_silent_mutation == '':
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext])

							else:
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext])

			# Create ngRNAs targeting (-) strand for (+) pegRNAs
			if pegid in pe_design[target_name]:
				for ng_minus in target_design[target_name]['ngRNA']['-']:
					ng_nick_ref_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
					nick_distance = ng_nick_ref_idx - pe_nick_ref_idx
					if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
						pe_design[target_name][pegid][1].append([ng_nick_ref_idx, reverse_complement(ng_spacer_sequence_edit), reverse_complement(ng_pam_edit), ng_annotate, '-', nick_distance])

	# Design pegRNAs targeting the (-) strand
	for peg_minus in target_design[target_name]['pegRNA']['-']:

		pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
		pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

		pe_annotate_constant = pe_annotate

		# See if pegRNA spacer can introduce all edits
		nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
		if nick2edit_length >= 0:

			# Loop through RTT lengths
			for rtt_length in rtt_length_list:

				# See if RT length can reach entire edit
				nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
				if nick2lastedit_length < rtt_length:

					# Loop through PBS lengths
					for pbs_length in pbs_length_list:
						pe_pam_ref_silent_mutation = ''

						# Construct pegRNA extension to encode intended edit(s)
						# pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
						# Patch for NGG PAMs - may need to build something more generalizable in the future
						if (silent_mutation.upper() == 'Y') and (pe_format == 'NNNNNNNNNNNNNNNNN/NNN[NGG]'):
						    
							if pe_annotate_constant == 'PAM_intact':

								nick_aa_index = int(pe_nick_edit_idx)%3
						        
								if nick_aa_index == 0:
									original_codon = edit_sequence[pe_nick_edit_idx - 6:pe_nick_edit_idx - 3].upper()

									if len(codon_swap_2[original_codon.upper()]) > 1:
										new_codon = codon_swap_2[original_codon][0][0].lower()
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 6] + new_codon + edit_sequence[pe_nick_edit_idx - 3:pe_nick_edit_idx + pbs_length]
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon)
										pe_annotate = 'PAM_mutated_silent_mutation'

									else:
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

								elif nick_aa_index == 1:
									original_codon = edit_sequence[pe_nick_edit_idx - 7:pe_nick_edit_idx - 4].upper()

									if len(codon_swap_0[original_codon.upper()]) > 1:
										new_codon = codon_swap_0[original_codon][0][0].lower()
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 7] + new_codon + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx + pbs_length]
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon[1:] + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx - 3].lower())
										pe_annotate = 'PAM_mutated_silent_mutation'

									else:
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

								elif nick_aa_index == 2:
									original_codon_1 = edit_sequence[pe_nick_edit_idx - 8:pe_nick_edit_idx - 5].upper()
									original_codon_2 = edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx - 2].upper()

									if len(codon_swap_1_1[original_codon_1.upper()]) > 1:
										new_codon = codon_swap_1_1[original_codon_1][0][0].lower()
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 8] + new_codon + edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx + pbs_length]
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon[2:] + original_codon_2[:2].lower())
										pe_annotate = 'PAM_mutated_silent_mutation'

									elif len(codon_swap_1_2[original_codon_2.upper()]) > 1:
										new_codon = codon_swap_1_2[original_codon_2][0][0].lower()
										pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:pe_nick_edit_idx + rtt_length])
										pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(original_codon_1[2:].lower() + new_codon[:2])
										pe_annotate = 'PAM_mutated_silent_mutation'

									else:
										pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

							else:
								pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

						else:
							pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length] ########

						# Check to see if pegRNA extension is within input sequence
						if len(pegRNA_ext) == (pbs_length + rtt_length):

							# Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
							if pegid not in pe_design[target_name]:

								# First list is for peg extension, second list is for nicking guide
								pe_design[target_name][pegid] = [[],[]]

							if pe_pam_ref_silent_mutation == '':
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext])
							
							else:
								pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext])

			# Create ngRNAs targeting (+) strand for (-) pegRNAs
			if pegid in pe_design[target_name]:
				for ng_plus in target_design[target_name]['ngRNA']['+']:
					ng_nick_ref_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
					nick_distance = ng_nick_ref_idx - pe_nick_ref_idx
					if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):
						pe_design[target_name][pegid][1].append([ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, '+', nick_distance])

	if counter%1000 == 0:
		logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter, total_regions))
	counter += 1

logger.info('Completed pegRNA and ngRNA search for %s out of %s sites ...' % (counter - 1, total_regions))

# Output pegRNAs
pegRNAs_summary_f = '%s_PrimeDesign.csv' % str(time.strftime("%Y%m%d_%I.%M.%S", time.localtime()))
logger.info('Writing pegRNA and ngRNA designs into output file %s ...' % pegRNAs_summary_f)

counter = 1
with open(out_dir + '/%s' % pegRNAs_summary_f, 'w') as f:
	
	f.write(','.join(map(str, ['Target_name', 'Target_sequence', 'pegRNA_number', 'gRNA_type', 'Spacer_sequence', 'PAM_sequence', 'Extension_sequence', 'Strand', 'Annotation', 'Nick_index', 'Distance_from_PE_nick', 'PBS_length', 'rtt_length','First_extension_nucleotide', 'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM', 'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM'])) + '\n')
	for target_name in pe_design:
		
		for pegid in pe_design[target_name]:

			# Write pegRNAs
			for pegRNA_entry in pe_design[target_name][pegid][0]:
				pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext = pegRNA_entry

				pegRNA_ext_first_base = pegRNA_ext[0]

				f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, pe_nick_ref_idx, '', pbs_length, rtt_length, pegRNA_ext_first_base, 'caccG' + pe_spacer_sequence[1:] + 'gtttt', 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence[1:]), 'gtgc' + pegRNA_ext, 'aaaa' + reverse_complement(pegRNA_ext)])) + '\n')

			# Write ngRNAs
			for ngRNA_entry in pe_design[target_name][pegid][1]:
				ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

				f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, ng_pam_edit, '', ng_strand, ng_annotate, ng_nick_ref_idx, nick_distance, '', '', '', 'caccG' + ng_spacer_sequence_edit[1:], 'aaac' + reverse_complement('G' + ng_spacer_sequence_edit[1:]), '', ''])) + '\n')

			counter += 1