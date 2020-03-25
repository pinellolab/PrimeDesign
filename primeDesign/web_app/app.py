# Import required libraries
import re
import os
import base64
import urllib.parse
import dash
import dash_table
import dash_bio as dashbio
import dash_html_components as html
import dash_bootstrap_components as dbc
from flask import Flask
import pandas as pd
from dash.dependencies import Input, Output, State, ClientsideFunction
import dash_core_components as dcc
import dash_html_components as html

'https://codepen.io/chriddyp/pen/bWLwgP.css'

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets = external_stylesheets)
server = app.server

peg_design_tmp = {'pegRNA group':[],'type':[], 'spacer sequence':[],'PAM':[],'strand':[],'peg-to-edit':[],'nick-to-peg':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'RTT length':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
df_tmp = pd.DataFrame.from_dict(peg_design_tmp)

app.layout = html.Div([

    html.Img(src=app.get_asset_url('primedesign_logo.png'), width = '20%', style = {'margin-bottom': '0px'}),
    # html.H6('Design tool for prime editing', style = {'color':'grey', 'margin-top': '0px'}),

    html.Div([

        dcc.Textarea(
            id='pe-sequence-input',
            placeholder='Enter sequence to prime edit ...\n\nEdit formatting: Substitutions (ATGC/CGTA)  |  Insertions (+ATGC)  |  Deletions (-ATGC)',
            value='',
            style={'width': '99.2%'}
        ),

        html.Label(id = 'input-check', children = '', style = {'font-weight':'bold'})
        ]),

    # html.Br(),

    html.Div([

        html.Div([

            html.H5('Visualize sequence'),
            dcc.RadioItems(
                id = 'sequence-option',
                options=[
                    {'label': 'Reference', 'value': 'ref'},
                    {'label': 'Edited', 'value': 'edit'},
                ],
                value='ref',
                labelStyle={'display': 'inline-block'}
            ),

            dashbio.SequenceViewer(
                id = 'reference-edit-sequence',
                sequence = '...',
                badge =False,
                charsPerLine = 150,
                sequenceMaxHeight = '10000px',
                search = True,
                coverage = [],
                legend = [{'name':'Substitution', 'color':'#1E90FF', 'underscore':False}, {'name':'Insertion', 'color':'#3CB371', 'underscore':False}, {'name':'Deletion', 'color':'#DC143C', 'underscore':False}, {'name':'Selected pegRNA spacer', 'color':'#d6d6d6', 'underscore':False}]
            ),

            html.Div(id='store-sequence', style={'display': 'none'}),

            ], style={'width': '97%','display': 'inline-block','border-radius': '5px','box-shadow': '2px 2px 2px lightgrey','background-color': '#fafafa','padding': '15px','margin': '10px'}),

        ], className = 'row'),
    
    html.Br(),

    html.Div([

        html.Div([

            html.H5('Prime editing parameters'),

            html.Div([

                html.Label(id = 'pbs-title', children = 'PBS length', style = {'font-weight':'bold'}),
                html.Span('?',
                      id = 'pbs-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Recommendation: ~13 nt',
                       target = 'pbs-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'pbs-info', children = 'Primer binding site', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'pbs-range',
                min=5,
                max=17,
                value=[11, 15],
                allowCross=False
            ),

            html.Div([

                html.Label(id = 'rtt-title', children = 'RTT length', style = {'font-weight':'bold'}),
                html.Span('?',
                      id = 'rtt-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Recommendation: ~10-16 nt',
                       target = 'rtt-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'rtt-info', children = 'Reverse transcription template', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'rtt-range',
                min=5,
                max=80,
                value=[10, 16],
                allowCross=False
            ),

            
            html.Div([
                html.Label(id = 'nick-dist-title', children = 'ngRNA distance', style = {'font-weight':'bold'}),
                html.Span('?',
                      id = 'nick-dist-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Recommendation: ~50 bp',
                       target = 'nick-dist-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'nick-dist-info', children = 'ngRNA to pegRNA distance', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'nick-dist-range',
                min=0,
                max=120,
                value=[0, 100],
                allowCross=False
            ),

            html.Label('Remove pegRNA extensions with C first base', style = {'font-weight':'bold'}),

            dcc.RadioItems(
                id = 'extfirstbase-option',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='yes',
                labelStyle={'display': 'inline-block'}
            ),

            ], className = 'three columns', style={'display': 'inline-block','border-radius': '5px','box-shadow': '2px 2px 2px lightgrey','background-color': '#fafafa','padding': '15px','margin-top': '10px', 'margin-bottom': '10px', 'margin-left': '10px', 'margin-right': '10px'}
            ),

        html.Div([

            html.H6('pegRNA spacer table'),

            # ['pegRNA group','type', 'spacer sequence','PAM','strand','peg-to-edit','nick-to-peg','pegRNA extension','PBS length','RTT length','annotation']
            dash_table.DataTable(
                id = 'peg-table',
                columns = [{'name': i, 'id': i} for i in ['spacer sequence','PAM','strand','peg-to-edit','annotation']],
                data = df_tmp.to_dict('records'),
                style_cell={'textAlign': 'left', 'padding': '5px',},
                # style_as_list_view=True,
                style_header={
                    'backgroundColor': 'white',
                    'fontWeight': 'bold'
                },
                sort_action = 'native',
                sort_mode = 'multi',
                # filter_action = 'native',
                row_selectable = 'multi',
            ),

            html.H6('pegRNA extension table'),

            dash_table.DataTable(
                id = 'pegext-table',
                columns = [{'name': i, 'id': i} for i in ['PBS length','RTT length','pegRNA extension']],
                data = df_tmp.to_dict('records'),
                style_cell={'textAlign': 'left', 'padding': '5px'},
                # style_as_list_view=True,
                style_header={
                    'backgroundColor': 'white',
                    'fontWeight': 'bold'
                },
                sort_action = 'native',
                sort_mode = 'multi',
                # filter_action = 'native',
            ),

            html.H6('ngRNA spacer table'),

            dash_table.DataTable(
                id = 'ng-table',
                columns = [{'name': i, 'id': i} for i in ['spacer sequence','PAM','strand','nick-to-peg','annotation']],
                data = df_tmp.to_dict('records'),
                style_cell={'textAlign': 'left', 'padding': '5px'},
                # style_as_list_view=True,
                style_header={
                    'backgroundColor': 'white',
                    'fontWeight': 'bold'
                },
                sort_action = 'native',
                sort_mode = 'multi',
                # filter_action = 'native',
            ),

            html.A(
                'Download all designs',
                id='download-link',
                download="PrimeDesign.csv",
                href="",
                target="_blank"
            ),

            html.Div(id='store-peg-table-total', style={'display': 'none'}),
            html.Div(id='store-peg-table', style={'display': 'none'}),

            ], className = 'nine columns', style={'display': 'inline-block','border-radius': '5px','box-shadow': '2px 2px 2px lightgrey','background-color': '#fafafa','padding': '15px','margin-top': '10px', 'margin-bottom': '10px', 'margin-left': '10px', 'margin-right': '10px'}),

        ], className = 'row'),
    ], style={})

@app.callback(Output('input-check', 'children'),
    [Input('pe-sequence-input','value')]
)

def update_input_check(input_sequence):

    if len(input_sequence) > 0:

        input_sequence = ''.join(input_sequence.split())
        if len(input_sequence) <= 10000:
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
                sequence_check = 'Error: Input sequence contains a character not in the following list: A,T,C,G,(,),+,-,/ ...'

            else:

                # Check formatting
                if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
                    if '((' not in format_check: # Checks both directions for nested parantheses
                        if '()' not in format_check: # Checks for empty annotations
                            if sum([1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']]) == 0:
                                sequence_check = 'Success: Input sequence has correct formatting'
                            else:
                                sequence_check = 'Error: Input sequence has more than one edit annotation per parantheses set or annotation outside of parantheses'
                        else:
                            sequence_check = 'Error: Input sequence has empty parantheses without an edit annotation (i.e. /,  + , -)'
                    else:
                        sequence_check = 'Error: Input sequence has nested parantheses which is not allowed'
                else:
                    sequence_check = 'Error: Input sequence does not have full sets of parantheses'

        else:
            sequence_check = 'Error: Input sequence has exceeded maximum length of 10kb'

    else:
        sequence_check = 'No input sequence with desired edits has been provided'

    return(sequence_check)

@app.callback([Output('reference-edit-sequence', 'sequence'), Output('reference-edit-sequence', 'coverage')],
    [Input('input-check','children'), Input('sequence-option', 'value')],
    state = [State('pe-sequence-input','value'), State('peg-table','selected_rows'), State('store-peg-table', 'children')]
)

def update_reference_sequence(input_check, sequence_option, input_sequence, selected_rows, store_peg_table):

    input_sequence = ''.join(input_sequence.split())
    editformat2sequence = {}
    editformat2sequencestore = {}
    index_shift = 0
    annotations = []
    if 'Success' in input_check:

        edit_idxs = [[m.start(), m.end()] for m in re.finditer('\(.*?\)', input_sequence)]
        for edit_idx in edit_idxs:

            edit = input_sequence[edit_idx[0]:edit_idx[1]]
            edit_length = edit_idx[1] - edit_idx[0]

            if sequence_option == 'ref':

                # Create edit format and number to sequence map
                if '/' in edit:
                    editformat2sequence[edit] = edit.split('/')[0].replace('(','')
                    annotations.append({'start':edit_idx[0] - index_shift, 'end':edit_idx[0] - index_shift + len(edit.split('/')[0].replace('(','')), 'color':'#1E90FF', 'bgcolor':'#e8f3ff', 'underscore':True})
                    index_shift += edit_length - len(edit.split('/')[0].replace('(',''))

                elif '+' in edit:
                    editformat2sequence[edit] = ''
                    annotations.append({'start':edit_idx[0] - index_shift, 'end':edit_idx[0] - index_shift, 'color':'#3CB371', 'bgcolor':'#ebf7f0', 'underscore':True})
                    index_shift += edit_length

                elif '-' in edit:
                    editformat2sequence[edit] = edit.split('-')[1].replace(')','')
                    annotations.append({'start':edit_idx[0] - index_shift, 'end':edit_idx[0] - index_shift + len(edit.split('-')[1].replace(')','')), 'color':'#DC143C', 'bgcolor':'#fbe7eb', 'underscore':True})
                    index_shift += edit_length - len(edit.split('-')[1].replace(')',''))

            elif sequence_option == 'edit':
                # Create edit format and number to sequence map
                if '/' in edit:
                    editformat2sequence[edit] = edit.split('/')[1].replace(')','')
                    annotations.append({'start':edit_idx[0] - index_shift, 'end':edit_idx[0] - index_shift + len(edit.split('/')[1].replace(')','')), 'color':'#1E90FF', 'bgcolor':'#e8f3ff', 'underscore':True})
                    index_shift += edit_length - len(edit.split('/')[1].replace(')',''))

                elif '+' in edit:
                    editformat2sequence[edit] = edit.split('+')[1].replace(')','')
                    annotations.append({'start':edit_idx[0] - index_shift, 'end':edit_idx[0] - index_shift + len(edit.split('+')[1].replace(')','')), 'color':'#3CB371', 'bgcolor':'#ebf7f0', 'underscore':True})
                    index_shift += edit_length -len(edit.split('+')[1].replace(')',''))

                elif '-' in edit:
                    editformat2sequence[edit] = ''
                    annotations.append({'start':edit_idx[0] - index_shift, 'end':edit_idx[0] - index_shift, 'color':'#DC143C', 'bgcolor':'#fbe7eb', 'underscore':True})
                    index_shift += edit_length

        for edit in editformat2sequence:
            input_sequence = input_sequence.replace(edit, editformat2sequence[edit])

        try:
            current_annotation_ranges = []
            for annotation in annotations:
                current_annotation_ranges.append([annotation['start'], annotation['end']])

            df_peg = pd.read_json(store_peg_table, orient='split')
            spacer_sequences = list(df_peg.loc[selected_rows, 'spacer sequence'].values)

            # Annotate pegRNA spacer sequences
            for spacer_sequence in spacer_sequences:

                try:
                    start_idx = re.search(spacer_sequence, input_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(spacer_sequence)
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                            annotations.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                            current_annotation_ranges.append([i, i + 1])

                except:
                    start_idx = re.search(reverse_complement(spacer_sequence), input_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(spacer_sequence)
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                            annotations.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                            current_annotation_ranges.append([i, i + 1])

        except:
            pass

    else:
        input_sequence = '...'
        reference_sequence = ''
        edit_sequence = ''

    return(input_sequence, annotations)

@app.callback(Output('pbs-title', 'children'),
    [Input('pbs-range','value')]
)

def update_pbs_title(pbs_range):
    return('PBS length: %s - %s nt' % (pbs_range[0], pbs_range[1]))

@app.callback(Output('rtt-title', 'children'),
    [Input('rtt-range','value')]
)

def update_pbs_title(rtt_range):
    return('RTT length: %s - %s nt' % (rtt_range[0], rtt_range[1]))

@app.callback(Output('nick-dist-title', 'children'),
    [Input('nick-dist-range','value')]
)

def update_pbs_title(nick_dist_range):
    return('Nicking distance: %s - %s bp' % (nick_dist_range[0], nick_dist_range[1]))

### Section to run pegDesigner code

# Helper functions
# IUPAC code map
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

# Extract reference and edited sequence information
def process_sequence(input_sequence):

    input_sequence = ''.join(input_sequence.split())

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


@app.callback([Output('peg-table', 'data'), Output('store-peg-table-total', 'children'), Output('store-peg-table', 'children')],
    [Input('input-check','children'), Input('pbs-range','value'), Input('rtt-range','value'), Input('nick-dist-range','value'), Input('extfirstbase-option','value')],
    state = [State('pe-sequence-input','value')]
)

def run_pegDesigner(input_check, pbs_range, rtt_range, nicking_distance_range, extfirstbase_filter, input_sequence):

    target_design = {}
    peg_design = {'pegRNA group':[],'type':[], 'spacer sequence':[],'PAM':[],'strand':[],'peg-to-edit':[],'nick-to-peg':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'RTT length':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}

    if 'Success' in input_check:

        input_sequence = ''.join(input_sequence.split())
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGG]'
        nicking_distance_minimum = nicking_distance_range[0]
        nicking_distance_maximum = nicking_distance_range[1]
        pbs_length_list = list(range(pbs_range[0], pbs_range[1] + 1))
        rtt_length_list = list(range(rtt_range[0], rtt_range[1] + 1))
        target_sequence = input_sequence#.upper()
        target_name = 'user-input'

        target_sequence = target_sequence.upper()
        editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(target_sequence)

        # Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
        target_design[target_name] = {'target_sequence':target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}

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

        # Remove annotations and convert into regex
        pe_format_rm_annotation = pe_format.replace('/', '').replace('[', '').replace(']', '')

        # Create PE format and PAM search sequences
        pe_format_search_plus = ''
        for base in pe_format_rm_annotation:
            pe_format_search_plus += iupac2bases(base)
        pe_format_search_minus = reverse_complement(pe_format_search_plus)

        pam_search = ''
        pam_sequence = pe_format_rm_annotation[pam_start_idx:pam_end_idx]
        for base in pam_sequence:
            pam_search += iupac2bases(base)

        ##### Initialize data storage for output
        counter = 1
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
                        if spacer_sequence_edit.upper() == spacer_sequence_ref.upper():
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
                        if spacer_sequence_edit.upper() == spacer_sequence_ref.upper():
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

            # Design pegRNAs targeting the (+) strand
            counter = 1
            counted = []
            for peg_plus in target_design[target_name]['pegRNA']['+']:

                pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus
                pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))

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

                                # Construct pegRNA extension to encode intended edit(s)
                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                                # Check to see if pegRNA extension is within input sequence
                                if len(pegRNA_ext) == (pbs_length + rtt_length):

                                    peg_design['pegRNA group'].append(counter)
                                    peg_design['type'].append('pegRNA')
                                    peg_design['spacer sequence'].append(pe_spacer_sequence)
                                    peg_design['PAM'].append(pe_pam_ref)
                                    peg_design['strand'].append('+')
                                    peg_design['peg-to-edit'].append(nick2lastedit_length)
                                    peg_design['nick-to-peg'].append('')
                                    peg_design['pegRNA extension'].append(pegRNA_ext)
                                    peg_design['extension first base'].append(pegRNA_ext[0])
                                    peg_design['PBS length'].append(pbs_length)
                                    peg_design['RTT length'].append(rtt_length)
                                    peg_design['annotation'].append(pe_annotate)
                                    peg_design['spacer top strand oligo'].append('caccG' + pe_spacer_sequence[1:] + 'gtttt')
                                    peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement('G' + pe_spacer_sequence[1:]))
                                    peg_design['pegRNA extension top strand oligo'].append('gtgc' + pegRNA_ext)
                                    peg_design['pegRNA extension bottom strand oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                                    counted.append(counter)

                    # Create ngRNAs targeting (-) strand for (+) pegRNAs
                    if counter in counted:
                        for ng_minus in target_design[target_name]['ngRNA']['-']:
                            ng_nick_ref_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
                            nick_distance = ng_nick_ref_idx - pe_nick_ref_idx
                            if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                                peg_design['pegRNA group'].append(counter)
                                peg_design['type'].append('ngRNA')
                                peg_design['spacer sequence'].append(reverse_complement(ng_spacer_sequence_edit))
                                peg_design['PAM'].append(reverse_complement(ng_pam_edit))
                                peg_design['strand'].append('-')
                                peg_design['peg-to-edit'].append('')
                                peg_design['nick-to-peg'].append(nick_distance)
                                peg_design['pegRNA extension'].append('')
                                peg_design['extension first base'].append('')
                                peg_design['PBS length'].append('')
                                peg_design['RTT length'].append('')
                                peg_design['annotation'].append(ng_annotate)
                                peg_design['spacer top strand oligo'].append('caccG' + reverse_complement(ng_spacer_sequence_edit)[1:])
                                peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement('G' + reverse_complement(ng_spacer_sequence_edit)[1:]))
                                peg_design['pegRNA extension top strand oligo'].append('')
                                peg_design['pegRNA extension bottom strand oligo'].append('')

                        counter += 1

            # Design pegRNAs targeting the (-) strand
            for peg_minus in target_design[target_name]['pegRNA']['-']:

                pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
                pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

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

                                # Construct pegRNA extension to encode intended edit(s)
                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

                                # Check to see if pegRNA extension is within input sequence
                                if len(pegRNA_ext) == (pbs_length + rtt_length):

                                    peg_design['pegRNA group'].append(counter)
                                    peg_design['type'].append('pegRNA')
                                    peg_design['spacer sequence'].append(reverse_complement(pe_spacer_sequence))
                                    peg_design['PAM'].append(reverse_complement(pe_pam_ref))
                                    peg_design['strand'].append('-')
                                    peg_design['peg-to-edit'].append(nick2lastedit_length)
                                    peg_design['nick-to-peg'].append('')
                                    peg_design['pegRNA extension'].append(pegRNA_ext)
                                    peg_design['extension first base'].append(pegRNA_ext[0])
                                    peg_design['PBS length'].append(pbs_length)
                                    peg_design['RTT length'].append(rtt_length)
                                    peg_design['annotation'].append(pe_annotate)
                                    peg_design['spacer top strand oligo'].append('caccG' + reverse_complement(pe_spacer_sequence)[1:] + 'gtttt')
                                    peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement('G' + reverse_complement(pe_spacer_sequence)[1:]))
                                    peg_design['pegRNA extension top strand oligo'].append('gtgc' + pegRNA_ext)
                                    peg_design['pegRNA extension bottom strand oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                                    counted.append(counter)

                    # Create ngRNAs targeting (+) strand for (-) pegRNAs
                    if counter in counted:
                        for ng_plus in target_design[target_name]['ngRNA']['+']:
                            ng_nick_ref_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
                            nick_distance = ng_nick_ref_idx - pe_nick_ref_idx
                            if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                                peg_design['pegRNA group'].append(counter)
                                peg_design['type'].append('ngRNA')
                                peg_design['spacer sequence'].append(ng_spacer_sequence_edit)
                                peg_design['PAM'].append(ng_pam_edit)
                                peg_design['strand'].append('+')
                                peg_design['peg-to-edit'].append('')
                                peg_design['nick-to-peg'].append(nick_distance)
                                peg_design['pegRNA extension'].append('')
                                peg_design['extension first base'].append('')
                                peg_design['PBS length'].append('')
                                peg_design['RTT length'].append('')
                                peg_design['annotation'].append(ng_annotate)
                                peg_design['spacer top strand oligo'].append('caccG' + ng_spacer_sequence_edit[1:])
                                peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement('G' + ng_spacer_sequence_edit[1:]))
                                peg_design['pegRNA extension top strand oligo'].append('')
                                peg_design['pegRNA extension bottom strand oligo'].append('')

                        counter += 1

        df = pd.DataFrame.from_dict(peg_design)

    else:
        df = {'pegRNA group':[],'type':[], 'spacer sequence':[],'PAM':[],'strand':[],'peg-to-edit':[],'nick-to-peg':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'RTT length':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df = pd.DataFrame.from_dict(peg_design)

    if extfirstbase_filter == 'yes':
        df = df[df['extension first base'] != 'C']
        df.reset_index(drop=True, inplace=True)

    df_pegs = df[df['type'] == 'pegRNA']
    df_pegs = df_pegs[['pegRNA group','spacer sequence','PAM','strand','peg-to-edit','annotation']].drop_duplicates()
    df_pegs.reset_index(drop=True, inplace=True)

    return(df_pegs.to_dict('records'), df.to_json(date_format='iso', orient='split'), df_pegs.to_json(date_format='iso', orient='split'))

# Trigger pegRNA extension and ngRNA tables with pegRNA spacer selection
@app.callback(Output('pegext-table', 'data'),
    [Input('peg-table','selected_rows')],
    [State('store-peg-table-total', 'children'), State('store-peg-table', 'children')]
)

def update_pegext_table(selected_row, store_peg_table_total, store_peg_table):

    try:
        # Open up stored peg table
        df_peg = pd.read_json(store_peg_table, orient='split')
        df_peg_total = pd.read_json(store_peg_table_total, orient='split')

        spacer_sequence = list(df_peg.loc[selected_row, 'spacer sequence'].values)
        df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(spacer_sequence)]
        df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
        df_pegext = df_pegext[['PBS length','RTT length','pegRNA extension']].drop_duplicates()
        df_pegext.reset_index(drop=True, inplace=True)

    except:
        df_pegext = {'pegRNA group':[],'type':[], 'spacer sequence':[],'PAM':[],'strand':[],'peg-to-edit':[],'nick-to-peg':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'RTT length':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df_pegext = pd.DataFrame.from_dict(df_pegext)

    return(df_pegext.to_dict('records'))

@app.callback(Output('ng-table', 'data'),
    [Input('peg-table','selected_rows')],
    [State('store-peg-table-total', 'children'), State('store-peg-table', 'children')]
)

def update_ng_table(selected_row, store_peg_table_total, store_peg_table):

    try:
        # Open up stored peg table
        df_peg = pd.read_json(store_peg_table, orient='split')
        df_peg_total = pd.read_json(store_peg_table_total, orient='split')

        peg_group = list(df_peg.loc[selected_row, 'pegRNA group'].values)
        df_ng = df_peg_total[df_peg_total['pegRNA group'].isin(peg_group)]
        df_ng = df_ng[df_ng['type'] == 'ngRNA']
        df_ng = df_ng[['spacer sequence','PAM','strand','nick-to-peg','annotation']].drop_duplicates()
        df_ng.reset_index(drop=True, inplace=True)

    except:
        df_ng = {'pegRNA group':[],'type':[], 'spacer sequence':[],'PAM':[],'strand':[],'peg-to-edit':[],'nick-to-peg':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'RTT length':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df_ng = pd.DataFrame.from_dict(df_ng)

    return(df_ng.to_dict('records'))

# Interface between table selection and sequenceviewer
@app.callback(Output('sequence-option', 'value'),
    [Input('peg-table','selected_rows')]
)

def see_active_cell(active_cell):#, store_peg_table, sequence_option_value):

    return('ref')

@app.callback(Output('download-link', 'href'),
    [Input('store-peg-table-total', 'children')]
)
def update_download_link(store_peg_table_total):

    try:
        # Open up stored peg table
        df_out = pd.read_json(store_peg_table_total, orient='split')

    except:
        df_out = {'pegRNA group':[],'type':[], 'spacer sequence':[],'PAM':[],'strand':[],'peg-to-edit':[],'nick-to-peg':[],'pegRNA extension':[],'PBS length':[],'RTT length':[],'annotation':[]}
        df_out = pd.DataFrame.from_dict(df_out)

    csv_string = df_out.to_csv(index = False, encoding='utf-8')
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_string)
    return csv_string

if __name__ == '__main__':
    app.run_server(debug = True, port = 9994, host = '0.0.0.0')
    # app.run_server(debug=True)