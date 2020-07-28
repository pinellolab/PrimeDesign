# Import required libraries
import re
import os
import io
import base64
import urllib.parse
from zipfile import ZipFile
import uuid
import time
import math
import glob
import dash
import subprocess
import dash_table
import dash_bio as dashbio
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_daq as daq
from flask import Flask, send_from_directory
import pandas as pd
from dash.dependencies import Input, Output, State, ClientsideFunction
import dash_core_components as dcc
import dash_html_components as html

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__)#, external_stylesheets = external_stylesheets)
app.config.suppress_callback_exceptions = True
server = app.server
server.secret_key = '\xfd\x00R\xb5\xbd\x83_t\xed\xdf\xc4\na\x08\xf7K\xc4\xfd\xa2do3\xa5\xdd'

UPLOAD_DIRECTORY = '/PrimeDesign/reports'
if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)

# Primevar mapping
primevar_map = {'clinvar':{'forward':{}, 'reverse':{},},'rs':{'forward':{}, 'reverse':{},}}
try:
    with open('/PrimeDesign/PrimeVar/PrimeVar_mapping.csv', 'r') as f:
        for line in f:

            try:
                variationid, rsid, direction, path = line.strip('\n').split(',')
                primevar_map['clinvar'][direction][variationid] = path
                primevar_map['rs'][direction][rsid] = path

            except:
                pass
except:
    pass

@server.route('/download/<path:path>')
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)

peg_design_tmp = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
df_tmp = pd.DataFrame.from_dict(peg_design_tmp)
 
def serve_layout():

    # session_id = str(uuid.uuid4())
    session_id = str(time.strftime("%Y%m%d_%I.%M.%S.", time.localtime())) + str(int(round(time.time() * 1000)))[-2:]

    # add genome wide example file
    genome_wide_examples = {'InputID': ['example_substitution', 'example_insertion', 'example_deletion',], 
    'InputSequence': ['CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGG(G/T)AGAGACCCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG',
    'CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGGGA(+GTAA)GAGACCCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG',
    'CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGGGAG(-AGAC)CCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG']}
    
    df_example = pd.DataFrame(data=genome_wide_examples)
    df_example.to_csv('/PrimeDesign/reports/PrimeDesign_genome_wide_example_file.csv', index=False)

    # add saturation mutagenesis example file
    saturation_mutagenesis_examples = {'InputID': ['example_saturation_mutagenesis',], 
    'InputSequence': ['CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCA(AATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGGGAGAGACCCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACC)GATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG']}
    
    df_example = pd.DataFrame(data=saturation_mutagenesis_examples)
    df_example.to_csv('/PrimeDesign/reports/PrimeDesign_saturation_mutagenesis_example_file.csv', index=False)

    return html.Div([

        dcc.Location(id='url', refresh=False),
        html.Div(session_id, id='session-id', style={'display': 'none'}),

        html.Div([
            html.Img(src=app.get_asset_url('primedesign_logo.png'), width = '350px', style = {'margin-bottom': '0px', 'margin-right': '20px', 'padding-left': '15px'}),

            dcc.Link('Design', href='/', style = {'color':'#6cb7ff', 'text-decoration':'none', 'margin-right':'25px', 'font-size':'18px'}),
            dcc.Link('PooledDesign', href='/pooled', style = {'color':'#6cb7ff', 'text-decoration':'none', 'margin-right':'25px', 'font-size':'18px'}),
            dcc.Link('PrimeVar', href='/primevar', style = {'color':'#6cb7ff', 'text-decoration':'none', 'margin-right':'25px', 'font-size':'18px'}),
            dcc.Link('About', href='/about', style = {'color':'#6cb7ff', 'text-decoration':'none', 'margin-right':'25px', 'font-size':'18px'}),
            dcc.Link('Help', href='/help', style = {'color':'#6cb7ff', 'text-decoration':'none', 'margin-right':'25px', 'font-size':'18px'}),

            dbc.Button(children = 'Watch a short tutorial!', outline = False, id='open', style = {'color':'#6cb7ff'}),
            dbc.Modal(
                [
                    dbc.ModalHeader("Navigating PrimeDesign"),
                    dbc.ModalBody(html.Img(src=app.get_asset_url('primedesign_design_demo.gif'), style = {'display':'block', 'margin-left':'auto', 'margin-right':'auto', 'width':'80%'}),),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="close", className="ml-auto")
                    ),
                ],
                id="modal",
            ),

            ]),

        html.Div(id='page-content')
    ])

app.layout = serve_layout

about_page = html.Div([

    html.Br(),

    html.H3('What is PrimeDesign?'),
    html.Div([
        '''PrimeDesign is a flexible and comprehensive design tool for prime editing.
        PrimeDesign can be utilized for the installation of 
        '''
        ], style = { 'display':'inline', 'color':'#6a6a6a'}),

    html.Span('substitution', style = {'color':'#1E90FF', 'display':'inline'}),
    html.Span(', ', style = {'display':'inline', 'color':'#6a6a6a'}),
    html.Span('insertion', style = {'color':'#3CB371', 'display':'inline'}),
    html.Span(', and ', style = {'display':'inline', 'color':'#6a6a6a'}),
    html.Span('deletion', style = {'color':'#DC143C', 'display':'inline'}),

    html.Div([
        ''' edits, and is generalizable for both single and combinatorial edits.
        Given an edit of interest, PrimeDesign identifies all possible prime editing guide RNAs (pegRNAs) and nicking guide RNAs (ngRNAs) within a specified parameter range for the optimization of prime editing.
        PrimeDesign also offers the functionality to design genome-wide and saturation mutagenesis libraries (PooledDesign) and provides a comprehensive database of pegRNA and ngRNA designs to install or correct ClinVar pathogenic variants (PrimeVar).
        In addition to the web application, PrimeDesign is also available as a stand-alone command line tool for more flexible and higher-throughput PrimeDesign functions.
        The command line tool is available here: 
        '''
        ], style = {'display':'inline', 'color':'#6a6a6a'}),

    html.A(
            'https://github.com/pinellolab/PrimeDesign',
            id='github-link',
            href="https://github.com/pinellolab/PrimeDesign",
            target='_blank',
            style = {'text-decoration':'none', 'display':'inline'}
        ),

    html.H3('Reference'),

    html.Div([

        '''Our manuscript is available on bioRxiv here: 
        '''

        ], style = {'color':'#6a6a6a', 'display':'inline'}),

    html.A(
            ' https://www.biorxiv.org/content/10.1101/2020.05.04.077750v1',
            id='biorxiv-link',
            href="https://www.biorxiv.org/content/10.1101/2020.05.04.077750v1",
            target='_blank',
            style = {'text-decoration':'none', 'display':'inline'}
        ),

    html.H3('Contact'),

    html.Div([

        '''If you have any questions or concerns, please don't hesitate to contact us at jyhsu (at) mit.edu.
        '''

        ], style = {'color':'#6a6a6a'}),

    html.H3('Labs'),

    html.A(
            'Pinello Lab',
            id='pinellolab-link',
            href="http://pinellolab.org/",
            target='_blank',
            style = {'text-decoration':'none', 'display':'inline', 'font-size':'20px', 'color':'#6a6a6a', 'margin-right':'15px'}
        ),

    html.Label('|', style = {'font-size':'20px', 'display':'inline', 'color':'#6a6a6a'}),

    html.A(
            'Joung Lab',
            id='jounglab-link',
            href="http://www.jounglab.org/",
            target='_blank',
            style = {'text-decoration':'none', 'display':'inline', 'font-size':'20px', 'color':'#6a6a6a', 'margin-right':'15px', 'margin-left':'15px'}
        ),

    html.Label('|', style = {'font-size':'20px', 'display':'inline', 'color':'#6a6a6a'}),

    html.A(
            'Liu Lab',
            id='liulab-link',
            href="https://liugroup.us/",
            target='_blank',
            style = {'text-decoration':'none', 'display':'inline', 'font-size':'20px', 'color':'#6a6a6a', 'margin-left':'15px'}
        ),

    ], style = {'padding': '15px','margin': '0px'}),

help_page = html.Div([

    html.Br(),

    html.H3('What is PrimeDesign?'),
    html.Div([
        '''PrimeDesign is a flexible and comprehensive design tool for prime editing.
        PrimeDesign can be utilized for the installation of 
        '''
        ], style = { 'display':'inline', 'color':'#6a6a6a'}),

    html.Span('substitution', style = {'color':'#1E90FF', 'display':'inline'}),
    html.Span(', ', style = {'display':'inline', 'color':'#6a6a6a'}),
    html.Span('insertion', style = {'color':'#3CB371', 'display':'inline'}),
    html.Span(', and ', style = {'display':'inline', 'color':'#6a6a6a'}),
    html.Span('deletion', style = {'color':'#DC143C', 'display':'inline'}),

    html.Div([
        ''' edits, and is generalizable for both single and combinatorial edits.
        Given an edit of interest, PrimeDesign identifies all possible prime editing guide RNAs (pegRNAs) and nicking guide RNAs (ngRNAs) within a specified parameter range for the optimization of prime editing.
        PrimeDesign also offers the functionality to design genome-wide and saturation mutagenesis libraries (PooledDesign) and provides a comprehensive database of pegRNA and ngRNA designs to install or correct ClinVar pathogenic variants (PrimeVar).
        In addition to the web application, PrimeDesign is also available as a stand-alone command line tool for more flexible and higher-throughput PrimeDesign functions.
        The command line tool is available here: 
        '''
        ], style = {'display':'inline', 'color':'#6a6a6a'}),

    html.A(
            'https://github.com/pinellolab/PrimeDesign',
            id='github-link',
            href="https://github.com/pinellolab/PrimeDesign",
            target='_blank',
            style = {'text-decoration':'none', 'display':'inline'}
        ),

    html.H3('How do you use PrimeDesign?'),

    html.H5('Input sequence'),
    html.Div([
        '''PrimeDesign only requires a single input that encodes both the reference and edit sequences. The edit encoding format is below:
        '''
        ], style = {'color':'#6a6a6a'}),

    html.Div([

        html.Br(),

        ], style = {'display':'block','line-height':'100%'}),

    html.Div([

        html.Span('Substitution', style = {'color':'#1E90FF', 'font-size':'20px'}),
        html.Span(': (reference/edit)', style = {'font-size':'20px'}),
        html.Br(),

        html.Span('Insertion', style = {'color':'#3CB371', 'font-size':'20px'}),
        html.Span(': (+insertion) or (/insertion)', style = {'font-size':'20px'}),
        html.Br(),

        html.Span('Deletion', style = {'color':'#DC143C', 'font-size':'20px'}),
        html.Span(': (-deletion) or (deletion/)', style = {'font-size':'20px'}),


        ], style = {'text-align':'center'}),
    
    html.Div([

        html.Br(),

        ], style = {'display':'block','line-height':'100%'}),

    html.Div([
        '''The input sequence can incorporate single or combinatorial edits by simply including the desired number of edit encodings. For example:
        '''
        ], style = {'color':'#6a6a6a'}),

    html.Div([

        html.Br(),

        ], style = {'display':'block','line-height':'100%'}),

    html.Div([

        html.Div([html.Span(['Reference sequence: '], style = {'font-weight':'bold'})], className = 'two columns'),

        html.Div([html.Span(['... CCTGCTTTCGCTGGGATCCAAGATTGGCAGCTGA'], style = {'font-family':'courier'}),
            html.Span(['A'], style = {'color':'#1E90FF', 'font-family':'courier'}),
            html.Span(['GCCG'], style = {'font-family':'courier'}),
            html.Span(['---'], style = {'color':'#3CB371', 'font-family':'courier'}),
            html.Span(['TTCC'], style = {'font-family':'courier'}),
            html.Span(['ATAG'], style = {'color':'#DC143C', 'font-family':'courier'}),
            html.Span(['TGAGTCCTTCGTCTGTGACTAACTGTGCCAAATCGTCTAGC ...'], style = {'font-family':'courier'}),
            ], className = 'ten columns'),

        ], className = 'row'),

    html.Div([

        html.Div([html.Span(['Edit sequence: '], style = {'font-weight':'bold'}),], className = 'two columns'),

        html.Div([html.Span(['... CCTGCTTTCGCTGGGATCCAAGATTGGCAGCTGA'], style = {'font-family':'courier'}),
            html.Span(['C'], style = {'color':'#1E90FF', 'font-family':'courier'}),
            html.Span(['GCCG'], style = {'font-family':'courier'}),
            html.Span(['CTT'], style = {'color':'#3CB371', 'font-family':'courier'}),
            html.Span(['TTCC'], style = {'font-family':'courier'}),
            html.Span(['----'], style = {'color':'#DC143C', 'font-family':'courier'}),
            html.Span(['TGAGTCCTTCGTCTGTGACTAACTGTGCCAAATCGTCTAGC ...'], style = {'font-family':'courier'}),
            ], className = 'ten columns'),

        ], className = 'row'),

    html.Div([

        html.Div([html.Span(['PrimeDesign input: '], style = {'font-weight':'bold'}),], className = 'two columns'),

        html.Div([html.Span(['... CCTGCTTTCGCTGGGATCCAAGATTGGCAGCTGA'], style = {'font-family':'courier'}),
            html.Span(['(A/C)'], style = {'color':'#1E90FF', 'font-family':'courier'}),
            html.Span(['GCCG'], style = {'font-family':'courier'}),
            html.Span(['(+CTT)'], style = {'color':'#3CB371', 'font-family':'courier'}),
            html.Span(['TTCC'], style = {'font-family':'courier'}),
            html.Span(['(-ATAG)'], style = {'color':'#DC143C', 'font-family':'courier'}),
            html.Span(['TGAGTCCTTCGTCTGTGACTAACTGTGCCAAATCGTCTAGC ...'], style = {'font-family':'courier'}),
            ], className = 'ten columns'),

        ], className = 'row'),

    html.Div([

        html.Br(),

        ], style = {'display':'block','line-height':'100%'}),

    html.Div([

        '''We recommend an input sequence length of >300 bp centered around the the edit(s) of interest to ensure that the complete set of pegRNAs and ngRNAs are designed.
        '''
        ], style = {'color':'#6a6a6a'}),

    html.Div([

        html.Br(),

        ], style = {'display':'block','line-height':'100%'}),

    html.H5('Navigating PrimeDesign'),

    html.Div([

        '''After you successfully construct a desired input sequence, check out the GIF below to learn how to navigate the PrimeDesign web application:
        '''

        ], style = {'color':'#6a6a6a'}),

    html.Div([

        html.Br(),

        ], style = {'display':'block','line-height':'100%'}),

    html.Img(src=app.get_asset_url('primedesign_design_demo.gif'), style = {'display':'block', 'margin-left':'auto', 'margin-right':'auto', 'width':'80%'}),

    ], style = {'padding': '15px','margin': '0px'}),

error_page = html.Div([

    html.Br(),
    html.Br(),
    html.Br(),

    html.H1('404 error: Page not found', style = {'text-align':'center'}),

    ]),

pooled_page = html.Div([

    html.Br(),

    html.Div([

        html.Div([html.H4('Step 1: Set design parameters', style = {'display':'inline', 'margin-right':'5px',}), html.Span('?', id = 'parameters-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'six columns', style = {'padding-top':'15px', 'padding-bottom':'10px'}),
        dbc.Tooltip('Set prime editing parameters to be used for pooled pegRNA and ngRNA design',
               target = 'parameters-tooltip',
               placement = 'right',
               style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
        ),

        html.Div([html.H4('Step 2: Upload file', style = {'display':'inline', 'margin-right':'5px',}),
            html.Span('?', id = 'design-tables-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'}, className = 'dot'),
            html.Span(html.A(children = 'Download genome wide example file', id='download-example-pool', download="PrimeDesign_example_file.csv", href="/download/PrimeDesign_genome_wide_example_file.csv", target="_blank", style = {'font-size':'15px', 'color':'#6cb7ff', 'text-decoration':'underline','float':'right','padding-bottom':'0px','margin-bottom':'0px'}))], className = 'six columns', style = {'padding-top':'15px', 'padding-bottom':'10px'}),

        dbc.Tooltip('Upload .csv file containing input sequences for pooled pegRNA and ngRNA design',
               target = 'design-tables-tooltip',
               placement = 'right',
               style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
         ),

        # html.Div([html.H4('Step 2b: Saturating mutagenesis', style = {'display':'inline', 'margin-right':'5px',}),
        #     html.Span('?', id = 'sm-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'}, className = 'dot')], className = 'four columns', style = {'padding-top':'15px', 'padding-bottom':'10px'}),

        # dbc.Tooltip('Input sequence for saturating mutagenesis. Assumes sequence is in-frame for coding sequences.',
        #        target = 'sm-tooltip',
        #        placement = 'right',
        #        style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
        #  ),

        ], className = 'row', style = {'padding-right': '15px', 'padding-left': '15px','margin': '0px'}),

    html.Div([

        html.Div([

            html.Div([

                html.Div([

                    html.Div([

                        html.Label(id = 'design-title-pool', children = 'Pooled design type', style = {'font-weight':'bold', 'margin-right':'5px'}),
                        html.Span('?',
                              id = 'design-tooltip',
                              style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                              className = 'dot'),

                        dbc.Tooltip('The genome-wide option designs pegRNAs independently per input sequence, while the saturation mutagenesis option first constructs tiling edits across the input sequence and then designs pegRNAs to install these generated edits',
                            target = 'design-tooltip',
                            placement = 'right',
                            style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                        ),

                    ], className='row', style={'display' : 'flex'}),

                    dcc.RadioItems(
                        id = 'design-option-pool',
                        options=[
                            {'label': 'Genome-wide', 'value': 'genome_wide'},
                            {'label': 'Saturation mutagenesis', 'value': 'saturation_mutagenesis'}
                        ],
                        value='genome_wide'
                    ),

                ], className = 'six columns'),

                html.Div(
                    id = 'satmut-type-container',
                    children = [

                        html.Div([

                            html.Label(children = 'Saturation mutagenesis type', style = {'font-weight':'bold', 'margin-right':'5px'}),
                            html.Span('?',
                                  id = 'satmut-design-tooltip',
                                  style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                                  className = 'dot'),

                            dbc.Tooltip('The base option introduces all single base substitutions and the amino acid option introduces all amino acid substitutions across specified portion of the input sequence',
                                target = 'satmut-design-tooltip',
                                placement = 'right',
                                style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                            ),

                        ], className='row', style={'display' : 'flex'}),

                        dcc.RadioItems(
                            id = 'satmut-type',
                            options=[
                                {'label': 'Base substitutions', 'value': 'base'},
                                {'label': 'Amino acid substitutions', 'value': 'aa'}
                            ],
                            value='base'
                        ),

                ], className = 'six columns', style = {'display':'none'}),
            
            ], className = 'row'),

            html.Br(),

            html.Div([

                html.Label(id = 'npegs-title-pool', children = 'Number of pegRNAs per edit', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'npegs-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('The number of pegRNAs to design per edit (if possible), ranked by pegRNA annotation and then nick-to-edit distance',
                       target = 'npegs-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'pbs-info', children = 'Number of pegRNAs to design per edit', style = {'color':'grey'}),
            dcc.Slider(
                id = 'npegs-pool',
                min=1,
                max=10,
                value=3,
            ),

            html.Div([

                html.Label(id = 'homology-downstream-title-pool', children = 'Homology downstream', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'homology-downstream-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Length of 5-30 nt recommended for short edits (<10 nt). Length of 30+ nt recommended for longer edits (>10 nt)',
                       target = 'homology-downstream-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'homology-downstream-info', children = 'Length of homology downstream of an edit for pegRNA designs', style = {'color':'grey'}),
            dcc.Slider(
                id = 'homology-downstream-pool',
                min=1,
                max=50,
                value=10,
            ),

            html.Div([

                html.Label(id = 'pbs-title-pool', children = 'PBS length', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'pbs-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 12-14 nt',
                       target = 'pbs-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'pbs-info', children = 'Primer binding site length to use for pegRNA designs', style = {'color':'grey'}),
            dcc.Slider(
                id = 'pbs-pool',
                min=5,
                max=17,
                value=14,
            ),

            html.Div([

                html.Label(id = 'rtt-title-pool', children = 'Maximum RTT length', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'rtt-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 50+ nt',
                       target = 'rtt-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'rtt-info', children = 'Maximum reverse transcription template length for pegRNA designs', style = {'color':'grey'}),
            dcc.Slider(
                id = 'rtt-pool',
                min=5,
                max=80,
                value=50,
            ),

            html.Div([

                html.Label(id = 'nngs-title-pool', children = 'Number of ngRNAs per pegRNA', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'nngs-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('The number of ngRNAs to design per pegRNA (if possible), ranked by ngRNA annotation and then ngRNA-to-pegRNA distance',
                       target = 'nngs-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'nngs-info', children = 'Number of ngRNAs to design per pegRNA', style = {'color':'grey'}),
            dcc.Slider(
                id = 'nngs-pool',
                min=0,
                max=10,
                value=3,
            ),
            
            html.Div([
                html.Label(id = 'nick-dist-title-pool', children = 'ngRNA distance', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'nick-dist-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 50+ bp (unless PE3b option available)',
                       target = 'nick-dist-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'nick-dist-info', children = 'ngRNA to pegRNA distance', style = {'color':'grey'}),
            dcc.Slider(
                id = 'nick-dist-pool',
                min=0,
                max=120,
                value=75,
            ),

            html.Div([
                html.Label(id = 'remove-first-c-base', children = 'Remove extensions with C first base', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'remove-first-c-base-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('pegRNA extensions that start with a C base may exhibit lower editing efficiencies',
                       target = 'remove-first-c-base-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            dcc.RadioItems(
                id = 'filter-c1-extension-option-pool',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='yes',
                labelStyle={'display': 'inline-block'}
            ),

            html.Div([
                html.Label(id = 'silent-mutation', children = 'Disrupt PAM with silent PAM mutation', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'silent-mutation-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip(children = 'Disrupting the PAM sequence via a silent mutation may improve prime editing efficiencies for coding sequence edits',
                       target = 'silent-mutation-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            dcc.RadioItems(
                id = 'silentmutation-option-pool',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='no',
                labelStyle={'display': 'inline-block'}
            ),

            ], className = 'six columns', style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin':'0px',}),# style={'display':'inline-block','border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin':'0px','float':'left', 'width':'47%'}),

        html.Div([

            dcc.Upload(
                id='upload-data',
                children=html.Div([

                    html.Div([

                        'Drag and Drop or ',
                        html.A('Select File')

                        ], style = {'padding-top':'40px', 'font-size':'25px'}),

                    html.Br(),

                    html.Div([

                        html.Div([

                            'Column names: ',

                            ], style = {'display':'inline','font-size':'15px', 'font-weight':'bold', 'color':'#6a6a6a'}),

                        html.Span('InputID,InputSequence', style = {'font-size':'15px', 'color':'#6a6a6a'}),

                        ]),

                    html.Br(),

                    html.Div(
                        id = 'genome-wide-format-container',
                        children = [

                        html.Div(['InputSequence format:'], style = {'font-size':'15px','font-weight':'bold', 'color':'#6a6a6a'}),
                        html.Span('Substitution', style = {'color':'#1E90FF', 'font-size':'15px'}),
                        html.Span(': (reference/edit)', style = {'font-size':'15px', 'color':'#6a6a6a'}),
                        html.Br(),

                        html.Span('Insertion', style = {'color':'#3CB371', 'font-size':'15px'}),
                        html.Span(': (+insertion) or (/insertion)', style = {'font-size':'15px', 'color':'#6a6a6a'}),
                        html.Br(),

                        html.Span('Deletion', style = {'color':'#DC143C', 'font-size':'15px'}),
                        html.Span(': (-deletion) or (deletion/)', style = {'font-size':'15px', 'color':'#6a6a6a'}),

                        ], style = {'display':'block'}),

                    html.Div(
                        id = 'saturation-mutation-format-container',
                        children = [

                        html.Div(['InputSequence format:'], style = {'font-size':'15px','font-weight':'bold', 'color':'#6a6a6a'}),
                        html.Span('Mutagenesis range', style = {'color':'#1E90FF', 'font-size':'15px'}),
                        html.Span(': (sequence)', style = {'font-size':'15px', 'color':'#6a6a6a'}),
                        html.Br(),

                        html.Span('Base subtitutions', style = {'font-size':'15px','font-weight':'bold', 'color':'#6a6a6a'}),
                        html.Span(': A-to-T,C,G, T-to-A,C,G, etc.', style = {'font-size':'15px', 'color':'#6a6a6a'}),
                        html.Br(),

                        html.Span('Amino acid substitutions', style = {'font-size':'15px','font-weight':'bold', 'color':'#6a6a6a'}),
                        html.Span(': AA-to-all other AAs', style = {'font-size':'15px', 'color':'#6a6a6a'}),

                        ], style = {'display':'none'}),

                    html.Br(),

                    html.Div([

                        html.Div([

                            'File format: ',

                            ], style = {'display':'inline','font-size':'15px', 'font-weight':'bold', 'color':'#6a6a6a'}),

                        html.Span('.csv (comma-separated)', style = {'font-size':'15px', 'color':'#6a6a6a'}),

                        ]),

                ]),
                style={
                    'font-size': '20px',
                    'width': '100%',
                    'height': '400px',
                    # 'lineHeight': '60px',
                    'borderWidth': '2px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '0px',
                    'color':'grey'
                },
            ),

            html.Div(id = 'input-check-pool', children = ['Update: No file uploaded'], style = {'color':'#ff4d4d', 'font-size':'20px', 'margin':'0px', 'padding':'0px'}),

            html.Br(),

            html.H4('Step 3: Download PrimeDesign summary'),

            dcc.Loading(id='loading-7', children=[

                html.H5(id='update-design-pool' , children = 'Design incomplete', style = {'color':'#6a6a6a', 'font-size':'25px'}),

                html.A(
                    children = ' ',
                    id='download-link-pool',
                    download="PrimeDesign_Pooled.csv",
                    href="",
                    target="_blank",
                    style = {'font-size':'25px', 'color':'#6cb7ff', 'text-decoration':'underline'}
                ),

            html.Div(id = 'design-pool-warning', children = []),

            ], type='default'),

            ], className = 'six columns'),

        ], className = 'row', style = {'padding-right': '15px', 'padding-left': '15px','margin': '0px'}),

    # html.Div([

    #     html.H4('Visualize pooled design'),




    #     ])

    ]),

database_page = html.Div([

    html.Div([

        html.Div([

            html.H4('Search PrimeVar'),

            dcc.Dropdown(
                id = 'editing-direction',
                options=[{'label':'Install the pathogenic variant', 'value':'forward'}, {'label':'Correct the pathogenic variant', 'value':'reverse'}],
                value = 'forward',
                style = {'width':'500px', 'min-width': '300px', 'margin-bottom':'5px'}
                ),

            dcc.Dropdown(
                id = 'primevar-id-search-type',
                options=[{'label':'dsbSNP rs#', 'value':'rs'}, {'label':'ClinVar VariationID', 'value':'clinvar'}],
                value = 'rs',
                style = {'width':'500px', 'min-width': '300px'}
                ),

            dcc.Input(
                id = 'primevar-id-search',
                placeholder='Enter variant # here (e.g. 113993960)',
                type='text',
                value='',
                style = {'width':'500px', 'min-width': '300px',}
                # size = '30',
            ),

            # dcc.Dropdown(
            #     id = 'primevar-id-search-type',
            #     options=[{'label':'dsbSNP rs#', 'value':'rs'}, {'label':'ClinVar VariationID', 'value':'clinvar'}],
            #     value = 'rs',
            #     style = {'width':'15px'}
            #     ),

            # dcc.Input(
            #     id = 'primevar-id-search',
            #     placeholder='Enter variant #',
            #     type='text',
            #     value='',
            #     size = '30'
            # ),

            html.Label(id = 'primevar-input-check', children = '', style = {'font-weight':'bold',}),

            ], className = 'six columns'),

        html.Div([

            html.H4('Recommended Designs', style = {'margin-bottom':'0px'}),

            html.Div([

                html.Div([

                    html.H6('pegRNA design', style = {'margin-top':'0px'}),

                    html.Div([

                        html.Label('Coming soon!'),

                        ], style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin': '0px'}),

                    ], className = 'six columns'),

                html.Div([

                    html.H6('ngRNA design', style = {'margin-top':'0px'}),

                    html.Div([

                        html.Label('Coming soon!'),

                        ], style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin': '0px'}),

                    ], className = 'six columns'),

                ], className = 'row', style = {'padding-right': '0px', 'padding-left': '0px','margin': '0px'}),

            ], className = 'six columns',),

        ], className = 'row', style = {'padding-right': '0px','padding-left': '0px','padding-top': '15px','padding-bottom': '15px','margin': '0px'}),

    html.Hr(style = {'margin-top':'5px', 'margin-bottom':'10px'}),

    html.Div([

        html.H5('Visualize sequence', style = {'padding-bottom':'0px'}),

        # dcc.Checklist(
        #         id = 'protein-option-db',
        #         options=[
        #             {'label': 'Visualize amino acid sequence (assumes sequence is in-frame)', 'value': 'protein'},
        #         ],
        #         value=[]
        #     ), 

        html.Div([

            html.H6('Reference DNA', style = {'margin':'0px', 'padding-bottom':'0px'}),
            html.Label('Select pegRNA spacer(s) in design table to visualize', style = {'color':'grey', 'margin-top':'0px'}),
            dashbio.SequenceViewer(
                id = 'reference-sequence-db',
                sequence = ' ',
                badge =False,
                charsPerLine = 90,
                sequenceMaxHeight = '10000px',
                search = False,
                coverage = [],
                # legend = [{'name':'Substitution', 'color':'#1E90FF', 'underscore':False}, {'name':'Insertion', 'color':'#3CB371', 'underscore':False}, {'name':'Deletion', 'color':'#DC143C', 'underscore':False}, {'name':'Selected pegRNA spacer', 'color':'#d6d6d6', 'underscore':False}]
            ),

            # html.Div(id = 'reference-protein-display-db', children = [

            #     html.H6('Reference Protein', style = {'margin':'0px', 'padding-bottom':'0px'}),
            #     dashbio.SequenceViewer(
            #         id = 'reference-protein-sequence-db',
            #         sequence = ' ',
            #         badge =False,
            #         charsPerLine = 80,
            #         sequenceMaxHeight = '10000px',
            #         search = False,
            #         coverage = [],
            #     ),

            #     ], style = {'display':'none'}),

            html.Div(id='store-sequence-db', style={'display': 'none'}),

            html.Span('Substitution', style = {'color':'#1E90FF'}),
            html.Span(' | '),
            html.Span('Deletion', style = {'color':'#DC143C'}),
            html.Span(' | '),
            html.Span('pegRNA spacer', style = {'color':'#3d3d3d', 'text-decoration':'underline'}),
            html.Span(' | '),
            html.Span('ngRNA spacer', style = {'color':'#808080'}),

            html.Hr(),

            html.H6('Edited DNA', style = {'margin':'0px', 'padding-bottom':'0px'}),
            html.Label('Select pegRNA extension(s) and ngRNA(s) in design tables to visualize', style = {'color':'grey', 'margin-top':'0px'}),
            dashbio.SequenceViewer(
                id = 'edit-sequence-db',
                sequence = ' ',
                badge =False,
                charsPerLine = 90,
                sequenceMaxHeight = '10000px',
                search = False,
                coverage = [],
                # legend = [{'name':'Substitution', 'color':'#1E90FF', 'underscore':False}, {'name':'Insertion', 'color':'#3CB371', 'underscore':False}, {'name':'Deletion', 'color':'#DC143C', 'underscore':False}, {'name':'Selected pegRNA spacer', 'color':'#d6d6d6', 'underscore':False}]
            ),

            html.Div(id='store-sequence2-db', style={'display': 'none'}),

            html.Span('Substitution', style = {'color':'#1E90FF'}),
            html.Span(' | '),
            html.Span('Insertion', style = {'color':'#3CB371'}),
            html.Span(' | '),
            html.Span('pegRNA spacer 1-17nt', style = {'color':'#3d3d3d', 'text-decoration':'underline'}),
            html.Span(' | '),
            html.Span('PBS', style = {'color':'#9f7fdf'}),
            html.Span(' | '),
            html.Span('RTT', style = {'color':'#ffa500'}),
            html.Span(' | '),
            html.Span('ngRNA spacer', style = {'color':'#808080'}),

            ], className = 'eight columns', style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin': '0px', 'height':'500px', 'overflow':'auto'}),
            
            html.Div([

                html.Div([

                    html.H6('pegRNA extension secondary structure', style = {'margin':'0px', 'padding-bottom':'0px'}),
                    html.Label('Select a pegRNA spacer and extension to visualize predicted secondary structure', style = {'color':'grey', 'margin-top':'0px'}),

                    html.Div([

                        html.Label(id = 'forna-option-db-left-text', children = 'Extension only', style = {'display':'inline-block','color':'grey', 'margin-right':'10px'}),
                        daq.ToggleSwitch(
                            id = 'forna-option-db',
                            labelPosition='bottom',
                            style = {'display':'inline-block'}
                            ),
                        html.Label(id = 'forna-option-db-right-text', children = 'Full pegRNA', style = {'display':'inline-block','color':'grey', 'margin-left':'10px', 'margin-right':'30px'}),

                        daq.NumericInput(
                            id = 'forna-temp-db',
                            min = -50,
                            max = 100,
                            value=37,
                            size=60,
                            style = {'display':'inline-block'}
                            ),

                        html.Label('°C', style = {'display':'inline-block','color':'grey','margin-left':'5px'}),

                        ], style = {'display':'inline-block', 'text-align': 'center'}),

                    dcc.Loading(
                        id = 'loading-8',
                        type = 'circle',
                        children = [

                            dashbio.FornaContainer(
                                id='forna-pegext-db',
                                # allowPanningAndZooming=False,
                                height=430
                            ),

                        ]
                        ),

                ]),

            ], className = 'four columns', style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px', 'margin': '0px', 'float':'right', 'height':'500px', 'overflow':'auto'}),

        ], className = 'row', style = {'padding-right': '0px', 'padding-left': '0px','margin': '0px'}),
    
    html.Br(),

    html.Div([

        html.Div([html.H5('Prime editing parameters', style = {'display':'inline', 'margin-right':'5px',}), html.Span('?', id = 'parameters-tooltip-db', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'three columns', style = {'padding-top':'15px'}),
        dbc.Tooltip('Interactively design pegRNAs with the parameter slides below',
               target = 'parameters-tooltip-db',
               placement = 'right',
               style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
         ),

        html.Div([html.H5('Design tables', style = {'display':'inline', 'margin-right':'5px',}), html.Span('?', id = 'design-tables-tooltip-db', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'nine columns', style = {'padding-top':'15px'}),
        dbc.Tooltip('Please select pegRNA spacer(s) to proceed with design of pegRNA extensions and ngRNAs',
               target = 'design-tables-tooltip-db',
               placement = 'right',
               style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
         ),

        ], className = 'row', style = {'padding-right': '0px', 'padding-left': '0px','margin': '0px'}),

    html.Div([

        html.Div([

            html.Div([

                html.Label(id = 'pbs-title-db', children = 'PBS length', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'pbs-tooltip-db',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 12-14 nt',
                       target = 'pbs-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'pbs-info-db', children = 'Primer binding site', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'pbs-range-db',
                min=10,
                max=17,
                value=[12, 16],
                allowCross=False
            ),

            html.Div([

                html.Label(id = 'rtt-title-db', children = 'RTT length', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'rtt-tooltip-db',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 10-20 nt',
                       target = 'rtt-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'rtt-info-db', children = 'Reverse transcription template', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'rtt-range-db',
                min=10,
                max=80,
                value=[10, 50],
                allowCross=False
            ),
            
            html.Div([
                html.Label(id = 'nick-dist-title-db', children = 'ngRNA distance', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'nick-dist-tooltip-db',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 50+ bp (unless PE3b option available)',
                       target = 'nick-dist-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'nick-dist-info-db', children = 'ngRNA to pegRNA distance', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'nick-dist-range-db',
                min=0,
                max=120,
                value=[0, 100],
                allowCross=False
            ),

            html.Div([
                html.Label(id = 'remove-first-c-base-db', children = 'Remove extensions with C first base', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'remove-first-c-base-tooltip-db',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('pegRNA extensions that start with a C base may exhibit lower editing efficiencies',
                       target = 'remove-first-c-base-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            dcc.RadioItems(
                id = 'filter-c1-extension-option-db',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='yes',
                labelStyle={'display': 'inline-block'}
            ),

            html.Div([
                html.Label(id = 'silent-mutation-db', children = 'Disrupt PAM with silent PAM mutation', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'silent-mutation-tooltip-db',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip(children = 'Disrupting the PAM sequence via a silent mutation may improve prime editing efficiencies for coding sequence edits',
                       target = 'silent-mutation-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            dcc.RadioItems(
                id = 'silentmutation-option-db',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='no',
                labelStyle={'display': 'inline-block'}
            ),

            ], className = 'three columns', style={'display': 'inline-block','border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin':'0px',}), #'float':'left','width':'25%'

        html.Div([

            html.Div([

                html.Div([html.H6('pegRNA spacers', style = {'display': 'inline', 'margin':'0px', 'margin-right':'5px'}), html.Span('?', id = 'pegspacer-tooltip-db', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'six columns'),

                dbc.Tooltip('Table of all possible pegRNA spacer designs given parameter ranges - Please select pegRNA spacer(s) to proceed with design',
                       target = 'pegspacer-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

                html.Div([

                    html.A(
                        children = 'Download designs',
                        id='download-link-db',
                        download="PrimeDesign_PrimeVar.csv",
                        href="",
                        target="_blank",
                        style = {'font-size':'20px', 'color':'#6cb7ff', 'text-decoration':'none'}
                    ),

                    ], className = 'six columns', style = {'text-align':'right', 'padding-bottom':'0px'}),

                ], className = 'row', style = {'display':'inline', 'margin':'0px'}),

            html.Label('Increase RTT length if no pegRNA spacer designs are available', style = {'color':'grey', 'margin-top':'0px'}),

            dcc.Loading(id='loading-4', children=[

                dash_table.DataTable(
                    id = 'peg-table-db',
                    columns = [{'name': i, 'id': i} for i in ['spacer sequence','PAM','strand','peg-to-edit distance','spacer GC content','annotation']],
                    data = df_tmp.to_dict('records'),
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    # style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'white',
                        # 'fontWeight': 'bold',
                        'font-family':'HelveticaNeue',
                        'font-size':'14px'
                    },
                    style_table={
                        'maxHeight': '300px',
                        'overflowY': 'scroll'
                    },
                    sort_action = 'native',
                    sort_mode = 'multi',
                    # filter_action = 'native',
                    row_selectable = 'single',
                    style_data_conditional=[{
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PAM_disrupted'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    },
                    {
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PAM_disrupted_silent_mutation'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    }]
                ),
            
            ], type='default'),

            html.H6('pegRNA extensions', style = {'display': 'inline', 'margin':'0px', 'margin-right':'5px'}),
            html.Span('?', id = 'pegext-tooltip-db', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot'),

            dbc.Tooltip('Table of all possible pegRNA extensions given parameter ranges - Please select pegRNA spacer(s) to proceed with design',
                       target = 'pegext-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

            html.Label('Please select pegRNA spacer(s) above to see associated extensions', style = {'color':'grey', 'margin-top':'0px'}),

            dcc.Loading(id='loading-5', children=[

                dash_table.DataTable(
                    id = 'pegext-table-db',
                    columns = [{'name': i, 'id': i} for i in ['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']],
                    data = df_tmp.to_dict('records'),
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    # style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'white',
                        # 'fontWeight': 'bold',
                        'font-family':'HelveticaNeue',
                        'font-size':'14px'
                    },
                    style_table={
                        'maxHeight': '300px',
                        'overflowY': 'scroll'
                    },
                    sort_action = 'native',
                    sort_mode = 'multi',
                    # filter_action = 'native',
                    row_selectable = 'single'
                    ),

                ], type='default'),

            html.H6('ngRNA spacers', style = {'display': 'inline', 'margin':'0px', 'margin-right':'5px'}),
            html.Span('?', id = 'ngspacer-tooltip-db', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot'),

            dbc.Tooltip('Table of all possible ngRNAs given parameter ranges - Please select pegRNA spacer(s) to proceed with design',
                       target = 'ngspacer-tooltip-db',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

            html.Label('Please select pegRNA spacer(s) above to see associated ngRNAs', style = {'color':'grey', 'margin-top':'0px'}),

            dcc.Loading(id='loading-6', children=[

                dash_table.DataTable(
                    id = 'ng-table-db',
                    columns = [{'name': i, 'id': i} for i in ['spacer sequence','PAM','strand','nick-to-peg distance','spacer GC content','annotation']],
                    data = df_tmp.to_dict('records'),
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    # style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'white',
                        # 'fontWeight': 'bold',
                        'font-family':'HelveticaNeue','font-size':'14px'

                    },
                    style_table={
                        'maxHeight': '300px',
                        'overflowY': 'scroll'
                    },
                    sort_action = 'native',
                    sort_mode = 'multi',
                    row_selectable = 'single',
                    # filter_action = 'native',
                    style_data_conditional=[{
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PE3b-seed'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    },
                    {
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PE3b-nonseed'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    },
                    ]
                    ),
                ], type='default'),

            html.Div(id='store-peg-table-total-db', style={'display': 'none'}),
            html.Div(id='store-peg-table-db', style={'display': 'none'}),


            ], className = 'nine columns', style={'display': 'inline-block','border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px',}), #'float':'right','width':'70%'

        ], className = 'row'), ####### END OF INSERT

    ], className = 'row', style = {'padding-right': '15px', 'padding-left': '15px','margin': '0px'}),

design_page = html.Div([

    html.Div([

        html.Div([

            html.H4('Input sequence', style = {'margin-right':'5px','display':'inline'}),
            html.Span('?', id = 'input-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot'),

            dcc.Checklist(
                id = 'example-option',
                options=[
                    {'label': 'Substitution', 'value': 'substitution'},
                    {'label': 'Insertion', 'value': 'insertion'},
                    {'label': 'Deletion examples', 'value': 'deletion'},
                ],
                value=[],
                labelStyle={'display': 'inline'}
            ),

            dbc.Tooltip('Edit formatting examples: Substitution (ATGC/CGTA)  |  Insertion (+ATGC) or (/ATGC)  |  Deletion (-ATGC) or (ATGC/)',
                       target = 'input-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

            dcc.Textarea(
                id='pe-sequence-input',
                placeholder='Enter sequence to prime edit or load example input sequence above ...\n\nEdit formatting examples: Substitution (ATGC/CGTA)  |  Insertion (+ATGC) or (/ATGC)  |  Deletion (-ATGC) or (ATGC/)',
                value='',
                style = {'width': '100%', 'margin': '0px', 'white-space':'hard'},
                className = 'textarea',
            ),

            html.Label(id = 'input-check', children = '', style = {'font-weight':'bold'}),
            
            ], className = 'six columns'),

        html.Div([

            html.H4('Recommended Designs', style = {'margin-bottom':'0px','display':'inline'}),

            html.Div([

                html.Div([

                    html.H6('pegRNA design', style = {'margin-top':'0px'}),

                    html.Div([

                        html.Div([

                            html.Div([html.Label('Annotation:', style = {'display':'inline-block', 'font-weight':'bold','padding-right':'5px'}), html.Label(id = 'pegrna-annotation-recommend', style = {'display':'inline-block'}),]),

                            ], className = 'row'),

                        html.Div([

                            html.Div([html.Label('PBS length:', style = {'display':'inline-block', 'font-weight':'bold','padding-right':'5px'}), html.Label(id = 'pegrna-pbs-recommend', style = {'display':'inline-block'})], className = 'six columns'),
                            html.Div([html.Label('RTT length:', style = {'display':'inline-block', 'font-weight':'bold','padding-right':'5px'}), html.Label(id = 'pegrna-rtt-recommend', style = {'display':'inline-block'})], className = 'six columns', style = {'padding-bottom':'0px'}),

                            ], className = 'row'),

                        html.Hr(style = {'margin-top':'7px', 'margin-bottom':'7px',}),

                        html.Label('Spacer oligo top:', style = {'font-weight':'bold'}),
                        html.Label(id = 'pegrna-spacer-recommend-top', style = {'overflow':'auto',}),
                        html.Label('Spacer oligo bottom:', style = {'font-weight':'bold'}),
                        html.Label(id = 'pegrna-spacer-recommend-bottom', style = {'overflow':'auto',}),
                        html.Label('Extension oligo top:', style = {'font-weight':'bold'}),
                        html.Label(id = 'pegrna-ext-recommend-top', style = {'overflow':'auto',}),
                        html.Label('Extension oligo bottom:', style = {'font-weight':'bold'}),
                        html.Label(id = 'pegrna-ext-recommend-bottom', style = {'overflow':'auto',}),


                        ], style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin': '0px', 'height':'260px', 'overflow':'auto'}),

                    ], className = 'six columns'),

                html.Div([

                    html.H6('ngRNA design', style = {'margin-top':'0px'}),

                    html.Div([

                        html.Div([

                            html.Div([html.Label('Annotation: ', style = {'display':'inline-block', 'font-weight':'bold','padding-right':'5px'}), html.Label(id = 'ngrna-annotation-recommend', style = {'display':'inline-block'}),]),

                            ], className = 'row'),

                        html.Div([

                            html.Div([html.Label('Nicking distance: ', style = {'display':'inline-block', 'font-weight':'bold','padding-right':'5px'}), html.Label(id = 'ngrna-distance-recommend', style = {'display':'inline-block'}),]),

                            ], className = 'row'),

                        html.Hr(style = {'margin-top':'7px', 'margin-bottom':'7px',}),

                        html.Label('Spacer oligo top:', style = {'font-weight':'bold'}),
                        html.Label(id = 'ngrna-spacer-recommend-top', style = {'overflow':'auto',}),
                        html.Label('Spacer oligo bottom:', style = {'font-weight':'bold'}),
                        html.Label(id = 'ngrna-spacer-recommend-bottom', style = {'overflow':'auto',}),

                        ], style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin': '0px', 'height':'260px', 'overflow':'auto'}),

                    ], className = 'six columns'),

                ], className = 'row', style = {'padding-right': '0px', 'padding-left': '0px','margin': '0px'}),

            ], className = 'six columns',),

        ], className = 'row', style = {'padding': '15px','margin': '0px'}),

    html.Div([

        html.Hr(style = {'margin-top':'5px', 'margin-bottom':'10px',}),

        ], style = {'padding-left':'15px', 'padding-right':'15px'}),

    html.Div([

        html.H5('Visualize sequence'),

        dcc.Checklist(
                id = 'protein-option',
                options=[
                    {'label': 'Visualize amino acid sequence (assumes sequence is in-frame)', 'value': 'protein'},
                ],
                value=[]
            ), 

        html.Div([

            html.H6('Reference DNA', style = {'margin':'0px', 'padding-bottom':'0px'}),
            html.Label('Select pegRNA spacer(s) in design table to visualize', style = {'color':'grey', 'margin-top':'0px'}),
            dashbio.SequenceViewer(
                id = 'reference-sequence',
                sequence = ' ',
                badge =False,
                charsPerLine = 90,
                sequenceMaxHeight = '10000px',
                search = False,
                coverage = [],
                # legend = [{'name':'Substitution', 'color':'#1E90FF', 'underscore':False}, {'name':'Insertion', 'color':'#3CB371', 'underscore':False}, {'name':'Deletion', 'color':'#DC143C', 'underscore':False}, {'name':'Selected pegRNA spacer', 'color':'#d6d6d6', 'underscore':False}]
            ),

            html.Div(id = 'reference-protein-display', children = [

                html.H6('Reference Protein', style = {'margin':'0px', 'padding-bottom':'0px'}),
                dashbio.SequenceViewer(
                    id = 'reference-protein-sequence',
                    sequence = ' ',
                    badge =False,
                    charsPerLine = 90,
                    sequenceMaxHeight = '10000px',
                    search = False,
                    coverage = [],
                ),

                ], style = {'display':'none'}),

            html.Div(id='store-sequence', style={'display': 'none'}),

            html.Span('Substitution', style = {'color':'#1E90FF'}),
            html.Span(' | '),
            html.Span('Deletion', style = {'color':'#DC143C'}),
            html.Span(' | '),
            html.Span('pegRNA spacer', style = {'color':'#3d3d3d', 'text-decoration':'underline'}),
            html.Span(' | '),
            html.Span('ngRNA spacer', style = {'color':'#808080'}),

            html.Hr(),

            #####
            html.H6('Edited DNA', style = {'margin':'0px', 'padding-bottom':'0px'}),
            html.Label('Select pegRNA extension(s) and ngRNA(s) in design tables to visualize', style = {'color':'grey', 'margin-top':'0px'}),
            dashbio.SequenceViewer(
                id = 'edit-sequence',
                sequence = ' ',
                badge =False,
                charsPerLine = 90,
                sequenceMaxHeight = '10000px',
                search = False,
                coverage = [],
                # legend = [{'name':'Substitution', 'color':'#1E90FF', 'underscore':False}, {'name':'Insertion', 'color':'#3CB371', 'underscore':False}, {'name':'Deletion', 'color':'#DC143C', 'underscore':False}, {'name':'Selected pegRNA spacer', 'color':'#d6d6d6', 'underscore':False}]
            ),

            html.Div(id = 'edit-protein-display', children = [

                html.H6('Edited Protein', style = {'margin':'0px', 'padding-bottom':'0px'}),
                dashbio.SequenceViewer(
                    id = 'edit-protein-sequence',
                    sequence = ' ',
                    badge =False,
                    charsPerLine = 90,
                    sequenceMaxHeight = '10000px',
                    search = False,
                    coverage = [],
                ),

                ], style = {'display':'none'}),

            html.Div(id='store-sequence2', style={'display': 'none'}),

            html.Span('Substitution', style = {'color':'#1E90FF'}),
            html.Span(' | '),
            html.Span('Insertion', style = {'color':'#3CB371'}),
            html.Span(' | '),
            html.Span('pegRNA spacer 1-17nt', style = {'color':'#3d3d3d', 'text-decoration':'underline'}),
            html.Span(' | '),
            html.Span('PBS', style = {'color':'#9f7fdf'}),
            html.Span(' | '),
            html.Span('RTT', style = {'color':'#ffa500'}),
            html.Span(' | '),
            html.Span('ngRNA spacer', style = {'color':'#808080'}),
            ######

            ], className = 'eight columns', style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin': '0px', 'height':'600px', 'overflow':'auto'}),
            
            html.Div([

                html.Div([

                    html.H6('pegRNA secondary structure', style = {'margin':'0px', 'padding-bottom':'0px'}),
                    html.Label('Select a pegRNA spacer and extension to visualize predicted secondary structure', style = {'color':'grey', 'margin-top':'0px'}),

                    html.Div([

                        html.Label(id = 'forna-option-left-text', children = 'Extension only', style = {'display':'inline-block','color':'grey', 'margin-right':'10px'}),
                        daq.ToggleSwitch(
                            id = 'forna-option',
                            labelPosition='bottom',
                            style = {'display':'inline-block'}
                            ),
                        html.Label(id = 'forna-option-right-text', children = 'Full pegRNA', style = {'display':'inline-block','color':'grey', 'margin-left':'10px', 'margin-right':'30px'}),

                        daq.NumericInput(
                            id = 'forna-temp',
                            min = -50,
                            max = 100,
                            value=37,
                            size=60,
                            style = {'display':'inline-block'}
                            ),

                        html.Label('°C', style = {'display':'inline-block','color':'grey','margin-left':'5px'}),

                        ], style = {'display':'inline-block', 'text-align': 'center'}),

                    dcc.Loading(
                        id = 'loading-9',
                        type = 'circle',
                        children = [

                            dashbio.FornaContainer(
                                id='forna-pegext',
                                # allowPanningAndZooming=False,
                                height=430
                            ),

                        ]
                        ),

                ]),

            ], className = 'four columns', style={'border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px', 'margin': '0px', 'float':'right', 'height':'600px', 'overflow':'auto'}),

        ], className = 'row', style = {'padding-right': '15px', 'padding-left': '15px','margin': '0px'}),
    
    html.Br(),

    html.Div([

        html.Div([html.H5('Prime editing parameters', style = {'display':'inline', 'margin-right':'5px',}), html.Span('?', id = 'parameters-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'three columns', style = {'padding-top':'15px'}),
        dbc.Tooltip('Interactively design pegRNAs with the parameter slides below',
               target = 'parameters-tooltip',
               placement = 'right',
               style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
         ),

        html.Div([html.H5('Design tables', style = {'display':'inline', 'margin-right':'5px',}), html.Span('?', id = 'design-tables-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'nine columns', style = {'padding-top':'15px'}),
        dbc.Tooltip('Please select pegRNA spacer(s) to proceed with design of pegRNA extensions and ngRNAs',
               target = 'design-tables-tooltip',
               placement = 'right',
               style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
         ),

        ], className = 'row', style = {'padding-right': '15px', 'padding-left': '15px','margin': '0px'}),

    html.Div([

        html.Div([

            html.Div([

                html.Label(id = 'pbs-title', children = 'PBS length', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'pbs-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 12-14 nt',
                       target = 'pbs-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'pbs-info', children = 'Primer binding site', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'pbs-range',
                min=7,
                max=17,
                value=[12, 14],
                allowCross=False
            ),

            html.Div([

                html.Label(id = 'rtt-title', children = 'RTT length', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'rtt-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 10-20 nt',
                       target = 'rtt-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            html.Label(id = 'rtt-info', children = 'Reverse transcription template', style = {'color':'grey'}),
            dcc.RangeSlider(
                id = 'rtt-range',
                min=10,
                max=80,
                value=[10, 20],
                allowCross=False
            ),
            
            html.Div([
                html.Label(id = 'nick-dist-title', children = 'ngRNA distance', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'nick-dist-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('Initial recommendation: 50+ bp (unless PE3b option available)',
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

            html.Div([
                html.Label(id = 'remove-first-c-base', children = 'Remove extensions with C first base', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'remove-first-c-base-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip('pegRNA extensions that start with a C base may exhibit lower editing efficiencies',
                       target = 'remove-first-c-base-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            dcc.RadioItems(
                id = 'filter-c1-extension-option',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='yes',
                labelStyle={'display': 'inline-block'}
            ),

            html.Div([
                html.Label(id = 'silent-mutation', children = 'Disrupt PAM with silent PAM mutation', style = {'font-weight':'bold', 'margin-right':'5px'}),
                html.Span('?',
                      id = 'silent-mutation-tooltip',
                      style={'font-size':'11px', 'textAlign': 'center', 'color': 'white'},
                      className = 'dot'),

                 dbc.Tooltip(children = 'Disrupting the PAM sequence via a silent mutation may improve prime editing efficiencies for coding sequence edits',
                       target = 'silent-mutation-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                 ),

            ], className='row', style={'display' : 'flex'}),

            dcc.RadioItems(
                id = 'silentmutation-option',
                options=[
                    {'label': 'Yes', 'value': 'yes'},
                    {'label': 'No', 'value': 'no'},
                ],
                value='no',
                labelStyle={'display': 'inline-block'}
            ),

            ], className = 'three columns', style={'display': 'inline-block','border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px','margin':'0px',}), #'float':'left','width':'25%'

        html.Div([

            html.Div([

                html.Div([html.H6('pegRNA spacers', style = {'display': 'inline', 'margin':'0px', 'margin-right':'5px'}), html.Span('?', id = 'pegspacer-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot')], className = 'six columns'),

                dbc.Tooltip('Table of all possible pegRNA spacer designs given parameter ranges - Please select pegRNA spacer(s) to proceed with design',
                       target = 'pegspacer-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

                html.Div([

                    html.A(
                        children = 'Download designs',
                        id='download-link',
                        download="PrimeDesign.csv",
                        href="",
                        target="_blank",
                        style = {'font-size':'20px', 'color':'#6cb7ff', 'text-decoration':'none'}
                    ),

                    ], className = 'six columns', style = {'text-align':'right', 'padding-bottom':'0px'}),

                ], className = 'row', style = {'display':'inline', 'margin':'0px'}),

            html.Label('Increase RTT length if no pegRNA spacer designs are available', style = {'color':'grey', 'margin-top':'0px'}),

            dcc.Loading(id='loading-1', children=[

                dash_table.DataTable(
                    id = 'peg-table',
                    columns = [{'name': i, 'id': i} for i in ['spacer sequence','PAM','strand','peg-to-edit distance','spacer GC content','annotation']],
                    data = df_tmp.to_dict('records'),
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    # style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'white',
                        # 'fontWeight': 'bold',
                        'font-family':'HelveticaNeue',
                        'font-size':'14px'
                    },
                    style_table={
                        'maxHeight': '300px',
                        'overflowY': 'scroll'
                    },
                    sort_action = 'native',
                    sort_mode = 'multi',
                    # filter_action = 'native',
                    row_selectable = 'single',
                    style_data_conditional=[{
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PAM_disrupted'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    },
                    {
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PAM_disrupted_silent_mutation'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    }]
                ),
            ], type='default'),

            html.H6('pegRNA extensions', style = {'display': 'inline', 'margin':'0px', 'margin-right':'5px'}),
            html.Span('?', id = 'pegext-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot'),

            dbc.Tooltip('Table of all possible pegRNA extensions given parameter ranges - Please select pegRNA spacer(s) to proceed with design',
                       target = 'pegext-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

            html.Label('Please select pegRNA spacer(s) above to see associated extensions', style = {'color':'grey', 'margin-top':'0px'}),

            dcc.Loading(id='loading-2', children=[

                dash_table.DataTable(
                    id = 'pegext-table',
                    columns = [{'name': i, 'id': i} for i in ['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']],
                    data = df_tmp.to_dict('records'),
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    # style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'white',
                        # 'fontWeight': 'bold',
                        'font-family':'HelveticaNeue',
                        'font-size':'14px'
                    },
                    style_table={
                        'maxHeight': '300px',
                        'overflowY': 'scroll'
                    },
                    sort_action = 'native',
                    sort_mode = 'multi',
                    # filter_action = 'native',
                    row_selectable = 'single'
                    ),
                ], type='default'),

            html.H6('ngRNA spacers', style = {'display': 'inline', 'margin':'0px', 'margin-right':'5px'}),
            html.Span('?', id = 'ngspacer-tooltip', style={'font-size':'11px', 'textAlign': 'center', 'color': 'white',}, className = 'dot'),

            dbc.Tooltip('Table of all possible ngRNAs given parameter ranges - Please select pegRNA spacer(s) to proceed with design',
                       target = 'ngspacer-tooltip',
                       placement = 'right',
                       style = {'background-color': '#C0C0C0', 'color': '#fff','border-radius': '6px',  'padding': '1px'}
                ),

            html.Label('Please select pegRNA spacer(s) above to see associated ngRNAs', style = {'color':'grey', 'margin-top':'0px'}),

            dcc.Loading(id='loading-3', children=[

                dash_table.DataTable(
                    id = 'ng-table',
                    columns = [{'name': i, 'id': i} for i in ['spacer sequence','PAM','strand','nick-to-peg distance','spacer GC content','annotation']],
                    data = df_tmp.to_dict('records'),
                    style_cell={'textAlign': 'left', 'padding': '5px'},
                    # style_as_list_view=True,
                    style_header={
                        'backgroundColor': 'white',
                        # 'fontWeight': 'bold',
                        'font-family':'HelveticaNeue','font-size':'14px'

                    },
                    style_table={
                        'maxHeight': '300px',
                        'overflowY': 'scroll'
                    },
                    sort_action = 'native',
                    sort_mode = 'multi',
                    row_selectable = 'single',
                    # filter_action = 'native',
                    style_data_conditional=[{
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PE3b-seed'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    },
                    {
                        'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PE3b-nonseed'},
                        'backgroundColor': "#62c096",
                        'color': 'white'
                    },
                    ]
                ),
            ], type='default'),
            
            html.Div(id='store-peg-table-total', style={'display': 'none'}),
            html.Div(id='store-peg-table', style={'display': 'none'}),


            ], className = 'nine columns', style={'display': 'inline-block','border-radius': '5px','box-shadow': '3px 3px 3px lightgrey','background-color': '#fafafa','padding': '15px',}), #'float':'right','width':'70%'

        ], className = 'row', style = {'padding-right': '15px', 'padding-left': '15px','margin': '0px'}), #'margin': '0px'
    
    html.Hr(),

])

# Modal
@app.callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

# Download file
def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    return html.A(filename, href=location)

# Multi page set up
# Update the index
@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):

    if pathname == '/':
        return(design_page)

    elif pathname == '/pooled':
        return(pooled_page)

    elif pathname == '/primevar':
        return(database_page)

    elif pathname == '/about':
        return(about_page)

    elif pathname == '/help':
        return(help_page)

    else:
        return(error_page)

# Load example data
@app.callback(Output('pe-sequence-input','value'),
    [Input('example-option', 'value')]
)

def update_input_check(example_values):
    
    if 'substitution' in example_values:
        if 'insertion' in example_values:
            if 'deletion' in example_values:
                return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGG(G/T)A(+GTAA)G(-AGAC)CCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with substitution, insertion, and deletion edits

            else:
                return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGG(G/T)A(+GTAA)GAGACCCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with substitution and insertion edits

        elif 'deletion' in example_values:
            return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGG(G/T)AG(-AGAC)CCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with substitution and deletion edits

        else:
            return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGG(G/T)AGAGACCCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with substitution edit

    elif 'insertion' in example_values:
        if 'deletion' in example_values:
            return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGGGA(+GTAA)G(-AGAC)CCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with insertion and deletion edits

        else:
            return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGGGA(+GTAA)GAGACCCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with insertion edit

    elif 'deletion' in example_values:
        return('CACACCTACACTGCTCGAAGTAAATATGCGAAGCGCGCGGCCTGGCCGGAGGCGTTCCGCGCCGCCACGTGTTCGTTAACTGTTGATTGGTGGCACATAAGCAATCGTAGTCCGTCAAATTCAGCTCTGTTATCCCGGGCGTTATGTGTCAAATGGCGTAGAACGGGATTGACTGTTTGACGGTAGCTGCTGAGGCGGGAG(-AGAC)CCTCCGTCGGGCTATGTCACTAATACTTTCCAAACGCCCCGTACCGATGCTGAACAAGTCGATGCAGGCTCCCGTCTTTGAAAAGGGGTAAACATACAAGTGGATAGATGATGGGTAGGGGCCTCCAATACATCCAACACTCTACGCCCTCTCCAAGAGCTAGAAGGGCACCCTGCAGTTGGAAAGGG') # Return example input with deletion edit

    else:
        return(None)

@app.callback([Output('input-check', 'children'), Output('input-check', 'style'),],
    [Input('pe-sequence-input','value')]
)

def update_input_check(input_sequence):

    if input_sequence is not None:

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
                sequence_check_style = {'color':'#ff4d4d'}

            else:

                # Check formatting
                if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
                    if '((' not in format_check: # Checks both directions for nested parantheses
                        if '()' not in format_check: # Checks for empty annotations
                            if sum([1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']]) == 0:
                                sequence_check = 'Success: Input sequence has correct formatting'
                                sequence_check_style = {'color':'#6bb6ff'}
                            else:
                                sequence_check = 'Error: Input sequence has more than one edit annotation per parantheses set or annotation outside of parantheses'
                                sequence_check_style = {'color':'#ff4d4d'}
                        else:
                            sequence_check = 'Error: Input sequence has empty parantheses without an edit annotation (i.e. /,  + , -)'
                            sequence_check_style = {'color':'#ff4d4d'}
                    else:
                        sequence_check = 'Error: Input sequence has nested parantheses which is not allowed'
                        sequence_check_style = {'color':'#ff4d4d'}
                else:
                    sequence_check = 'Error: Input sequence does not have full sets of parantheses'
                    sequence_check_style = {'color':'#ff4d4d'}

        else:
            sequence_check = 'Error: Input sequence has exceeded maximum length of 10kb'
            sequence_check_style = {'color':'#ff4d4d'}

    else:
        sequence_check = 'No input sequence with desired edits has been provided'
        sequence_check_style = {'color':'#ff4d4d'}

    return(sequence_check, sequence_check_style)

@app.callback([Output('reference-sequence', 'sequence'), Output('reference-sequence', 'coverage'), Output('reference-protein-sequence', 'sequence'), Output('reference-protein-sequence', 'coverage'), Output('edit-sequence', 'sequence'), Output('edit-sequence', 'coverage'), Output('edit-protein-sequence', 'sequence'), Output('edit-protein-sequence', 'coverage')],
    [Input('input-check','children'), Input('peg-table','selected_rows'), Input('pegext-table','selected_rows'), Input('ng-table','selected_rows')],
    state = [State('pe-sequence-input','value'), State('pbs-range','value'), State('rtt-range','value'), State('nick-dist-range','value'), State('store-peg-table', 'children'), State('store-peg-table-total', 'children')]
)

def update_reference_sequence(input_check, selected_rows_peg, selected_rows_pegext, selected_rows_ng, input_sequence, pbs_range, rtt_range, nicking_distance_range, store_peg_table, store_peg_table_total):

    annotations_ref = []
    annotations_edit = []
    annotations_aa_ref = []
    annotations_aa_edit = []

    if input_sequence is not None:
        input_sequence = ''.join(input_sequence.split())
        reference_sequence = input_sequence
        edit_sequence = input_sequence
        editformat2sequence_ref = {}
        editformat2sequence_edit = {}
        index_shift_ref = 0
        index_shift_edit = 0

        if 'Success' in input_check:

            edit_idxs = [[m.start(), m.end()] for m in re.finditer('\(.*?\)', input_sequence)]
            for edit_idx in edit_idxs:

                edit = input_sequence[edit_idx[0]:edit_idx[1]]
                edit_length = edit_idx[1] - edit_idx[0]

                # Create edit format and number to sequence map
                if '/' in edit:
                    editformat2sequence_ref[edit] = edit.split('/')[0].replace('(','')

                    if len(edit.split('/')[1].replace(')','')) == 0:
                        annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref + len(edit.split('/')[0].replace('(','')), 'color':'#DC143C', 'bgcolor':'#fbe7eb',})# 'underscore':True})
                    else:
                        annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref + len(edit.split('/')[0].replace('(','')), 'color':'#1E90FF', 'bgcolor':'#e8f3ff',})# 'underscore':True})
                    
                    aa_start = math.floor(float(edit_idx[0] - index_shift_ref)/3.0)
                    aa_stop = math.ceil(float(edit_idx[0] - index_shift_ref + len(edit.split('/')[0].replace('(','')))/3.0)

                    if len(edit.split('/')[1].replace(')','')) == 0:
                        annotation_entry = {'start':aa_start, 'end':aa_stop, 'color':'#DC143C', 'bgcolor':'#fbe7eb',}# 'underscore':True}
                    else:
                        annotation_entry = {'start':aa_start, 'end':aa_stop, 'color':'#1E90FF', 'bgcolor':'#e8f3ff',}# 'underscore':True}

                    if annotation_entry not in annotations_aa_ref:
                        annotations_aa_ref.append(annotation_entry)
                    
                    index_shift_ref += edit_length - len(edit.split('/')[0].replace('(',''))

                elif '+' in edit:
                    editformat2sequence_ref[edit] = ''
                    annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref, 'color':'#3CB371', 'bgcolor':'#ebf7f0',})# 'underscore':True})

                    # aa_start = math.floor(float(edit_idx[0] - index_shift)/3.0)
                    # aa_stop = math.ceil(float(edit_idx[0] - index_shift)/3.0)
                    # annotations_aa.append({'start':aa_start, 'end':aa_stop, 'color':'#3CB371', 'bgcolor':'#ebf7f0', 'underscore':True})

                    index_shift_ref += edit_length

                elif '-' in edit:
                    editformat2sequence_ref[edit] = edit.split('-')[1].replace(')','')
                    annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref + len(edit.split('-')[1].replace(')','')), 'color':'#DC143C', 'bgcolor':'#fbe7eb',})# 'underscore':True})

                    aa_start = math.floor(float(edit_idx[0] - index_shift_ref)/3.0)
                    aa_stop = math.ceil(float(edit_idx[0] - index_shift_ref + len(edit.split('-')[1].replace(')','')))/3.0)

                    annotation_entry = {'start':aa_start, 'end':aa_stop, 'color':'#DC143C', 'bgcolor':'#fbe7eb',}# 'underscore':True}
                    if annotation_entry not in annotations_aa_ref:
                        annotations_aa_ref.append(annotation_entry)

                    index_shift_ref += edit_length - len(edit.split('-')[1].replace(')',''))

                # Create edit format and number to sequence map
                if '/' in edit:
                    editformat2sequence_edit[edit] = edit.split('/')[1].replace(')','')

                    if len(edit.split('/')[0].replace('(','')) == 0:
                        annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit + len(edit.split('/')[1].replace(')','')), 'color':'#3CB371', 'bgcolor':'#ebf7f0',})# 'underscore':True})
                    else:
                        annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit + len(edit.split('/')[1].replace(')','')), 'color':'#1E90FF', 'bgcolor':'#e8f3ff',})# 'underscore':True})

                    aa_start = math.floor(float(edit_idx[0] - index_shift_edit)/3.0)
                    aa_stop = math.ceil(float(edit_idx[0] - index_shift_edit + len(edit.split('/')[1].replace(')','')))/3.0)

                    if len(edit.split('/')[0].replace('(','')) == 0:
                        annotation_entry = {'start':aa_start, 'end':aa_stop, 'color':'#3CB371', 'bgcolor':'#ebf7f0',}# 'underscore':True}
                    else:
                        annotation_entry = {'start':aa_start, 'end':aa_stop, 'color':'#1E90FF', 'bgcolor':'#e8f3ff',}# 'underscore':True}

                    if annotation_entry not in annotations_aa_edit:
                        annotations_aa_edit.append(annotation_entry)

                    index_shift_edit += edit_length - len(edit.split('/')[1].replace(')',''))

                elif '+' in edit:
                    editformat2sequence_edit[edit] = edit.split('+')[1].replace(')','')
                    annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit + len(edit.split('+')[1].replace(')','')), 'color':'#3CB371', 'bgcolor':'#ebf7f0',})# 'underscore':True})

                    aa_start = math.floor(float(edit_idx[0] - index_shift_edit)/3.0)
                    aa_stop = math.ceil(float(edit_idx[0] - index_shift_edit + len(edit.split('+')[1].replace(')','')))/3.0)

                    annotation_entry = {'start':aa_start, 'end':aa_stop, 'color':'#3CB371', 'bgcolor':'#ebf7f0',}# 'underscore':True}
                    if annotation_entry not in annotations_aa_edit:
                        annotations_aa_edit.append(annotation_entry)

                    index_shift_edit += edit_length -len(edit.split('+')[1].replace(')',''))

                elif '-' in edit:
                    editformat2sequence_edit[edit] = ''
                    annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit, 'color':'#DC143C', 'bgcolor':'#fbe7eb',})# 'underscore':True})

                    # aa_start = math.floor(float(edit_idx[0] - index_shift)/3.0)
                    # aa_stop = math.ceil(float(edit_idx[0] - index_shift)/3.0)
                    # annotations_aa.append({'start':aa_start, 'end':aa_stop, 'color':'#DC143C', 'bgcolor':'#fbe7eb', 'underscore':True})

                    index_shift_edit += edit_length

            for edit in editformat2sequence_ref:
                reference_sequence = reference_sequence.replace(edit, editformat2sequence_ref[edit])

            for edit in editformat2sequence_edit:
                edit_sequence = edit_sequence.replace(edit, editformat2sequence_edit[edit])

            aa_sequence_ref = sequencetoaa(reference_sequence)
            aa_sequence_edit = sequencetoaa(edit_sequence)

            # print(aa_sequence_ref)
            # print(aa_sequence_edit)

            # Visualizing pegRNA spacer in reference sequence
            try:
                current_annotation_ranges_ref = []
                current_annotation_ranges_edit = []
                for annotation in annotations_ref:
                    current_annotation_ranges_ref.append([annotation['start'], annotation['end']])

                for annotation in annotations_edit:
                    current_annotation_ranges_edit.append([annotation['start'], annotation['end']])

                df_peg = pd.read_json(store_peg_table, orient='split')
                spacer_sequences = list(df_peg.loc[selected_rows_peg, 'spacer sequence'].values)

                # Annotate pegRNA spacer sequences
                for spacer_sequence in spacer_sequences:

                    try:
                        start_idx = re.search(spacer_sequence, reference_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(spacer_sequence)
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                                annotations_ref.append({'start':i, 'end':i + 1, 'underscore':True})
                                current_annotation_ranges_ref.append([i, i + 1])

                        start_idx = re.search(spacer_sequence[:17], edit_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(spacer_sequence) - 3
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                                annotations_edit.append({'start':i, 'end':i + 1, 'underscore':True})
                                current_annotation_ranges_edit.append([i, i + 1])

                    except:
                        start_idx = re.search(reverse_complement(spacer_sequence), reference_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(spacer_sequence)
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                                annotations_ref.append({'start':i, 'end':i + 1, 'underscore':True})
                                current_annotation_ranges_ref.append([i, i + 1])

                        start_idx = re.search(reverse_complement(spacer_sequence[:17]), edit_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(spacer_sequence) - 3
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                                annotations_edit.append({'start':i, 'end':i + 1, 'underscore':True})
                                current_annotation_ranges_edit.append([i, i + 1])

            except:
                pass

            # Visualizing pegRNA extension in edit sequence
            try:
                current_annotation_ranges = []
                for annotation in annotations_edit:
                    current_annotation_ranges.append([annotation['start'], annotation['end']])

                df_peg = pd.read_json(store_peg_table, orient='split')
                df_peg_total = pd.read_json(store_peg_table_total, orient='split')

                peg_group = list(df_peg.loc[selected_rows_peg, 'spacer sequence'].values)
                df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(peg_group)]
                df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
                df_pegext = df_pegext[['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']].drop_duplicates()

                # df_pegext = df_pegext[(df_pegext['PBS length'] >= pbs_range[0]) & (df_pegext['PBS length'] <= pbs_range[1])]
                # df_pegext = df_pegext[(df_pegext['RTT length'] >= rtt_range[0]) & (df_pegext['RTT length'] <= rtt_range[1])]
                df_pegext = df_pegext.sort_values(['PBS length', 'RTT length'])

                df_pegext.reset_index(drop=True, inplace=True)
                pegext_sequences = list(df_pegext.loc[selected_rows_pegext, 'pegRNA extension'].values)
                pbs_lengths = list(df_pegext.loc[selected_rows_pegext, 'PBS length'].values)

                # Annotate pegRNA spacer sequences
                for pegext_sequence, pbs_length in zip(pegext_sequences, pbs_lengths):

                    try:

                        start_idx = re.search(pegext_sequence, edit_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(pegext_sequence)

                        for entry in annotations_edit:
                            start_entry = entry['start']
                            if start_entry in range(stop_idx - pbs_length, stop_idx):
                                entry['bgcolor'] = '#cebeef'

                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                                annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#ffdb99', })
                                current_annotation_ranges.append([i, i + 1])

                    except:

                        start_idx = re.search(reverse_complement(pegext_sequence), edit_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(pegext_sequence)

                        for entry in annotations_edit:
                            start_entry = entry['start']
                            if start_entry in range(start_idx, start_idx + pbs_length):
                                entry['bgcolor'] = '#cebeef'

                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                                annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#ffdb99', })
                                current_annotation_ranges.append([i, i + 1])

            except:
                pass

            # Visualizing ngRNA spacer in edit sequence
            try:
                current_annotation_ranges_ref = []
                current_annotation_ranges_edit = []
                for annotation in annotations_ref:
                    current_annotation_ranges_ref.append([annotation['start'], annotation['end']])

                for annotation in annotations_edit:
                    current_annotation_ranges_edit.append([annotation['start'], annotation['end']])

                df_peg = pd.read_json(store_peg_table, orient='split')
                df_peg_total = pd.read_json(store_peg_table_total, orient='split')

                peg_group = list(df_peg.loc[selected_rows_peg, 'pegRNA group'].values)
                df_ng = df_peg_total[df_peg_total['pegRNA group'].isin(peg_group)]
                df_ng = df_ng[df_ng['type'] == 'ngRNA']
                df_ng = df_ng[['spacer sequence','PAM','strand','nick-to-peg distance','spacer GC content','annotation']].drop_duplicates()

                df_ng = df_ng[(abs(df_ng['nick-to-peg distance']) >= nicking_distance_range[0]) & (abs(df_ng['nick-to-peg distance']) <= nicking_distance_range[1])]
                df_ng = df_ng.sort_values(['nick-to-peg distance'])

                df_ng.reset_index(drop=True, inplace=True)
                ngRNA_sequences = list(df_ng.loc[selected_rows_ng, 'spacer sequence'].values)

                # Annotate pegRNA spacer sequences
                for ngRNA_sequence in ngRNA_sequences:

                    try:

                        start_idx = re.search(ngRNA_sequence, edit_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(ngRNA_sequence)
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                                annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                                current_annotation_ranges_edit.append([i, i + 1])

                        try:
                            start_idx = re.search(ngRNA_sequence, reference_sequence, re.IGNORECASE).start()
                            stop_idx = start_idx + len(ngRNA_sequence)
                            for i in range(start_idx, stop_idx):
                                if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                                    annotations_ref.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                                    current_annotation_ranges_ref.append([i, i + 1])

                        except:
                            pass

                    except:

                        start_idx = re.search(reverse_complement(ngRNA_sequence), edit_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(ngRNA_sequence)
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                                annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                                current_annotation_ranges_edit.append([i, i + 1])

                        try:
                            start_idx = re.search(reverse_complement(ngRNA_sequence), reference_sequence, re.IGNORECASE).start()
                            stop_idx = start_idx + len(ngRNA_sequence)
                            for i in range(start_idx, stop_idx):
                                if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                                    annotations_ref.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                                    current_annotation_ranges_ref.append([i, i + 1])

                        except:
                            pass

            except:
                pass

        else:
            reference_sequence = ' '
            edit_sequence = ' '
            aa_sequence_ref = ' '
            aa_sequence_edit = ' '

    else:
        reference_sequence = ' '
        edit_sequence = ' '
        aa_sequence_ref = ' '
        aa_sequence_edit = ' '

    return(reference_sequence, annotations_ref, aa_sequence_ref, annotations_aa_ref, edit_sequence, annotations_edit, aa_sequence_edit, annotations_aa_edit)

# Visualize protein sequence
@app.callback([Output('reference-protein-display', 'style'), Output('edit-protein-display', 'style')],
    [Input('protein-option','value')]
)

def protein_display(protein_option):

    if 'protein' in protein_option:
        return({'display':'block'}, {'display':'block'})
    else:
        return({'display':'none'}, {'display':'none'})

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


### PooledDesign
@app.callback(Output('npegs-title-pool', 'children'),
    [Input('npegs-pool','value')]
)

def update_npegs_title(n_pegs):
    return('Number of pegRNAs per edit: %s' % (n_pegs))

@app.callback(Output('homology-downstream-title-pool', 'children'),
    [Input('homology-downstream-pool','value')]
)

def update_homology_downstream_title(homology_downstream):
    return('Length of homology downstream: %s nt' % (homology_downstream))

@app.callback(Output('pbs-title-pool', 'children'),
    [Input('pbs-pool','value')]
)

def update_pbs_title(pbs_range):
    return('PBS length: %s nt' % (pbs_range))

@app.callback(Output('rtt-title-pool', 'children'),
    [Input('rtt-pool','value')]
)

def update_rtt_title(rtt_range):
    return('Maximum RTT length: %s nt' % (rtt_range))

@app.callback(Output('nngs-title-pool', 'children'),
    [Input('nngs-pool','value')]
)

def update_ngrnas_title(n_ngRNAs):
    return('Number of ngRNAs per pegRNA: %s' % (n_ngRNAs))

@app.callback(Output('nick-dist-title-pool', 'children'),
    [Input('nick-dist-pool','value')]
)

def update_nickdist_title(nick_dist_range):
    return('Nicking distance: %s bp' % (nick_dist_range))

### Section to run PrimeDesign code

# Helper functions
def gc_content(sequence):
    sequence = sequence.upper()
    GC_count = sequence.count('G') + sequence.count('C')
    GC_content = float(GC_count)/float(len(sequence))

    return("%.2f" % GC_content)

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
        elif base == 'a':
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
    'TGG':['Trp','W', 1],'TGA':['End','*', 0.52],'TGT':['Cys','C', 0.45],'TGC':['Cys','C', 0.55],
    'TAG':['End','*', 0.2],'TAA':['End','*', 0.28],'TAT':['Tyr','Y', 0.43],'TAC':['Tyr','Y', 0.57],
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

for codon in aa2codon:
    aa2codon[codon] = sorted(aa2codon[codon], key = lambda x: x[1], reverse = True)

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

def sequencetoaa(sequence):
    sequence = sequence.upper()
    codon_list = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    if len(codon_list[-1]) != 3:
        codon_list = codon_list[:-1]
    aa_sequence = ''.join(map(str, [codon_dict[x][1] for x in codon_list]))
    return(aa_sequence)

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

# Process sequence for saturating mutagenesis
def saturating_mutagenesis_input_sequences(target_name, target_sequence, sm_type):

    # Check formatting is correct
    format_check = ''
    for i in target_sequence:
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

    # Check for correct formatting of saturating mutagenesis input
    if len(target_sequence) != sum([1 if x in ['A','T','C','G', '(',')'] else 0 for x in target_sequence.upper()]):
        logger.error('Input sequence %s contains a character not in the following list: A,T,C,G,(,) ...' % str(target_sequence))
        sys.exit(1)

    # Check formatting
    if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
        if format_check.count('(') == 1:
            pass
        else:
            logger.error('Input sequence %s has more than one set of parantheses ...' % str(target_sequence))
            sys.exit(1)
    else:
        logger.error('Input sequence %s does not have full sets of parantheses ...' % str(target_sequence))
        sys.exit(1)

    parantheses_start = target_sequence.find('(')
    parantheses_stop = target_sequence.find(')')

    sequence_left = target_sequence[:parantheses_start]
    sequence_right = target_sequence[parantheses_stop + 1:]
    sequence_to_edit = target_sequence[parantheses_start + 1:parantheses_stop]

    sm_target_sequence_list = []
    sm_target_name_list = []

    if sm_type == 'aa':

        for base_index in range(0, len(sequence_to_edit), 3):

            codon_ref = sequence_to_edit[base_index:base_index + 3]

            if len(codon_ref) == 3:

                aa_ref = codon_dict[codon_ref][1]
                inner_sequence_left = sequence_to_edit[:base_index]
                inner_sequence_right = sequence_to_edit[base_index + 3:]

                aa_edit_list = [x for x in aa2codon if x != aa_ref]
                for aa_edit in aa_edit_list:

                    codon_edit = aa2codon[aa_edit][0][0]
                    sm_target_name_list.append('%s_%s_%sto%s' % (target_name, str(int(base_index/3 + 1)), aa_ref, aa_edit))
                    sm_target_sequence_list.append(sequence_left + inner_sequence_left + '(%s/%s)' % (codon_ref, codon_edit) + inner_sequence_right + sequence_right)

    elif sm_type == 'base':

        base_list = ['A','T','C','G']
        for base_index in range(len(sequence_to_edit)):

            base_ref = sequence_to_edit[base_index]
            inner_sequence_left = sequence_to_edit[:base_index]
            inner_sequence_right = sequence_to_edit[base_index + 1:]

            base_edit_list = [x for x in base_list if x != base_ref.upper()]
            for base_edit in base_edit_list:

                sm_target_name_list.append('%s_%s_%sto%s' % (target_name, str(base_index), base_ref, base_edit))
                sm_target_sequence_list.append(sequence_left + inner_sequence_left + '(%s/%s)' % (base_ref, base_edit) + inner_sequence_right + sequence_right)

    return(sm_target_name_list, sm_target_sequence_list)

@app.callback([Output('peg-table', 'data'), Output('store-peg-table-total', 'children'), Output('store-peg-table', 'children'), Output('pegrna-spacer-recommend-top', 'children'), Output('pegrna-spacer-recommend-bottom', 'children'), Output('pegrna-ext-recommend-top', 'children'), Output('pegrna-ext-recommend-bottom', 'children'), Output('pegrna-annotation-recommend', 'children'), Output('pegrna-pbs-recommend', 'children'), Output('pegrna-rtt-recommend', 'children'), Output('ngrna-spacer-recommend-top', 'children'), Output('ngrna-spacer-recommend-bottom', 'children'), Output('ngrna-annotation-recommend', 'children'), Output('ngrna-distance-recommend', 'children')],
    [Input('input-check','children'), Input('pbs-range','value'), Input('rtt-range','value'), Input('nick-dist-range','value'), Input('filter-c1-extension-option','value'), Input('silentmutation-option','value')],
    state = [State('pe-sequence-input','value'), State('session-id', 'children')]
)

def run_primedesign(input_check, pbs_range, rtt_range, nicking_distance_range, filter_c1_extension, silent_mutation, input_sequence, session_id):


    target_design = {}
    peg_design = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[],}

    if 'Success' in input_check:

        input_sequence = ''.join(input_sequence.split())
        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGG]'
        # nicking_distance_minimum = nicking_distance_range[0]
        # nicking_distance_maximum = nicking_distance_range[1]
        pbs_length_list = list(range(pbs_range[0], pbs_range[1] + 1))
        rtt_length_list = list(range(rtt_range[0], rtt_range[1] + 1))

        if 14 not in pbs_length_list:
            pbs_length_list.append(14)
            pbs_length_list = sorted(pbs_length_list)

        if 80 not in rtt_length_list:
            rtt_length_list.append(80)
            rtt_length_list = sorted(rtt_length_list)

        nicking_distance_minimum = 0
        nicking_distance_maximum = 120
        # pbs_length_list = list(range(5, 18))
        # rtt_length_list = list(range(5, 81))

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
                                pe_annotate = 'PAM_disrupted'

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
                                pe_annotate = 'PAM_disrupted'

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
                            if spacer_sequence_edit.upper()[-10:] == spacer_sequence_ref.upper()[-10:]:
                                ng_annotate = 'PE3b-nonseed'
                            else:
                                ng_annotate = 'PE3b-seed'

                        # Store ngRNA spacer
                        nick_ref_idx = re.search(full_search_ref, reference_sequence).end() - (pe_format_length - cut_idx)
                        nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
                        nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
                        target_design[target_name]['ngRNA']['+'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

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
                            if spacer_sequence_edit.upper()[:10] == spacer_sequence_ref.upper()[:10]:
                                ng_annotate = 'PE3b-nonseed'
                            else:
                                ng_annotate = 'PE3b-seed'

                        # Store ngRNA spacer
                        nick_ref_idx = re.search(full_search_ref, reference_sequence).start() + (pe_format_length - cut_idx)
                        nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
                        nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
                        target_design[target_name]['ngRNA']['-'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

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

                pe_annotate_constant = pe_annotate

                # See if pegRNA spacer can introduce all edits
                nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
                if nick2edit_length >= 0:

                    # Loop through RTT lengths
                    silent_mutation_edit = ''
                    for rtt_length in rtt_length_list:

                        # See if RT length can reach entire edit
                        nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
                        if nick2lastedit_length < rtt_length:

                            # Loop through PBS lengths
                            for pbs_length in pbs_length_list:
                                pe_pam_ref_silent_mutation = ''

                                # Construct pegRNA extension to encode intended edit(s)
                                if (pe_nick_edit_idx - rtt_length) < 0:
                                    rtt_length = pe_nick_edit_idx

                                # Patch for NGG PAMs - may need to build something more generalizable in the future
                                if silent_mutation == 'yes':
                                    
                                    if pe_annotate_constant == 'PAM_intact':

                                        nick_aa_index = int(pe_nick_edit_idx)%3
                                        
                                        if nick_aa_index == 0:
                                            original_codon = edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 6].upper()

                                            if len(codon_swap_0[original_codon.upper()]) > 0:
                                                new_codon = codon_swap_0[original_codon][0][0].lower()
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 3] + new_codon + edit_sequence[pe_nick_edit_idx + 6:pe_nick_edit_idx + rtt_length])
                                                pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 3] + new_codon + edit_sequence[pe_nick_edit_idx + 6:]

                                            else:
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                                        elif nick_aa_index == 1:

                                            original_codon_1 = edit_sequence[pe_nick_edit_idx + 2:pe_nick_edit_idx + 5].upper()
                                            original_codon_2 = edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + 8].upper()

                                            if len(codon_swap_1_1[original_codon_1.upper()]) > 0:

                                                new_codon = codon_swap_1_1[original_codon_1][0][0].lower()
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 2] + new_codon + edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + rtt_length])
                                                pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon[1:] + original_codon_2[:1].lower()
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 2] + new_codon + edit_sequence[pe_nick_edit_idx + 5:]

                                            elif len(codon_swap_1_2[original_codon_2.upper()]) > 0:

                                                new_codon = codon_swap_1_2[original_codon_2][0][0].lower()
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:pe_nick_edit_idx + rtt_length])
                                                pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + original_codon_1[1:].lower() + new_codon[:1]
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:]

                                            else:
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                                        elif nick_aa_index == 2:
                                            original_codon = edit_sequence[pe_nick_edit_idx + 4:pe_nick_edit_idx + 7].upper()

                                            if len(codon_swap_2[original_codon.upper()]) > 0:
                                                new_codon = codon_swap_2[original_codon][0][0].lower()
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 4] + new_codon + edit_sequence[pe_nick_edit_idx + 7:pe_nick_edit_idx + rtt_length])
                                                pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 4].lower() + new_codon[:2]
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 4] + new_codon + edit_sequence[pe_nick_edit_idx + 7:]

                                            else:
                                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                                    else:
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                                else:
                                    pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                                # Check to see if pegRNA extension is within input sequence
                                if len(pegRNA_ext) == (pbs_length + rtt_length):

                                    peg_design['pegRNA group'].append(counter)
                                    peg_design['type'].append('pegRNA')
                                    peg_design['spacer sequence'].append(pe_spacer_sequence)
                                    peg_design['spacer GC content'].append(gc_content(pe_spacer_sequence))

                                    if pe_pam_ref_silent_mutation == '':
                                        peg_design['PAM'].append(pe_pam_ref)
                                    else:
                                        peg_design['PAM'].append(pe_pam_ref_silent_mutation)

                                    peg_design['strand'].append('+')
                                    peg_design['peg-to-edit distance'].append(nick2lastedit_length)
                                    peg_design['nick-to-peg distance'].append('')
                                    peg_design['pegRNA extension'].append(pegRNA_ext)
                                    peg_design['extension first base'].append(pegRNA_ext[0])
                                    peg_design['PBS length'].append(pbs_length)
                                    peg_design['PBS GC content'].append(gc_content(pegRNA_ext[rtt_length:]))
                                    peg_design['RTT length'].append(rtt_length)
                                    peg_design['RTT GC content'].append(gc_content(pegRNA_ext[:rtt_length]))
                                    peg_design['annotation'].append(pe_annotate)

                                    if pe_spacer_sequence[0] == 'G':
                                        peg_design['spacer top strand oligo'].append('cacc' + pe_spacer_sequence + 'gtttt')
                                        peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement(pe_spacer_sequence))

                                    else:
                                        peg_design['spacer top strand oligo'].append('caccG' + pe_spacer_sequence + 'gtttt')
                                        peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement('G' + pe_spacer_sequence))

                                    peg_design['pegRNA extension top strand oligo'].append('gtgc' + pegRNA_ext)
                                    peg_design['pegRNA extension bottom strand oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                                    counted.append(counter)

                    # Create ngRNAs targeting (-) strand for (+) pegRNAs
                    if counter in counted:
                        for ng_minus in target_design[target_name]['ngRNA']['-']:
                            ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
                            nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

                            # change ngRNAs that overlap silent mutation
                            if (silent_mutation == 'yes') and (len(silent_mutation_edit) > 0):
                                ng_spacer_sequence_edit = silent_mutation_edit[ng_edit_start_idx:ng_edit_end_idx]

                                mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
                                if len(mutation_indices) > 0:
                                    if len([1 for x in mutation_indices if x >= 10]) > 0:
                                        ng_annotate = 'PE3b-seed'

                                    else:
                                        ng_annotate = 'PE3b-nonseed'
                                else:
                                    ng_annotate = 'PE3'

                            if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                                peg_design['pegRNA group'].append(counter)
                                peg_design['type'].append('ngRNA')
                                peg_design['spacer sequence'].append(reverse_complement(ng_spacer_sequence_edit))
                                peg_design['spacer GC content'].append(gc_content(reverse_complement(ng_spacer_sequence_edit)))
                                peg_design['PAM'].append(reverse_complement(ng_pam_edit))
                                peg_design['strand'].append('-')
                                peg_design['peg-to-edit distance'].append('')
                                peg_design['nick-to-peg distance'].append(nick_distance)
                                peg_design['pegRNA extension'].append('')
                                peg_design['extension first base'].append('')
                                peg_design['PBS length'].append('')
                                peg_design['PBS GC content'].append('')
                                peg_design['RTT length'].append('')
                                peg_design['RTT GC content'].append('')
                                peg_design['annotation'].append(ng_annotate)

                                if reverse_complement(ng_spacer_sequence_edit)[0] == 'G':
                                    peg_design['spacer top strand oligo'].append('cacc' + reverse_complement(ng_spacer_sequence_edit))
                                    peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement(reverse_complement(ng_spacer_sequence_edit)))

                                else:
                                    peg_design['spacer top strand oligo'].append('caccG' + reverse_complement(ng_spacer_sequence_edit))
                                    peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement('G' + reverse_complement(ng_spacer_sequence_edit)))

                                peg_design['pegRNA extension top strand oligo'].append('')
                                peg_design['pegRNA extension bottom strand oligo'].append('')

                        counter += 1

            # Design pegRNAs targeting the (-) strand
            for peg_minus in target_design[target_name]['pegRNA']['-']:

                pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
                pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

                pe_annotate_constant = pe_annotate

                # See if pegRNA spacer can introduce all edits
                nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
                if nick2edit_length >= 0:

                    # Loop through RTT lengths
                    silent_mutation_edit = ''
                    for rtt_length in rtt_length_list:

                        # See if RT length can reach entire edit
                        nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
                        if nick2lastedit_length < rtt_length:

                            # Loop through PBS lengths
                            for pbs_length in pbs_length_list:
                                pe_pam_ref_silent_mutation = ''

                                # Construct pegRNA extension to encode intended edit(s)
                                # pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                                if (pe_nick_edit_idx - rtt_length) < 0:
                                    rtt_length = pe_nick_edit_idx

                                # Patch for NGG PAMs - may need to build something more generalizable in the future
                                if silent_mutation == 'yes':
                                    
                                    if pe_annotate_constant == 'PAM_intact':

                                        nick_aa_index = int(pe_nick_edit_idx)%3
                                        
                                        if nick_aa_index == 0:
                                            original_codon = edit_sequence[pe_nick_edit_idx - 6:pe_nick_edit_idx - 3].upper()

                                            if len(codon_swap_2[original_codon.upper()]) > 0:
                                                new_codon = codon_swap_2[original_codon][0][0].lower()
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 6] + new_codon + edit_sequence[pe_nick_edit_idx - 3:pe_nick_edit_idx + pbs_length]
                                                pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon)
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 6] + new_codon + edit_sequence[pe_nick_edit_idx - 3:]

                                            else:
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

                                        elif nick_aa_index == 1:
                                            original_codon = edit_sequence[pe_nick_edit_idx - 7:pe_nick_edit_idx - 4].upper()

                                            if len(codon_swap_0[original_codon.upper()]) > 0:
                                                new_codon = codon_swap_0[original_codon][0][0].lower()
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 7] + new_codon + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx + pbs_length]
                                                pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon[1:] + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx - 3].lower())
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 7] + new_codon + edit_sequence[pe_nick_edit_idx - 4:]

                                            else:
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

                                        elif nick_aa_index == 2:
                                            original_codon_1 = edit_sequence[pe_nick_edit_idx - 8:pe_nick_edit_idx - 5].upper()
                                            original_codon_2 = edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx - 2].upper()

                                            if len(codon_swap_1_1[original_codon_1.upper()]) > 0:
                                                new_codon = codon_swap_1_1[original_codon_1][0][0].lower()
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 8] + new_codon + edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx + pbs_length]
                                                pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon[2:] + original_codon_2[:2].lower())
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 8] + new_codon + edit_sequence[pe_nick_edit_idx - 5:]

                                            elif len(codon_swap_1_2[original_codon_2.upper()]) > 0:
                                                new_codon = codon_swap_1_2[original_codon_2][0][0].lower()
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 5] + new_codon + edit_sequence[pe_nick_edit_idx - 2:pe_nick_edit_idx + pbs_length]
                                                pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(original_codon_1[2:].lower() + new_codon[:2])
                                                pe_annotate = 'PAM_disrupted_silent_mutation'
                                                silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 5] + new_codon + edit_sequence[pe_nick_edit_idx - 2:]

                                            else:
                                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

                                    else:
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

                                else:
                                    pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length] ########

                                # Check to see if pegRNA extension is within input sequence
                                if len(pegRNA_ext) == (pbs_length + rtt_length):

                                    peg_design['pegRNA group'].append(counter)
                                    peg_design['type'].append('pegRNA')
                                    peg_design['spacer sequence'].append(reverse_complement(pe_spacer_sequence))
                                    peg_design['spacer GC content'].append(gc_content(reverse_complement(pe_spacer_sequence)))

                                    if pe_pam_ref_silent_mutation == '':
                                        peg_design['PAM'].append(reverse_complement(pe_pam_ref))
                                    else:
                                        peg_design['PAM'].append(pe_pam_ref_silent_mutation)

                                    peg_design['strand'].append('-')
                                    peg_design['peg-to-edit distance'].append(nick2lastedit_length)
                                    peg_design['nick-to-peg distance'].append('')
                                    peg_design['pegRNA extension'].append(pegRNA_ext)
                                    peg_design['extension first base'].append(pegRNA_ext[0])
                                    peg_design['PBS length'].append(pbs_length)
                                    peg_design['PBS GC content'].append(gc_content(pegRNA_ext[rtt_length:]))
                                    peg_design['RTT length'].append(rtt_length)
                                    peg_design['RTT GC content'].append(gc_content(pegRNA_ext[:rtt_length]))
                                    peg_design['annotation'].append(pe_annotate)

                                    if reverse_complement(pe_spacer_sequence)[0] == 'G':
                                        peg_design['spacer top strand oligo'].append('cacc' + reverse_complement(pe_spacer_sequence) + 'gtttt')
                                        peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement(reverse_complement(pe_spacer_sequence)))

                                    else:
                                        peg_design['spacer top strand oligo'].append('caccG' + reverse_complement(pe_spacer_sequence) + 'gtttt')
                                        peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement('G' + reverse_complement(pe_spacer_sequence)))

                                    peg_design['pegRNA extension top strand oligo'].append('gtgc' + pegRNA_ext)
                                    peg_design['pegRNA extension bottom strand oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                                    counted.append(counter)

                    # Create ngRNAs targeting (+) strand for (-) pegRNAs
                    if counter in counted:
                        for ng_plus in target_design[target_name]['ngRNA']['+']:
                            ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
                            nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

                            if (silent_mutation == 'yes') and (len(silent_mutation_edit) > 0):
                                ng_spacer_sequence_edit = silent_mutation_edit[ng_edit_start_idx:ng_edit_end_idx]

                                mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
                                if len(mutation_indices) > 0:
                                    if len([1 for x in mutation_indices if x >= 10]) > 0:
                                        ng_annotate = 'PE3b-seed'

                                    else:
                                        ng_annotate = 'PE3b-nonseed'
                                else:
                                    ng_annotate = 'PE3'

                            if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                                peg_design['pegRNA group'].append(counter)
                                peg_design['type'].append('ngRNA')
                                peg_design['spacer sequence'].append(ng_spacer_sequence_edit)
                                peg_design['spacer GC content'].append(gc_content(ng_spacer_sequence_edit))
                                peg_design['PAM'].append(ng_pam_edit)
                                peg_design['strand'].append('+')
                                peg_design['peg-to-edit distance'].append('')
                                peg_design['nick-to-peg distance'].append(nick_distance)
                                peg_design['pegRNA extension'].append('')
                                peg_design['extension first base'].append('')
                                peg_design['PBS length'].append('')
                                peg_design['PBS GC content'].append('')
                                peg_design['RTT length'].append('')
                                peg_design['RTT GC content'].append('')
                                peg_design['annotation'].append(ng_annotate)

                                if ng_spacer_sequence_edit[0] == 'G':
                                    peg_design['spacer top strand oligo'].append('cacc' + ng_spacer_sequence_edit)
                                    peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement(ng_spacer_sequence_edit))

                                else:
                                    peg_design['spacer top strand oligo'].append('caccG' + ng_spacer_sequence_edit)
                                    peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement('G' + ng_spacer_sequence_edit))

                                peg_design['pegRNA extension top strand oligo'].append('')
                                peg_design['pegRNA extension bottom strand oligo'].append('')

                        counter += 1

        df = pd.DataFrame.from_dict(peg_design)

    else:
        df = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df = pd.DataFrame.from_dict(peg_design)

    df_pegs = df[df['type'] == 'pegRNA']

    # Find recommended pegRNA
    if len(df_pegs.sort_values(['annotation', 'peg-to-edit distance'])['pegRNA group']) > 0:

        edit_effective_length = max([edit_span_length_w_ref, edit_span_length_w_edit])
        if edit_effective_length <= 1:
            homology_downstream_recommended = 9
        elif edit_effective_length <= 5:
            homology_downstream_recommended = 14
        elif edit_effective_length <= 10:
            homology_downstream_recommended = 19
        elif edit_effective_length <= 15:
            homology_downstream_recommended = 24
        else:
            homology_downstream_recommended = 34

        pegrna_group = df_pegs.sort_values(['annotation', 'peg-to-edit distance'])['pegRNA group'].values[0]
        rtt_length_optimal = min(df_pegs[df_pegs['pegRNA group'] == pegrna_group]['peg-to-edit distance']) + homology_downstream_recommended
        rtt_max = max(df_pegs[(df_pegs['pegRNA group'] == pegrna_group)]['RTT length'].values)

        # find_peg = 0
        extension_first_base = 'C'
        while (extension_first_base == 'C') and (rtt_length_optimal < rtt_max):
            rtt_length_optimal += 1
            extension_first_base = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == rtt_max)]['pegRNA extension'].values[0][rtt_max-int(rtt_length_optimal):rtt_max][0]

        # if len(df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == (rtt_length_optimal - 1))]) > 0:
        if extension_first_base != 'C':

            pbs_extension_recommended = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == rtt_max)]['pegRNA extension'].values[0][rtt_max:]
            rtt_extension_max = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == rtt_max)]['pegRNA extension'].values[0][:rtt_max]
            extension_recommended = rtt_extension_max[-rtt_length_optimal:] + pbs_extension_recommended

            peg_spacer_top_recommended = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == rtt_max)]['spacer top strand oligo'].values[0]
            peg_spacer_bottom_recommended = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == rtt_max)]['spacer bottom strand oligo'].values[0]
            peg_ext_top_recommended = 'gtgc' + extension_recommended
            peg_ext_bottom_recommended = 'aaaa' + reverse_complement(extension_recommended)

            peg_annotation_recommended = ' %s' % str(df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == 14) & (df_pegs['RTT length'] == rtt_max)]['annotation'].values[0]).replace('_', ' ')
            peg_pbs_recommended = '14 nt'
            peg_rtt_recommended = '%s nt' % str(rtt_length_optimal)

            # Find recommended ngRNA
            df_ngs = df[(df['type'] == 'ngRNA') & (df['pegRNA group'] == pegrna_group)]
            df_ngs['optimal_distance'] = abs(abs(df_ngs['nick-to-peg distance']) - 75)
            df_ngs = df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])

            if len(df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['spacer top strand oligo']) > 0:

                ng_spacer_top_recommended = df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['spacer top strand oligo'].values[0]
                ng_spacer_bottom_recommended = df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['spacer bottom strand oligo'].values[0]
                ng_annotation_recommended = ' %s' % str(df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['annotation'].values[0]).replace('_', ' ')
                ng_distance_recommended = ' %s bp' % str(df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['nick-to-peg distance'].values[0])

            else:

                ng_spacer_top_recommended = 'n/a'
                ng_spacer_bottom_recommended = 'n/a'
                ng_annotation_recommended = 'n/a'
                ng_distance_recommended = 'n/a'

        else:

            peg_spacer_top_recommended = 'n/a'
            peg_spacer_bottom_recommended = 'n/a'
            peg_ext_top_recommended = 'n/a'
            peg_ext_bottom_recommended = 'n/a'

            peg_annotation_recommended = ' n/a'
            peg_pbs_recommended = ' n/a'
            peg_rtt_recommended = ' n/a'

            ng_spacer_top_recommended = 'n/a'
            ng_spacer_bottom_recommended = 'n/a'
            ng_annotation_recommended = ' n/a'
            ng_distance_recommended = ' n/a'

    else:

        peg_spacer_top_recommended = ''
        peg_spacer_bottom_recommended = ''
        peg_ext_top_recommended = ''
        peg_ext_bottom_recommended = ''

        peg_annotation_recommended = ''
        peg_pbs_recommended = ''
        peg_rtt_recommended = ''

        ng_spacer_top_recommended = ''
        ng_spacer_bottom_recommended = ''
        ng_annotation_recommended = ''
        ng_distance_recommended = ''

    # Filter dataframes
    if filter_c1_extension == 'yes':
        df = df[(df['extension first base'] != 'C')]

    pbs_range_list = list(range(pbs_range[0], pbs_range[1] + 1)) + ['']
    rtt_range_list = list(range(rtt_range[0], rtt_range[1] + 1)) + ['']

    pegrna_groups_to_keep = list(set(df_pegs[df_pegs['RTT length'].isin(rtt_range_list)]['pegRNA group'].values))

    df = df[df['PBS length'].isin(pbs_range_list)]
    df_pegs = df_pegs[df_pegs['PBS length'].isin(pbs_range_list)]

    df = df[df['RTT length'].isin(rtt_range_list)]
    df_pegs = df_pegs[df_pegs['RTT length'].isin(rtt_range_list)]

    df = df[df['pegRNA group'].isin(pegrna_groups_to_keep)]
    df_pegs = df_pegs[df_pegs['pegRNA group'].isin(pegrna_groups_to_keep)]

    df.reset_index(drop=True, inplace=True)
    df.to_csv('/PrimeDesign/reports/PrimeDesign_%s.csv' % session_id)

    # # Filter pegrna dataframe
    # df = df[((df['type'] == 'pegRNA') & (df['RTT length'] >= rtt_range[0]) & (df['RTT length'] <= rtt_range[1]) & (df['PBS length'] >= pbs_range[0]) & (df['PBS length'] <= pbs_range[1])) | ((df['type'] == 'ngRNA') & (abs(df['nick-to-peg distance']) >= nicking_distance_range[0]) & (abs(df['nick-to-peg distance']) <= nicking_distance_range[1]))]
    # df_pegs = df_pegs[(df_pegs['RTT length'] >= rtt_range[0]) & (df_pegs['RTT length'] <= rtt_range[1])]

    df_pegs = df_pegs[['pegRNA group','spacer sequence','PAM','strand','peg-to-edit distance','spacer GC content','annotation']].drop_duplicates()
    df_pegs = df_pegs.sort_values('peg-to-edit distance')
    df_pegs.reset_index(drop=True, inplace=True)

    return(df_pegs.to_dict('records'), df.to_json(date_format='iso', orient='split'), df_pegs.to_json(date_format='iso', orient='split'), peg_spacer_top_recommended, peg_spacer_bottom_recommended, peg_ext_top_recommended, peg_ext_bottom_recommended, peg_annotation_recommended, peg_pbs_recommended, peg_rtt_recommended, ng_spacer_top_recommended, ng_spacer_bottom_recommended, ng_annotation_recommended, ng_distance_recommended)

@app.callback(Output('pegext-table', 'data'),
    [Input('peg-table','selected_rows'), Input('pbs-range','value'), Input('rtt-range','value'), Input('filter-c1-extension-option','value'), Input('store-peg-table-total', 'children'), Input('store-peg-table', 'children')]
)

def update_pegext_table(selected_row, pbs_range, rtt_range, filter_c1_extension, store_peg_table_total, store_peg_table):

    if selected_row:

        try:
            # Open up stored peg table
            df_peg = pd.read_json(store_peg_table, orient='split')
            df_peg_total = pd.read_json(store_peg_table_total, orient='split')

            spacer_sequence = list(df_peg.loc[selected_row, 'spacer sequence'].values)
            df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(spacer_sequence)]
            df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
            df_pegext = df_pegext[['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']].drop_duplicates()

            # df_pegext = df_pegext[(df_pegext['PBS length'] >= pbs_range[0]) & (df_pegext['PBS length'] <= pbs_range[1])]
            # df_pegext = df_pegext[(df_pegext['RTT length'] >= rtt_range[0]) & (df_pegext['RTT length'] <= rtt_range[1])]
            df_pegext = df_pegext.sort_values(['PBS length', 'RTT length'])

            df_pegext.reset_index(drop=True, inplace=True)

        except:
            df_pegext = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
            df_pegext = pd.DataFrame.from_dict(df_pegext)

    else:
        df_pegext = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df_pegext = pd.DataFrame.from_dict(df_pegext)

    return(df_pegext.to_dict('records'))

@app.callback(Output('ng-table', 'data'),
    [Input('peg-table','selected_rows'), Input('nick-dist-range', 'value'), Input('store-peg-table-total', 'children'), Input('store-peg-table', 'children')]
)

def update_ng_table(selected_row, nicking_distance_range, store_peg_table_total, store_peg_table):

    if selected_row:

        try:
            # Open up stored peg table
            df_peg = pd.read_json(store_peg_table, orient='split')
            df_peg_total = pd.read_json(store_peg_table_total, orient='split')

            peg_group = list(df_peg.loc[selected_row, 'pegRNA group'].values)
            df_ng = df_peg_total[df_peg_total['pegRNA group'].isin(peg_group)]
            df_ng = df_ng[df_ng['type'] == 'ngRNA']
            df_ng = df_ng[['spacer sequence','PAM','strand','nick-to-peg distance','spacer GC content','annotation']].drop_duplicates()

            df_ng = df_ng[(abs(df_ng['nick-to-peg distance']) >= nicking_distance_range[0]) & (abs(df_ng['nick-to-peg distance']) <= nicking_distance_range[1])]
            df_ng = df_ng.sort_values(['nick-to-peg distance'])

            df_ng.reset_index(drop=True, inplace=True)

        except:
            df_ng = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
            df_ng = pd.DataFrame.from_dict(df_ng)

    else:
        df_ng = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df_ng = pd.DataFrame.from_dict(df_ng)

    return(df_ng.to_dict('records'))

@app.callback(Output('download-link', 'href'),
    [Input('input-check','children')],
    state = [State('session-id', 'children')]
)
def update_download_link(input_check, session_id):
    return('/download/PrimeDesign_%s.csv' % session_id)

##### PooledDesign
@app.callback(Output('download-link-pool', 'href'),
    [Input('input-check-pool','children')],
    state = [State('session-id', 'children')]
)
def update_download_link(input_check, session_id):
    return('/download/PrimeDesign_Pooled_%s.csv' % session_id)

@app.callback([Output('input-check-pool', 'children'), Output('input-check-pool', 'style'),],
    [Input('upload-data','contents')],
    [State('upload-data', 'filename'),
    State('upload-data', 'last_modified'),
    State('design-option-pool', 'value')]
)

def update_input_check(contents, filename, list_of_dates, pooled_design_type):

    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        return('There was an error processing this file.', {'color':'#ff4d4d', 'font-size':'20px'})

    if len(df.index) > 10000:
        sequence_check = 'Error: Input file exceeds the 10k input sequence limit'
        sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
        return(sequence_check, sequence_check_style)

    if len(df.index) == 0:
        sequence_check = 'Error: Input file does not contain any input sequences. Make sure to include column names InputID and InputSequence'
        sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
        return(sequence_check, sequence_check_style)

    for index, row in df.iterrows():

        if len(row) != 2:
            sequence_check = 'Error: Input file does not contain 2 columns (SequencID, InputSequence)'
            sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
            return(sequence_check, sequence_check_style)

        input_sequence = ''.join(row[1])
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

            if pooled_design_type == 'saturation_mutagenesis':

                # Check for correct formatting of saturating mutagenesis input
                if len(input_sequence) != sum([1 if x in ['A','T','C','G', '(',')'] else 0 for x in input_sequence.upper()]):
                    sequence_check = 'Error: Input sequence contains a character not in the following list: A,T,C,G,(,) ...'
                    sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                    return(sequence_check, sequence_check_style)

                else:
                    # Check formatting
                    if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
                        if format_check.count('(') == 1:
                            pass

                        else:
                            sequence_check = 'Error: Input sequence has more than one set of parantheses'
                            sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                            return(sequence_check, sequence_check_style)
                    else:
                        sequence_check = 'Error: Input sequence does not have a full set of parantheses'
                        sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                        return(sequence_check, sequence_check_style)

            else:

                # Check composition of input sequence
                if len(input_sequence) != sum([1 if x in ['A','T','C','G','(',')','+','-','/'] else 0 for x in input_sequence.upper()]):
                    sequence_check = 'Error: Input sequence contains a character not in the following list: A,T,C,G,(,),+,-,/ ...'
                    sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                    return(sequence_check, sequence_check_style)

                else:

                    # Check formatting
                    if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
                        if '((' not in format_check: # Checks both directions for nested parantheses
                            if '()' not in format_check: # Checks for empty annotations
                                if sum([1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']]) == 0:
                                    pass

                                else:
                                    sequence_check = 'Error: Input sequence has more than one edit annotation per parantheses set or annotation outside of parantheses'
                                    sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                                    return(sequence_check, sequence_check_style)
                            else:
                                sequence_check = 'Error: Input sequence has empty parantheses without an edit annotation (i.e. /,  + , -)'
                                sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                                return(sequence_check, sequence_check_style)
                        else:
                            sequence_check = 'Error: Input sequence has nested parantheses which is not allowed'
                            sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                            return(sequence_check, sequence_check_style)
                    else:
                        sequence_check = 'Error: Input sequence does not have full sets of parantheses'
                        sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
                        return(sequence_check, sequence_check_style)

        else:
            sequence_check = 'Error: Input sequence has exceeded maximum length of 10kb'
            sequence_check_style = {'color':'#ff4d4d', 'font-size':'20px'}
            return(sequence_check, sequence_check_style)

        # sequence_check = 'No input sequence with desired edits has been provided'
        # sequence_check_style = {'color':'#ff4d4d'}

    sequence_check = 'Successfully uploaded input file'
    sequence_check_style = {'color':'#6bb6ff', 'font-size':'20px'}

    return(sequence_check, sequence_check_style)

@app.callback([Output('update-design-pool', 'children'), Output('design-pool-warning', 'children'), Output('download-link-pool', 'children')],
    [Input('input-check-pool','children')],
    state = [State('upload-data','contents'), State('upload-data', 'filename'), State('design-option-pool','value'), State('satmut-type','value'), State('npegs-pool','value'), State('homology-downstream-pool','value'), State('pbs-pool','value'), State('rtt-pool','value'), State('nngs-pool','value'), State('nick-dist-pool','value'), State('filter-c1-extension-option-pool','value'), State('silentmutation-option-pool','value'), State('session-id', 'children')]
)

def run_primedesign_pooled(input_check, contents, filename, pool_type, satmut_type, number_of_pegrnas, homology_downstream, pbs_length_pooled, rtt_max_length_pooled, number_of_ngrnas, nicking_distance_pooled, filter_c1_extension, silent_mutation, session_id):

    target_design = {}
    peg_count_dict = {}
    warning_list = []
    peg_design = {'Target_name':[], 'Target_sequence':[], 'pegRNA_number':[],'gRNA_type':[], 'Spacer_sequence':[],'Spacer_GC_content':[],'PAM_sequence':[],'Strand':[],'pegRNA-to-edit_distance':[],'ngRNA-to-pegRNA_distance':[],'pegRNA_extension':[], 'First_extension_nucleotide':[],'PBS_length':[],'PBS_GC_content':[],'RTT_length':[],'RTT_GC_content':[],'Annotation':[],'Spacer_top_strand_oligo':[], 'Spacer_bottom_strand_oligo':[], 'pegRNA_extension_top_strand_oligo':[], 'pegRNA_extension_bottom_strand_oligo':[]}
    if 'Success' in input_check:

        content_type, content_string = contents.split(',')

        decoded = base64.b64decode(content_string)
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df_in = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df_in = pd.read_excel(io.BytesIO(decoded))

        pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGG]'
        nicking_distance_minimum = 0
        nicking_distance_maximum = nicking_distance_pooled + 100
        # pbs_length_list = list(range(pbs_range[0], pbs_range[1] + 1))
        # rtt_length_list = list(range(rtt_range[0], rtt_range[1] + 1))

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

        for index, row in df_in.iterrows():

            target_sequence = str(''.join(row[1])).upper()
            target_name = str(row[0])

            if pool_type == 'saturation_mutagenesis':

                sm_target_name_list, sm_target_sequence_list = saturating_mutagenesis_input_sequences(target_name, target_sequence, satmut_type)

                for sm_target_name, sm_target_sequence in zip(sm_target_name_list, sm_target_sequence_list):

                    editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(sm_target_sequence)

                    # Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
                    target_design[sm_target_name] = {'target_sequence':sm_target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}
                    peg_count_dict[sm_target_name] = 0

            else:

                editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(target_sequence)

                # Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
                target_design[target_name] = {'target_sequence':target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}
                peg_count_dict[target_name] = 0

        ##### Initialize data storage for output
        pe_design = {}
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
            edit_sequence = target_design[target_name]['edit_sequence']

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
                                pe_annotate = 'PAM_disrupted'

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
                                pe_annotate = 'PAM_disrupted'

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
                            if spacer_sequence_edit.upper()[-10:] == spacer_sequence_ref.upper()[-10:]:
                                ng_annotate = 'PE3b-nonseed'
                            else:
                                ng_annotate = 'PE3b-seed'

                        # Store ngRNA spacer
                        nick_ref_idx = re.search(full_search_ref, reference_sequence).end() - (pe_format_length - cut_idx)
                        nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
                        nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
                        target_design[target_name]['ngRNA']['+'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

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
                            if spacer_sequence_edit.upper()[:10] == spacer_sequence_ref.upper()[:10]:
                                ng_annotate = 'PE3b-nonseed'
                            else:
                                ng_annotate = 'PE3b-seed'

                        # Store ngRNA spacer
                        nick_ref_idx = re.search(full_search_ref, reference_sequence).start() + (pe_format_length - cut_idx)
                        nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
                        nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
                        target_design[target_name]['ngRNA']['-'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

            # Grab index information of edits to introduce to target sequence
            edit_start_in_ref = int(target_design[target_name]['edit_start_in_ref'])
            edit_stop_in_ref_rev = int(target_design[target_name]['edit_stop_in_ref_rev'])
            edit_span_length_w_ref = int(target_design[target_name]['edit_span_length'][0])
            edit_span_length_w_edit = int(target_design[target_name]['edit_span_length'][1])

            # Initialize pegRNA and ngRNA design dictionary
            pe_design[target_name] = {}

            # # Design for genome-wide or saturation mutagenesis screening applications
            # if genome_wide_design or saturation_mutagenesis:

            # Design pegRNAs targeting the (+) strand
            for peg_plus in target_design[target_name]['pegRNA']['+']:

                pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus
                # pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))

                pe_annotate_constant = pe_annotate

                # See if pegRNA spacer can introduce all edits
                nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
                if nick2edit_length >= 0:

                    # See if RTT length can reach entire edit with homology downstream constraint
                    silent_mutation_edit = ''
                    nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
                    rtt_length = nick2lastedit_length + homology_downstream
                    if rtt_length < rtt_max_length_pooled:

                        pbs_length = pbs_length_pooled
                        pe_pam_ref_silent_mutation = ''

                        # Construct pegRNA extension to encode intended edit(s)

                        # Patch for NGG PAMs - may need to build something more generalizable in the future
                        if silent_mutation == 'yes':
                            
                            if pe_annotate_constant == 'PAM_intact':

                                nick_aa_index = int(pe_nick_edit_idx)%3
                                
                                if nick_aa_index == 0:
                                    original_codon = edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 6].upper()

                                    if len(codon_swap_0[original_codon.upper()]) > 0:
                                        new_codon = codon_swap_0[original_codon][0][0].lower()
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 3] + new_codon + edit_sequence[pe_nick_edit_idx + 6:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 3] + new_codon + edit_sequence[pe_nick_edit_idx + 6:pe_nick_edit_idx + rtt_max_length_pooled])
                                        pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 3] + new_codon + edit_sequence[pe_nick_edit_idx + 6:]

                                    else:
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])

                                elif nick_aa_index == 1:
                                    original_codon_1 = edit_sequence[pe_nick_edit_idx + 2:pe_nick_edit_idx + 5].upper()
                                    original_codon_2 = edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + 8].upper()

                                    if len(codon_swap_1_1[original_codon_1.upper()]) > 0:
                                        new_codon = codon_swap_1_1[original_codon_1][0][0].lower()
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 2] + new_codon + edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 2] + new_codon + edit_sequence[pe_nick_edit_idx + 5:pe_nick_edit_idx + rtt_max_length_pooled])
                                        pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + new_codon[1:] + original_codon_2[:1].lower()
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 2] + new_codon + edit_sequence[pe_nick_edit_idx + 5:]

                                    elif len(codon_swap_1_2[original_codon_2.upper()]) > 0:
                                        new_codon = codon_swap_1_2[original_codon_2][0][0].lower()
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:pe_nick_edit_idx + rtt_max_length_pooled])
                                        pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + original_codon_1[1:].lower() + new_codon[:1]
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:]

                                    else:
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])

                                elif nick_aa_index == 2:
                                    original_codon = edit_sequence[pe_nick_edit_idx + 4:pe_nick_edit_idx + 7].upper()

                                    if len(codon_swap_2[original_codon.upper()]) > 0:
                                        new_codon = codon_swap_2[original_codon][0][0].lower()
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 4] + new_codon + edit_sequence[pe_nick_edit_idx + 7:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 4] + new_codon + edit_sequence[pe_nick_edit_idx + 7:pe_nick_edit_idx + rtt_max_length_pooled])
                                        pe_pam_ref_silent_mutation = pe_pam_ref + '-to-' + edit_sequence[pe_nick_edit_idx + 3:pe_nick_edit_idx + 4].lower() + new_codon[:2]
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx + 4] + new_codon + edit_sequence[pe_nick_edit_idx + 7:]

                                    else:
                                        pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])

                            else:
                                pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
                                pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])

                        else:
                            pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])
                            pegRNA_ext_max = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_max_length_pooled])

                        # Check to see if pegRNA extension is within input sequence
                        if 'PAM_disrupted' in pe_annotate:
                            pe_annotate_code = 0
                        else:
                            pe_annotate_code = 1

                        pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))
                        if len(pegRNA_ext) == (pbs_length + rtt_length):

                            # Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
                            if pegid not in pe_design[target_name]:

                                # First list is for peg extension, second list is for nicking guide
                                pe_design[target_name][pegid] = [[],[]]

                            if pe_pam_ref_silent_mutation == '':
                                pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length])

                            else:
                                pe_design[target_name][pegid][0].append([pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref_silent_mutation, pe_annotate, '+', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length])

                    # Create pegID
                    if 'PAM_disrupted' in pe_annotate:
                        pe_annotate_code = 0
                    else:
                        pe_annotate_code = 1

                    pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))
                    # Create ngRNAs targeting (-) strand for (+) pegRNAs
                    if pegid in pe_design[target_name]:
                    # if counter in counted:
                        for ng_minus in target_design[target_name]['ngRNA']['-']:
                            ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
                            nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

                            if (silent_mutation == 'yes') and (len(silent_mutation_edit) > 0):
                                ng_spacer_sequence_edit = silent_mutation_edit[ng_edit_start_idx:ng_edit_end_idx]

                                mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
                                if len(mutation_indices) > 0:
                                    if len([1 for x in mutation_indices if x >= 10]) > 0:
                                        ng_annotate = 'PE3b-seed'

                                    else:
                                        ng_annotate = 'PE3b-nonseed'
                                else:
                                    ng_annotate = 'PE3'

                            if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                                if ng_annotate == 'PE3b-seed':
                                    ng_code = 0
                                elif ng_annotate == 'PE3b-nonseed':
                                    ng_code = 1
                                else:
                                    ng_code = 2

                                pe_design[target_name][pegid][1].append([str(ng_code) + '0'*(3 - len(str(abs(abs(nick_distance) - nicking_distance_pooled)))) + str(abs(abs(nick_distance) - nicking_distance_pooled)), ng_nick_ref_idx, reverse_complement(ng_spacer_sequence_edit), reverse_complement(ng_pam_edit), ng_annotate, '-', nick_distance])

                        pe_design[target_name][pegid][1] = sorted(pe_design[target_name][pegid][1])

            # Design pegRNAs targeting the (-) strand
            for peg_minus in target_design[target_name]['pegRNA']['-']:

                pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
                # pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

                pe_annotate_constant = pe_annotate

                # See if pegRNA spacer can introduce all edits
                nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
                if nick2edit_length >= 0:

                    # See if RT length can reach entire edit
                    silent_mutation_edit = ''
                    nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
                    rtt_length = nick2lastedit_length + homology_downstream
                    if rtt_length < rtt_max_length_pooled:

                        pbs_length = pbs_length_pooled
                        pe_pam_ref_silent_mutation = ''

                        # Construct pegRNA extension to encode intended edit(s)

                        # Patch for NGG PAMs - may need to build something more generalizable in the future
                        if silent_mutation == 'yes':
                            
                            if pe_annotate_constant == 'PAM_intact':

                                nick_aa_index = int(pe_nick_edit_idx)%3
                                
                                if nick_aa_index == 0:
                                    original_codon = edit_sequence[pe_nick_edit_idx - 6:pe_nick_edit_idx - 3].upper()

                                    if len(codon_swap_2[original_codon.upper()]) > 0:
                                        new_codon = codon_swap_2[original_codon][0][0].lower()
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 6] + new_codon + edit_sequence[pe_nick_edit_idx - 3:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx - 6] + new_codon + edit_sequence[pe_nick_edit_idx - 3:pe_nick_edit_idx + pbs_length]
                                        pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon)
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 6] + new_codon + edit_sequence[pe_nick_edit_idx - 3:]

                                    else:
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]

                                elif nick_aa_index == 1:
                                    original_codon = edit_sequence[pe_nick_edit_idx - 7:pe_nick_edit_idx - 4].upper()

                                    if len(codon_swap_0[original_codon.upper()]) > 0:
                                        new_codon = codon_swap_0[original_codon][0][0].lower()
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 7] + new_codon + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx - 7] + new_codon + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx + pbs_length]
                                        pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon[1:] + edit_sequence[pe_nick_edit_idx - 4:pe_nick_edit_idx - 3].lower())
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 7] + new_codon + edit_sequence[pe_nick_edit_idx - 4:]

                                    else:
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]

                                elif nick_aa_index == 2:
                                    original_codon_1 = edit_sequence[pe_nick_edit_idx - 8:pe_nick_edit_idx - 5].upper()
                                    original_codon_2 = edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx - 2].upper()

                                    if len(codon_swap_1_1[original_codon_1.upper()]) > 0:
                                        new_codon = codon_swap_1_1[original_codon_1][0][0].lower()
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 8] + new_codon + edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx - 8] + new_codon + edit_sequence[pe_nick_edit_idx - 5:pe_nick_edit_idx + pbs_length]
                                        pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(new_codon[2:] + original_codon_2[:2].lower())
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 8] + new_codon + edit_sequence[pe_nick_edit_idx - 5:]

                                    elif len(codon_swap_1_2[original_codon_2.upper()]) > 0:
                                        new_codon = codon_swap_1_2[original_codon_2][0][0].lower()
                                        # pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + 5] + new_codon + edit_sequence[pe_nick_edit_idx + 8:pe_nick_edit_idx + rtt_length])
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx - 5] + new_codon + edit_sequence[pe_nick_edit_idx - 2:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx - 5] + new_codon + edit_sequence[pe_nick_edit_idx - 2:pe_nick_edit_idx + pbs_length]
                                        pe_pam_ref_silent_mutation = reverse_complement(pe_pam_ref) + '-to-' + reverse_complement(original_codon_1[2:].lower() + new_codon[:2])
                                        pe_annotate = 'PAM_disrupted_silent_mutation'
                                        silent_mutation_edit = edit_sequence[:pe_nick_edit_idx - 5] + new_codon + edit_sequence[pe_nick_edit_idx - 2:]

                                    else:
                                        pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                                        pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]

                            else:
                                pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                                pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]

                        else:
                            pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                            pegRNA_ext_max = edit_sequence[pe_nick_edit_idx - rtt_max_length_pooled:pe_nick_edit_idx + pbs_length]

                        # Check to see if pegRNA extension is within input sequence
                        if 'PAM_disrupted' in pe_annotate:
                            pe_annotate_code = 0
                        else:
                            pe_annotate_code = 1

                        pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))
                        # Check to see if pegRNA extension is within input sequence
                        if len(pegRNA_ext) == (pbs_length + rtt_length):

                            # Initiate entry for new pegRNA spacers that are close enough to edit window based on RTT length parameter list
                            if pegid not in pe_design[target_name]:

                                # First list is for peg extension, second list is for nicking guide
                                pe_design[target_name][pegid] = [[],[]]

                            if pe_pam_ref_silent_mutation == '':
                                pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), reverse_complement(pe_pam_ref), pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length])
                            
                            else:
                                pe_design[target_name][pegid][0].append([pe_nick_ref_idx, reverse_complement(pe_spacer_sequence), pe_pam_ref_silent_mutation, pe_annotate, '-', pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length])

                    # Create pegID
                    if 'PAM_disrupted' in pe_annotate:
                        pe_annotate_code = 0
                    else:
                        pe_annotate_code = 1

                    pegid = '_'.join(map(str, [str(pe_annotate_code) + '0'*(3 - len(str(abs(nick2lastedit_length)))) + str(abs(nick2lastedit_length)), pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

                    # Create ngRNAs targeting (+) strand for (-) pegRNAs
                    if pegid in pe_design[target_name]:
                    # if counter in counted:
                        for ng_plus in target_design[target_name]['ngRNA']['+']:
                            ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
                            nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

                            if (silent_mutation == 'yes') and (len(silent_mutation_edit) > 0):
                                ng_spacer_sequence_edit = silent_mutation_edit[ng_edit_start_idx:ng_edit_end_idx]

                                mutation_indices = [i for i, a in enumerate(ng_spacer_sequence_edit) if a.islower()]
                                if len(mutation_indices) > 0:
                                    if len([1 for x in mutation_indices if x >= 10]) > 0:
                                        ng_annotate = 'PE3b-seed'

                                    else:
                                        ng_annotate = 'PE3b-nonseed'
                                else:
                                    ng_annotate = 'PE3'

                            if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                                if ng_annotate == 'PE3b-seed':
                                    ng_code = 0
                                elif ng_annotate == 'PE3b-nonseed':
                                    ng_code = 1
                                else:
                                    ng_code = 2

                                pe_design[target_name][pegid][1].append([str(ng_code) + '0'*(3 - len(str(abs(abs(nick_distance) - nicking_distance_pooled)))) + str(abs(abs(nick_distance) - nicking_distance_pooled)), ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, '+', nick_distance])

                        pe_design[target_name][pegid][1] = sorted(pe_design[target_name][pegid][1])

            # Sort pegRNAs and ngRNAs and filter for top designs
            pe_design[target_name] = dict(sorted(pe_design[target_name].items(), key=lambda v: int(v[0].split('_')[0])))

        counter = 1
        for target_name in pe_design:

            if filter_c1_extension == 'yes':
                
                peg_count = 0
                for pegid in list(pe_design[target_name].keys()):

                    # Write pegRNAs
                    pegRNA_entry = pe_design[target_name][pegid][0][0]
                    pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length = pegRNA_entry

                    pegRNA_ext_first_base = pegRNA_ext[0]
                    spacer_gc_content = gc_content(pe_spacer_sequence)
                    pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
                    rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

                    if pegRNA_ext_first_base.upper() == 'C':

                        # Find minimum non-C index extending the RTT template
                        shift_rtt_index_list = [pegRNA_ext_max[:-len(pegRNA_ext)][::-1].upper().find('A'), pegRNA_ext_max[:-len(pegRNA_ext)][::-1].upper().find('G'), pegRNA_ext_max[:-len(pegRNA_ext)][::-1].upper().find('T')]
                        shift_rtt_index_list = [x for x in shift_rtt_index_list if x != -1]

                        # Make sure there are non-C indices
                        if len(shift_rtt_index_list) > 0:

                            shift_rtt_index = min(shift_rtt_index_list) + 1
                            rtt_length += shift_rtt_index

                            pegRNA_ext = pegRNA_ext_max[-len(pegRNA_ext) - shift_rtt_index:]
                            pegRNA_ext_first_base = pegRNA_ext[0]
                            spacer_gc_content = gc_content(pe_spacer_sequence)
                            pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
                            rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

                            peg_design['Target_name'].append(target_name)
                            peg_design['Target_sequence'].append(target_design[target_name]['target_sequence'])
                            peg_design['pegRNA_number'].append(counter)
                            peg_design['gRNA_type'].append('pegRNA')
                            peg_design['Spacer_sequence'].append(pe_spacer_sequence)
                            peg_design['Spacer_GC_content'].append(spacer_gc_content)
                            peg_design['PAM_sequence'].append(pe_pam_ref)
                            peg_design['Strand'].append(pe_strand)
                            peg_design['pegRNA-to-edit_distance'].append(nick2lastedit_length)
                            peg_design['ngRNA-to-pegRNA_distance'].append('')
                            peg_design['pegRNA_extension'].append(pegRNA_ext)
                            peg_design['First_extension_nucleotide'].append(pegRNA_ext_first_base)
                            peg_design['PBS_length'].append(pbs_length)
                            peg_design['PBS_GC_content'].append(pbs_gc_content)
                            peg_design['RTT_length'].append(rtt_length)
                            peg_design['RTT_GC_content'].append(rtt_gc_content)
                            peg_design['Annotation'].append(pe_annotate)

                            if pe_spacer_sequence[0].upper() == 'G':
                                peg_design['Spacer_top_strand_oligo'].append('cacc' + pe_spacer_sequence + 'gtttt')
                                peg_design['Spacer_bottom_strand_oligo'].append('ctctaaaac' + reverse_complement(pe_spacer_sequence))

                            else:
                                peg_design['Spacer_top_strand_oligo'].append('caccG' + pe_spacer_sequence + 'gtttt')
                                peg_design['Spacer_bottom_strand_oligo'].append('ctctaaaac' + reverse_complement('G' + pe_spacer_sequence))

                            peg_design['pegRNA_extension_top_strand_oligo'].append('gtgc' + pegRNA_ext)
                            peg_design['pegRNA_extension_bottom_strand_oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                            # f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, 'caccG' + pe_spacer_sequence[1:] + 'gtttt', 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence[1:]), 'gtgc' + pegRNA_ext, 'aaaa' + reverse_complement(pegRNA_ext)])) + '\n')
                            peg_count += 1
                            peg_count_dict[target_name] += 1

                    else:

                        peg_design['Target_name'].append(target_name)
                        peg_design['Target_sequence'].append(target_design[target_name]['target_sequence'])
                        peg_design['pegRNA_number'].append(counter)
                        peg_design['gRNA_type'].append('pegRNA')
                        peg_design['Spacer_sequence'].append(pe_spacer_sequence)
                        peg_design['Spacer_GC_content'].append(spacer_gc_content)
                        peg_design['PAM_sequence'].append(pe_pam_ref)
                        peg_design['Strand'].append(pe_strand)
                        peg_design['pegRNA-to-edit_distance'].append(nick2lastedit_length)
                        peg_design['ngRNA-to-pegRNA_distance'].append('')
                        peg_design['pegRNA_extension'].append(pegRNA_ext)
                        peg_design['First_extension_nucleotide'].append(pegRNA_ext_first_base)
                        peg_design['PBS_length'].append(pbs_length)
                        peg_design['PBS_GC_content'].append(pbs_gc_content)
                        peg_design['RTT_length'].append(rtt_length)
                        peg_design['RTT_GC_content'].append(rtt_gc_content)
                        peg_design['Annotation'].append(pe_annotate)

                        if pe_spacer_sequence[0].upper() == 'G':
                            peg_design['Spacer_top_strand_oligo'].append('cacc' + pe_spacer_sequence + 'gtttt')
                            peg_design['Spacer_bottom_strand_oligo'].append('ctctaaaac' + reverse_complement(pe_spacer_sequence))

                        else:
                            peg_design['Spacer_top_strand_oligo'].append('caccG' + pe_spacer_sequence + 'gtttt')
                            peg_design['Spacer_bottom_strand_oligo'].append('ctctaaaac' + reverse_complement('G' + pe_spacer_sequence))
                        
                        peg_design['pegRNA_extension_top_strand_oligo'].append('gtgc' + pegRNA_ext)
                        peg_design['pegRNA_extension_bottom_strand_oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                        # f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, 'caccG' + pe_spacer_sequence[1:] + 'gtttt', 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence[1:]), 'gtgc' + pegRNA_ext, 'aaaa' + reverse_complement(pegRNA_ext)])) + '\n')
                        peg_count += 1
                        peg_count_dict[target_name] += 1

                    # Write ngRNAs
                    for ngRNA_entry in pe_design[target_name][pegid][1][:number_of_ngrnas]:
                        ng_code, ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

                        spacer_gc_content = gc_content(ng_spacer_sequence_edit)

                        peg_design['Target_name'].append(target_name)
                        peg_design['Target_sequence'].append(target_design[target_name]['target_sequence'])
                        peg_design['pegRNA_number'].append(counter)
                        peg_design['gRNA_type'].append('ngRNA')
                        peg_design['Spacer_sequence'].append(ng_spacer_sequence_edit)
                        peg_design['Spacer_GC_content'].append(spacer_gc_content)
                        peg_design['PAM_sequence'].append(ng_pam_edit)
                        peg_design['Strand'].append(ng_strand)
                        peg_design['pegRNA-to-edit_distance'].append('')
                        peg_design['ngRNA-to-pegRNA_distance'].append(nick_distance)
                        peg_design['pegRNA_extension'].append('')
                        peg_design['First_extension_nucleotide'].append('')
                        peg_design['PBS_length'].append('')
                        peg_design['PBS_GC_content'].append('')
                        peg_design['RTT_length'].append('')
                        peg_design['RTT_GC_content'].append('')
                        peg_design['Annotation'].append(ng_annotate)

                        if ng_spacer_sequence_edit[0].upper() == 'G':
                            peg_design['Spacer_top_strand_oligo'].append('cacc' + ng_spacer_sequence_edit)
                            peg_design['Spacer_bottom_strand_oligo'].append('aaac' + reverse_complement(ng_spacer_sequence_edit))

                        else:
                            peg_design['Spacer_top_strand_oligo'].append('caccG' + ng_spacer_sequence_edit)
                            peg_design['Spacer_bottom_strand_oligo'].append('aaac' + reverse_complement('G' + ng_spacer_sequence_edit))

                        peg_design['pegRNA_extension_top_strand_oligo'].append('')
                        peg_design['pegRNA_extension_bottom_strand_oligo'].append('')

                        # f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', 'caccG' + ng_spacer_sequence_edit[1:], 'aaac' + reverse_complement('G' + ng_spacer_sequence_edit[1:]), '', ''])) + '\n')

                    counter += 1

                    if peg_count == number_of_pegrnas:
                        break

            else:

                for pegid in list(pe_design[target_name].keys())[:number_of_pegrnas]:

                    # Write pegRNAs
                    pegRNA_entry = pe_design[target_name][pegid][0][0]
                    pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, pe_strand, pbs_length, rtt_length, pegRNA_ext, pegRNA_ext_max, nick2lastedit_length = pegRNA_entry

                    pegRNA_ext_first_base = pegRNA_ext[0]
                    spacer_gc_content = gc_content(pe_spacer_sequence)
                    pbs_gc_content = gc_content(pegRNA_ext[rtt_length:])
                    rtt_gc_content = gc_content(pegRNA_ext[:rtt_length])

                    peg_design['Target_name'].append(target_name)
                    peg_design['Target_sequence'].append(target_design[target_name]['target_sequence'])
                    peg_design['pegRNA_number'].append(counter)
                    peg_design['gRNA_type'].append('pegRNA')
                    peg_design['Spacer_sequence'].append(pe_spacer_sequence)
                    peg_design['Spacer_GC_content'].append(spacer_gc_content)
                    peg_design['PAM_sequence'].append(pe_pam_ref)
                    peg_design['Strand'].append(pe_strand)
                    peg_design['pegRNA-to-edit_distance'].append(nick2lastedit_length)
                    peg_design['ngRNA-to-pegRNA_distance'].append('')
                    peg_design['pegRNA_extension'].append(pegRNA_ext)
                    peg_design['First_extension_nucleotide'].append(pegRNA_ext_first_base)
                    peg_design['PBS_length'].append(pbs_length)
                    peg_design['PBS_GC_content'].append(pbs_gc_content)
                    peg_design['RTT_length'].append(rtt_length)
                    peg_design['RTT_GC_content'].append(rtt_gc_content)
                    peg_design['Annotation'].append(pe_annotate)

                    if pe_spacer_sequence[0].upper() == 'G':
                        peg_design['Spacer_top_strand_oligo'].append('cacc' + pe_spacer_sequence + 'gtttt')
                        peg_design['Spacer_bottom_strand_oligo'].append('ctctaaaac' + reverse_complement(pe_spacer_sequence))

                    else:
                        peg_design['Spacer_top_strand_oligo'].append('caccG' + pe_spacer_sequence + 'gtttt')
                        peg_design['Spacer_bottom_strand_oligo'].append('ctctaaaac' + reverse_complement('G' + pe_spacer_sequence))
                    
                    peg_design['pegRNA_extension_top_strand_oligo'].append('gtgc' + pegRNA_ext)
                    peg_design['pegRNA_extension_bottom_strand_oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                    peg_count_dict[target_name] += 1

                    # f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'pegRNA', pe_spacer_sequence, spacer_gc_content, pe_pam_ref, pegRNA_ext, pe_strand, pe_annotate, nick2lastedit_length, pe_nick_ref_idx, '', pbs_length, pbs_gc_content, rtt_length, rtt_gc_content, pegRNA_ext_first_base, 'caccG' + pe_spacer_sequence[1:] + 'gtttt', 'ctctaaaac' + reverse_complement('G' + pe_spacer_sequence[1:]), 'gtgc' + pegRNA_ext, 'aaaa' + reverse_complement(pegRNA_ext)])) + '\n')

                    # Sort ngRNAs

                    # Write ngRNAs
                    for ngRNA_entry in pe_design[target_name][pegid][1][:number_of_ngrnas]:
                        ng_code, ng_nick_ref_idx, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate, ng_strand, nick_distance = ngRNA_entry

                        spacer_gc_content = gc_content(ng_spacer_sequence_edit)

                        peg_design['Target_name'].append(target_name)
                        peg_design['Target_sequence'].append(target_design[target_name]['target_sequence'])
                        peg_design['pegRNA_number'].append(counter)
                        peg_design['gRNA_type'].append('ngRNA')
                        peg_design['Spacer_sequence'].append(ng_spacer_sequence_edit)
                        peg_design['Spacer_GC_content'].append(spacer_gc_content)
                        peg_design['PAM_sequence'].append(ng_pam_edit)
                        peg_design['Strand'].append(ng_strand)
                        peg_design['pegRNA-to-edit_distance'].append('')
                        peg_design['ngRNA-to-pegRNA_distance'].append(nick_distance)
                        peg_design['pegRNA_extension'].append('')
                        peg_design['First_extension_nucleotide'].append('')
                        peg_design['PBS_length'].append('')
                        peg_design['PBS_GC_content'].append('')
                        peg_design['RTT_length'].append('')
                        peg_design['RTT_GC_content'].append('')
                        peg_design['Annotation'].append(ng_annotate)

                        if ng_spacer_sequence_edit[0].upper() == 'G':
                            peg_design['Spacer_top_strand_oligo'].append('cacc' + ng_spacer_sequence_edit)
                            peg_design['Spacer_bottom_strand_oligo'].append('aaac' + reverse_complement(ng_spacer_sequence_edit))

                        else:
                            peg_design['Spacer_top_strand_oligo'].append('caccG' + ng_spacer_sequence_edit)
                            peg_design['Spacer_bottom_strand_oligo'].append('aaac' + reverse_complement('G' + ng_spacer_sequence_edit))
                        
                        peg_design['pegRNA_extension_top_strand_oligo'].append('')
                        peg_design['pegRNA_extension_bottom_strand_oligo'].append('')

                        # f.write(','.join(map(str, [target_name, target_design[target_name]['target_sequence'], counter, 'ngRNA', ng_spacer_sequence_edit, spacer_gc_content, ng_pam_edit, '', ng_strand, ng_annotate, '', ng_nick_ref_idx, nick_distance, '', '', '', '', '', 'caccG' + ng_spacer_sequence_edit[1:], 'aaac' + reverse_complement('G' + ng_spacer_sequence_edit[1:]), '', ''])) + '\n')

                    counter += 1

        df = pd.DataFrame.from_dict(peg_design)

        # if extfirstbase_filter == 'yes':
        #     df = df[df['First_extension_nucleotide'] != 'C']
        #     df.reset_index(drop=True, inplace=True)

        df.to_csv('/PrimeDesign/reports/PrimeDesign_Pooled_%s.csv' % session_id)

        for x in peg_count_dict:
            if peg_count_dict[x] != number_of_pegrnas:
                warning_list.append(html.H5('%s: %s designs' % (str(x), str(peg_count_dict[x])), style = {'color':'#ff4d4d', 'font-size':'15px'}))

        if len(warning_list) > 0:
            warning_list = [html.H5('Inputs with less than the desired number of pegRNA designs', style = {'color':'#ff4d4d', 'font-size':'15px', 'font-weight':'bold'})] + warning_list

        return('Design completed', warning_list, 'Download designs')
    else:

        return('Design incomplete', ' ',' ')

# Download pooled example file update
# html.A(children = 'Download genome wide example file', id='download-example-pool', download="PrimeDesign_genome_wide_example_file.csv", href="/download/PrimeDesign_genome_wide_example_file.csv"

@app.callback([Output('download-example-pool', 'children'),Output('download-example-pool', 'download'),Output('download-example-pool', 'href'),],
    [Input('design-option-pool', 'value')],
)

def update_example_file(pooled_design_type):

    if pooled_design_type == 'saturation_mutagenesis':
        return('Download saturation mutagenesis example file', 'PrimeDesign_saturation_mutagenesis_example_file.csv', '/download/PrimeDesign_saturation_mutagenesis_example_file.csv')
    else:
        return('Download genome wide example file', 'PrimeDesign_genome_wide_example_file.csv', '/download/PrimeDesign_genome_wide_example_file.csv')

# PooledDesign logic
@app.callback(Output('satmut-type-container', 'style'),
    [Input('design-option-pool','value')],
)

def update_input_check(pooled_design_type):

    if pooled_design_type == 'saturation_mutagenesis':
        return({'display':'block'})
    else:
        return({'display':'none'})

@app.callback([Output('genome-wide-format-container', 'style'), Output('saturation-mutation-format-container', 'style'),],
    [Input('design-option-pool','value')],
)

def update_formatting_container(pooled_design_type):

    if pooled_design_type == 'saturation_mutagenesis':
        return({'display':'none'}, {'display':'block'})
    else:
        return({'display':'block'}, {'display':'none'})


### Logic for PrimeVar

@app.callback([Output('primevar-input-check', 'children'), Output('primevar-input-check', 'style'),],
    [Input('primevar-id-search-type','value'), Input('primevar-id-search','value'), Input('editing-direction', 'value')]
)

def update_input_check(primevar_id_search_type, primevar_id_search, editing_direction):

    if (primevar_id_search is not None):

        try:
            primevar_id_search = primevar_id_search.replace('rs', '')
            tmp = primevar_map[primevar_id_search_type][editing_direction][primevar_id_search]
            sequence_check = 'Success'
            sequence_check_style = {'color':'#6bb6ff'}

        except:
            sequence_check = 'Designs are not available'
            sequence_check_style = {'color':'#ff4d4d'}


    else:
        sequence_check = 'No variant ID has been provided'
        sequence_check_style = {'color':'#ff4d4d'}

    return(sequence_check, sequence_check_style)


@app.callback([Output('reference-sequence-db', 'sequence'), Output('reference-sequence-db', 'coverage'), Output('edit-sequence-db', 'sequence'), Output('edit-sequence-db', 'coverage')],
    [Input('primevar-input-check','children'), Input('peg-table-db','selected_rows'), Input('pegext-table-db','selected_rows'), Input('ng-table-db','selected_rows')],
    state = [State('pbs-range-db','value'), State('rtt-range-db','value'), State('nick-dist-range-db','value'), State('primevar-id-search-type','value'), State('primevar-id-search','value'), State('editing-direction','value'), State('store-peg-table-db', 'children'), State('store-peg-table-total-db', 'children')]
)

def update_reference_sequence(primevar_input, selected_rows_peg, selected_rows_pegext, selected_rows_ng, pbs_range, rtt_range, nicking_distance_range, primevar_id_search_type, primevar_id_search, editing_direction, store_peg_table, store_peg_table_total):

    annotations_ref = []
    annotations_edit = []

    if 'Success' in primevar_input:

        primevar_id_search = primevar_id_search.replace('rs', '')
        zip_file = ZipFile('/PrimeDesign/PrimeVar/' + primevar_map[primevar_id_search_type][editing_direction][primevar_id_search])
        df = pd.read_csv(zip_file.open(primevar_map[primevar_id_search_type][editing_direction][primevar_id_search].replace('.zip', '')), names = ['Target_name', 'Target_sequence', 'pegRNA group', 'type', 'spacer sequence', 'spacer GC content', 'PAM', 'pegRNA extension', 'strand', 'annotation', 'peg-to-edit distance', 'Nick_index', 'nick-to-peg distance', 'PBS length', 'PBS GC content', 'RTT length', 'RTT GC content', 'extension first base', 'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM', 'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM'], header = None)

        input_sequence = df.Target_sequence[0]

        reference_sequence = input_sequence
        edit_sequence = input_sequence
        editformat2sequence_ref = {}
        editformat2sequence_edit = {}
        index_shift_ref = 0
        index_shift_edit = 0

        edit_idxs = [[m.start(), m.end()] for m in re.finditer('\(.*?\)', input_sequence)]
        for edit_idx in edit_idxs:

            edit = input_sequence[edit_idx[0]:edit_idx[1]]
            edit_length = edit_idx[1] - edit_idx[0]

            # Create edit format and number to sequence map
            if '/' in edit:
                editformat2sequence_ref[edit] = edit.split('/')[0].replace('(','')

                if len(edit.split('/')[1].replace(')','')) == 0:
                    annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref + len(edit.split('/')[0].replace('(','')), 'color':'#DC143C', 'bgcolor':'#fbe7eb', 'underscore':True})
                else:
                    annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref + len(edit.split('/')[0].replace('(','')), 'color':'#1E90FF', 'bgcolor':'#e8f3ff', 'underscore':True})
                
                index_shift_ref += edit_length - len(edit.split('/')[0].replace('(',''))

            elif '+' in edit:
                editformat2sequence_ref[edit] = ''
                annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref, 'color':'#3CB371', 'bgcolor':'#ebf7f0', 'underscore':True})

                index_shift_ref += edit_length

            elif '-' in edit:
                editformat2sequence_ref[edit] = edit.split('-')[1].replace(')','')
                annotations_ref.append({'start':edit_idx[0] - index_shift_ref, 'end':edit_idx[0] - index_shift_ref + len(edit.split('-')[1].replace(')','')), 'color':'#DC143C', 'bgcolor':'#fbe7eb', 'underscore':True})

                index_shift_ref += edit_length - len(edit.split('-')[1].replace(')',''))

            # Create edit format and number to sequence map
            if '/' in edit:
                editformat2sequence_edit[edit] = edit.split('/')[1].replace(')','')

                if len(edit.split('/')[0].replace('(','')) == 0:
                    annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit + len(edit.split('/')[1].replace(')','')), 'color':'#3CB371', 'bgcolor':'#ebf7f0', 'underscore':True})
                else:
                    annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit + len(edit.split('/')[1].replace(')','')), 'color':'#1E90FF', 'bgcolor':'#e8f3ff', 'underscore':True})

                index_shift_edit += edit_length - len(edit.split('/')[1].replace(')',''))

            elif '+' in edit:
                editformat2sequence_edit[edit] = edit.split('+')[1].replace(')','')
                annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit + len(edit.split('+')[1].replace(')','')), 'color':'#3CB371', 'bgcolor':'#ebf7f0', 'underscore':True})

                index_shift_edit += edit_length -len(edit.split('+')[1].replace(')',''))

            elif '-' in edit:
                editformat2sequence_edit[edit] = ''
                annotations_edit.append({'start':edit_idx[0] - index_shift_edit, 'end':edit_idx[0] - index_shift_edit, 'color':'#DC143C', 'bgcolor':'#fbe7eb', 'underscore':True})

                index_shift_edit += edit_length

        for edit in editformat2sequence_ref:
            reference_sequence = reference_sequence.replace(edit, editformat2sequence_ref[edit])

        for edit in editformat2sequence_edit:
            edit_sequence = edit_sequence.replace(edit, editformat2sequence_edit[edit])

        # Visualizing pegRNA spacer in reference sequence
        try:
            current_annotation_ranges_ref = []
            current_annotation_ranges_edit = []
            for annotation in annotations_ref:
                current_annotation_ranges_ref.append([annotation['start'], annotation['end']])

            for annotation in annotations_edit:
                current_annotation_ranges_edit.append([annotation['start'], annotation['end']])

            df_peg = pd.read_json(store_peg_table, orient='split')
            spacer_sequences = list(df_peg.loc[selected_rows_peg, 'spacer sequence'].values)

            # Annotate pegRNA spacer sequences
            for spacer_sequence in spacer_sequences:

                try:
                    start_idx = re.search(spacer_sequence, reference_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(spacer_sequence)
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                            annotations_ref.append({'start':i, 'end':i + 1, 'underscore':True})
                            current_annotation_ranges_ref.append([i, i + 1])

                    start_idx = re.search(spacer_sequence[:17], edit_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(spacer_sequence) - 3
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                            annotations_edit.append({'start':i, 'end':i + 1, 'underscore':True})
                            current_annotation_ranges_edit.append([i, i + 1])

                except:
                    start_idx = re.search(reverse_complement(spacer_sequence), reference_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(spacer_sequence)
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                            annotations_ref.append({'start':i, 'end':i + 1, 'underscore':True})
                            current_annotation_ranges_ref.append([i, i + 1])

                    start_idx = re.search(reverse_complement(spacer_sequence[:17]), edit_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(spacer_sequence) - 3
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                            annotations_edit.append({'start':i, 'end':i + 1, 'underscore':True})
                            current_annotation_ranges_edit.append([i, i + 1])

                # try:
                #     start_idx = re.search(spacer_sequence, reference_sequence, re.IGNORECASE).start()
                #     stop_idx = start_idx + len(spacer_sequence)
                #     for i in range(start_idx, stop_idx):
                #         if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                #             annotations_ref.append({'start':i, 'end':i + 1, 'bgcolor':'#dedede'})
                #             current_annotation_ranges.append([i, i + 1])

                # except:
                #     start_idx = re.search(reverse_complement(spacer_sequence), reference_sequence, re.IGNORECASE).start()
                #     stop_idx = start_idx + len(spacer_sequence)
                #     for i in range(start_idx, stop_idx):
                #         if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                #             annotations_ref.append({'start':i, 'end':i + 1, 'bgcolor':'#dedede'})
                #             current_annotation_ranges.append([i, i + 1])

        except:
            pass

        # Visualizing pegRNA extension in edit sequence
        try:
            current_annotation_ranges = []
            for annotation in annotations_edit:
                current_annotation_ranges.append([annotation['start'], annotation['end']])

            df_peg = pd.read_json(store_peg_table, orient='split')
            df_peg_total = pd.read_json(store_peg_table_total, orient='split')

            peg_group = list(df_peg.loc[selected_rows_peg, 'spacer sequence'].values)
            df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(peg_group)]
            df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
            df_pegext = df_pegext[['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']].drop_duplicates()

            df_pegext = df_pegext[(df_pegext['PBS length'] >= pbs_range[0]) & (df_pegext['PBS length'] <= pbs_range[1])]
            df_pegext = df_pegext[(df_pegext['RTT length'] >= rtt_range[0]) & (df_pegext['RTT length'] <= rtt_range[1])]

            df_pegext.reset_index(drop=True, inplace=True)
            pegext_sequences = list(df_pegext.loc[selected_rows_pegext, 'pegRNA extension'].values)
            pbs_lengths = list(df_pegext.loc[selected_rows_pegext, 'PBS length'].values)

            # Annotate pegRNA spacer sequences
            for pegext_sequence, pbs_length in zip(pegext_sequences, pbs_lengths):
                pbs_length = int(pbs_length)

                try:

                    start_idx = re.search(pegext_sequence, edit_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(pegext_sequence)

                    for entry in annotations_edit:
                        start_entry = entry['start']
                        if start_entry in range(stop_idx - pbs_length, stop_idx):
                            entry['bgcolor'] = '#cebeef'

                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                            annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#ffdb99', })
                            current_annotation_ranges.append([i, i + 1])

                except:

                    start_idx = re.search(reverse_complement(pegext_sequence), edit_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(pegext_sequence)

                    for entry in annotations_edit:
                        start_entry = entry['start']
                        if start_entry in range(start_idx, start_idx + pbs_length):
                            entry['bgcolor'] = '#cebeef'

                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges]) == 0:
                            annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#ffdb99', })
                            current_annotation_ranges.append([i, i + 1])

        except:
            pass

        # Visualizing ngRNA spacer in edit sequence
        try:
            current_annotation_ranges_ref = []
            current_annotation_ranges_edit = []

            for annotation in annotations_ref:
                current_annotation_ranges_ref.append([annotation['start'], annotation['end']])

            for annotation in annotations_edit:
                current_annotation_ranges_edit.append([annotation['start'], annotation['end']])

            df_peg = pd.read_json(store_peg_table, orient='split')
            df_peg_total = pd.read_json(store_peg_table_total, orient='split')

            peg_group = list(df_peg.loc[selected_rows_peg, 'pegRNA group'].values)
            df_ng = df_peg_total[df_peg_total['pegRNA group'].isin(peg_group)]
            df_ng = df_ng[df_ng['type'] == 'ngRNA']
            df_ng = df_ng[['spacer sequence','PAM','strand','nick-to-peg distance','spacer GC content','annotation']].drop_duplicates()

            df_ng = df_ng[(abs(df_ng['nick-to-peg distance']) >= nicking_distance_range[0]) & (abs(df_ng['nick-to-peg distance']) <= nicking_distance_range[1])]

            df_ng.reset_index(drop=True, inplace=True)
            ngRNA_sequences = list(df_ng.loc[selected_rows_ng, 'spacer sequence'].values)

            # Annotate pegRNA spacer sequences
            for ngRNA_sequence in ngRNA_sequences:

                try:
                    start_idx = re.search(ngRNA_sequence, edit_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(ngRNA_sequence)
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                            annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                            current_annotation_ranges_edit.append([i, i + 1])

                    try:
                        start_idx = re.search(ngRNA_sequence, reference_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(ngRNA_sequence)
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                                annotations_ref.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                                current_annotation_ranges_ref.append([i, i + 1])

                    except:
                        pass

                except:
                    start_idx = re.search(reverse_complement(ngRNA_sequence), edit_sequence, re.IGNORECASE).start()
                    stop_idx = start_idx + len(ngRNA_sequence)
                    for i in range(start_idx, stop_idx):
                        if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_edit]) == 0:
                            annotations_edit.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                            current_annotation_ranges_edit.append([i, i + 1])

                    try:
                        start_idx = re.search(reverse_complement(ngRNA_sequence), reference_sequence, re.IGNORECASE).start()
                        stop_idx = start_idx + len(ngRNA_sequence)
                        for i in range(start_idx, stop_idx):
                            if sum([1 if (x[0] <= i < x[1]) else 0 for x in current_annotation_ranges_ref]) == 0:
                                annotations_ref.append({'start':i, 'end':i + 1, 'bgcolor':'#d6d6d6'})
                                current_annotation_ranges_ref.append([i, i + 1])

                    except:
                        pass

        except:
            pass

    else:
        reference_sequence = ' '
        edit_sequence = ' '

    return(reference_sequence, annotations_ref, edit_sequence, annotations_edit)

@app.callback([Output('peg-table-db', 'data'), Output('store-peg-table-total-db', 'children'), Output('store-peg-table-db', 'children')],
    [Input('primevar-input-check','children'), Input('pbs-range-db','value'), Input('rtt-range-db','value'), Input('nick-dist-range-db','value'), Input('filter-c1-extension-option-db','value'), Input('silentmutation-option-db','value')],
    state = [State('primevar-id-search-type','value'), State('primevar-id-search','value'), State('editing-direction', 'value'), State('session-id', 'children')]
)

def update_database_tables(primevar_input, pbs_range, rtt_range, nicking_distance_range, filter_c1_extension, silent_mutation, primevar_id_search_type, primevar_id_search, editing_direction, session_id):

    if 'Success' in primevar_input:

        primevar_id_search = primevar_id_search.replace('rs', '')
        zip_file = ZipFile('/PrimeDesign/PrimeVar/' + primevar_map[primevar_id_search_type][editing_direction][primevar_id_search])
        df = pd.read_csv(zip_file.open(primevar_map[primevar_id_search_type][editing_direction][primevar_id_search].replace('.zip', '')), names = ['Target_name', 'Target_sequence', 'pegRNA group', 'type', 'spacer sequence', 'spacer GC content', 'PAM', 'pegRNA extension', 'strand', 'annotation', 'peg-to-edit distance', 'Nick_index', 'nick-to-peg distance', 'PBS length', 'PBS GC content', 'RTT length', 'RTT GC content', 'extension first base', 'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM', 'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM'], header = None)

        df = df.replace('PAM_mutated','PAM_disrupted')

        if filter_c1_extension == 'yes':
            df = df[df['extension first base'] != 'C']
            df.reset_index(drop=True, inplace=True)

        df_pegs = df[df['type'] == 'pegRNA']
        df_pegs = df_pegs[(df_pegs['RTT length'] >= rtt_range[0]) & (df_pegs['RTT length'] <= rtt_range[1])]
        df_pegs = df_pegs[['pegRNA group','spacer sequence','PAM','strand','peg-to-edit distance','spacer GC content','annotation']].drop_duplicates()
        df_pegs['spacer GC content'] = df_pegs['spacer GC content'].round(2)
        df_pegs = df_pegs.sort_values('peg-to-edit distance')
        df_pegs.reset_index(drop=True, inplace=True)

        df.to_csv('/PrimeDesign/reports/PrimeDesign_PrimeVar_%s.csv' % session_id)

        # if primevar_id_search_type == 'rs':
        #     df.to_csv('/PrimeDesign/reports/PrimeVar_dbSNPrs%s.csv' % (str(primevar_id_search)))

        # else:
        #     df.to_csv('/PrimeDesign/reports/PrimeVar_ClinVarVariationID%s.csv' % (str(primevar_id_search)))

    else:
        peg_design = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df = pd.DataFrame.from_dict(peg_design)
        df_pegs = pd.DataFrame.from_dict(peg_design)

    return(df_pegs.to_dict('records'), df.to_json(date_format='iso', orient='split'), df_pegs.to_json(date_format='iso', orient='split'))

# Trigger pegRNA extension and ngRNA tables with pegRNA spacer selection for PrimeVar database
@app.callback(Output('pegext-table-db', 'data'),
    [Input('peg-table-db','selected_rows'), Input('store-peg-table-total-db', 'children'), Input('store-peg-table-db', 'children')],
    state = [State('pbs-range-db','value'), State('rtt-range-db','value'),]
)

def update_pegext_table(selected_row, store_peg_table_total, store_peg_table, pbs_range, rtt_range):

    if selected_row:

        try:
            # Open up stored peg table
            df_peg = pd.read_json(store_peg_table, orient='split')
            df_peg_total = pd.read_json(store_peg_table_total, orient='split')

            spacer_sequence = list(df_peg.loc[selected_row, 'spacer sequence'].values)
            df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(spacer_sequence)]
            df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
            df_pegext = df_pegext[['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']].drop_duplicates()
            df_pegext['PBS GC content'] = df_pegext['PBS GC content'].round(2)
            df_pegext['RTT GC content'] = df_pegext['RTT GC content'].round(2)
            df_pegext = df_pegext[(df_pegext['PBS length'] >= pbs_range[0]) & (df_pegext['PBS length'] <= pbs_range[1])]
            df_pegext = df_pegext[(df_pegext['RTT length'] >= rtt_range[0]) & (df_pegext['RTT length'] <= rtt_range[1])]
            df_pegext.reset_index(drop=True, inplace=True)

        except:
            df_pegext = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
            df_pegext = pd.DataFrame.from_dict(df_pegext)

    else:
        df_pegext = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df_pegext = pd.DataFrame.from_dict(df_pegext)

    return(df_pegext.to_dict('records'))

@app.callback(Output('ng-table-db', 'data'),
    [Input('peg-table-db','selected_rows'), Input('store-peg-table-total-db', 'children'), Input('store-peg-table-db', 'children')],
    state = [State('nick-dist-range-db','value')]
)

def update_ng_table(selected_row, store_peg_table_total, store_peg_table, nicking_distance_range):

    if selected_row:

        try:
            # Open up stored peg table
            df_peg = pd.read_json(store_peg_table, orient='split')
            df_peg_total = pd.read_json(store_peg_table_total, orient='split')

            peg_group = list(df_peg.loc[selected_row, 'pegRNA group'].values)
            df_ng = df_peg_total[df_peg_total['pegRNA group'].isin(peg_group)]
            df_ng = df_ng[df_ng['type'] == 'ngRNA']
            df_ng = df_ng[['spacer sequence','PAM','strand','nick-to-peg distance','spacer GC content','annotation']].drop_duplicates()
            df_ng['spacer GC content'] = df_ng['spacer GC content'].round(2)
            df_ng = df_ng[(abs(df_ng['nick-to-peg distance']) >= nicking_distance_range[0]) & (abs(df_ng['nick-to-peg distance']) <= nicking_distance_range[1])]
            df_ng.reset_index(drop=True, inplace=True)

        except:
            df_ng = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
            df_ng = pd.DataFrame.from_dict(df_ng)

    else:
        df_ng = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[],'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[]}
        df_ng = pd.DataFrame.from_dict(df_ng)

    return(df_ng.to_dict('records'))

@app.callback(Output('pbs-title-db', 'children'),
    [Input('pbs-range-db','value')]
)

def update_pbs_title(pbs_range):
    return('PBS length: %s - %s nt' % (pbs_range[0], pbs_range[1]))

@app.callback(Output('rtt-title-db', 'children'),
    [Input('rtt-range-db','value')]
)

def update_pbs_title(rtt_range):
    return('RTT length: %s - %s nt' % (rtt_range[0], rtt_range[1]))

@app.callback(Output('nick-dist-title-db', 'children'),
    [Input('nick-dist-range-db','value')]
)

def update_pbs_title(nick_dist_range):
    return('Nicking distance: %s - %s bp' % (nick_dist_range[0], nick_dist_range[1]))

@app.callback(Output('download-link-db', 'href'),
    [Input('primevar-input-check','children')],
    state = [State('session-id', 'children')]
    # state = [State('primevar-id-search-type','value'), State('primevar-id-search','value'), State('editing-direction', 'value')]
)
def update_download_link(input_check, session_id):
    return('/download/PrimeDesign_PrimeVar_%s.csv' % session_id)

# RNA folding design
@app.callback(Output('forna-pegext', 'sequences'),
    [Input('pegext-table','selected_rows'), Input('peg-table','selected_rows'), Input('forna-temp','value'), Input('forna-option','value'),],
    state = [State('store-peg-table', 'children'), State('store-peg-table-total', 'children')],
)
def show_selected_sequences(selected_rows_pegext, selected_rows_peg, temperature, fold_option, store_peg_table, store_peg_table_total):
    
    if selected_rows_pegext is not None:

        df_peg = pd.read_json(store_peg_table, orient='split')
        df_peg_total = pd.read_json(store_peg_table_total, orient='split')

        peg_group = list(df_peg.loc[selected_rows_peg, 'spacer sequence'].values)
        df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(peg_group)]
        df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
        df_pegext = df_pegext[['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']].drop_duplicates()

        df_pegext = df_pegext.sort_values(['PBS length', 'RTT length'])

        df_pegext.reset_index(drop=True, inplace=True)
        peg_spacer_sequence = peg_group[0]
        tracr_sequence = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
        pegext_sequence = list(df_pegext.loc[selected_rows_pegext, 'pegRNA extension'].values)[0]

        pegRNA_complete = peg_spacer_sequence + tracr_sequence + pegext_sequence

        if fold_option: # right -> full pegrna
            result = subprocess.run(['seqfold', pegRNA_complete.replace('T','U').replace('t','U'), '-v', '-t', str(temperature)], stdout=subprocess.PIPE)
            pegext_structure = result.stdout.split(b'\n')[1].decode("utf-8")
            return [{'sequence':pegRNA_complete, 'structure':pegext_structure}]
            
        else: # left -> peg extension
            result = subprocess.run(['seqfold', pegext_sequence.replace('T','U').replace('t','U'), '-v', '-t', str(temperature)], stdout=subprocess.PIPE)
            pegext_structure = result.stdout.split(b'\n')[1].decode("utf-8")
            return [{'sequence':pegext_sequence, 'structure':pegext_structure}]

    else:
        return [{'sequence':'', 'structure':''}]

# RNA folding primevar
@app.callback(Output('forna-pegext-db', 'sequences'),
    [Input('pegext-table-db','selected_rows'), Input('peg-table-db','selected_rows'), Input('forna-temp-db','value'), Input('forna-option-db','value'),],
    state = [State('pbs-range-db','value'), State('rtt-range-db','value'), State('store-peg-table-db', 'children'), State('store-peg-table-total-db', 'children')],
)
def show_selected_sequences(selected_rows_pegext, selected_rows_peg, temperature, fold_option, pbs_range, rtt_range, store_peg_table, store_peg_table_total):
    
    if selected_rows_pegext is not None:

        df_peg = pd.read_json(store_peg_table, orient='split')
        df_peg_total = pd.read_json(store_peg_table_total, orient='split')

        peg_group = list(df_peg.loc[selected_rows_peg, 'spacer sequence'].values)
        df_pegext = df_peg_total[df_peg_total['spacer sequence'].isin(peg_group)]
        df_pegext = df_pegext[df_pegext['type'] == 'pegRNA']
        df_pegext = df_pegext[['PBS length','PBS GC content','RTT length','RTT GC content','pegRNA extension']].drop_duplicates()

        df_pegext = df_pegext[(df_pegext['PBS length'] >= pbs_range[0]) & (df_pegext['PBS length'] <= pbs_range[1])]
        df_pegext = df_pegext[(df_pegext['RTT length'] >= rtt_range[0]) & (df_pegext['RTT length'] <= rtt_range[1])]

        df_pegext.reset_index(drop=True, inplace=True)
        peg_spacer_sequence = peg_group[0]
        tracr_sequence = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
        pegext_sequence = list(df_pegext.loc[selected_rows_pegext, 'pegRNA extension'].values)[0]

        pegRNA_complete = peg_spacer_sequence + tracr_sequence + pegext_sequence

        if fold_option: # right -> full pegrna
            result = subprocess.run(['seqfold', pegRNA_complete.replace('T','U').replace('t','U'), '-v', '-t', str(temperature)], stdout=subprocess.PIPE)
            pegext_structure = result.stdout.split(b'\n')[1].decode("utf-8")
            return [{'sequence':pegRNA_complete, 'structure':pegext_structure}]
            
        else: # left -> peg extension
            result = subprocess.run(['seqfold', pegext_sequence.replace('T','U').replace('t','U'), '-v', '-t', str(temperature)], stdout=subprocess.PIPE)
            pegext_structure = result.stdout.split(b'\n')[1].decode("utf-8")
            return [{'sequence':pegext_sequence, 'structure':pegext_structure}]

    else:
        return [{'sequence':'', 'structure':''}]

if __name__ == '__main__':
    app.run_server(debug = True, port = 9994, host = '0.0.0.0')
    # app.run_server(debug=True)