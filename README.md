# PrimeDesign

PrimeDesign is a software tool for the flexible and comprehensive design of prime editing. PrimeDesign is an edit-centric design tool for the installation of substitution, insertion, and deletion edits, and is generalizable for both single and combinatorial edits. Taking a single input that encodes both the reference and edit sequence, PrimeDesign enumerates all possible pegRNA protospacers, pegRNA extensions, and ngRNAs within set parameter ranges to install the desired edit(s).

## Installation with Docker

With Docker, no installation is required - the only dependence is Docker itself. Users will not need to deal with installation and configuration issues.

Docker can be downloaded freely here: [https://store.docker.com/search?offering=community&type=edition](https://store.docker.com/search?offering=community&type=edition)

Following the installation of Docker, simply execute the following command in the terminal for a local version of PrimeDesign:
* ```docker pull pinellolab/primedesign```

## PrimeDesign web application

Run the PrimeDesign web application locally with the command:

```
docker run -p 9994:9994 pinellolab/primedesign:latest primedesign_webapp
```
After execution of the command, the user will have a local instance of the website accessible at the URL: http://localhost:9994

## PrimeDesign command line tool

Run the PrimeDesign command line interface (CLI) with the command in the terminal:

```
docker run -v ${PWD}/:/DATA -w /DATA pinellolab/primedesign primedesign_cli [options]
```

Users can specify the following options:
```
-f, --file
      Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (Required)

-pe_format, --pe_format
      Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM] (Default: NNNNNNNNNNNNNNNNN/NNN[NGG]). Examples: NNNNNNNNNNNNNNNNN/NNN[NGG], NNNNNNNNNNNNNNNNN/NNN[NG]
-pbs, --pbs_length_list
      List of primer binding site (PBS) lengths for the pegRNA extension (Default: 10 to 15 nt). Example: 12 13 14 15
-rtt, --rtt_length_list
      List of reverse transcription (RT) template lengths for the pegRNA extension (Default: 10 to 30 nt). Example: 10 15 20
-nick_dist_min, --nicking_distance_minimum
      Minimum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 0 bp).
-nick_dist_max, --nicking_distance_maximum
      Maximum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 100 bp).
-filter_c1, --filter_c1_extension
      Option to filter against pegRNA extensions that start with a C base.
-silent_mut, --silent_mutation
      Introduce silent mutation into PAM assuming sequence is in-frame. Currently only available with SpCas9.
-genome_wide, --genome_wide_design
      Whether or not this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters.
-sat_mut, --saturation_mutagenesis
      Saturation mutagenesis design with prime editing. The 'aa' option makes amino acid changes and the 'base' option makes DNA base changes. (Options: aa, base)
-n_pegrnas, --number_of_pegrnas
      The maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 3)
n_ngrnas, --number_of_ngrnas
      The maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3)
-nick_dist_pooled, --nicking_distance_pooled
      The nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority (PE3b seed -> PE3b non-seed), followed by nicking distance closest to this parameter. (Default: 75 bp)
-homology_downstream, --homology_downstream
      For pooled designs (genome_wide or saturation_mutagenesis needs to be indicated), this parameter determines the RT extension length downstream of an edit for pegRNA designs. (Default: 10)
-pbs_pooled, --pbs_length_pooled
      The PBS length to design pegRNAs for pooled design applications. (Default: 14 nt)
rtt_pooled, --rtt_max_length_pooled
      The maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt)
-out, --out_dir
      Name of output directory. (Default: ./DATETIMESTAMP_PrimeDesign)
```
## PrimeDesign input sequence format

PrimeDesign takes in a single input sequence to design pegRNAs and ngRNAs for prime editing. The input sequence encodes both the reference and edited sequence with the following formatting:

* Substitution:     (reference/edit)
* Insertion:        (+insertion) or (/insertion)
* Deletion:         (-deletion) or (deletion/)

### Examples
**Reference sequence:** GCCTGTGACTAACTGCGCCAAAACGTCTTCCAATCCCCTTATCCAATTTA

**Substitution edit sequence:** GCCTGTGACTAACTGCGCCAAAACG**A**CTTCCAATCCCCTTATCCAATTTA<br/>
**PrimeDesign input sequence:** GCCTGTGACTAACTGCGCCAAAACG(T/A)CTTCCAATCCCCTTATCCAATTTA

**Insertion edit sequence:** GCCTGTGACTAACTGCGCCAAAACGT**CTT**CTTCCAATCCCCTTATCCAATTTA<br/>
**PrimeDesign input sequence:** GCCTGTGACTAACTGCGCCAAAACGT(+CTT)CTTCCAATCCCCTTATCCAATTTA

**Deletion edit sequence:** GCCTGTGACTAACTGCGCCAAAAC----TCCAATCCCCTTATCCAATTTA<br/>
**PrimeDesign input sequence:** GCCTGTGACTAACTGCGCCAAAAC(-GTCT)TCCAATCCCCTTATCCAATTTA

**Combinatorial edit sequence:** GCCTGTGACTAACTGC**T**CCA**ATCG**AAACGTC----AATCCCCTTATCCAATTTA<br/>
**PrimeDesign input sequence:** GCCTGTGACTAACTGC(G/T)CCA(+ATCG)AAACGTC(-TTCC)AATCCCCTTATCCAATTTA
