# PrimeDesign

PrimeDesign is a software tool for the flexible design of prime editing.

## Installation with Docker

With Docker, no installation is required - the only dependence is Docker itself. Users will not need to deal with installation and configuration issues. Docker will do all the dirty work for you!

Docker can be downloaded freely here: [https://store.docker.com/search?offering=community&type=edition](https://store.docker.com/search?offering=community&type=edition)

To get a local copy of CRISPR-SURF, simply execute the following command:
* ```docker pull pinellolab/primedesign```

## PrimeDesign web application

Run the PrimeDesign web application locally with the command:

```
docker run -p 9994:9994 pinellolab/primedesign:latest primedesign_webapp
```

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
      List of primer binding site (PBS) lengths for the pegRNA extension (Default: 10 to 16 nt). Example: 12 13 14 15
-rt, --rt_length_list
      List of reverse transcription (RT) template lengths for the pegRNA extension (Default: 10 to 16 nt). Example: 10 15 20
-nick_dist_min, --nicking_distance_minimum
      Minimum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 0).
-nick_dist_max, --nicking_distance_maximum
      Maximum nicking distance for designing ngRNAs upstream and downstream of a pegRNA (Default: 100).
-out, --out_dir
      Name of output directory. (Default: ./DATETIMESTAMP_PrimeDesign)
```
