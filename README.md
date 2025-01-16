# <img src="https://cs.wellesley.edu/~btjaden/Popcorn/img/Popcorn_trans.png" width=100> [Popcorn](https://cs.wellesley.edu/~btjaden/Popcorn) <img src="https://cs.wellesley.edu/~btjaden/Popcorn/img/Popcorn_trans.png" width=100>
==========

### [Popcorn](https://cs.wellesley.edu/~btjaden/Popcorn) is a tool for predicting whether ***prokaryotic*** genomic sequences are ***coding*** or ***noncoding***.

#### POPCORN: PrOkaryotic Prediction of Coding OR Noncoding<BR>

EXAMPLE USAGE: &nbsp;&nbsp;&nbsp;&nbsp;`python Popcorn.py -s ACGTACGTACGT`<BR>
EXAMPLE USAGE: &nbsp;&nbsp;&nbsp;&nbsp;`python Popcorn.py -f *.fa`<BR>
EXAMPLE USAGE: &nbsp;&nbsp;&nbsp;&nbsp;`python Popcorn.py -s ACGTACGTACGT -g genome_dir`<BR>
EXAMPLE USAGE: &nbsp;&nbsp;&nbsp;&nbsp;`python Popcorn.py -f *.fa -g genome_dir`<BR>

*****   Required argument   *****

	-s STRING	Genomic sequence, e.g., ACGTACGTACGT
				Either -s or -f flag is required but not both
				Use -s to make a prediction for a single
				sequence provided on the command line
	-f STRING	File of genomic sequences either in FASTA format
				or with each sequence separated by a blank line
				Either -s or -f flag is required but not both
				Use -f to make predictions for one or more
				sequences in a provided file

*****   Optional argument (RECOMMENDED)  *****

	-g STRING	Path to directory containing two files:
				- *.fna (genome in FASTA format)
				- *.gff (list of genes in GFF format)
				Files may be gzipped or not

*****   Optional arguments  *****

	-o STRING	File to which results should be output
				(default is standard out)
	-m STRING	File containing the trained ML model
				(default is DataFiles/model.pickle with -g flag)
				(default is DataFiles/model.all.pickle without -g flag)
	-z STRING	File containing genome statistics when -g flag is not used
				(default is DataFiles/genome.all.pickle)
	-h		print USAGE and DESCRIPTION, ignore all other flags
	-help		print USAGE and DESCRIPTION, ignore all other flags
