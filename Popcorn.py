import sys, os, math, pickle, gzip, pathlib
import warnings; warnings.filterwarnings('ignore')


####################################
##########   PARAMETERS   ##########
####################################

SHORT_LENGTH = 150  # Definition of short ORF - CDS with at most this length in NTs
MIN_CODING_LENGTH = 6  # Sequences shorter than this are penalized
VALID_NT_SEQUENCE_THRESHOLD = 0.75  # Sequences must contain at least this percentage of NTs
CATEGORY = {0:'NONCODING', 1:'CODING'}
MODEL_FILE_GENOME = os.path.join('DataFiles', 'model.pickle')
MODEL_FILE = os.path.join('DataFiles', 'model.all.pickle')
STATS_FILE_GENOME = 'genome_stats.pickle'
STATS_FILE = os.path.join('DataFiles', 'genome.all.pickle')
CUSTOM_STATS_FILE = False
SEQUENCES, SEQUENCE_NAMES, GENOME_DIR, GENOME_FILE, GENES_FILE = [], [], '', '', ''
OUTPUT_FILE = sys.stdout


################################
##########   CODONS   ##########
################################

START_CODONS, STOP_CODONS = {'ATG':True, 'GTG':True}, {'TAA':True, 'TAG':True, 'TGA':True}
CODON_TO_AA = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
AA_TO_CODON = {'K': {'AAA': True, 'AAG': True}, 'N': {'AAC': True, 'AAT': True}, 'T': {'ACA': True, 'ACC': True, 'ACG': True, 'ACT': True}, 'R': {'AGA': True, 'AGG': True, 'CGA': True, 'CGC': True, 'CGG': True, 'CGT': True}, 'S': {'AGC': True, 'AGT': True, 'TCA': True, 'TCC': True, 'TCG': True, 'TCT': True}, 'I': {'ATA': True, 'ATC': True, 'ATT': True}, 'M': {'ATG': True}, 'Q': {'CAA': True, 'CAG': True}, 'H': {'CAC': True, 'CAT': True}, 'P': {'CCA': True, 'CCC': True, 'CCG': True, 'CCT': True}, 'L': {'CTA': True, 'CTC': True, 'CTG': True, 'CTT': True, 'TTA': True, 'TTG': True}, 'E': {'GAA': True, 'GAG': True}, 'D': {'GAC': True, 'GAT': True}, 'A': {'GCA': True, 'GCC': True, 'GCG': True, 'GCT': True}, 'G': {'GGA': True, 'GGC': True, 'GGG': True, 'GGT': True}, 'V': {'GTA': True, 'GTC': True, 'GTG': True, 'GTT': True}, '*': {'TAA': True, 'TAG': True, 'TGA': True}, 'Y': {'TAC': True, 'TAT': True}, 'C': {'TGC': True, 'TGT': True}, 'W': {'TGG': True}, 'F': {'TTC': True, 'TTT': True}}
AAs = dict.fromkeys(list('ACDEFGHIKLMNPQRSTVWY'), True)


###############################
##########   USAGE   ##########
###############################

def command():
        printUsage = False
        for arg in sys.argv:
                if (arg.lower() == '-h') or (arg.lower() == '-help') or (arg.lower() == '--help'): printUsage = True
        if (len(sys.argv) < 2) or (printUsage):
                sys.stderr.write("\nPOPCORN: PrOkaryotic Prediction of Coding OR Noncoding\n")
                sys.stderr.write("Version 1.0\n")
                sys.stderr.write("Popcorn predicts whether genomic sequences are coding or noncoding\n\n")
                sys.stderr.write("EXAMPLE USAGE:   python Popcorn.py -s ACGTACGTACGT\n")
                sys.stderr.write("EXAMPLE USAGE:   python Popcorn.py -f *.fa\n")
                sys.stderr.write("EXAMPLE USAGE:   python Popcorn.py -s ACGTACGTACGT -g genome_dir\n")
                sys.stderr.write("EXAMPLE USAGE:   python Popcorn.py -f *.fa -g genome_dir\n")
                sys.stderr.write("\n*****   Required argument   *****\n\n")
                sys.stderr.write("\t-s STRING\tGenomic sequence, e.g., ACGTACGTACGT\n")
                sys.stderr.write("\t\t\t\tEither -s or -f flag is required but not both\n")
                sys.stderr.write("\t\t\t\tUse -s to make a prediction for a single\n")
                sys.stderr.write("\t\t\t\tsequence provided on the command line\n")
                sys.stderr.write("\t-f STRING\tFile of genomic sequences either in FASTA format\n")
                sys.stderr.write("\t\t\t\tor with each sequence separated by a blank line\n")
                sys.stderr.write("\t\t\t\tEither -s or -f flag is required but not both\n")
                sys.stderr.write("\t\t\t\tUse -f to make predictions for one or more\n")
                sys.stderr.write("\t\t\t\tsequences in a provided file\n")
                sys.stderr.write("\n*****   Optional argument (RECOMMENDED)  *****\n\n")
                sys.stderr.write("\t-g STRING\tPath to directory containing two files:\n")
                sys.stderr.write("\t\t\t\t- *.fna (genome in FASTA format)\n")
                sys.stderr.write("\t\t\t\t- *.gff (list of genes in GFF format)\n")
                sys.stderr.write("\t\t\t\tFiles may be gzipped or not\n")
                sys.stderr.write("\n*****   Optional arguments  *****\n\n")
                sys.stderr.write("\t-o STRING\tFile to which results should be output\n")
                sys.stderr.write("\t\t\t\t(default is standard out)\n")
                sys.stderr.write("\t-m STRING\tFile containing the trained ML model\n")
                sys.stderr.write("\t\t\t\t(default is DataFiles" + os.path.sep + "model.pickle with -g flag)\n")
                sys.stderr.write("\t\t\t\t(default is DataFiles" + os.path.sep + "model.all.pickle without -g flag)\n")
                sys.stderr.write("\t-z STRING\tFile containing genome statistics when -g flag is not used\n")
                sys.stderr.write("\t\t\t\t(default is DataFiles" + os.path.sep + "genome.all.pickle)\n")
                sys.stderr.write("\t-h\t\tprint USAGE and DESCRIPTION, ignore all other flags\n")
                sys.stderr.write("\t-help\t\tprint USAGE and DESCRIPTION, ignore all other flags\n")
                sys.stderr.write("\n")
                sys.exit(1)


def arguments():
        global SEQUENCES, SEQUENCE_NAMES, GENOME_DIR, STATS_FILE, CUSTOM_STATS_FILE, GENOME_FILE, GENES_FILE, MODEL_FILE, OUTPUT_FILE
        sequence, sequence_file = '', ''
        model_file_flag = False
        for i in range(1, len(sys.argv)):
                if (sys.argv[i] == '-s'): sequence = sys.argv[i+1]
                elif (sys.argv[i] == '-f'): sequence_file = sys.argv[i+1]
                elif (sys.argv[i] == '-g'): GENOME_DIR = sys.argv[i+1]
                elif (sys.argv[i] == '-o'): OUTPUT_FILE = open(sys.argv[i+1], 'w')
                elif (sys.argv[i] == '-m'): MODEL_FILE = sys.argv[i+1]; model_file_flag = True
                elif (sys.argv[i] == '-z'): STATS_FILE = sys.argv[i+1]; CUSTOM_STATS_FILE = True
        if (len(sequence) > 0): SEQUENCES, SEQUENCE_NAMES = read_in_sequence(sequence)
        elif (len(sequence_file) > 0): SEQUENCES, SEQUENCE_NAMES = read_in_sequence_file(sequence_file)
        else:
                sys.stderr.write('\n' + 'Error - as a command line argument, either the -s flag is required followed by a genomic sequence or the -f flag is required followed by the name of a file containing one or more genomic sequences' + '\n')
                sys.stderr.write('Please try again with a different command line argument. Thanks!' + '\n\n')
                sys.exit(1)
        if (len(GENOME_DIR) > 0):  # Specific genome
                if (not pathlib.Path(GENOME_DIR).is_dir()):
                        sys.stderr.write('\n' + 'Error - the flag -g should be followed by the path to a directory containing genome files but this does not appear to be a valid directory: ' + GENOME_DIR + '\n')
                        sys.stderr.write('Please try again with a different command line argument. Thanks!' + '\n\n')
                        sys.exit(1)
                if (not CUSTOM_STATS_FILE) and (pathlib.Path(os.path.join(GENOME_DIR, STATS_FILE_GENOME)).is_file()):  # Already calculated stats
                        STATS_FILE = os.path.join(GENOME_DIR, STATS_FILE_GENOME)
                else:  # Need genome/genes files to calculate stats
                        filelist = os.listdir(GENOME_DIR)
                        for f in filelist:
                                if (f.lower().endswith('.fna') or f.lower().endswith('.fna.gz')) and ('_rna_from_' not in f.lower()): GENOME_FILE = f
                                if (f.lower().endswith('.gff') or f.lower().endswith('.gff.gz')): GENES_FILE = f
                        if (len(GENOME_FILE) == 0):
                                sys.stderr.write('\n' + 'Error - unable to locate a FASTA formatted genome file *.fna or *.fna.gz in the provided directory: ' + GENOME_DIR + '\n\n')
                                sys.exit(1)
                        if (len(GENES_FILE) == 0):
                                sys.stderr.write('\n' + 'Error - unable to locate a GFF formatted file *.gff or *.gff.gz containing gene information in the provided directory: ' + GENOME_DIR + '\n\n')
                                sys.exit(1)
                if (not model_file_flag): MODEL_FILE = MODEL_FILE_GENOME  # Use genome specific ML model
        if (not pathlib.Path(MODEL_FILE).is_file()):
                sys.stderr.write('\n' + 'Error - could not locate the pickle file containing the trained ML model: ' + MODEL_FILE + '\n')
                sys.stderr.write('The pickled model file can be downloaded from:  https://github.com/btjaden/Popcorn' + '\n\n')
                sys.exit(1)


##############################
#####   OUTPUT RESULTS   #####
##############################

def pad(s):
        if (len(s) < 11): return s + ' '*(11-len(s))
        return s


def output_results(out_file, seq_names, probs, preds):
        out_file.write('Sequence ID' + '\t' + 'Coding Probability' + '\t' + 'Prediction' + '\n')
        for i in range(len(probs)):
                out_file.write(pad(seq_names[i]) + '\t' + str(probs[i]) + '\t' + preds[i] + '\n')
        out_file.close()


###################################################
##########   SEQUENCES AND GENOME INFO   ##########
###################################################

def is_valid_sequence(s):
        NT_count = s.count('A') + s.count('C') + s.count('G') + s.count('T')
        if (float(NT_count) / len(s) < VALID_NT_SEQUENCE_THRESHOLD): return False
        return True


def read_in_sequence(sequence):
        s = sequence.upper().replace('U', 'T')
        if (not is_valid_sequence(s)):
                sys.stderr.write('\n' + 'Error - the flag -s should be followed by a valid genomic sequence but this does not appear to be valid: ' + sequence + '\n')
                sys.stderr.write('Please try again with a different command line argument. Thanks!' + '\n\n')
                sys.exit(1)
        return [s], ['Seq_1']


def read_in_sequence_file(sequence_file):
        sequences, sequence_names = [], []
        if (pathlib.Path(sequence_file).is_file()):  # Input is a filename
                with open(sequence_file, 'r') as in_file:
                        line = in_file.readline()
                        while (line == ''): line = in_file.readline()  # Ignore blank header lines
                        if (line.startswith('>')):  # FASTA file
                                s, s_name = '', ''
                                while (line != ''):
                                        if (line.startswith('>')):
                                                if (len(s) > 0):
                                                        sequences.append(s)
                                                        if (s_name == ''): s_name = 'Seq_' + str(len(sequences))
                                                        sequence_names.append(s_name)
                                                        s = ''
                                                s_name = line[1:].strip()
                                        else: s += line.strip()
                                        line = in_file.readline()
                                if (len(s) > 0):
                                        sequences.append(s)
                                        if (s_name == ''): s_name = 'Seq_' + str(len(sequences))
                                        sequence_names.append(s_name)
                        else:  # File where sequences are separated by blank lines
                                s = ''
                                while (line != ''):
                                        if (line.strip() == ''):
                                                if (len(s) > 0):
                                                        sequences.append(s)
                                                        sequence_names.append('Seq_' + str(len(sequences)))
                                                        s = ''
                                        else: s += line.strip()
                                        line = in_file.readline()
                                if (len(s) > 0): sequences.append(s); sequence_names.append('Seq_' + str(len(sequences)))
        else:  # Not a valid file
                sys.stderr.write('\n' + 'Error - the flag -f should be followed by a file containing genomic sequences but this does not appear to be a valid file: ' + sequence_file + '\n')
                sys.stderr.write('Please try again with a different command line argument. Thanks!' + '\n\n')
                sys.exit(1)
        for i in range(len(sequences)):
                s = sequences[i].upper().replace('U', 'T')
                sequences[i] = s
                if (not is_valid_sequence(s)): sys.stderr.write('Warning - in the provided file (' + sequence_file + ') the following sequence does not appear to be a valid genomic sequence:' + '\n' + s + '\n\n')
        return sequences, sequence_names


def read_in_genome(filename):
        genome = []
        name, sequence = '', ''
        if (filename.lower().endswith('.gz')): in_file = gzip.open(filename, 'rt')
        else: in_file = open(filename, 'r')
        line = in_file.readline()
        while (line != ''):
                if (line.startswith('>')):
                        if (len(name) > 0): genome.append((name, sequence))
                        name = line.strip()[1:]
                        sequence = ''
                else: sequence += line.strip()
                line = in_file.readline()
        if (len(name) > 0): genome.append((name, sequence))
        in_file.close()
        return genome


def read_in_genes(filename):
        genes, CDS, sORFs, ncRNAs = [], [], [], []
        if (filename.lower().endswith('.gz')): in_file = gzip.open(filename, 'rt')
        else: in_file = open(filename, 'r')
        line = in_file.readline()
        while (line != ''):
                if (not line.startswith('#')):  # Ignore comments
                        parse_line = line.strip().split('\t')
                        if (parse_line[2] != 'region') and (parse_line[2] != 'gene'):
                                replicon_ID, gene_type, start, stop, strand, gene_info = parse_line[0], parse_line[2], int(parse_line[3]), int(parse_line[4]), parse_line[6], parse_line[-1]
                                length = stop - start + 1
                                parse_gene_info = gene_info.split(';')
                                gene_name, gene_ID = '', ''
                                for p in parse_gene_info:
                                        if (p.startswith('gene=')): gene_name = p[5:]
                                        if (p.startswith('locus_tag=')): gene_ID = p[10:]
                                genes.append((replicon_ID, gene_type, start, stop, strand, gene_ID + ':::' + gene_name))
                                if (gene_type == 'CDS'):  # CDS and sORF
                                        if (length <= SHORT_LENGTH): sORFs.append((replicon_ID, 'sORF', start, stop, strand, gene_ID + ':::' + gene_name))
                                        else: CDS.append((replicon_ID, gene_type, start, stop, strand, gene_ID + ':::' + gene_name))
                                elif (gene_type == 'ncRNA'):  # ncRNA
                                        ncRNAs.append((replicon_ID, 'ncRNA', start, stop, strand, gene_ID + ':::' + gene_name))
                line = in_file.readline()
        in_file.close()
        return genes, CDS, sORFs, ncRNAs


def get_IGs(genome_dict, genes):
        IGs, coords = [], {}
        for ID, sequence in genome_dict.items(): coords[ID] = [''] * len(sequence)
        for g in genes:
                replicon_ID, gene_type, start, stop, strand, name = g
                for i in range(start-1, min(stop, len(coords[replicon_ID]))): coords[replicon_ID][i] = name
        for ID in coords:
                i, start, stop, previous_gene = 0, -1, -1, '?'
                while (i < len(coords[ID])):
                        if (coords[ID][i] == ''):  # An IG - not annotated gene
                                if (start == -1): start = i
                                stop = i
                        else:  # An annotated gene - not an IG
                                if (start >= 0): IGs.append((ID, 'IG', start+1, stop+1, '?', 'IG...' + previous_gene + '...' + coords[ID][i]))
                                start, previous_gene = -1, coords[ID][i]
                        i += 1
                if (start >= 0): IGs.append((ID, 'IG', start+1, stop+1, '?', 'IG...' + previous_gene + '...?'))
        return IGs


def get_genome_info():
        genome = read_in_genome(os.path.join(GENOME_DIR, GENOME_FILE))
        genes, CDS, sORFs, ncRNAs = read_in_genes(os.path.join(GENOME_DIR, GENES_FILE))

        # Map replicon ID (key) to its genomic sequence (value)
        genome_dict = {}
        for g in genome:
                ID = g[0].split()[0]
                genome_dict[ID] = g[1]

        IGs = get_IGs(genome_dict, genes)  # Determine Intergenic regions
        return genome_dict, genes, CDS, sORFs, ncRNAs, IGs


def reverse(s):
        return s[::-1]


def complement(s):
        s = s.replace('C', '?')
        s = s.replace('G', 'C')
        s = s.replace('?', 'G')
        s = s.replace('T', '?')
        s = s.replace('A', 'T')
        s = s.replace('?', 'A')
        return s


def reverse_complement(s):
        return reverse(complement(s))


def find_longest_ORF(seq, check_both_strands=True):  # Returns the longest ORF in a sequence
        ORFs = []
        strands = ['+', '-']
        if (not check_both_strands): strands = ['+']
        for strand in strands:
                s = seq
                if (strand == '-'): s = reverse_complement(seq)
                for frame in range(0, 3):  # For each of the three reading frames
                        # Find next stop codon
                        j, old_j = frame, frame-3
                        while (j < len(s)-2):
                                if (s[j:j+3] in STOP_CODONS):
                                        # Find start codon for the current stop codon
                                        i = old_j + 3
                                        while (s[i:i+3] not in START_CODONS) and (i < j): i += 3
                                        if (i < j): ORFs.append((i, j+2, strand))
                                        old_j = j
                                j += 3

        # Get longest ORF
        longest_length, longest_index = -1, -1
        for idx, ORF in enumerate(ORFs):
                length = ORF[1] - ORF[0] +1
                if (length > longest_length): longest_length, longest_index = length, idx
        if (longest_index >= 0):
                start, stop, strand = ORFs[longest_index][0], ORFs[longest_index][1], ORFs[longest_index][2]
                if (strand == '+'): return start, stop, '+', seq[start:stop+1]
                else: return start, stop, '-', reverse_complement(seq)[start:stop+1]
        else: return -1, -1, '+', ''


###################################################################
##########   STATS USED FOR CALCULATING FEATURE VALUES   ##########
###################################################################

def get_GC_content(genome_dict):
        GCs, total = 0, 0
        for ID in genome_dict:
                GCs += genome_dict[ID].count('G') + genome_dict[ID].count('C')
                total += len(genome_dict[ID])
        return float(GCs) / total


def get_max_synonymous_codon_count(codon, codon_counts):
        aa = CODON_TO_AA[codon]
        synonymous_codons = AA_TO_CODON[aa]
        max_codon, max_count = '', 0
        for c in synonymous_codons:
                if (c in codon_counts) and (codon_counts[c] > max_count):
                        max_count = codon_counts[c]
                        max_codon = c
        return max_count


def get_codon_weights(genome_dict, genes):
        codon_counts = {}
        for g in genes:
                replicon_ID, gene_type, start, stop, strand, name = g
                if (gene_type == 'CDS'):
                        sequence = genome_dict[replicon_ID][start-1:stop]
                        if (strand == '-'): sequence = reverse_complement(sequence)
                        for i in range(0, len(sequence)-2, 3):
                                codon = sequence[i:i+3]
                                if (codon not in codon_counts): codon_counts[codon] = 0
                                codon_counts[codon] += 1
        codon_weights = {}
        for codon in codon_counts:
                if (codon in CODON_TO_AA): codon_weights[codon] = float(codon_counts[codon]) / get_max_synonymous_codon_count(codon, codon_counts)
        return codon_weights


def get_hexamer_frequencies(genome_dict, genes, IGs):
        # Get hexamer usage for coding regions
        hexamers = {'coding':{}, 'noncoding':{}}
        for g in genes:
                replicon_ID, gene_type, start, stop, strand, name = g
                if (gene_type == 'CDS'):
                        sequence = genome_dict[replicon_ID][start-1:stop]
                        if (strand == '-'): sequence = reverse_complement(sequence)
                        for i in range(0, len(sequence)-5, 3):
                                hexamer = sequence[i:i+6]
                                if (hexamer not in hexamers['coding']): hexamers['coding'][hexamer] = 0
                                hexamers['coding'][hexamer] += 1
        total_hexamers = 0
        for k,v in hexamers['coding'].items(): total_hexamers += v
        for k,v in hexamers['coding'].items(): hexamers['coding'][k] /= float(total_hexamers)

        # Get hexamer usage for noncoding regions
        for ig in IGs:
                replicon_ID, gene_type, start, stop, strand, name = ig
                if (gene_type == 'IG'):  # Should all be IG
                        sequence = genome_dict[replicon_ID][start-1:stop]
                        ORF_start, ORF_stop, ORF_strand, sequence = find_longest_ORF(sequence)
                        for i in range(0, len(sequence)-5, 3):
                                hexamer = sequence[i:i+6]
                                if (hexamer not in hexamers['noncoding']): hexamers['noncoding'][hexamer] = 0
                                hexamers['noncoding'][hexamer] += 1
        total_hexamers = 0
        for k,v in hexamers['noncoding'].items(): total_hexamers += v
        for k,v in hexamers['noncoding'].items(): hexamers['noncoding'][k] /= float(total_hexamers)
        return hexamers


# Used when STATS file already exists
def get_stats_existing(stats_filename):
        if (stats_filename.endswith('.gz')): f = gzip.open(stats_filename, 'rb')
        else: f = open(stats_filename, 'rb')
        try: GC_content, codon_weights, hexamers = pickle.load(f)
        except:
                f.close()
                sys.stderr.write('\n' + 'Error - could not read in the stats file ' + stats_filename + '\n')
                sys.stderr.write('The pickled stats file can be downloaded from:  https://github.com/btjaden/Popcorn' + '\n')
                sys.stderr.write('If the problem persists, one possible cause may be that the version of Python and its packages that you are using are not compatible with this file' + '\n\n')
                sys.exit(1)
        f.close()
        return GC_content, codon_weights, hexamers


# Used when STATS may need to be calculated
def get_stats():
        # Need to calculate stats file if genome is specified and
        # if stats file hasn't been created previously
        if (not CUSTOM_STATS_FILE) and (len(GENOME_DIR) > 0) and (not STATS_FILE.endswith(STATS_FILE_GENOME)):
                genome_dict, genes, CDS, sORFs, ncRNAs, IGs = get_genome_info()
                GC_content = get_GC_content(genome_dict)
                codon_weights = get_codon_weights(genome_dict, genes)
                hexamers = get_hexamer_frequencies(genome_dict, genes, IGs)
                try:  # Try writing the genome stats file for future uses
                        with open(os.path.join(GENOME_DIR, STATS_FILE_GENOME), 'wb') as out_file:
                                pickle.dump((GC_content, codon_weights, hexamers), out_file)
                except: None  # Unable to write stats file for genome
        else:  # Read in existing stats file
                if (not pathlib.Path(STATS_FILE).is_file()):
                        sys.stderr.write('\n' + 'Error - could not locate the pickle file containing genome statistics: ' + STATS_FILE + '\n')
                        sys.stderr.write('The pickled stats file can be downloaded from:  https://github.com/btjaden/Popcorn' + '\n\n')
                        sys.exit(1)
                GC_content, codon_weights, hexamers = get_stats_existing(STATS_FILE)
        return GC_content, codon_weights, hexamers


########################################
##########   FEATURE VALUES   ##########
########################################

def GC_score(s, GC_content):
        if (len(s) == 0): return 0.0
        return ((s.count('C') + s.count('G')) / float(len(s))) - GC_content


def CAI(s, codon_weights):
        if (len(s) < 3): return 0.0
        cai = 0.0
        for i in range(0, len(s)-2, 3):
                if (s[i:i+3] in codon_weights): cai += math.log(codon_weights[s[i:i+3]])
        cai /= len(s)//3
        cai = math.exp(cai)
        return cai


# FICKETT parameters
WEIGHTS_POSITION = {'A':0.26, 'C':0.18, 'G':0.31, 'T':0.33}
WEIGHTS_CONTENT = {'A':0.11, 'C':0.12, 'G':0.15, 'T':0.14}
POS_PROB_A = [0.22] * 11 + [0.20, 0.34, 0.45, 0.68, 0.58, 0.93, 0.84, 0.68, 0.94]
POS_PROB_C = [0.23] * 11 + [0.30, 0.33, 0.51, 0.48, 0.66, 0.81, 0.70, 0.70, 0.80]
POS_PROB_G = [0.08] * 11 + [0.08, 0.16, 0.27, 0.48, 0.53, 0.64, 0.74, 0.88, 0.90]
POS_PROB_T = [0.09] * 11 + [0.09, 0.20, 0.54, 0.44, 0.69, 0.68, 0.91, 0.97, 0.97]
CON_PROB_A = [0.21] * 17 + [0.81, 0.81, 0.65, 0.65, 0.67, 0.67, 0.49, 0.49, 0.62, 0.62, 0.55, 0.55, 0.44, 0.44, 0.49, 0.49, 0.28, 0.28]
CON_PROB_C = [0.31] * 17 + [0.39, 0.39, 0.44, 0.44, 0.43, 0.43, 0.59, 0.59, 0.59, 0.59, 0.64, 0.64, 0.51, 0.51, 0.64, 0.64, 0.82, 0.82]
CON_PROB_G = [0.29] * 17 + [0.33, 0.33, 0.41, 0.41, 0.41, 0.41, 0.73, 0.73, 0.64, 0.64, 0.64, 0.64, 0.47, 0.47, 0.54, 0.54, 0.40, 0.40]
CON_PROB_T = [0.58] * 17 + [0.51, 0.51, 0.69, 0.69, 0.56, 0.56, 0.75, 0.75, 0.55, 0.55, 0.40, 0.40, 0.39, 0.39, 0.24, 0.24, 0.28, 0.28]
PROB_POSITION = {'A':POS_PROB_A, 'C':POS_PROB_C, 'G':POS_PROB_G, 'T':POS_PROB_T}
PROB_CONTENT = {'A':CON_PROB_A, 'C':CON_PROB_C, 'G':CON_PROB_G, 'T':CON_PROB_T}

def fickett(s):
        counts = {'A':[0,0,0], 'C':[0,0,0], 'G':[0,0,0], 'T':[0,0,0]}
        for start in range(3):
                for i in range(start, len(s), 3):
                        if (s[i] in counts): counts[s[i]][start] += 1
        positions, contents = {}, {}
        for NT in counts:
                positions[NT] = max(counts[NT]) / float(min(counts[NT]) + 1)
                if (len(s) == 0): contents[NT] = 0.0
                else: contents[NT] = sum(counts[NT]) / float(len(s))
        TESTCODE = 0.0
        for NT in counts:
                weight1, weight2 = WEIGHTS_POSITION[NT], WEIGHTS_CONTENT[NT]
                if (positions[NT] >= 1.9): prob1 = PROB_POSITION[NT][-1]
                else: prob1 = PROB_POSITION[NT][int(positions[NT] * 10)]
                if (contents[NT] >= 0.33): prob2 = PROB_CONTENT[NT][-1]
                else: prob2 = PROB_CONTENT[NT][int(contents[NT] * 100)]
                TESTCODE += (prob1 * weight1) + (prob2 * weight2)
        return TESTCODE


def hexamer_score(s, hexamers):
        hex_score = 0.0
        count = 0
        for i in range(0, len(s)-5, 3):
                hexamer = s[i:i+6]
                if (hexamer in hexamers['coding']) and (hexamer in hexamers['noncoding']):
                        hex_score += hexamers['coding'][hexamer] / hexamers['noncoding'][hexamer]
                        count += 1
        if (count == 0): return 0.0
        hex_score /= count
        return hex_score


# Isoelectric Point parameters
positive_pKs = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
negative_pKs = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
pKcterminal = {"D": 4.55, "E": 4.75}
pKnterminal = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7}
charged_aas = ("K", "R", "H", "D", "E", "C", "Y")

def translate(s):
        aa = ''
        for i in range(0, len(s)-2, 3):
                if (s[i:i+3] in CODON_TO_AA): aa += CODON_TO_AA[s[i:i+3]]
                else: aa += '?'
        return aa

def pH_charge(pos_pKs, neg_pKs, charged, pH):
        pos_charge = 0.0
        for aa, pK in pos_pKs.items():
                partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
                pos_charge += charged[aa] * partial_charge
        neg_charge = 0.0
        for aa, pK in neg_pKs.items():
                partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
                neg_charge += charged[aa] * partial_charge
        return pos_charge - neg_charge

def iso_point(pos_pKs, neg_pKs, charged, pH, min_value, max_value):
        charge = pH_charge(pos_pKs, neg_pKs, charged, pH)
        if ((max_value - min_value) > 0.0001):
                if (charge > 0.0): min_value = pH
                else: max_value = pH
                return iso_point(pos_pKs, neg_pKs, charged, (min_value + max_value) / 2, min_value, max_value)
        return pH

def IP(s):
        if (len(s) < 3): return 12.0
        s = translate(s)
        aa_content = {}
        for aa in AAs: aa_content[aa] = 0
        for aa in s:
                if (aa in AAs): aa_content[aa] += 1
        charged = {'Nterm':1.0, 'Cterm':1.0}
        for aa in charged_aas: charged[aa] = float(aa_content[aa])
        pos_pKs = positive_pKs.copy()
        neg_pKs = negative_pKs.copy()
        nterm, cterm = s[0], s[-1]
        if nterm in pKnterminal: pos_pKs["Nterm"] = pKnterminal[nterm]
        if cterm in pKcterminal: neg_pKs["Cterm"] = pKcterminal[cterm]
        return iso_point(pos_pKs, neg_pKs, charged, 7.775, 4.05, 12)


def get_feature_values(s, GC_content, codon_weights, hexamers):
        start, stop, strand, ORF_seq = find_longest_ORF(s, False)
        x = []
        x.append(GC_score(s, GC_content))
        x.append(CAI(s, codon_weights))
        x.append(fickett(s))
        x.append(hexamer_score(s, hexamers))
        x.append(IP(s))
        x.append(len(ORF_seq))
        if (len(s) == 0): x.append(0.0)
        else: x.append(float(len(ORF_seq)) / len(s))
        x.append(GC_score(ORF_seq, GC_content))
        x.append(CAI(ORF_seq, codon_weights))
        x.append(fickett(ORF_seq))
        x.append(hexamer_score(ORF_seq, hexamers))
        x.append(IP(ORF_seq))
        if (len(s) < MIN_CODING_LENGTH):
                for i in range(4): x[i] /= 2
                x[5] *= 1.5
        return [x]


def run_Popcorn(seqs, model_filename=None, stats_filename=None):
        probabilities, predictions = [], []
        if (model_filename is None): model_filename = MODEL_FILE
        if (stats_filename is not None): GC_content, codon_weights, hexamers = get_stats_existing(stats_filename)
        else: GC_content, codon_weights, hexamers = get_stats()
        try:  # Load model
                with open(model_filename, 'rb') as f: model, scaler = pickle.load(f)
        except:
                sys.stderr.write('\n' + 'Error - could not read in the pickle file containing the trained ML model: ' + model_filename + '\n')
                sys.stderr.write('The pickled model file can be downloaded from:  https://github.com/btjaden/Popcorn' + '\n\n')
                sys.stderr.write('If the problem persists, one possible cause may be that the version of Python and its packages that you are using are not compatible with this file' + '\n\n')
                sys.exit(1)
        for s in seqs:  # Make coding/noncoding prediction for each sequence
                x = get_feature_values(s, GC_content, codon_weights, hexamers)
                x_scaled = scaler.transform(x)
                y_pred = model.predict(x_scaled)
                y_pred_proba = model.predict_proba(x_scaled)
                probabilities.append(y_pred_proba[0][1])
                predictions.append(CATEGORY[y_pred[0]])
        return probabilities, predictions


##############################
##########   MAIN   ##########
##############################

if __name__ == "__main__":
        command()
        arguments()
        probabilities, predictions = run_Popcorn(SEQUENCES)
        output_results(OUTPUT_FILE, SEQUENCE_NAMES, probabilities, predictions)

