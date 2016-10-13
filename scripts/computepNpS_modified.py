import numpy as np
import pandas as pd
import pysam
import sys

# Functions to compute piNpiS based on reads mapped onto a genome.


# The genetic code dictionary
genetic_code_U = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

genetic_code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
       "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}



# Reading a fasta file
def readFasta (file):
    try:
        f=open(file, 'r')
    except IOError:
        print ("Unknown file: ",file)
        sys.exit()
    seqs = list()
    nams = list()
    seq=""
    for l in f:
        if l.startswith('>'):
            if seq != "" and seqname != "":
                seqs.append(seq)
                nams.append(seqname)
            seq=""
            seqname = l.replace(">","")
        else:
            seq=seq+l.strip()
    # last sequence
    if seq != "" and seqname != "":
        seqs.append(seq)
        nams.append(seqname)
    f.close()
    return (nams, seqs )


class Annotation (object):
    name = ""
    start = 0
    end = 0

    # The class "constructor" - It's actually an initializer
    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end
    def __repr__ (self):
        return self.name + ": " + str (self.start) + " to " + str (self.end)
    def __str__ (self):
        return self.name + ": " + str (self.start) + " to " + str (self.end)


def make_annotation (name, start, end, correction):
    annotation = Annotation (name, (int(start)-1) + correction, (int (end)-1) + correction)
    return annotation


def determin_codon_posit_coordinates(annots_list):
    positionall = list()
    position1 = list()
    position2 = list()
    position3 = list()
    for x in annots_list:
        if x.start>0 and 'UTR' not in x.name:
            start = x.start
            end = x.end
            positionall += range (start, end)
            position1 += range (start, end, 3)
            position2 += range (start+1, end, 3)
            position3 += range (start+2, end, 3)
    return positionall, position1, position2, position3


def read_annotations (file):
    try:
        f= open (file, 'r')
    except IOError:
        print ("Unknown file " + file)
        sys.exit()
    line = ""
    annots_list = list ()
    line_list = list ()
    correction = - 1846 # Correction to remove the promoter and other sequence from the plasmid
    for l in f:
        if ('CDS' in l or 'UTR' in l) and 'DEN2' not in l and 'Beta-lactamase' not in l:
            line_list = l.split()
            annots_list.append (make_annotation (line_list [12], line_list [3], line_list [4], correction))
    f.close()
    annots_list.sort (key=lambda x: x.start)
    print ('The annotations and their positions are: \n',annots_list)
    return annots_list


def update_coverage_piN_piS(coverage, piN, piS, piStop, read_codon, read_codon_quality, reference_codon, genetic_code, position, codon_quality_threshold):
    if len(read_codon) == len(reference_codon) == 3 :
        average_codon_quality = (read_codon_quality[0] + read_codon_quality[1] + read_codon_quality[2]) /3
        #print(average_codon_quality)
        if (average_codon_quality > codon_quality_threshold):
            coverage[position] += 1
            if read_codon != reference_codon:
                if (genetic_code[read_codon] == genetic_code[reference_codon]):
                    piS[position] += 1
                    #print("SYN: "+read_codon + "<>"+ reference_codon)
                else:
                    if genetic_code[read_codon] == "STOP":
                        piStop[position] += 1
                    else:
                        piN[position] += 1
                        #print("NONSYN: "+read_codon + "<>"+ reference_codon)
    return piN, piS, piStop


def compute_pN_pS(bamfile, csv_outputfile, consensus_sequence, position1s, position2s, position3s, codon_quality_threshold):
    samfile = pysam.AlignmentFile(bamfile, "rb" ) # open a file handle for the bam file under the name samfile
    coverage = list () # creating an empty list that will take coverage by position (column)
    position = list ()
    expected_number_of_errors = list ()
    As = list ()
    Cs = list ()
    Ts = list ()
    Gs = list ()
    Ns = list ()
    number_of_codons = int(len(consensus_sequence) / 3)
    piN = [0] * number_of_codons
    piS = [0] * number_of_codons
    piStop = [0] * number_of_codons
    coverage = [0] * number_of_codons
    for read in samfile.fetch('NGC_virus', 0, samfile.lengths[0]):
        if not read.is_unmapped and read.mapping_quality > 30: # by default we filter on mapping quality
            start = read.reference_start
            end = read.reference_end
            length = read.reference_length
            length_codons = int(read.reference_length/3)
            read_sequence = read.query_alignment_sequence
            read_qualities = read.query_alignment_qualities
            reference_sequence = consensus_sequence[start:end]
            if (start in position1s):
                for i in range(length_codons):
                    read_codon = read_sequence[i*3:i*3+3]
                    read_codon_quality = read_qualities[i*3:i*3+3]
                    reference_codon = reference_sequence[i*3:i*3+3]
                    piN, piS, piStop = update_coverage_piN_piS(coverage, piN, piS, piStop, read_codon, read_codon_quality, reference_codon, genetic_code, int(start/3)+i, codon_quality_threshold)
            elif (start in position2s):
                for i in range(length_codons):
                    read_codon = read_sequence[i*3+2:i*3+2+3]
                    read_codon_quality = read_qualities[i*3+2:i*3+2+3]
                    reference_codon = reference_sequence[i*3+2:i*3+2+3]
                    piN, piS, piStop = update_coverage_piN_piS(coverage, piN, piS, piStop, read_codon, read_codon_quality, reference_codon, genetic_code, int(start/3)+i, codon_quality_threshold)
            else:
                for i in range(length_codons):
                    read_codon = read_sequence[i*3+1:i*3+1+3]
                    read_codon_quality = read_qualities[i*3+1:i*3+1+3]
                    reference_codon = reference_sequence[i*3+1:i*3+1+3]
                    piN, piS, piStop = update_coverage_piN_piS(coverage, piN, piS, piStop, read_codon, read_codon_quality, reference_codon, genetic_code, int(start/3)+i, codon_quality_threshold)
    piN_piS_dataframe = pd.DataFrame({'piN':piN, 'piS':piS, "piStop":piStop, 'coverage':coverage})
    piN_piS_dataframe.to_csv( csv_outputfile )
    return(piN_piS_dataframe)



def compute_positionAndConsensusCodon(csv_outputfile, consensus_sequence, position1s, position2s, position3s):
    positions = list ()
    codonPositions = list ()
    start = min(position1s)
    end = max(position3s)
    number_of_codons = int(len(consensus_sequence) / 3)
    reference_codons=list()
    for i in range(start, end, 3):
    #for i in range(number_of_codons):
        reference_codons.append( consensus_sequence[i:i+3])
        positions.append(i)
        codonPositions.append(int(i/3))
    pos_cons_dataframe = pd.DataFrame({'position':positions, 'codonPosition':codonPositions, 'sequence':reference_codons})
    pos_cons_dataframe.to_csv( csv_outputfile, index=False )
    return(pos_cons_dataframe)


if __name__ == '__main__':
    bamfile = "reads_subset.bam" # Bam file
    csv_outputfile = "piNpiS.csv" # output CSV file
    fastaFile = "NGC_virus_reference.fastaOK.txt" # Fasta file containing the reference viral genome
    annotations = 'ref_geno_annotation_relabel.gff' # Annotations, including the protein coding genes
    codon_quality_threshold = 39 # Codoon quality threshold: a codon must have a quality at least this good to be included in the analysis
    name, consensus_sequence = readFasta(fastaFile)
    consensus_sequence[0] = consensus_sequence[0].upper()
    annots_list = read_annotations(annotations)
    positionall, position1, position2, position3 = determin_codon_posit_coordinates(annots_list)
    compute_pN_pS(bamfile, csv_outputfile, consensus_sequence[0], position1, position2, position3, codon_quality_threshold)
    csv_outputfile = "position_and_sequence.csv" # output CSV file
    compute_positionAndConsensusCodon(csv_outputfile, consensus_sequence[0], position1, position2, position3)
