import numpy as np
import pandas as pd
import pysam 

def extracts_seqs_from_genbankformat (filename):
	try:
		genbankfile = open (filename, 'r')
	except IOError:
		print ("Unknown file " + filename)
		sys.exit()
	interestingline = False
	sequence = ""
	for l in genbankfile:
		if "ORIGIN" in l:
			interestingline = True
		elif "//" in l:
			interestingline = False
		elif interestingline:
			linelist = l.split ()
			linelist.pop (0)
			for e in linelist:	
				sequence = sequence + e
	return sequence


def extracts_seqs_from_fasta (filename):
	try:
		fastafile = open (filename, 'r')
	except IOError:
		print ("Unknown file " + filename)
		sys.exit()
	sequence = ""
	for l in fastafile:
		if ">" in l:
			pass
		else:
			sequence = sequence + l.strip()
	return sequence


def determin_codon_posit(annotation, site_value):
    start = annotation.start
    end = annotation.end
    positionall=range (start, end)
    position1=range (start, end, 3)
    position2=range (start+1, end, 3)
    position3=range (start+2, end, 3)
    return [site_value [i] for i in positionall], [site_value [i] for i in position1], [site_value [i] for i in position2], [site_value [i] for i in position3]

def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]
    
    
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
    
    

def getting_cover_ntfreqs(bamfile, csv_outputfile, sequencingquality_on, alignquality_on, pairedread_on):

	samfile = pysam.AlignmentFile(bamfile, "rb" ) # open a file handle for the bam file under the name samfile
	indels_by_column_list = list ()
	coverage = list () # creating an empty list that will take coverage by position (column)
	position = list ()
	expected_number_of_errors = list ()
	As = list ()
	Cs = list ()
	Ts = list ()
	Gs = list ()
	Ns = list ()
	majorbases = list ()
	majorsequence = ""
	secondbase = list ()
	number_of_indels_dictionary = dict()
	number_of_indels_by_length_dict = dict()
	
	for read in samfile.fetch('NGC_virus', 0, samfile.lengths[0]):
		length = len(read.get_blocks())
		if number_of_indels_dictionary.__contains__(length):
			number_of_indels_dictionary [length] += 1
		else:
			number_of_indels_dictionary [length] = 1

	number_of_indels_dataframe = pd.DataFrame({'number_of_indels_per_read': list (number_of_indels_dictionary.keys()), 'number_of_reads_with_given_indel_number': list (number_of_indels_dictionary.values())})
	#pd.DataFrame.from_dict(number_of_indels_dictionary)
	number_of_indels_dataframe.to_csv (csv_outputfile + "_indels_per_read")
#	print (indels_dataframe.describe ())

# ITERATES COLUMN BY COLUMN	
	for pileupcolumn in samfile.pileup('NGC_virus', 0, samfile.lengths[0], max_depth = 2000000): # iterates over alignment columns
#		coverage.append (pileupcolumn.n) # pileupcolumn.n gives the number of reads mapping to that column. here these numbers are appended to the coverage list as the program loops column by column.
		position.append (pileupcolumn.pos)
		A=C=T=G=N=indels_in_single_column=0 # shortcut for A=0, C=0 etc.
		number_of_errors = 0
# WITHIN COLUMN, ITERATES ROW BY ROW
		for pileupread in pileupcolumn.pileups:
			if (pairedread_on and pileupread.alignment.is_proper_pair) or (not pairedread_on) :
				if (alignquality_on and pileupread.alignment.mapping_quality>=30) or (not alignquality_on):
					if not pileupread.is_del and not pileupread.is_refskip :  # query position is None if is_del or is_refskip is set.
						if (sequencingquality_on and pileupread.alignment.query_qualities[pileupread.query_position]>=30) or (not sequencingquality_on): # skips lower hierarchy loops if sequencing filter is on AND seq quality is low.
							number_of_errors = number_of_errors + pow (10.0, (-float (pileupread.alignment.mapping_quality)/10.0))
							if pileupread.indel != 0:
								indels_in_single_column += 1
								if number_of_indels_by_length_dict.__contains__(pileupread.indel): #checks if we already have a key equal to the current 'pileupread.indel' value, which provides the size of the indel starting at that position.
									number_of_indels_by_length_dict [pileupread.indel] += 1 # if the indel size key above was present, then increments that indel size by one
								else:
									number_of_indels_by_length_dict [pileupread.indel] = 1 # if the indel size key above was absent, then creates it as new key and sets the value to one
							if pileupread.alignment.query_sequence[pileupread.query_position] == "A":
								A=A+1 # counts the number of As at a site
							elif pileupread.alignment.query_sequence[pileupread.query_position] == "C":
								C=C+1
							elif pileupread.alignment.query_sequence[pileupread.query_position] == "T":
								T=T+1
							elif pileupread.alignment.query_sequence[pileupread.query_position] == "G":
								G=G+1
							else:
								N=N+1
		if A+C+T+G+N == 0:	# this if else bit is to avoid 'Division by Zero' error in case coverage is = 0
			coverage.append (1)
		else:
			coverage.append (A+C+T+G+N)
		expected_number_of_errors.append(number_of_errors)
		indels_by_column_list.append(indels_in_single_column)
		As.append(A) # populates the empty As list with the NUMBER of As by genome site.
		Cs.append(C)
		Ts.append(T)
		Gs.append(G)
		Ns.append(N)
		l= [A, C, T, G]
		l.sort()
		majorbases.append (l[3]) # Populates the majorbases list with the numbers of the most frequent nuke by genome site
		secondbase.append (l[2]) # Populates the secondbases list with the numbers of the 2nd most frequent nuke by genome site
		if A == max (A, C, G, T):
			majorbaseid = "a"
		elif C == max (A, C, G, T):
			majorbaseid = "c"
		elif G == max (A, C, G, T):
			majorbaseid = "g"
		else:
			majorbaseid = "t"
		majorsequence = majorsequence + majorbaseid # creates the CONSENSUS SEQUENCE (populates a STRING with ACTGs to identify the most common nuke
	
	def ratio(x,y):
		return x/y

	majorbase_ratio = list(map(ratio, majorbases, coverage))
	secondbase_ratio = list (map (ratio, secondbase, coverage))
	
	probability_of_seq_error = list (map (ratio, expected_number_of_errors, coverage)) # Creates a per-column probability of alignment error.
	
#	print (majorsequence[19:38])

#	print (majorsequence[0:38])

	counts_dataframe = pd.DataFrame ({'position': position, 'coverage':coverage, 'As':As, 'Cs':Cs, 'Ts':Ts, 'Gs':Gs, 'Ns':Ns, 'majorbases':majorbases, 'majorbase_ratio': majorbase_ratio, 'secondbase': secondbase, 'secondbase_ratio': secondbase_ratio, 'majorsequence': list (majorsequence), 'expected_number_of_errors':expected_number_of_errors, 'probability_of_seq_error': probability_of_seq_error, 'number_of_indels_by_position' : indels_by_column_list} )

	counts_dataframe.to_csv (csv_outputfile)
	print (counts_dataframe.describe ())
	samfile.close()

#	indels_dataframe = pd.DataFrame({'number_of_indels_per_read': list (number_of_indels_dictionary.keys()), 'number_of_reads_with_gien_indel_number': list (number_of_indels_dictionary.values())})

	indels_sizes_dataframe = pd.DataFrame ({'indel_size': list (number_of_indels_by_length_dict.keys()), 'number_of_indels_by_size': list (number_of_indels_by_length_dict.values())})
	#.from_dict(number_of_indels_by_length_dict)
	indels_sizes_dataframe.to_csv (csv_outputfile + "_indels_lengths")
	print (number_of_indels_dataframe.describe ())
	print (indels_sizes_dataframe.describe ())




def getting_position_correction (reference_genome, majorsequence):
	majorseq_sample20 = majorsequence [19:39]
	refgenome = extracts_seqs_from_genbankformat (reference_genome)
	print (len (refgenome))
	return refgenome.find (majorseq_sample20) - 19
	
	
