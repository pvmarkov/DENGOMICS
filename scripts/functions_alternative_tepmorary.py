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

def determin_codon_posit(annotation, site_value):
    start = annotation.start
    end = annotation.end
    position1=range (start, end, 3)
    position2=range (start+1, end, 3)
    position3=range (start+2, end, 3)
    return [site_value [i] for i in position1], [site_value [i] for i in position2], [site_value [i] for i in position3]

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
    annotation = Annotation (name, int(start) + correction, int (end) + correction)
    return annotation
    
    

def getting_cover_ntfreqs(bamfile, csv_outputfile):

	samfile = pysam.AlignmentFile(bamfile, "rb" ) # open a file handle for the bam file under the name samfile
	coverage = list ()
	position = list ()
	As = list ()
	Cs = list ()
	Ts = list ()
	Gs = list ()
	Ns = list ()
	majorbases = list ()
	majorsequence = ""
	for pileupcolumn in samfile.pileup('NGC_virus', 0, samfile.lengths[0], max_depth = 2000000):
		coverage.append (pileupcolumn.n)
		position.append (pileupcolumn.pos)
		A=C=T=G=N=0 # shortcut for A=0, C=0 etc.
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.mapping_quality>=30:  # query position is None if is_del or is_refskip is set.
				if pileupread.alignment.query_sequence[pileupread.query_position] == "A":
					A=A+1
				elif pileupread.alignment.query_sequence[pileupread.query_position] == "C":
					C=C+1
				elif pileupread.alignment.query_sequence[pileupread.query_position] == "T":
					T=T+1
				elif pileupread.alignment.query_sequence[pileupread.query_position] == "G":
					G=G+1
				else:
					N=N+1
		As.append (A)
		Cs.append (C)
		Ts.append (T)
		Gs.append (G)
		Ns.append (N)
		majorbases.append (max(A, C, T, G))
		if A == max (A, C, G, T):
			majorbaseid = "a"
		elif C == max (A, C, G, T):
			majorbaseid = "c"
		elif G == max (A, C, G, T):
			majorbaseid = "g"
		else:
			majorbaseid = "t"
		majorsequence = majorsequence + majorbaseid


	print (majorsequence[19:38])

	print (majorsequence[0:38])

	counts_dataframe = pd.DataFrame ({'position': position, 'coverage':coverage, 'As':As, 'Cs':Cs, 'Ts':Ts, 'Gs':Gs, 'Ns':Ns, 'majorbases':majorbases, 'majorsequence': list (majorsequence)})

	counts_dataframe.to_csv (csv_outputfile)
	print (counts_dataframe.describe ())
	samfile.close()



def getting_position_correction (reference_genome, majorsequence):
	majorseq_sample20 = majorsequence [19:39]
	refgenome = extracts_seqs_from_genbankformat (reference_genome)
	print (len (refgenome))
	return refgenome.find (majorseq_sample20) - 19
	
	
