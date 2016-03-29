import numpy as np
import pandas as pd
import pysam 


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