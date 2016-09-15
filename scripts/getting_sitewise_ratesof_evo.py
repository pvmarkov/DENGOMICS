import numpy as np
import pandas as pd
import pysam 


codon_number = 0
synonymous_substitutions_by_codon = list ()
non_synonymous_substitutions_by_codon = list ()
non_synonymous_substitutions_by_position_in_codon = list ()
codon_numbers = list ()
codon_numbers_second_file = list ()

# open a file for reading
try:
	codon_substititons_file = open ('/Users/pvmarkov/dengue/data/DENV4_complete_human_alignment_viatranslation_renamed_mapping_counts_per_branch_per_site_per_type_1.count', 'r')
except IOError:
	print ("Unknown file " + filename)
	sys.exit()
# parse line by line
for l in codon_substititons_file:
	branch_counter = 0
	substitutions_sum = 0
	if l.find ("sites") > -1: # skipping headings line
		pass
	else:
		branches_list = l.split("\t") # creates a list from each line in turn, using the tab as a separator
		#print (branches_list)
		for branch in branches_list:
			branch_counter = branch_counter + 1
			if branch_counter == 1:
				if int (branch) + 1 < 2299:
					codon_number = int (branch) + 1 - 5 # '+1' is to convert 'pythonian' counting of positions to normal numbers. the '-5' bit is an adjustment of the codon positions in the database alignment in comparison to the reference genome. '-5' is kept separate from '+1' for clarity.
				else:
					codon_number = int (branch) + 1 - 2 # the '-2' now is a change in the adjustment after aa position 2299 where codon positions in the database alignment is now misaligned in comparison to the reference genome by just two positions. '-2' is kept separate from '+1' for clarity.
			else:
				if branch == str ("nan"):
					branch = 0
				substitutions_sum = substitutions_sum + float (branch)
			#print (branch)
		#print (l [0:19])
		#print (codon_number)
		#print (substitutions_sum)
		codon_numbers.append(codon_number)
		synonymous_substitutions_by_codon.append(substitutions_sum) # populates the empty As list with the NUMBER of As by genome site.
#print (codon_positions)				
#print (synonymous_substitutions_by_codon)


# THIS SECOND PART IS IDENTICAL TO ABOVE ONLY USING THE NONGYSNONYMOUS SUBSTITUTIONS SOURCE FILE.
codon_number = 0

# open a file for reading
try:
	codon_substititons_file = open ('/Users/pvmarkov/dengue/data/DENV4_complete_human_alignment_viatranslation_renamed_mapping_counts_per_branch_per_site_per_type_2.count', 'r')
except IOError:
	print ("Unknown file " + filename)
	sys.exit()
# parse line by line
for l in codon_substititons_file:
	branch_counter = 0
	substitutions_sum = 0
	if l.find ("sites") > -1: # skipping headings line
		pass
	else:
		branches_list = l.split("\t")
		#print (branches_list)
		for branch in branches_list:
			branch_counter = branch_counter + 1
			if branch_counter == 1:
				codon_number = int (branch) + 1  - 5 # the '-5' bit is an adjustment of the codon positions in the database alignment in comparison to the reference genome. '-5' is kept separate from '+1' for clarity.
			else:
				if branch == str ("nan"):
					branch = 0
				substitutions_sum = substitutions_sum + float (branch)
			#print (branch)
		#print (l [0:19])
#		print (codon_number)
#		print (substitutions_sum)
		codon_numbers_second_file.append(codon_number)
		non_synonymous_substitutions_by_codon.append(substitutions_sum) # populates the empty As list with the NUMBER of As by genome site.
		non_synonymous_substitutions_by_position_in_codon.append(substitutions_sum/2)
def ratio(x,y):
	return x/y


dn_ds = list(map(ratio, non_synonymous_substitutions_by_codon, synonymous_substitutions_by_codon))
dn_ds_per_position_in_codon = list(map(ratio, non_synonymous_substitutions_by_position_in_codon, synonymous_substitutions_by_codon))


dn_ds_by_codon_dataframe = pd.DataFrame ({'codon_numbers': codon_numbers, 'dn_ds': dn_ds, 'dn_ds_per_position_in_codon': dn_ds_per_position_in_codon, 'synonymous_substitutions':synonymous_substitutions_by_codon, 'non_synonymous_substitutions' : non_synonymous_substitutions_by_codon})
dn_ds_by_codon_dataframe.to_csv ('/Users/pvmarkov/dengue/data/dn_ds_dataframe_refined.csv')

print (len(codon_numbers))
print (len(dn_ds))
print (len(codon_numbers_second_file))
print (len(synonymous_substitutions_by_codon))
print (len(non_synonymous_substitutions_by_codon))