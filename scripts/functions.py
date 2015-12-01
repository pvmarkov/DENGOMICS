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
	