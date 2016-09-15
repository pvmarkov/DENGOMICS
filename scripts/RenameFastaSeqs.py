#! /usr/bin/env python
#renames fasta seqs

import math
import sys
#print sys.argv

argv=sys.argv[1:]
file=argv[0]
out=argv[1]
ref=argv[2]

try:
    f=open(file, 'r')
except IOError:
    print ("Unknown file: ",file)
    sys.exit()

try:
    fout=open(out, 'w')
except IOError:
    print ("Unknown file: ",out)
    sys.exit()

try:
    fref=open(ref, 'w')
except IOError:
    print ("Unknown file: ",ref)
    sys.exit()

for l in f:
    if '>' in l:
        liste=l.split()
        num = liste[0].split("|")[1]
        fout.write(">"+ liste[5].replace('/', '_').replace('-', '_').replace(',', '') + "_" + num + '\n')
        fref.write(l.strip()+'\t'+ liste[5].replace('/', '_').replace('-', '_').replace(',', '') + "_" + num + '\n')
    else:
        fout.write(l)
f.close()
fout.close()
fref.close()
