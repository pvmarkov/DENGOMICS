#! /usr/bin/env python

import sys

argv=sys.argv[1:]
file=argv[0]
out=argv[0].split('.')[0]+'.phy'

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

hash=dict()
l_old=""
seq=""
for l in f:
     if '>' in l:
         hash[l_old]=seq
         seq=""
         l_old=l
     else:
         l=l.strip().replace(' ','')
         seq=seq+l
hash[l_old]=seq

print (hash.__len__())
print (hash.keys())
f.close()


fout.write(str(hash.__len__()-1)+" "+str(len(hash[l_old]))+"\n")
for eachseq in hash:
    if (eachseq.strip()=="") :
    	pass
    else :
    	fout.write (eachseq.replace('>','').replace('\n','  '))
    	fout.write (hash[eachseq]+'\n')
fout.close()
