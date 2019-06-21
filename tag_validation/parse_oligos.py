#!/usr/bin/env python

import difflib
import sys
import tre
import gzip

# define barcode format; build regex objects for approximate string matching
linker1 = "CCTAGTCGCGTAGAC"
l1reg = tre.compile(linker1)
linker1Length = len(linker1)

# define Fuzzyness for tre matching
fz = tre.Fuzzyness(maxins=0,maxdel=0,maxsub=1)

# pull in read for parsing
def readread(s):
    return [s.readline().rstrip('\n'), s.readline().rstrip('\n'), s.readline().rstrip('\n'), s.readline().rstrip('\n')]

def diff_letters(a, b):
    return sum ( a[i] != b[i] for i in range(len(a)) )

def parseRead(s, o):
    fastqRead=readread(s)
    while (fastqRead[0]):
        err=[]
        #find linker1 sequence
        l1 = l1reg.search(fastqRead[1],fz)
        
        #if not found, write error
        if hasattr(l1,'groups') == False:
            err.append('NoLinker1\n')
            
        #if found, ensure sufficient bases downstream, then clip
        if hasattr(l1,'groups') == True:
            linker1Pos = l1.groups()[0][0]
            
            #check if enough downstream bases
            if linker1Pos+linker1Length+14 > len(fastqRead[1]):
                err.append('WrongLinker1Pos\n')
        
        if len(err) == 0:
            o.write(fastqRead[0].split(" ",2)[0] + "\t" + fastqRead[1][linker1Pos+linker1Length:linker1Pos+linker1Length+14] + "\n")
        else:
            o.write(fastqRead[0].split(" ",2)[0] + "\t" + err[0])
        fastqRead = readread(s)

if __name__ == '__main__':
    print sys.argv
    o = gzip.open(sys.argv[2], "w")
    s = gzip.open(sys.argv[1], "r")
    parseRead(s, o)
