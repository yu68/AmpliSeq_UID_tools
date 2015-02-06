import sys, os, argparse
import itertools
from Bio import SeqIO
from time import time

print >>sys.stdout, "Checking Biopython installed or not... "
try:
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
except:
    print >>sys.stdout, "# Installing biopython... "
    os.system("easy_install --user -f http://biopython.org/DIST/ biopython")
print >>sys.stdout, "Biopython installed and loaded"

'''
dir(hsp):
['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']
'''

def ParseArg():
    p=argparse.ArgumentParser(description="Extract UID sequences and put it in fastq/a name", epilog="require: Bio")
    group=p.add_mutually_exclusive_group()
    group.add_argument("-f","--fasta",action='store_true',help='add this option for fasta input file')
    group.add_argument("-q","--fastq",action='store_true',help='add this option for fastq input file')
    p.add_argument('input',type=str,help='input fastq/fasta file with UID and identifiers')
    p.add_argument('-I','--identifier',dest='identifier',type=str,help='identifier sequence file')
    p.add_argument('-l','--uidLen',type=int,default=10,help='length of UID, will search for +-1 basepair, default=10')
    p.add_argument('-H','--head',action='store_true',help='set this if uid is at the beginning of the reads, otherwise, uid is at the end of reads')
    p.add_argument('-m','--max_mis',dest='max_mis',type=int,default=2, help="max(mismatch+indel) allowed for identifier match, otherwise move reads into 'unassigned' file. default: 2")
    p.add_argument('-o','--output',type=str, help="output fastq/fasta file name")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def fuzzy_substring(needle, haystack):
    """Calculates the fuzzy match of needle in haystack,
    using a modified version of the Levenshtein distance
    algorithm.
    The function is modified from the levenshtein function
    in the bktree module by Adam Hupp
    http://ginstrom.com/scribbles/2007/12/01/fuzzy-substring-
    matching-with-levenshtein-distance-in-python/"""
    m, n = len(needle), len(haystack)

    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m

    row1 = [0] * (n+1)
    minS=m
    for i in range(0,m):
        row2 = [i+1]
        for j in range(0,n):
            cost = ( needle[i] != haystack[j] )

            row2.append(   min(row1[j+1]+1, # deletion
                               row2[j]+1, #insertion
                               row1[j]+cost) #substitution
                          )
            if i == m-1:
                if row2[j+1] <= minS:
                    minS=row2[j+1]
                    end=j+1
        row1 = row2
    return minS, end
'''
#TEST:
print (fuzzy_substring('ACTC', 'C_ ATCG'))
print (fuzzy_substring('ACTC', 'C_ ACTGG'))
print (fuzzy_substring("ACTAAC", "ACTAACTAGCCATGCAATGGCTAG"))
'''

def Main():
    args=ParseArg()
    uidLen=args.uidLen

    if args.fastq:
       type="fastq"
    elif args.fasta:
       type="fasta"

    output=open("%s.%s"%(args.output,type),'w')
    output_unassign=open("unassign_%s.%s"%(args.output,type),'w')
    uid_info = open("uid_info.txt",'w')    
 
    #----------- read identifiers ----------
    identifiers=[]
    for i in open(args.identifier,'r'):
        i=i.strip()
        identifiers.append(i)
        identi_len=len(i)

    print >>sys.stdout, "\nStart to extract UID"
    for record in SeqIO.parse(args.input, type):
        miScore = identi_len
        record.description = record.description.split(" ")[1]
        if args.head:
            read_seq = record.seq[(uidLen-1):(uidLen+identi_len+1)]
        else:
            read_seq = record.seq[(-uidLen-identi_len-1):(-uidLen+2)]
        for i in identifiers:
            score,j=fuzzy_substring(i,read_seq)
            if score<miScore:
                barcode=i
                end=j
                miScore=score
        if miScore>args.max_mis:
            SeqIO.write(record, output_unassign, type)
        else:
            if args.head:
                uidSeq = record.seq[:(uidLen-1+end-identi_len)]
                record.id += '_%s'%(uidSeq)
                SeqIO.write(record[(uidLen-1+end-identi_len):],output,type)
            else:
                uidSeq = record.seq[(-uidLen-identi_len-1+end):]
                record.id += '_%s'%(uidSeq)
                SeqIO.write(record[:(-uidLen-identi_len-1+end)],output,type)
            print >>uid_info, "%s\t%d\t%s"%(uidSeq,miScore,record.id)
    output.close()
    output_unassign.close()
    uid_info.close()

if __name__=="__main__":
    Main()
