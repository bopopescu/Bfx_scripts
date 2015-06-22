import pysam
import os,sys,getopt
from types import *
from Bio.Seq import Seq


def showDetail(read):
	print "aend:	",read.aend
	print "alen:	",read.alen
	print "bin:		",read.bin
	print "cigar:	",read.cigar
	print "compare:	",read.compare
	print "fancy_str:",read.fancy_str
	print "flag:	",read.flag
	print "is_duplicate:",read.is_duplicate
	print "is_paired:",read.is_paired
	print "is_proper_pair:",read.is_proper_pair
	print "is_qcfail:",read.is_qcfail
	print "is_read1:",read.is_read1
	print "is_read2:",read.is_read2
	print "is_reverse:",read.is_reverse
	print "is_secondary:",read.is_secondary
	print "is_unmapped:",read.is_unmapped
	print "isize:	",read.isize
	print "mapq:	",read.mapq
	print "mate_is_reverse:",read.mate_is_reverse
	print "mate_is_unmapped:",read.mate_is_unmapped
	print "mpos:	",read.mpos
	print "mrnm:	",read.mrnm
	print "opt:		",read.opt
	print "overlap:	",read.overlap
	print "pnext:	",read.pnext
	print "pos:		",read.pos
	print "positions:",read.positions
	print "qstart:	",read.qstart
	print "qend:	",read.qend
	print "qlen:	",read.qlen
	print "qname:	",read.qname
	print "qqual:	",read.qqual
	print "query:	",read.query
	print "rlen:	",read.rlen
	print "rname:	",read.rname
	print "rnext:	",read.rnext
	print "seq:		",read.seq
	print "tags:	",read.tags
	print "tid:		",read.tid
	print "tlen:	",read.tlen
	return

def reverseString(a):
    r = list(a)
    r.reverse()
    b = ''.join(r)
    return b

################
# processFiles #
################
def processFiles(seqfile):
    bamfile = pysam.Samfile(seqfile,'rb')
    ref = bamfile.references[0]
    refName = [ref]
    refLength = [int(bamfile.lengths[0])]
    gene_pos = {'1b_Con1_full_reference_seq':{'ns5b':{'nterm':7599,'cterm':9371}},
                '1a_H77_full_reference_seq':{'ns5b':{'nterm':7602,'cterm':9374}}}

    read_pool = bamfile.fetch(bamfile.references[0], gene_pos[ref]['ns5b']['nterm'],gene_pos[ref]['ns5b']['cterm'])
 #   outBAMfile = pysam.Samfile(seqfile+".extracted.bam","wb",referencenames=refName,referencelengths=refLength)
    outFASTQfile = open(seqfile+".extracted.fastq",'w')

    for read in read_pool:
        seqlen = len(read.seq)
        if read.pos >= gene_pos[ref]['ns5b']['nterm'] and read.pos+seqlen <= gene_pos[ref]['ns5b']['cterm']:
            #outBAMfile.write(read)
            if read.is_reverse == True:
                seq = Seq(read.query)
                rc = seq.reverse_complement().tostring()
                rq = reverseString(read.qqual)
                outFASTQfile.write("@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n")
            else:
                outFASTQfile.write("@"+read.qname+"\n"+read.query+"\n+\n"+read.qqual+"\n")
#            outFASTQfile.write("@"+read.qname+"\n"+read.query+"\n+\n"+read.qqual+"\n")

        # Longer than gene N-term
        elif read.pos < gene_pos[ref]['ns5b']['nterm']:
          #  c = gene_pos[ref]['ns5b']['nterm'] - (read.pos - read.qstart + 1) - 1
          #  read.seq = read.seq[c:]
            q = gene_pos[ref]['ns5b']['nterm'] - read.pos - 1
            #if type(read.qqual) != NoneType: read.qqual = read.qqual[q:]
         #   read.query = read.query[q:]
         #   outBAMfile.write(read)
            if read.is_reverse == True:
                seq = Seq(read.query[q:])
                rc = seq.reverse_complement().tostring()
                rq = reverseString(read.qqual[q:])
                outFASTQfile.write("@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n")
            else:
                outFASTQfile.write("@"+read.qname+"\n"+read.query[q:]+"\n+\n"+read.qqual[q:]+"\n")
#            outFASTQfile.write("@"+read.qname+"\n"+read.query[q:]+"\n+\n"+read.qqual[q:]+"\n")

        # Longer than C-term
        elif ((read.pos-read.qstart) + len(read.seq)) > gene_pos[ref]['ns5b']['cterm']:
            s = gene_pos[ref]['ns5b']['cterm']
            if read.pos <= s:
              #  c = s - (read.pos + read.qstart)
              #  read.seq = read.seq[:c]
                q = s - read.pos
                if read.is_reverse == True:
                    seq = Seq(read.query[:q])
                    rc = seq.reverse_complement().tostring()
                    rq = reverseString(read.qqual[:q])
                    outFASTQfile.write("@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n")
                else:

                #if type(read.qqual) != NoneType: read.qqual = read.qqual[:q]
              #  read.query = read.query[:q]
              #  outBAMfile.write(read)
                    outFASTQfile.write("@"+read.qname+"\n"+read.query[:q]+"\n+\n"+read.qqual[:q]+"\n")
#            outFASTQfile.write("@"+read.qname+"\n"+read.query[:q]+"\n+\n"+read.qqual[:q]+"\n")

    #outBAMfile.close()
    outFASTQfile.close()
    print "Finished extracting trimmed reads\n\n"


# iter = samfile.fetch('1b_Con1_full_reference_seq',7599,7600)
# iter = samfile.fetch('1b_Con1_full_reference_seq',9371,9372)








###########
# Main    #
###########
def main(argv):
    seqfile=''
    
    # Read command line arguments
    try:
       opts, args = getopt.getopt(argv,"f:h:",["file,help"])
    except getopt.GetoptError:
       usage()
       sys.exit(2)
    for opt, arg in opts:
       if opt in ("-h","--help"):
          usage()
          sys.exit()
       elif opt in ("-f","--file"): seqfile = arg
   
    if seqfile !='':
       processFiles(seqfile)
    else:
       print opts
       print sys.argv
       usage()
       sys.exit()



if __name__ == '__main__':
    main(sys.argv[1:])