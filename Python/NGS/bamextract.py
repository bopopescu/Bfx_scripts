import pysam
import os,sys,getopt,glob
from types import *
from Bio.Seq import Seq


#############
# Usage     #
#############
def usage():
    code =  "\n\n#############\nProgram parses a .bam file and extracts reads matching user defined start and stop positions \n"
    code += "Output format is required, can be either fasta or fastq\n"
    code += "Usage: $ python bamextract.py [-f] bam_filename [-b] begin_pos [-e] end_pos [-o] output_type\n\n"
    code += "-f [filename] name of bam file to process\n"
    code += "-b 5 prime start index in sequence to extract\n"
    code += "-e 3 prime stop index in sequence to extract\n"
    code += "-o Output format, must be either 'fasta' or 'fastq'\n"
    code += "-h help\n\n#############\n\n"
    print code
    return


##############################################################################
# showDetail                                                                 #
# A utility function that dumps all the values in a PySam AlignedRead object #
##############################################################################
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


####################################################################
# reverseString                                                    #
# Converts string to list, reverses order, converts back to string #
####################################################################
def reverseString(a):
    r = list(a)
    r.reverse()
    b = ''.join(r)
    return b


##############################################################################################################
# extractRegion                                                                                              #
# Extracts all reads from bam file within start and stop indices. Is greedy, takes those reads that          #
# touch indices. If read touches indecies, trims reads and qual scores to that index. Reverse complements    #
# reads that are reversed. Writes output to a fastq file so it can be properly imported to any other program #
##############################################################################################################
def extractRegion(bamfile,start,stop,output):
    pysam.index(bamfile)                # must create a .bai index for any bam file to be read or fetch won't work
    bam = pysam.Samfile(bamfile,'rb')   # and must be done before bamfile is opened
    ref = bam.references[0]             # Get name of reference reads aligned to in bam
    outfile = open(bamfile+".extracted."+output,'w')

    # Get the reads in region of interest
    read_pool = bam.fetch(bam.references[0], start,stop)
    
    # Process reads
    for read in read_pool:
        if read.is_reverse == True:                     # all reverse reads in a bam file have been reverse 
            seq = Seq(read.query)                       # complemented already so they need to be reverse 
            rc = seq.reverse_complement().tostring()    # complemented again, along with the quality scores
            rq = reverseString(read.qqual)              # to write correctly to the fastq
            if output == 'fastq':
                outfile.write("@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n")
            elif output == 'fasta':
                outfile.write('>'+read.qname+'\n'+rc+'\n')
        else:
            if output == 'fastq':
                outfile.write("@"+read.qname+"\n"+read.query+"\n+\n"+read.qqual+"\n")
            elif output == 'fasta':
                outfile.write('>'+read.qname+'\n'+read.query+'\n')

    outfile.close()
    return


################
# processFiles #
################
def processFiles(seqfile,start,stop,output):
#    if seqfile == 'all':
#        path =""
#        for infile in glob.glob (os.path.join(path,"*")):
#            if infile[-4:] == '.bam':
#                extractRegion(infile)
#    else:

    extractRegion(seqfile,start,stop,output)

    print "Finished extracting trimmed reads\n\n"



###########
# Main    #
###########
def main(argv):
    seqfile = start = stop = output = ''
    
    # Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:b:e:o:h:",["file,begin,end,output,help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--file"): seqfile = arg
        elif opt in ("-b","--begin"): start = int(arg)
        elif opt in ("-e","--end"): stop = int(arg)
        elif opt in ("-o","--output"): output = arg
    
    if output.lower() not in ['fasta','fastq']:
        print "Output type is not fasta or fastq\n"
        print opts
        print sys.argv
        usage()
        sys.exit()
    
    if seqfile !='' and start != '' and stop != '' and output !='':
        processFiles(seqfile,start,stop,output)
    else:
        print opts
        print sys.argv
        usage()
        sys.exit()



if __name__ == '__main__':
    main(sys.argv[1:])