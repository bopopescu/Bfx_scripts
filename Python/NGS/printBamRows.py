import pysam
import os,sys,getopt,glob
from types import *
from Bio.Seq import Seq


#############
# Usage     #
#############
def usage():
    code =  "\n\n#############\nProgram parses .bam files and trims reads to NS5B seq only, writing .fastq file\n"
    code += "Usage: $ python bamtrim2.py [-f] bam_filename \n\n"
    code += "-f [filename] name of bam file to process\n"
    code += "-a [all] must add a value after -a, will process all bam files in current directory\n"
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
def extractRecords(bamfile,records):
    pysam.index(bamfile)                # must create a .bai index for any bam file to be read or fetch won't work
    bam = pysam.Samfile(bamfile,'rb')   # and must be done before bamfile is opened
    ref = bam.references[0]             # Get name of reference reads aligned to in bam
    outfile = open(bamfile+".records",'w')

    # Need to keep this dictionary up-to-date with references you expect to see 
    gene_pos = {'1b_Con1_full_reference_seq':{'ns5b':{'nterm':7599,'cterm':9371}},
                '1a_H77_full_reference_seq':{'ns5b':{'nterm':7602,'cterm':9374}},
                'H77_genome':{'ns5b':{'nterm':7602,'cterm':9374}},
                '12_146520_1GC':{'ns5b':{'nterm':1,'cterm':1773}},
                'JFH-1_genome':{'ns5b':{'nterm':7666,'cterm':9443}}}

    # Get the reads in region of interest
    read_pool = bam.fetch(bam.references[0], gene_pos[ref]['ns5b']['nterm'],gene_pos[ref]['ns5b']['cterm'])
    
    # Process reads
    count = 0
    outfile.write("qname\taend\talen\tbin\tcigar\tflag\tis_duplicate\tis_paired\tis_proper_pair\tis_qcfail\tis_read1\tis_read2\tis_reverse\tis_secondary\tisize\tmapq\tmate_is_reverse\tmate_is_unmapped\tmpos\tmrnm\tpnext\tpos\tpositions\tqstart\tqend\tqlen\tqqual\trlen\trname\trnext\tseq\ttags\ttid\ttlen\n")

    for read in read_pool:
        if records > count:
            outfile.write(  str(read.qname)+'\t'+
                            str(read.aend)+'\t'+
                            str(read.alen)+'\t'+
                            str(read.bin)+'\t'+
                            str(read.cigar)+'\t'+
                            str(read.flag)+'\t'+
                            str(read.is_duplicate)+'\t'+
                            str(read.is_paired)+'\t'+
                            str(read.is_proper_pair)+'\t'+
                            str(read.is_qcfail)+'\t'+
                            str(read.is_read1)+'\t'+
                            str(read.is_read2)+'\t'+
                            str(read.is_reverse)+'\t'+
                            str(read.is_secondary)+'\t'+
                            str(read.isize)+'\t'+
                            str(read.mapq)+'\t'+
                            str(read.mate_is_reverse)+'\t'+
                            str(read.mate_is_unmapped)+'\t'+
                            str(read.mpos)+'\t'+
                            str(read.mrnm)+'\t'+
                            str(read.pnext)+'\t'+
                            str(read.pos)+'\t'+
                            str(read.positions)+'\t'+
                            str(read.qstart)+'\t'+
                            str(read.qend)+'\t'+
                            str(read.qlen)+'\t'+
                            str(read.qqual)+'\t'+
                            str(read.query)+'\t'+
                            str(read.rlen)+'\t'+
                            str(read.rname)+'\t'+
                            str(read.rnext)+'\t'+
                            str(read.seq)+'\t'+
                            str(read.tags)+'\t'+
                            str(read.tid)+'\t'+
                            str(read.tlen)+'\n')
            count += 1


    outfile.close()
    return


################
# processFiles #
################
def processFiles(seqfile,records):

    extractRecords(seqfile,records)

    print "Finished extracting trimmed reads\n\n"



###########
# Main    #
###########
def main(argv):
    seqfile = records = ''
    
    # Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:r:h:",["file,records,help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--file"): seqfile = arg
        elif opt in ("-r","--records"): records = int(arg)
    
    if seqfile !='' and records !='':
        processFiles(seqfile,records)
    else:
        print opts
        print sys.argv
        usage()
        sys.exit()




if __name__ == '__main__':
    main(sys.argv[1:])