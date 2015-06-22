#!/usr/bin/env python
# encoding: utf-8
"""
indel_scan.py
Process Bam file to filter out indel errors

Created by Mark Evans on 2012-11-09.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys, os, getopt, subprocess, glob, random, time, re
import pysam as ps
from Bio.Seq import Seq
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter
from commands import *  # This should be replaced in Python 2.7+ with subprocess, but we need it for Bali

class glbl:
    jobs =[]
    indel_reads = {}    # {pos: [(read,indel_len,indel_type),],}
 #   multi_indel_reads   # ???


#############
# Usage     #
#############
def usage():
    code =  "\n\n#############\nProgram parses fastq files and trims reads to NS5B seq only \n\n"
    code += "Usage: $ python ns5b_fastq_trim.py [-s] fastq_filename \n\n"
    code += "-s [seq] name of fastq file to process\n"
    code += "-h help\n\n#############\n\n"
    print code
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


# Call Clustalw2
################################
def callClustal(fn,path):
    ref_aln = q_aln = ''
    subprocess.call(["clustalw2","-align","-infile="+path+fn+".seq", "-outfile="+path+fn+".aln"])
    
    # Parse the alignment file and get the aligned seq strings.  Should now be equal length
    align = AlignIO.read(path+fn+'.aln','clustal')
    for s in align:
        if   s.id == 'Ref': ref_aln = seq.seq.tostring()
        elif s.id == 'Query': q_aln = seq.seq.tostring()
    os.remove(fn+'.seq')
    os.remove(fn+'.aln')
    os.remove(fn+'.dnd')
    return (ref_aln,q_aln)


##################################################################################################
# getPolymerPositions                                                                            #
# find coordinates of all homopolymers >2 in sequence and convert coords to genome ref positions #
##################################################################################################
def getPolymerPositions(seq,start):
    # Note: positions are python string coords, so entering them directly into slice will give right result
    # ex. ([11,13],) when used thusly seq[11:13] will give 'TT' when seq is 'ATTG' and T is pos 11 and 12
    polyA = polyT = polyC = polyG = []
    
    for x in [(m.start(0), m.end(0)) for m in re.finditer(r'AA+', seq)]:
        polyA.append([x[0]+start,x[1]+start])
    
    for x in [(m.start(0), m.end(0)) for m in re.finditer(r'TT+', seq)]:
        polyT.append([x[0]+start,x[1]+start])
    
    for x in [(m.start(0), m.end(0)) for m in re.finditer(r'CC+', seq)]:
        polyC.append([x[0]+start,x[1]+start])
    
    for x in [(m.start(0), m.end(0)) for m in re.finditer(r'GG+', seq)]:
        polyG.append([x[0]+start,x[1]+start])
    
    return {'A':polyA,'T':polyT,'C':polyC,'G':polyG}


# assignWork
##########################
def assignWork(read,seq,path,poly_refs):
    # Use of cigars may obviate the need for using Clustal
    # Original cigar nomenclature looks like this, with Pos showing position on the reference:
    # POS: 5, CIGAR: 13M1D122M1I17M2I60M
    # Pysam represents cigars as tuples, so
    # cigar: [(0, 13), (2, 1), (0, 122), (1, 1), (0, 17), (1, 2), (0, 60)]
    # where the first index represents the operation (0=match,1=insertion,2=deletion)
    # and the second index represents the length of the event
    # pos info from pos, qstart, qend, qlen, rlen

    # Can get indices of all polymer locations in seq with this, repeat for each base:
    # [(m.start(0), m.end(0)) for m in re.finditer(r'TT+', d)] produces something like this: [(5, 8), (24, 27)]

    # ex qname=read_15002
    
    # create list of all reference positions that are also polymer positions (using polymer = 2+)
    temp = []
    for res in poly_refs:
        for pair in poly_refs[res]:
            for x in range(pair[0],pair[1]):
                temp.append(x)
    poly_pos = list(set(temp))
    poly_pos.sort()  # now have unique, sorted list of positions

    outFASTQfile = open(bamfile+".indel_corr.fastq",'w')

    if len(read.cigar) > 1:
        marker = 0
        idx = 0
        error_pos = []
        for cig in read.cigar:  # cigar:  [(0, 13), (2, 1), (0, 122), (1, 1), (0, 17), (1, 2), (0, 60)]
            
            # Match
            if cig[0] == 0:     
                marker = marker + cig[1]    # keep track of our overall length as we go through the cigar
            
            # Insertion
            elif cig[0] == 1:   
                idx = marker    # idx will be 1 less than the current pos (first pos of indel)
                # idx-1 because read.positions is zero-based array and cigar is 1-based count. insertions are present in read.positions
                # so testing before and first pos of indel for poly_pos membership
                if read.positions[idx-1] in poly_pos or read.positions[idx] in poly_pos:    
                    for z in range(1,cig[1]+1):     # make sure to add 'error' for each base of the indel
                        error_pos.append(idx+z)     # we are adding array index positions, not refseq positions
                marker = marker + cig[1]            # update marker to move counter to end of indel
            
            # Deletion
            elif cig[0] == 2:   
                idx = marker                        # idx will be 1 less than the current pos (first pos of indel)
                # deletions are not present in read.positions, so testing before, after and first pos of indel for poly_pos membership
                if read.positions[idx-1] in poly_pos or read.positions[idx] in poly_pos or read.positions[idx-1]+cig[1]+1 in poly_pos: 
                    for z in range(1,cig[1]+1):     # make sure to add 'error' for each base of the indel
                        error_pos.append(idx+z)     # we are adding array index positions, not refseq positions
                marker = marker + cig[1]            # update marker to move counter to end of indel
        
        # Edit read
        seq_list = list(read.query)     # convert the sequence to a zero-based array
        qual_list = list(read.qqual)    # convert the quality scores to a zero-based array, should be the same length as the sequence
        if len(error_pos != 0):
            error_pos.sort()                # make sure error indexes are sorted lo -> hi
            error_pos.reverse()             # then reverse them so now hi -> lo
            for p in error_pos:             # loop through array
                del(seq_list[p])            # remove index from seq
                del(qual_list[p])           # remove index from qual

        newseq = ''.join(seq_list)      # convert list back to a string
        newqual = ''.join(qual_list)    # convert list back to a string

        # write result
        if read.is_reverse == True:
            s = Seq(newseq)
            rc = s.reverse_complement().tostring()
            rq = reverseString(newqual)
            outFASTAfile.write("@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n")
        else:
            outFASTAfile.write("@"+read.qname+"\n"+newseq+"\n+\n"+newqual+"\n")
    else:
        outFASTAfile.write("@"+read.qname+"\n"+read.query+"\n+\n"+read.qqual+'\n')











################
# processFiles #
################
def processFiles(seqfile,threshold,width):
    # Need to keep this dictionary up-to-date with references you expect to see 
    gene_pos = {'1b_Con1_full_reference_seq':{'ns5b':{'nterm':7599,'cterm':9371,'seq':'TCGATGTCCTACACATGGACAGGCGCCCTGATCACGCCATGCGCTGCGGAGGAAACCAAGCTGCCCATCAATGCACTGAGCAACTCTTTGCTCCGTCACCACAACTTGGTCTATGCTACAACATCTCGCAGCGCAAGCCTGCGGCAGAAGAAGGTCACCTTTGACAGACTGCAGGTCCTGGACGACCACTACCGGGACGTGCTCAAGGAGATGAAGGCGAAGGCGTCCACAGTTAAGGCTAAACTTCTATCCGTGGAGGAAGCCTGTAAGCTGACGCCCCCACATTCGGCCAGATCTAAATTTGGCTATGGGGCAAAGGACGTCCGGAACCTATCCAGCAAGGCCGTTAACCACATCCGCTCCGTGTGGAAGGACTTGCTGGAAGACACTGAGACACCAATTGACACCACCATCATGGCAAAAAATGAGGTTTTCTGCGTCCAACCAGAGAAGGGGGGCCGCAAGCCAGCTCGCCTTATCGTATTCCCAGATTTGGGGGTTCGTGTGTGCGAGAAAATGGCCCTTTACGATGTGGTCTCCACCCTCCCTCAGGCCGTGATGGGCTCTTCATACGGATTCCAATACTCTCCTGGACAGCGGGTCGAGTTCCTGGTGAATGCCTGGAAAGCGAAGAAATGCCCTATGGGCTTCGCATATGACACCCGCTGTTTTGACTCAACGGTCACTGAGAATGACATCCGTGTTGAGGAGTCAATCTACCAATGTTGTGACTTGGCCCCCGAAGCCAGACAGGCCATAAGGTCGCTCACAGAGCGGCTTTACATCGGGGGCCCCCTGACTAATTCTAAAGGGCAGAACTGCGGCTATCGCCGGTGCCGCGCGAGCGGTGTACTGACGACCAGCTGCGGTAATACCCTCACATGTTACTTGAAGGCCGCTGCGGCCTGTCGAGCTGCGAAGCTCCAGGACTGCACGATGCTCGTATGCGGAGACGACCTTGTCGTTATCTGTGAAAGCGCGGGGACCCAAGAGGACGAGGCGAGCCTACGGGCCTTCACGGAGGCTATGACTAGATACTCTGCCCCCCCTGGGGACCCGCCCAAACCAGAATACGACTTGGAGTTGATAACATCATGCTCCTCCAATGTGTCAGTCGCGCACGATGCATCTGGCAAAAGGGTGTACTATCTCACCCGTGACCCCACCACCCCCCTTGCGCGGGCTGCGTGGGAGACAGCTAGACACACTCCAGTCAATTCCTGGCTAGGCAACATCATCATGTATGCGCCCACCTTGTGGGCAAGGATGATCCTGATGACTCATTTCTTCTCCATCCTTCTAGCTCAGGAACAACTTGAAAAAGCCCTAGATTGTCAGATCTACGGGGCCTGTTACTCCATTGAGCCACTTGACCTACCTCAGATCATTCAACGACTCCATGGCCTTAGCGCATTTTCACTCCATAGTTACTCTCCAGGTGAGATCAATAGGGTGGCTTCATGCCTCAGGAAACTTGGGGTACCGCCCTTGCGAGTCTGGAGACATCGGGCCAGAAGTGTCCGCGCTAGGCTACTGTCCCAGGGGGGGAGGGCTGCCACTTGTGGCAAGTACCTCTTCAACTGGGCAGTAAGGACCAAGCTCAAACTCACTCCAATCCCGGCTGCGTCCCAGTTGGATTTATCCAGCTGGTTCGTTGCTGGTTACAGCGGGGGAGACATATATCACAGCCTGTCTCGTGCCCGACCCCGCTGGTTCATGTGGTGCCTACTCCTACTTTCTGTAGGGGTAGGCATCTATCTACTCCCCAACCGA'}},
                '1a_H77_full_reference_seq':{'ns5b':{'nterm':7602,'cterm':9374, 'seq':'TCAATGTCTTATTCCTGGACAGGCGCACTCGTCACCCCGTGCGCTGCGGAAGAACAAAAACTGCCCATCAACGCACTGAGCAACTCGTTGCTACGCCATCACAATCTGGTGTATTCCACCACTTCACGCAGTGCTTGCCAAAGGCAGAAGAAAGTCACATTTGACAGACTGCAAGTTCTGGACAGCCATTACCAGGACGTGCTCAAGGAGGTCAAAGCAGCGGCGTCAAAAGTGAAGGCTAACTTGCTATCCGTAGAGGAAGCTTGCAGCCTGACGCCCCCACATTCAGCCAAATCCAAGTTTGGCTATGGGGCAAAAGACGTCCGTTGCCATGCCAGAAAGGCCGTAGCCCACATCAACTCCGTGTGGAAAGACCTTCTGGAAGACAGTGTAACACCAATAGACACTACCATCATGGCCAAGAACGAGGTTTTCTGCGTTCAGCCTGAGAAGGGGGGTCGTAAGCCAGCTCGTCTCATCGTGTTCCCCGACCTGGGCGTGCGCGTGTGCGAGAAGATGGCCCTGTACGACGTGGTTAGCAAGCTCCCCCTGGCCGTGATGGGAAGCTCCTACGGATTCCAATACTCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCTCGTATGATACCCGCTGTTTTGACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTGACCTGGACCCCCAAGCCCGCGTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGGGCCCTCTTACCAATTCAAGGGGGGAAAACTGCGGCTACCGCAGGTGCCGCGCGAGCGGCGTACTGACAACTAGCTGTGGTAACACCCTCACTTGCTACATCAAGGCCCGGGCAGCCTGTCGAGCCGCAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGACTTAGTCGTTATCTGTGAAAGTGCGGGGGTCCAGGAGGACGCGGCGAGCCTGAGAGCCTTCACGGAGGCTATGACCAGGTACTCCGCCCCCCCCGGGGACCCCCCACAACCAGAATACGACTTGGAGCTTATAACATCATGCTCCTCCAACGTGTCAGTCGCCCACGACGGCGCTGGAAAGAGGGTCTACTACCTTACCCGTGACCCTACAACCCCCCTCGCGAGAGCCGCGTGGGAGACAGCAAGACACACTCCAGTCAATTCCTGGCTAGGCAACATAATCATGTTTGCCCCCACACTGTGGGCGAGGATGATACTGATGACCCATTTCTTTAGCGTCCTCATAGCCAGGGATCAGCTTGAACAGGCTCTTAACTGTGAGATCTACGGAGCCTGCTACTCCATAGAACCACTGGATCTACCTCCAATCATTCAAAGACTCCATGGCCTCAGCGCATTTTCACTCCACAGTTACTCTCCAGGTGAAATCAATAGGGTGGCCGCATGCCTCAGAAAACTTGGGGTCCCGCCCTTGCGAGCTTGGAGACACCGGGCCCGGAGCGTCCGCGCTAGGCTTCTGTCCAGAGGAGGCAGGGCTGCCATATGTGGCAAGTACCTCTTCAACTGGGCAGTAAGAACAAAGCTCAAACTCACTCCAATAGCGGCCGCTGGCCGGCTGGACTTGTCCGGTTGGTTCACGGCTGGCTACAGCGGGGGAGACATTTATCACAGCGTGTCTCATGCCCGGCCCCGCTGGTTCTGGTTTTGCCTACTCCTGCTCGCTGCAGGGGTAGGCATCTACCTCCTCCCCAACCGA'}},
                'H77_genome':{'ns5b':{'nterm':7602,'cterm':9374,'seq':'TCAATGTCTTATTCCTGGACAGGCGCACTCGTCACCCCGTGCGCTGCGGAAGAACAAAAACTGCCCATCAACGCACTGAGCAACTCGTTGCTACGCCATCACAATCTGGTGTATTCCACCACTTCACGCAGTGCTTGCCAAAGGCAGAAGAAAGTCACATTTGACAGACTGCAAGTTCTGGACAGCCATTACCAGGACGTGCTCAAGGAGGTCAAAGCAGCGGCGTCAAAAGTGAAGGCTAACTTGCTATCCGTAGAGGAAGCTTGCAGCCTGACGCCCCCACATTCAGCCAAATCCAAGTTTGGCTATGGGGCAAAAGACGTCCGTTGCCATGCCAGAAAGGCCGTAGCCCACATCAACTCCGTGTGGAAAGACCTTCTGGAAGACAGTGTAACACCAATAGACACTACCATCATGGCCAAGAACGAGGTTTTCTGCGTTCAGCCTGAGAAGGGGGGTCGTAAGCCAGCTCGTCTCATCGTGTTCCCCGACCTGGGCGTGCGCGTGTGCGAGAAGATGGCCCTGTACGACGTGGTTAGCAAGCTCCCCCTGGCCGTGATGGGAAGCTCCTACGGATTCCAATACTCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCTCGTATGATACCCGCTGTTTTGACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTGACCTGGACCCCCAAGCCCGCGTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGGGCCCTCTTACCAATTCAAGGGGGGAAAACTGCGGCTACCGCAGGTGCCGCGCGAGCGGCGTACTGACAACTAGCTGTGGTAACACCCTCACTTGCTACATCAAGGCCCGGGCAGCCTGTCGAGCCGCAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGACTTAGTCGTTATCTGTGAAAGTGCGGGGGTCCAGGAGGACGCGGCGAGCCTGAGAGCCTTCACGGAGGCTATGACCAGGTACTCCGCCCCCCCCGGGGACCCCCCACAACCAGAATACGACTTGGAGCTTATAACATCATGCTCCTCCAACGTGTCAGTCGCCCACGACGGCGCTGGAAAGAGGGTCTACTACCTTACCCGTGACCCTACAACCCCCCTCGCGAGAGCCGCGTGGGAGACAGCAAGACACACTCCAGTCAATTCCTGGCTAGGCAACATAATCATGTTTGCCCCCACACTGTGGGCGAGGATGATACTGATGACCCATTTCTTTAGCGTCCTCATAGCCAGGGATCAGCTTGAACAGGCTCTTAACTGTGAGATCTACGGAGCCTGCTACTCCATAGAACCACTGGATCTACCTCCAATCATTCAAAGACTCCATGGCCTCAGCGCATTTTCACTCCACAGTTACTCTCCAGGTGAAATCAATAGGGTGGCCGCATGCCTCAGAAAACTTGGGGTCCCGCCCTTGCGAGCTTGGAGACACCGGGCCCGGAGCGTCCGCGCTAGGCTTCTGTCCAGAGGAGGCAGGGCTGCCATATGTGGCAAGTACCTCTTCAACTGGGCAGTAAGAACAAAGCTCAAACTCACTCCAATAGCGGCCGCTGGCCGGCTGGACTTGTCCGGTTGGTTCACGGCTGGCTACAGCGGGGGAGACATTTATCACAGCGTGTCTCATGCCCGGCCCCGCTGGTTCTGGTTTTGCCTACTCCTGCTCGCTGCAGGGGTAGGCATCTACCTCCTCCCCAACCGA'}},
                'JFH-1_genome':{'ns5b':{'nterm':7666,'cterm':9443,'seq':'CTCCATGTCATACTCCTGGACCGGGGCTCTAATAACTCCCTGTAGCCCCGAAGAGGAAAAGTTGCCAATCAACCCTTTGAGTAACTCGCTGTTGCGATACCATAACAAGGTGTACTGTACAACATCAAAGAGCGCCTCACAGAGGGCTAAAAAGGTAACTTTTGACAGGACGCAAGTGCTCGACGCCCATTATGACTCAGTCTTAAAGGACATCAAGCTAGCGGCTTCCAAGGTCAGCGCAAGGCTCCTCACCTTGGAGGAGGCGTGCCAGTTGACTCCACCCCATTCTGCAAGATCCAAGTATGGATTCGGGGCCAAGGAGGTCCGCAGCTTGTCCGGGAGGGCCGTTAACCACATCAAGTCCGTGTGGAAGGACCTCCTGGAAGACCCACAAACACCAATTCCCACAACCATCATGGCCAAAAATGAGGTGTTCTGCGTGGACCCCGCCAAGGGGGGTAAGAAACCAGCTCGCCTCATCGTTTACCCTGACCTCGGCGTCCGGGTCTGCGAGAAAATGGCCCTCTATGACATTACACAAAAGCTTCCTCAGGCGGTAATGGGAGCTTCCTATGGCTTCCAGTACTCCCCTGCCCAACGGGTGGAGTATCTCTTGAAAGCATGGGCGGAAAAGAAGGACCCCATGGGTTTTTCGTATGATACCCGATGCTTCGACTCAACCGTCACTGAGAGAGACATCAGGACCGAGGAGTCCATATACCAGGCCTGCTCCCTGCCCGAGGAGGCCCGCACTGCCATACACTCGCTGACTGAGAGACTTTACGTAGGAGGGCCCATGTTCAACAGCAAGGGTCAAACCTGCGGTTACAGACGTTGCCGCGCCAGCGGGGTGCTAACCACTAGCATGGGTAACACCATCACATGCTATGTGAAAGCCCTAGCGGCCTGCAAGGCTGCGGGGATAGTTGCGCCCACAATGCTGGTATGCGGCGATGACCTAGTAGTCATCTCAGAAAGCCAGGGGACTGAGGAGGACGAGCGGAACCTGAGAGCCTTCACGGAGGCCATGACCAGGTACTCTGCCCCTCCTGGTGATCCCCCCAGACCGGAATATGACCTGGAGCTAATAACATCCTGTTCCTCAAATGTGTCTGTGGCGTTGGGCCCGCGGGGCCGCCGCAGATACTACCTGACCAGAGACCCAACCACTCCACTCGCCCGGGCTGCCTGGGAAACAGTTAGACACTCCCCTATCAATTCATGGCTGGGAAACATCATCCAGTATGCTCCAACCATATGGGTTCGCATGGTCCTAATGACACACTTCTTCTCCATTCTCATGGTCCAAGACACCCTGGACCAGAACCTCAACTTTGAGATGTATGGATCAGTATACTCCGTGAATCCTTTGGACCTTCCAGCCATAATTGAGAGGTTACACGGGCTTGACGCCTTTTCTATGCACACATACTCTCACCACGAACTGACGCGGGTGGCTTCAGCCCTCAGAAAACTTGGGGCGCCACCCCTCAGGGTGTGGAAGAGTCGGGCTCGCGCAGTCAGGGCGTCCCTCATCTCCCGTGGAGGGAAAGCGGCCGTTTGCGGCCGATATCTCTTCAATTGGGCGGTGAAGACCAAGCTCAAACTCACTCCATTGCCGGAGGCGCGCCTACTGGACTTATCCAGTTGGTTCACCGTCGGCGCCGGCGGGGGCGACATTTTTCACAGCGTGTCGCGCGCCCGACCCCGCTCATTACTCTTCGGCCTACTCCTACTTTTCGTAGGGGTAGGCCTCTTCCTACTCCCCGCTCGGTAGA'}}}
    
    ps.index(seqfile)
    bam = ps.Samfile(seqfile,'rb')
    outFASTQfile = open(seqfile+".indel_corrected.fastq",'w')
    ref = bam.references[0]
    poly_refs = getPolymerPositions(gene_pos[ref]['ns5b']['seq'],gene_pos[ref]['ns5b']['nterm'])
    read_pool = bam.fetch(bam.references[0], gene_pos[ref]['ns5b']['nterm'],gene_pos[ref]['ns5b']['cterm'])


    for read in read_pool:
        # create list of all reference positions that are also polymer positions (using polymer = 2+)
        temp = []
        for res in poly_refs:
            for pair in poly_refs[res]:
                for x in range(pair[0],pair[1]):
                    temp.append(x)
        poly_pos = list(set(temp))
        poly_pos.sort()  # now have unique, sorted list of positions


        if len(read.cigar) > 1:
            marker = 0
            idx = 0
            error_pos = []
            for cig in read.cigar:  # cigar:  [(0, 13), (2, 1), (0, 122), (1, 1), (0, 17), (1, 2), (0, 60)]
            
                # Match
                if cig[0] == 0:     
                    marker = marker + cig[1]    # keep track of our overall length as we go through the cigar
                
                # Insertion
                elif cig[0] == 1:   
                    idx = marker    # idx will be 1 less than the current pos (first pos of indel)
                    # idx-1 because read.positions is zero-based array and cigar is 1-based count. insertions are present in read.positions
                    # so testing before and first pos of indel for poly_pos membership
                    if read.positions[idx-1] in poly_pos or read.positions[idx] in poly_pos:    
                        for z in range(1,cig[1]+1):     # make sure to add 'error' for each base of the indel
                            error_pos.append(idx+z)     # we are adding array index positions, not refseq positions
                    marker = marker + cig[1]            # update marker to move counter to end of indel
                
                # Deletion
                elif cig[0] == 2:   
                    idx = marker                        # idx will be 1 less than the current pos (first pos of indel)
                    # deletions are not present in read.positions, so testing before, after and first pos of indel for poly_pos membership
                    if read.positions[idx-1] in poly_pos or read.positions[idx] in poly_pos or read.positions[idx-1]+cig[1]+1 in poly_pos: 
                        for z in range(1,cig[1]+1):     # make sure to add 'error' for each base of the indel
                            error_pos.append(idx+z)     # we are adding array index positions, not refseq positions
                    marker = marker + cig[1]            # update marker to move counter to end of indel
                
            # Edit read
            seq_list = list(read.query)     # convert the sequence to a zero-based array
            qual_list = list(read.qqual)    # convert the quality scores to a zero-based array, should be the same length as the sequence

            if len(error_pos) != 0:
                error_pos.sort()                # make sure error indexes are sorted lo -> hi
                error_pos.reverse()             # then reverse them so now hi -> lo
                print "length seq_list = ",len(seq_list)
                print "error_pos = ",error_pos
                print read.cigar
                for p in error_pos:             # loop through array
                    print "p = ",p
                    del(seq_list[p])            # remove index from seq
                    del(qual_list[p])           # remove index from qual

            newseq = ''.join(seq_list)      # convert list back to a string
            newqual = ''.join(qual_list)    # convert list back to a string

            # write result
            if read.is_reverse == True:
                s = Seq(newseq)
                rc = s.reverse_complement().tostring()
                rq = reverseString(newqual)
                outFASTQfile.write("@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n")
            else:
                outFASTQfile.write("@"+read.qname+"\n"+newseq+"\n+\n"+newqual+"\n")
        else:
            outFASTQfile.write("@"+read.qname+"\n"+read.query+"\n+\n"+read.qqual+'\n')



    print "QC complete\n\n"





###########
# Main    #
###########
def main(argv):
    seqfile=''
    threshold = 1
    width = 0
    
    # Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:t:w:h:",["file,threshold,width,help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-f","--file"):    seqfile = arg
        elif opt in ("-t","--threshold"):   threshold = int(arg)
        elif opt in ("-w","--width"):   width = int(arg)
        
    if seqfile !='':
       processFiles(seqfile,threshold,width)
    else:
       print opts
       print sys.argv
       usage()
       sys.exit()



if __name__ == '__main__':
    main(sys.argv[1:])