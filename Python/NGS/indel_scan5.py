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
from operator import itemgetter
from commands import *  # This should be replaced in Python 2.7+ with subprocess, but we need it for Bali

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


##############################################################################
# showMiniDetail                                                             #
# A utility function that dumps all the values in a PySam AlignedRead object #
##############################################################################
def showMiniDetail(read):
    print "**************************\n"
    print "qname:   ",read.qname
    print "aend:    ",read.aend
    print "alen:    ",read.alen
    print "cigar:   ",read.cigar
    print "is_unmapped:",read.is_unmapped
    print "is_reverse:",read.is_reverse
    print "pos:     ",read.pos
    print "positions:",read.positions
    print "qstart:  ",read.qstart
    print "qend:    ",read.qend
    print "qlen:    ",read.qlen
    print "qqual:   ",read.qqual
    print "query:   ",read.query
    print "rlen:    ",read.rlen
    print "rname:   ",read.rname
    print "rnext:   ",read.rnext
    print "seq:     ",read.seq
    print "tlen:    ",read.tlen
    print "**************************\n"
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


#####################################################################################
# forceFastQ                                                                        #
# Formats read into fastq format string, performing reverse complement if necessary #
#####################################################################################
def forceFastQ(qname,seq,qual,rev):
    fastq = ''
    if rev == True:                     # all reverse reads in a bam file have been reverse 
        newseq = Seq(seq)                       # complemented already so they need to be reverse 
        rc = newseq.reverse_complement().tostring()    # complemented again, along with the quality scores
        rq = reverseString(qual)              # to write correctly to the fastq
        fastq = "@"+qname+"\n"+rc+"\n+\n"+rq+"\n"
    else:
        fastq = "@"+qname+"\n"+seq+"\n+\n"+qual+"\n"
    return fastq


def writeFastQ(read):
    fastq = ''
    if read.is_reverse == True:                     # all reverse reads in a bam file have been reverse 
        seq = Seq(read.query)                       # complemented already so they need to be reverse 
        rc = seq.reverse_complement().tostring()    # complemented again, along with the quality scores
        rq = reverseString(read.qqual)              # to write correctly to the fastq
        fastq = "@"+read.qname+"\n"+rc+"\n+\n"+rq+"\n"
    else:
        fastq = "@"+read.qname+"\n"+read.query+"\n+\n"+read.qqual+"\n"
    return fastq


################
# processFiles #
################
def processFiles(seqfile):
    # Need to keep this dictionary up-to-date with references you expect to see 
    gene_pos = {'1b_Con1_full_reference_seq':{'ns5b':{'nterm':7599,'cterm':9371,'seq':'TCGATGTCCTACACATGGACAGGCGCCCTGATCACGCCATGCGCTGCGGAGGAAACCAAGCTGCCCATCAATGCACTGAGCAACTCTTTGCTCCGTCACCACAACTTGGTCTATGCTACAACATCTCGCAGCGCAAGCCTGCGGCAGAAGAAGGTCACCTTTGACAGACTGCAGGTCCTGGACGACCACTACCGGGACGTGCTCAAGGAGATGAAGGCGAAGGCGTCCACAGTTAAGGCTAAACTTCTATCCGTGGAGGAAGCCTGTAAGCTGACGCCCCCACATTCGGCCAGATCTAAATTTGGCTATGGGGCAAAGGACGTCCGGAACCTATCCAGCAAGGCCGTTAACCACATCCGCTCCGTGTGGAAGGACTTGCTGGAAGACACTGAGACACCAATTGACACCACCATCATGGCAAAAAATGAGGTTTTCTGCGTCCAACCAGAGAAGGGGGGCCGCAAGCCAGCTCGCCTTATCGTATTCCCAGATTTGGGGGTTCGTGTGTGCGAGAAAATGGCCCTTTACGATGTGGTCTCCACCCTCCCTCAGGCCGTGATGGGCTCTTCATACGGATTCCAATACTCTCCTGGACAGCGGGTCGAGTTCCTGGTGAATGCCTGGAAAGCGAAGAAATGCCCTATGGGCTTCGCATATGACACCCGCTGTTTTGACTCAACGGTCACTGAGAATGACATCCGTGTTGAGGAGTCAATCTACCAATGTTGTGACTTGGCCCCCGAAGCCAGACAGGCCATAAGGTCGCTCACAGAGCGGCTTTACATCGGGGGCCCCCTGACTAATTCTAAAGGGCAGAACTGCGGCTATCGCCGGTGCCGCGCGAGCGGTGTACTGACGACCAGCTGCGGTAATACCCTCACATGTTACTTGAAGGCCGCTGCGGCCTGTCGAGCTGCGAAGCTCCAGGACTGCACGATGCTCGTATGCGGAGACGACCTTGTCGTTATCTGTGAAAGCGCGGGGACCCAAGAGGACGAGGCGAGCCTACGGGCCTTCACGGAGGCTATGACTAGATACTCTGCCCCCCCTGGGGACCCGCCCAAACCAGAATACGACTTGGAGTTGATAACATCATGCTCCTCCAATGTGTCAGTCGCGCACGATGCATCTGGCAAAAGGGTGTACTATCTCACCCGTGACCCCACCACCCCCCTTGCGCGGGCTGCGTGGGAGACAGCTAGACACACTCCAGTCAATTCCTGGCTAGGCAACATCATCATGTATGCGCCCACCTTGTGGGCAAGGATGATCCTGATGACTCATTTCTTCTCCATCCTTCTAGCTCAGGAACAACTTGAAAAAGCCCTAGATTGTCAGATCTACGGGGCCTGTTACTCCATTGAGCCACTTGACCTACCTCAGATCATTCAACGACTCCATGGCCTTAGCGCATTTTCACTCCATAGTTACTCTCCAGGTGAGATCAATAGGGTGGCTTCATGCCTCAGGAAACTTGGGGTACCGCCCTTGCGAGTCTGGAGACATCGGGCCAGAAGTGTCCGCGCTAGGCTACTGTCCCAGGGGGGGAGGGCTGCCACTTGTGGCAAGTACCTCTTCAACTGGGCAGTAAGGACCAAGCTCAAACTCACTCCAATCCCGGCTGCGTCCCAGTTGGATTTATCCAGCTGGTTCGTTGCTGGTTACAGCGGGGGAGACATATATCACAGCCTGTCTCGTGCCCGACCCCGCTGGTTCATGTGGTGCCTACTCCTACTTTCTGTAGGGGTAGGCATCTATCTACTCCCCAACCGA'}},
                '1a_H77_full_reference_seq':{'ns5b':{'nterm':7602,'cterm':9374, 'seq':'TCAATGTCTTATTCCTGGACAGGCGCACTCGTCACCCCGTGCGCTGCGGAAGAACAAAAACTGCCCATCAACGCACTGAGCAACTCGTTGCTACGCCATCACAATCTGGTGTATTCCACCACTTCACGCAGTGCTTGCCAAAGGCAGAAGAAAGTCACATTTGACAGACTGCAAGTTCTGGACAGCCATTACCAGGACGTGCTCAAGGAGGTCAAAGCAGCGGCGTCAAAAGTGAAGGCTAACTTGCTATCCGTAGAGGAAGCTTGCAGCCTGACGCCCCCACATTCAGCCAAATCCAAGTTTGGCTATGGGGCAAAAGACGTCCGTTGCCATGCCAGAAAGGCCGTAGCCCACATCAACTCCGTGTGGAAAGACCTTCTGGAAGACAGTGTAACACCAATAGACACTACCATCATGGCCAAGAACGAGGTTTTCTGCGTTCAGCCTGAGAAGGGGGGTCGTAAGCCAGCTCGTCTCATCGTGTTCCCCGACCTGGGCGTGCGCGTGTGCGAGAAGATGGCCCTGTACGACGTGGTTAGCAAGCTCCCCCTGGCCGTGATGGGAAGCTCCTACGGATTCCAATACTCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCTCGTATGATACCCGCTGTTTTGACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTGACCTGGACCCCCAAGCCCGCGTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGGGCCCTCTTACCAATTCAAGGGGGGAAAACTGCGGCTACCGCAGGTGCCGCGCGAGCGGCGTACTGACAACTAGCTGTGGTAACACCCTCACTTGCTACATCAAGGCCCGGGCAGCCTGTCGAGCCGCAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGACTTAGTCGTTATCTGTGAAAGTGCGGGGGTCCAGGAGGACGCGGCGAGCCTGAGAGCCTTCACGGAGGCTATGACCAGGTACTCCGCCCCCCCCGGGGACCCCCCACAACCAGAATACGACTTGGAGCTTATAACATCATGCTCCTCCAACGTGTCAGTCGCCCACGACGGCGCTGGAAAGAGGGTCTACTACCTTACCCGTGACCCTACAACCCCCCTCGCGAGAGCCGCGTGGGAGACAGCAAGACACACTCCAGTCAATTCCTGGCTAGGCAACATAATCATGTTTGCCCCCACACTGTGGGCGAGGATGATACTGATGACCCATTTCTTTAGCGTCCTCATAGCCAGGGATCAGCTTGAACAGGCTCTTAACTGTGAGATCTACGGAGCCTGCTACTCCATAGAACCACTGGATCTACCTCCAATCATTCAAAGACTCCATGGCCTCAGCGCATTTTCACTCCACAGTTACTCTCCAGGTGAAATCAATAGGGTGGCCGCATGCCTCAGAAAACTTGGGGTCCCGCCCTTGCGAGCTTGGAGACACCGGGCCCGGAGCGTCCGCGCTAGGCTTCTGTCCAGAGGAGGCAGGGCTGCCATATGTGGCAAGTACCTCTTCAACTGGGCAGTAAGAACAAAGCTCAAACTCACTCCAATAGCGGCCGCTGGCCGGCTGGACTTGTCCGGTTGGTTCACGGCTGGCTACAGCGGGGGAGACATTTATCACAGCGTGTCTCATGCCCGGCCCCGCTGGTTCTGGTTTTGCCTACTCCTGCTCGCTGCAGGGGTAGGCATCTACCTCCTCCCCAACCGA'}},
                'H77_genome':{'ns5b':{'nterm':7602,'cterm':9374,'seq':'TCAATGTCTTATTCCTGGACAGGCGCACTCGTCACCCCGTGCGCTGCGGAAGAACAAAAACTGCCCATCAACGCACTGAGCAACTCGTTGCTACGCCATCACAATCTGGTGTATTCCACCACTTCACGCAGTGCTTGCCAAAGGCAGAAGAAAGTCACATTTGACAGACTGCAAGTTCTGGACAGCCATTACCAGGACGTGCTCAAGGAGGTCAAAGCAGCGGCGTCAAAAGTGAAGGCTAACTTGCTATCCGTAGAGGAAGCTTGCAGCCTGACGCCCCCACATTCAGCCAAATCCAAGTTTGGCTATGGGGCAAAAGACGTCCGTTGCCATGCCAGAAAGGCCGTAGCCCACATCAACTCCGTGTGGAAAGACCTTCTGGAAGACAGTGTAACACCAATAGACACTACCATCATGGCCAAGAACGAGGTTTTCTGCGTTCAGCCTGAGAAGGGGGGTCGTAAGCCAGCTCGTCTCATCGTGTTCCCCGACCTGGGCGTGCGCGTGTGCGAGAAGATGGCCCTGTACGACGTGGTTAGCAAGCTCCCCCTGGCCGTGATGGGAAGCTCCTACGGATTCCAATACTCACCAGGACAGCGGGTTGAATTCCTCGTGCAAGCGTGGAAGTCCAAGAAGACCCCGATGGGGTTCTCGTATGATACCCGCTGTTTTGACTCCACAGTCACTGAGAGCGACATCCGTACGGAGGAGGCAATTTACCAATGTTGTGACCTGGACCCCCAAGCCCGCGTGGCCATCAAGTCCCTCACTGAGAGGCTTTATGTTGGGGGCCCTCTTACCAATTCAAGGGGGGAAAACTGCGGCTACCGCAGGTGCCGCGCGAGCGGCGTACTGACAACTAGCTGTGGTAACACCCTCACTTGCTACATCAAGGCCCGGGCAGCCTGTCGAGCCGCAGGGCTCCAGGACTGCACCATGCTCGTGTGTGGCGACGACTTAGTCGTTATCTGTGAAAGTGCGGGGGTCCAGGAGGACGCGGCGAGCCTGAGAGCCTTCACGGAGGCTATGACCAGGTACTCCGCCCCCCCCGGGGACCCCCCACAACCAGAATACGACTTGGAGCTTATAACATCATGCTCCTCCAACGTGTCAGTCGCCCACGACGGCGCTGGAAAGAGGGTCTACTACCTTACCCGTGACCCTACAACCCCCCTCGCGAGAGCCGCGTGGGAGACAGCAAGACACACTCCAGTCAATTCCTGGCTAGGCAACATAATCATGTTTGCCCCCACACTGTGGGCGAGGATGATACTGATGACCCATTTCTTTAGCGTCCTCATAGCCAGGGATCAGCTTGAACAGGCTCTTAACTGTGAGATCTACGGAGCCTGCTACTCCATAGAACCACTGGATCTACCTCCAATCATTCAAAGACTCCATGGCCTCAGCGCATTTTCACTCCACAGTTACTCTCCAGGTGAAATCAATAGGGTGGCCGCATGCCTCAGAAAACTTGGGGTCCCGCCCTTGCGAGCTTGGAGACACCGGGCCCGGAGCGTCCGCGCTAGGCTTCTGTCCAGAGGAGGCAGGGCTGCCATATGTGGCAAGTACCTCTTCAACTGGGCAGTAAGAACAAAGCTCAAACTCACTCCAATAGCGGCCGCTGGCCGGCTGGACTTGTCCGGTTGGTTCACGGCTGGCTACAGCGGGGGAGACATTTATCACAGCGTGTCTCATGCCCGGCCCCGCTGGTTCTGGTTTTGCCTACTCCTGCTCGCTGCAGGGGTAGGCATCTACCTCCTCCCCAACCGA'}},
                'JRCSF':{'ns5b':{'nterm':0,'cterm':2547,'seq':'ATGAGAGTGAAGGGGATCAGGAAGAATTATCAGCACTTGTGGAAAGGGGGCATCTTGCTCCTTGGGACATTAATGATCTGTAGTGCTGTAGAAAAGTTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAAGAAACAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGAAAATGTAACAGAAGATTTTAACATGTGGAAAAATAACATGGTAGAACAGATGCAGGAGGATGTAATCAATTTATGGGATCAAAGCTTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTACTTTAAATTGCAAAGATGTGAATGCTACTAATACCACTAGTAGTAGTGAGGGAATGATGGAGAGAGGAGAAATAAAAAACTGCTCTTTCAATATCACCAAAAGCATAAGAGATAAGGTGCAGAAAGAATATGCTCTTTTTTATAAACTGGATGTAGTACCAATAGATAATAAGAATAATACCAAATATAGGTTAATAAGTTGTAACACCTCAGTCATTACACAAGCCTGTCCAAAGGTATCCTTTGAACCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAGTGTAATAATAAGACATTCAATGGAAAAGGACAATGTAAAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGCTAAATGGCAGTCTAGCAGAAGAAAAGGTTGTAATTAGATCTGACAATTTTACGGACAATGCTAAAACCATAATAGTACAGCTGAATGAATCTGTAAAAATTAATTGTACAAGGCCCAGCAACAATACAAGAAAAAGTATACATATAGGACCAGGGAGAGCATTTTATACAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGTAACATTAGTAGAGCACAATGGAATAACACTTTAAAACAGATAGTTGAAAAATTAAGAGAACAATTTAATAATAAAACAATAGTCTTTACTCACTCCTCAGGAGGGGATCCAGAAATTGTAATGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGAATGATACTGAAAAGTCAAGTGGCACTGAAGGAAATGACACCATCATACTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGGAAGTGGGAAAAGCAATGTATGCTCCTCCCATTAAAGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGAGATGGTGGTAAAAATGAGAGTGAGATCGAGATCTTCAGACCTGGAGGAGGAGACATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAAAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGTCAATGACACTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAACAGCAAAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATATGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAGTCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCATGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATACTAGTTGGAGTAATAAATCTCTGGATAGTATTTGGAATAACATGACCTGGATGGAGTGGGAAAAAGAAATTGAGAATTACACAAACACAATATACACCCTAATTGAAGAATCGCAGATCCAACAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTGACATAACAAAATGGCTGTGGTATATAAAAATATTCATAATGATAGTAGGAGGCTTGATAGGTTTAAGAATAGTTTTTTCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATACTCACCCTTATCGTTTCAGACCCTCCTCCCAGCAACGAGGGGACCCGACAGGCCCGAAGGAATCGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCGGACAATTAGTGAACGGATTCTTAGCACTTATCTGGGTCGACCTGCGGAGCCTGTTCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGACTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAATCCTGAAATACTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTTAATGCCACAGCTATAGCAGTAGCTGAGGGGACAGATAGGATTATAGAAGTAGTACAAAGAGTTTATAGGGCTATTCTCCACATACCTACAAGAATAAGACAGGGCTTGGAAAGGGCTTTGCTATAA'}},
                'JFH-1_genome':{'ns5b':{'nterm':7666,'cterm':9443,'seq':'CTCCATGTCATACTCCTGGACCGGGGCTCTAATAACTCCCTGTAGCCCCGAAGAGGAAAAGTTGCCAATCAACCCTTTGAGTAACTCGCTGTTGCGATACCATAACAAGGTGTACTGTACAACATCAAAGAGCGCCTCACAGAGGGCTAAAAAGGTAACTTTTGACAGGACGCAAGTGCTCGACGCCCATTATGACTCAGTCTTAAAGGACATCAAGCTAGCGGCTTCCAAGGTCAGCGCAAGGCTCCTCACCTTGGAGGAGGCGTGCCAGTTGACTCCACCCCATTCTGCAAGATCCAAGTATGGATTCGGGGCCAAGGAGGTCCGCAGCTTGTCCGGGAGGGCCGTTAACCACATCAAGTCCGTGTGGAAGGACCTCCTGGAAGACCCACAAACACCAATTCCCACAACCATCATGGCCAAAAATGAGGTGTTCTGCGTGGACCCCGCCAAGGGGGGTAAGAAACCAGCTCGCCTCATCGTTTACCCTGACCTCGGCGTCCGGGTCTGCGAGAAAATGGCCCTCTATGACATTACACAAAAGCTTCCTCAGGCGGTAATGGGAGCTTCCTATGGCTTCCAGTACTCCCCTGCCCAACGGGTGGAGTATCTCTTGAAAGCATGGGCGGAAAAGAAGGACCCCATGGGTTTTTCGTATGATACCCGATGCTTCGACTCAACCGTCACTGAGAGAGACATCAGGACCGAGGAGTCCATATACCAGGCCTGCTCCCTGCCCGAGGAGGCCCGCACTGCCATACACTCGCTGACTGAGAGACTTTACGTAGGAGGGCCCATGTTCAACAGCAAGGGTCAAACCTGCGGTTACAGACGTTGCCGCGCCAGCGGGGTGCTAACCACTAGCATGGGTAACACCATCACATGCTATGTGAAAGCCCTAGCGGCCTGCAAGGCTGCGGGGATAGTTGCGCCCACAATGCTGGTATGCGGCGATGACCTAGTAGTCATCTCAGAAAGCCAGGGGACTGAGGAGGACGAGCGGAACCTGAGAGCCTTCACGGAGGCCATGACCAGGTACTCTGCCCCTCCTGGTGATCCCCCCAGACCGGAATATGACCTGGAGCTAATAACATCCTGTTCCTCAAATGTGTCTGTGGCGTTGGGCCCGCGGGGCCGCCGCAGATACTACCTGACCAGAGACCCAACCACTCCACTCGCCCGGGCTGCCTGGGAAACAGTTAGACACTCCCCTATCAATTCATGGCTGGGAAACATCATCCAGTATGCTCCAACCATATGGGTTCGCATGGTCCTAATGACACACTTCTTCTCCATTCTCATGGTCCAAGACACCCTGGACCAGAACCTCAACTTTGAGATGTATGGATCAGTATACTCCGTGAATCCTTTGGACCTTCCAGCCATAATTGAGAGGTTACACGGGCTTGACGCCTTTTCTATGCACACATACTCTCACCACGAACTGACGCGGGTGGCTTCAGCCCTCAGAAAACTTGGGGCGCCACCCCTCAGGGTGTGGAAGAGTCGGGCTCGCGCAGTCAGGGCGTCCCTCATCTCCCGTGGAGGGAAAGCGGCCGTTTGCGGCCGATATCTCTTCAATTGGGCGGTGAAGACCAAGCTCAAACTCACTCCATTGCCGGAGGCGCGCCTACTGGACTTATCCAGTTGGTTCACCGTCGGCGCCGGCGGGGGCGACATTTTTCACAGCGTGTCGCGCGCCCGACCCCGCTCATTACTCTTCGGCCTACTCCTACTTTTCGTAGGGGTAGGCCTCTTCCTACTCCCCGCTCGGTAGA'}}}
    
    ps.index(seqfile)
    bam = ps.Samfile(seqfile,'rb')
    outFASTQfile = open(seqfile+".indel_corrected.fastq",'w')
    ERRfile = open(seqfile+".multi_indel_err.fastq",'w')
    changeFile = open(seqfile+".changed",'w')
    ref = bam.references[0]
    poly_refs = getPolymerPositions(gene_pos[ref]['ns5b']['seq'],gene_pos[ref]['ns5b']['nterm'])
    read_pool = bam.fetch(bam.references[0], gene_pos[ref]['ns5b']['nterm'],gene_pos[ref]['ns5b']['cterm'])
    refseq = gene_pos[ref]['ns5b']['seq']
    offset = gene_pos[ref]['ns5b']['nterm']

    # create list of all reference positions that are also polymer positions (using polymer = 2+ aa)
    temp = []
    for res in poly_refs:
        for pair in poly_refs[res]:
            for x in range(pair[0],pair[1]):
                temp.append(x)
    poly_pos = list(set(temp))
    poly_pos.sort()  # now have unique, sorted list of positions
    

    for read in read_pool:
        # cigar:  [(0, 13), (2, 1), (0, 122), (1, 1), (0, 17), (1, 2), (0, 60)]
        if read.is_unmapped == False:
            if len(read.cigar) > 1: # len 1 should equal a complete match, no indels
                marker = 0
                idx = 0
                ins_error_pos = del_error_pos = []
                del_detail = {}
                ops = []    # operations types in the read
                # collect all operations for this cigar
                for op,oplen in read.cigar:  
                    ops.append(op)

                # if ops include 1 or 2, then an indel occurred
                if 1 in ops or 2 in ops:    # Has indel
                    # Is there one indel or more than one
                    if ops.count(1) + ops.count(2) == 1:
                        # Loop through cigar again to find the indel
                        for op,oplen in read.cigar:

                            # Insertion that is not in frame
                            if op == 1 and oplen % 3 != 0:
                    #            print "procesing insertion......"
                                idx = marker
                                # Assume insertion does not occur BEFORE homopolymer, only after or within
                                try:
                                    if read.positions[idx] in poly_pos or read.positions[idx+1] in poly_pos:
                                        for z in range(1,oplen+1):         # make sure to add 'error' for each base of the indel
                                            ins_error_pos.append(idx+z)     # we are adding array index positions, not refseq positions
                                        seq_list = list(read.query)     # convert the sequence to a zero-based array
                                        qual_list = list(read.qqual)    # convert the quality scores to a zero-based array, should be the same length as the sequence
                                    
                                        if len(ins_error_pos) != 0:
                                            ins_error_pos.sort()                # make sure error indexes are sorted lo -> hi
                                            ins_error_pos.reverse()             # then reverse them so now hi -> lo
                                      #  showMiniDetail(read)
                                            if read.qname == 'read_154':
                                                showMiniDetail(read)
                                                print "length seq_list = ",len(seq_list)
                                                print "ins_error_pos = ",ins_error_pos
                                            for p in ins_error_pos:             # loop through array
                                                del(seq_list[p-1])            # remove index from seq
                                                del(qual_list[p-1])           # remove index from qual
                                        newseq = ''.join(seq_list)      # convert list back to a string
                                        newqual = ''.join(qual_list)    # convert list back to a string
                                        changeFile.write(read.qname+'\t'+str(read.cigar)+'\t'+str(read.query)+'\t'+str(newseq)+'\n')
                                        outFASTQfile.write(forceFastQ(read.qname,newseq,newqual,read.is_reverse))
                                
                                # Indel is frameshift, but not in homopolymer so probably real
                                    else:
                                        outFASTQfile.write(writeFastQ(read))
                                except:
                                    print "########## Error ###############"
                                    print 'Error:',sys.exc_info()[0]
                                    print 'read_name: ',read.qname
                                    showMiniDetail(read)
                                    print 'idx = ',idx
                                    print 'read.positions = ',read.positions
                                    print 'poly_pos = ',poly_pos
                                    print "################################"
                                    sys.exit()
                            # Deletion that is not in frame
                            elif op == 2 and oplen % 3 != 0:
                                idx = marker                        # idx will be 1 less than the current pos (first pos of indel)
                                # deletions are not present in read.positions, so testing before, after and first pos of indel for poly_pos membership
                                if read.positions[idx] in poly_pos or read.positions[idx]+oplen+1 in poly_pos or read.positions[idx+1] in poly_pos: 
                                    del_error_pos.append(idx)     # only capture index of first pos for this array 
                                    if oplen == 1:
                                        del_detail[idx] = (oplen,read.positions[idx-1]+1)    # del len 1, get corresponding refpos
                                    else:
                                        del_detail[idx] =(oplen,read.positions[idx-1]+1,read.positions[idx-1]+1+oplen)  # del len > 1, get refpos start and stop vals
                                    seq_list = list(read.query)     # convert the sequence to a zero-based array
                                    qual_list = list(read.qqual)    # convert the quality scores to a zero-based array, should be the same length as the sequence
                                    if len(del_error_pos) != 0:
                                        if read.qname=='read_106':
                                            showMiniDetail(read)
                                            print "length seq_list = ",len(seq_list)
                                            print "del_detail = ",del_detail
                                        for p in del_error_pos:             # loop through array
                                            if del_detail.has_key(p):
                                                if del_detail[p][0] == 1:
                                                  #      print "seq_list.insert(",p+1,",refseq[",del_detail[p][1]-offset,"]"
                                                    seq_list.insert(p,refseq[del_detail[p][1]-offset])
                                                    qual_list.insert(p,qual_list[p])
                                                else:
                                                    seq_list.insert(p,refseq[del_detail[p][1]-offset:del_detail[p][2]-offset])
                                                    qual_list.insert(p,del_detail[p][0]*qual_list[p])
                                
                                    newseq = ''.join(seq_list)      # convert list back to a string
                                    newqual = ''.join(qual_list)    # convert list back to a string
                                    changeFile.write(read.qname+'\t'+str(read.cigar)+'\t'+str(read.query)+'\t'+str(newseq)+'\n')
                                    outFASTQfile.write(forceFastQ(read.qname,newseq,newqual,read.is_reverse))
                                
                                # Indel is frameshift, but not in homopolymer so probably real
                                else:
                                    outFASTQfile.write(writeFastQ(read))

                            # Insertion or deletion but preserves frame
                            elif (op == 1 or op == 2) and oplen % 3 == 0:
                                outFASTQfile.write(writeFastQ(read))
                            marker += oplen

                    # has more than one indel, not dealing with these in this version, so these are excluded
                    else:   
                        ERRfile.write(writeFastQ(read))
                # has some other kind of issue, but not indel
                else:   
                    outFASTQfile.write(writeFastQ(read))
            # No indels, just matches
            else:
                outFASTQfile.write(writeFastQ(read))

    outFASTQfile.close()
    ERRfile.close()
    changeFile.close()
    print "QC complete\n\n"





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
        elif opt in ("-f","--file"):    seqfile = arg

        
    if seqfile !='':
       processFiles(seqfile)
    else:
       print opts
       print sys.argv
       usage()
       sys.exit()



if __name__ == '__main__':
    main(sys.argv[1:])