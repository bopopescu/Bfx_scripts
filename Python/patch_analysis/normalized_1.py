#!/usr/bin/env python
# encoding: utf-8
"""
normalized.py

A python script to change the absolute value of ASA to RSA
or RSA to absolute ASA.

By Mark Evans ( evans.markc@gmail.com )
Created 11.30.2011

Based upon the normalized.pl perl script by Hemajit Singh found here 
http://hoa.netasa.org
BMC Structural Biology 2009, 9:25
doi:10.1186/1472-6807-9-25

"""

import sys, os, getopt
from decimal import *

# ALA-X-ALA values are assigned to a dictionary
ALA_X_ALA={"A":"110.2","D":"144.1","C":"140.4","E":"174.7","F":"200.7","G":"78.7","H":"181.9","I":"185.0","K":"205.7","L":"183.1","M":"200.1","N":"146.4","P":"141.9","Q":"178.6","R":"229.0","S":"117.2","T":"138.7","V":"153.7","W":"240.5","Y":"213.7","X":"10"}


##################################
# Define usage syntax for script #
##################################   
def usage():
    code =  "\n\n#############\nASA normalization Utility\n"
    code += "Usage: $ python normalized.py -o normalization_method -a ASA_filename [-h HOA_filename] \n"
    code += "-o [r,a] r = convert absolute ASA to relative ASA, a = convert relative ASA to absolute ASA\n"
    code += "-a [ASA_filename] path to ASA file\n"
    code += "-h [HOA_filename] optional path to HOA file. If not included, default ALA-X-ALA values will be used\n"
    code += "ASA files should be in two tab-seperated columns. First column for residue, second for ASA \n#############\n\n"
    print code
    return


################################################
# processFiles                                 #
# Open input files, read and process each line #
################################################
def processFiles(asa_file,option,mode,hoa_file=''):
    
    print "Beginning normalization...\n"
    RESULT = open(asa_file+".normalized","w")
    
    ASAFile = open(asa_file,'r')
    asalines = ASAFile.readlines()
    ASAFile.close()
    asa_result = 0

    if mode=='hoa':
        ter = 'X'
        HOAFile = open(hoa_file,'r')
        tripeptide = {}
        for line in HOAFile:
            tripeptide[line.split('\t')[0].replace('\s','').rstrip()] = line.split('\t')[1].replace('\s','').rstrip()
        HOAFile.close()

        # first residue
        aa1 = asalines[0][0]
        aa2 = asalines[1][0]
        asa = asalines[0].split('\t')[1].rstrip()
        if aa1 == 'X': asa = 0
        tri = '-'.join([ter,aa1,aa2])
        
        if not tripeptide.has_key(tri): # ALA-X-ALA value is used if no HOA tripeptide is available
            tri = aa1
            tripeptide[tri] = ALA_X_ALA[aa1]
        
        if   option == 'r': asa_result = ("%.5f" % (Decimal(asa)/Decimal(tripeptide[tri])))  # Relative ASA option
        elif option == 'a': asa_result = ("%1.0f" % (Decimal(asa)*Decimal(tripeptide[tri])))  # Absolute ASA option
        RESULT.write(aa1+'\t'+asa_result+'\n')

        # Loop for rest of ASA file
        for i in range(2,len(asalines)):
            aa1 = asalines[i-2][0]
            aa2 = asalines[i-1][0]
            aa3 = asalines[i][0]
            asa = asalines[i-1].split('\t')[1].rstrip()
            tri = '-'.join([aa1,aa2,aa3])
            if aa2 == 'X': asa = 0
            
            if not tripeptide.has_key(tri): # ALA-X-ALA value is used if no HOA tripeptide is available
                tri = aa2
                tripeptide[tri] = ALA_X_ALA[aa2]
            
            if   option == 'r': asa_result = ("%.5f" % (Decimal(asa)/Decimal(tripeptide[tri])))  # Relative ASA option
            elif option == 'a': asa_result = ("%1.0f" % (Decimal(asa)*Decimal(tripeptide[tri])))  # Absolute ASA option
            RESULT.write(aa2+'\t'+asa_result+'\n')
        
        # Last residue
        aa1 = asalines[len(asalines)-2][0]
        aa2 = asalines[len(asalines)-1][0]
        asa = asalines[len(asalines)-1].split('\t')[1].rstrip()
        tri = '-'.join([aa1,aa2,ter])
        if aa2 == 'X': asa = 0

        if not tripeptide.has_key(tri): # ALA-X-ALA value is used if no HOA tripeptide is available
            tri = aa2
            tripeptide[tri] = ALA_X_ALA[aa2]
        
        if   option == 'r': asa_result = ("%.5f" % (Decimal(asa)/Decimal(tripeptide[tri])))  # Relative ASA option
        elif option == 'a': asa_result = ("%1.0f" % (Decimal(asa)*Decimal(tripeptide[tri])))  # Absolute ASA option
        RESULT.write(aa2+'\t'+asa_result+'\n')

    else:
        # Mode is Ala-X-Ala
        for i in range(0,len(asalines)):
            aa1 = asalines[i][0]
            asa = asalines[i].split('\t')[1].rstrip()
            if aa1 == 'X': asa = 0
            
            if   option == 'r': asa_result = ("%.5f" % (Decimal(asa)/Decimal(ALA_X_ALA[aa1])))   # Relative ASA option
            elif option == 'a': asa_result = ("%1.0f" % (Decimal(asa)*Decimal(ALA_X_ALA[aa1])))   # Absolute ASA option
            RESULT.write(aa1+'\t'+asa_result+'\n')
    
    RESULT.close()
    print "Normalization process complete\n\n"
    sys.exit()




#################################################################
# Main                                                          #
# Get cmd line arguments and call the main processing function  #
#################################################################
def main(argv):
    asa_file =''
    hoa_file=''
    option = ''
    mode =''

    # Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"o:a:h:",["option,asa_file,hoa_file"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    if len(opts)== 3:
        mode = 'hoa'
        for opt, arg in opts:
            if opt.lower() in ("-h","--hoa_file"): hoa_file = arg
            elif opt.lower() in ("-a","--asa_file"): asa_file = arg
            elif opt.lower() in ("-o","--option"): option = arg
   
        if hoa_file !='' and asa_file !='':
            processFiles(asa_file,option,mode,hoa_file)
        else:
            usage()
            sys.exit()
    elif len(opts) == 2:
        mode = 'Ala'
        for opt, arg in opts:
            if opt.lower() in ("-a","--asa_file"): asa_file = arg
            elif opt.lower() in ("-o","--option"): option = arg
        if asa_file !='':
            processFiles(asa_file,option,mode)
    else:
        usage()
        sys.exit()


if __name__ == '__main__':
    main(sys.argv[1:])