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


##############################
# Declare 'global' variables #
##############################
ALA_X_ALA={"A":"110.2","D":"144.1","C":"140.4","E":"174.7","F":"200.7","G":"78.7","H":"181.9","I":"185.0","K":"205.7","L":"183.1","M":"200.1","N":"146.4","P":"141.9","Q":"178.6","R":"229.0","S":"117.2","T":"138.7","V":"153.7","W":"240.5","Y":"213.7","X":"10"}
tripeptide = {}


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


#########################################################
# processHOA                                            #
# Read HOA datafile and load into tripeptide dictionary #
#########################################################
def processHOA(hoa_file):
    HOAFile = open(hoa_file,'r')
    for line in HOAFile:
        tripeptide[line.split('\t')[0].replace('\s','').rstrip()] = line.split('\t')[1].replace('\s','').rstrip()
    HOAFile.close()    
    return


################################################################
# checkHOAKey                                                  #
# If tripeptide does not have the key, default to the          #
# ALA-X-ALA value for that amino acid and add it to dictionary #
################################################################
def checkHOAKey(tri,aa):
    if tripeptide.has_key(tri): return
    else: 
        tripeptide[aa] = ALA_X_ALA[aa]
        return


##########################################################
# getASAResult                                           #
# Performs proper conversion depending upon option value #
# and returns properly formatted value                   #
##########################################################
def getASAResult(num,denom,option):
    if   option == 'r': asa_result = ("%.5f" % (Decimal(num)/Decimal(denom)))   # Relative ASA option
    elif option == 'a': asa_result = ("%1.0f" % (Decimal(num)*Decimal(denom)))  # Absolute ASA option
    return asa_result


############################################################
# getNTermResidue                                          #
# Properly construct N-term residue for lookup in HOA data #
# Works with any N-term, including those that result from  #
# interchain breaks in PDB / DSSP files                    #
############################################################
def getNTermResidue(i,asalines):
    aa1 = asalines[i-1][0]
    aa2 = asalines[i][0]
    asa = asalines[i-1].split('\t')[1].rstrip()
    if aa1 == 'X': asa = 0
    tri = '-'.join(['X',aa1,aa2])
    return (tri,aa1,asa)


###########################################################
# getCTermResidue                                          #
# Properly construct C-term residue for lookup in HOA data #
# Works with any C-term, including those that result from  #
# interchain breaks in PDB / DSSP files                    #
############################################################
def getCTermResidue(i,asalines):
    aa1 = asalines[i-2][0]
    aa2 = asalines[i-1][0]
    asa = asalines[i-1].split('\t')[1].rstrip()
    tri = '-'.join([aa1,aa2,'X'])
    if aa2 == 'X': asa = 0
    return (tri,aa2,asa)


################################################
# processFiles                                 #
# Open input files, read and process each line #
################################################
def processFiles(asa_file,option,mode,hoa_file=''):

    print "Beginning normalization...\n"
    RESULT = open(asa_file+".normalized","w")

    # Read in all the ASA data into an array
    ASAFile = open(asa_file,'r')
    asalines = ASAFile.readlines()
    ASAFile.close()
    asa_result = 0

    if mode=='hoa':

        # Read in all of the HOA data into a dictionary
        processHOA(hoa_file)

        # Loop through ASA file
        for i in range(1,len(asalines)+1):
            if i == 1:  # First residue
                tri,aa,asa = getNTermResidue(i,asalines)
            elif i == len(asalines):   # Last residue
                tri,aa,asa = getCTermResidue(i,asalines)
            else:
                aa1 = asalines[i-2][0]
                aa2 = asalines[i-1][0]
                aa3 = asalines[i][0]
                asa = asalines[i-1].split('\t')[1].rstrip()
                tri = '-'.join([aa1,aa2,aa3])
                if aa2 == 'X': asa = 0

                # Handle a DSSP chain break character, which is a '!'
                if aa3 == '!':   
                    tri,aa,asa = getCTermResidue(i,asalines)
                elif aa2 == '!':
                    RESULT.write(aa2+'\tchain break\n')
                    continue
                elif aa1 == '!':
                    tri,aa,asa = getNTermResidue(i,asalines)
                else: aa = aa2

            checkHOAKey(tri,aa)
            asa_result = getASAResult(asa,tripeptide[tri],option)
            RESULT.write(aa+'\t'+asa_result+'\n')

    else:
        # Mode is Ala-X-Ala
        for i in range(0,len(asalines)):
            aa1 = asalines[i][0]
            asa = asalines[i].split('\t')[1].rstrip()
            if aa1 == 'X': asa = 0
            asa_result = getASAResult(asa,ALA_X_ALA[aa1],option)
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
            if opt == "-h": hoa_file = arg
            elif opt == "-a": asa_file = arg
            elif opt == "-o": option = arg.lower()
   
        if hoa_file !='' and asa_file !='':
            processFiles(asa_file,option,mode,hoa_file)
        else:
            usage()
            sys.exit()
    elif len(opts) == 2:
        mode = 'Ala'
        for opt, arg in opts:
            if opt == "-a": asa_file = arg
            elif opt == "-o": option = arg.lower()
        if asa_file !='':
            processFiles(asa_file,option,mode)
    else:
        usage()
        sys.exit()


if __name__ == '__main__':
    main(sys.argv[1:])