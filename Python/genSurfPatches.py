#!/usr/bin/env python
# encoding: utf-8
"""
genSurfPatches.py

Created by Mark Evans on 2011-11-29.
Copyright (c) 2011 __Monogram Biosciences__. All rights reserved.

Purpose is to take PDB model file, surface accessibility and atomic depth
data and generate non-redundant list of surface patches

Atomic depth data derived from DPX
http://hydra.icgeb.trieste.it/dpx/

Surface accessibility data derived from DSSP
http://swift.cmbi.ru.nl/gv/dssp/

DSSP Output file explained (from http://alumni.cs.ucr.edu/~vladimir/lib/formats/3DSSP_Output_File.html )
Columns Header  Definition
1 - 5       #           Line number
7 - 11      RESIDUE     Residue number, as in PDB file
12                      PDB chain
14          AA          Single letter code amino acid, ! for missing residues, X for non-standard
17          STRUCTURE   Secondary structure assignment, see the following legend table
19 - 25     STRUCTURE    
            BP1  
            BP2  
36 - 38     ACC         Accessible surface area
            N-H->O   
            O->H-N   
            N-H->O   
            O->H-N   
            TCO  
            KAPPA    
            ALPHA    
104 - 109   PHI         Angle
110 - 115   PSI         Angle
116 - 122   X-CA        C Alpha Coords
123 - 129   Y-CA        C Alpha Coords
130 - 136   Z-CA        C Alpha Coords

Need to create hybrid data file that merges DSSP and Dpx outputs that is tab delimited
which will be used for input here.  Format of datafile should be 5 columns, tab-delimited
Residue\tChain\tAA\tACC\tdpx_Avg without header row

"""

import sys, os, getopt
from decimal import *
from time import localtime, strftime
from Bio.PDB import *

##############################
# Declare 'global' variables #
##############################
asa_data={}
ALA_X_ALA={"A":"110.2","D":"144.1","C":"140.4","E":"174.7","F":"200.7","G":"78.7","H":"181.9","I":"185.0","K":"205.7","L":"183.1","M":"200.1","N":"146.4","P":"141.9","Q":"178.6","R":"229.0","S":"117.2","T":"138.7","V":"153.7","W":"240.5","Y":"213.7","X":"10"}
tripeptide = {}
aa_lookup = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

##################################
# Define usage syntax for script #
##################################   
def usage():
    code =  "\n\n#############\ngenSurfPatches Utility\n"
    code += "Usage: $ python genSurfPatches.py -p pdbmodel_filename -d data_filename -s patch_size \n"
    code += "-p [pdb] name of structural model file in PDB format\n"
    code += "-d [data] name of file containing surface accessibility and dpx data\n"
    code += "-s [size] of patch in angstroms. Antibody epitope size of 20 is default\n"
    code += "-h [help]\n#############\n\n"
    print code
    return


#########################
# Returns current time  #
#########################
def gt():
    return strftime("%d %b %Y %H:%M:%S",localtime())


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
def normalizeASA(asa_file,option,mode,hoa_file=''):

    print "Beginning normalization...\n"

    # Read in all the ASA data into an array
    ASAFile = open(asa_file,'r')
    asalines = ASAFile.readlines()
    ASAFile.close()
    os.remove(asa_file)
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
               #     RESULT.write(aa2+'\tchain break\n')
                    asa_data[i-1]['rsa'] = 'chain break'
                    continue
                elif aa1 == '!':
                    tri,aa,asa = getNTermResidue(i,asalines)
                else: aa = aa2

            checkHOAKey(tri,aa)
            asa_result = getASAResult(asa,tripeptide[tri],option)
  #          print asa,tripeptide[tri],option,asa_result
            asa_data[i-1]['rsa'] = asa_result

    # Only need this if want to use Ala-X-Ala values
    else:
        # Mode is Ala-X-Ala
        for i in range(0,len(asalines)):
            aa1 = asalines[i][0]
            asa = asalines[i].split('\t')[1].rstrip()
            if aa1 == 'X': asa = 0
            asa_result = getASAResult(asa,ALA_X_ALA[aa1],option)
            asa_data[i]['rsa'] = asa_result
    
    print "Normalization process complete\n\n"
    return


###############################################
# processFiles                                #
# Open input file, read and process each line #
###############################################
def processFiles(pdbfile,datafile,size):
    
    # Load all residue data (ASA, dpx)
    DATA = open(datafile,'r')
    c = 0
    for line in DATA:
        vals = line.rstrip().split('\t')
        asa_data[c] = {'pos':vals[0],'chain':vals[1],'aa':vals[2],'asa':vals[3],'dpx':vals[4]}
        c += 1

    # Write out ASA tmp file
    pos = asa_data.keys()
    pos.sort()
    TMP = open('tmp_asa_input.txt','w')
    for x in pos:
        if asa_data[x]['aa']!='!': TMP.write(asa_data[x]['aa']+'\t'+asa_data[x]['asa']+'\n') 
        else:  TMP.write('!\t0\n')
    TMP.close()

    # Normalize ASA values using HOA
    normalizeASA('tmp_asa_input.txt','r','hoa','HOA.txt')

    # Remove chain breaks from the residue data structure
    print "length of asa_data = "+str(len(asa_data))
    asa_data2 = {}
    c = 0
    for x in sorted(asa_data.keys()):
        if asa_data[x]['aa'] == '!': 
            print "chain break at "+str(x)+",removing..."
        else: 
            asa_data2[c] = asa_data[x]
            c += 1
    print "length of asa_data2 now = "+str(len(asa_data2))
#    print str(asa_data2)

    # Load PDB model
    print "Assuming for now that PDB only has one chain...\n"
    p = PDBParser()
    structure = p.get_structure('Antigen_model',pdbfile)

    # Begin going through structure to make patches
    final_patches = {}
    residues = Selection.unfold_entities(structure,'R')
    PATCHLIST = open('dev_patchlist.txt','w')
    PATCHLIST.write("patch_num\tpatch_readable\tpatch\tRSA\tDPX\tDist\n")
    pnum = 1
    for i in range(0,len(residues)):
        #
        # Add check here that the current residue meets "surface" critera
        #
        # if pdb file contains 'other' residues (Het, fusion protein aa, etc) then length of residues
        # does not correspond to length/pos of asa_data2, because we would have trimmed that out of the ASA file
        # easiest solution is to trim out 'bad' fusion aa from the pdb file so that the residue length matches the ASA file
        # (..I think..)
        #
        if residues[i].resname != 'HOH':  # make sure residue is not het water
            if (float(asa_data2[i]['rsa']) > 0.2) and (float(asa_data2[i]['dpx']) <= 2.5): # 0, 1.5 or 2.5 Angstrom depths are options
                current_residue_ca = residues[i]['CA']
                working_patch = [int(residues[i].id[1]),]
                working_readable = [aa_lookup[residues[i].resname]+str(residues[i].id[1]),]
                working_dpx = [asa_data2[i]['dpx'].rstrip(),]
                working_rsa = [asa_data2[i]['rsa'],]
                working_dist = []
                for j in range(0,len(residues)):
                    if i != j :
                        if residues[j].resname != 'HOH':  # make sure residue is not het water
#                            print 'i=',str(i),' and j=',str(j)

                            if current_residue_ca - residues[j]['CA'] <= int(size):
 #                               print "residue ",residues[j].id[1]
 #                               print ' distance:',str(current_residue_ca - residues[j]['CA'])
 #                               print ' rsa:',str(asa_data2[j]['rsa'])
 #                               print ' DPX: ',str(asa_data2[j]['dpx'])
                                if (float(asa_data2[j]['rsa']) > 0.2) and (float(asa_data2[j]['dpx']) <= 0):  # 0, 1.5 or 2.5 Angstrom depths are options
 #                                   print 'patch: ',residues[j].id[1],'RSA: ',str(asa_data2[j]['rsa']),' DPX: ',asa_data2[j]['dpx']
                                    working_patch.append(int(residues[j].id[1]))
                                    working_dpx.append(asa_data2[j]['dpx'].rstrip())
                                    working_rsa.append(asa_data2[j]['rsa'])
                                    working_dist.append(str(current_residue_ca - residues[j]['CA']))
                                    working_readable.append(aa_lookup[residues[j].resname]+str(residues[j].id[1]))
                        else: print '==========> residue j '+str(j)+' is HOH but asa_data2[i]='#+str(asa_data2[i])+'/n*****/n*****'
                if len(working_patch) != 1:
                    idx =''
                    for x in sorted(working_patch): idx = idx + str(x)

                    if not final_patches.has_key(idx): 
                        readable =''; wpatch=''; wrsa=''; wdpx=''; wdist='';
                        final_patches[idx] = str(working_patch.sort())
                        for index, val in enumerate(working_readable):  readable = readable + val + ','
                        for index, val in enumerate(working_patch): wpatch = wpatch + str(val) + ','
                        for index, val in enumerate(working_rsa): wrsa = wrsa + str(val) + ','
                        for index, val in enumerate(working_dpx): wdpx = wdpx + str(val) + ','
                        for index, val in enumerate(working_dist): wdist = wdist + str(val) + ','
                        PATCHLIST.write('patch '+str(pnum)+'\t'+readable[:-1]+'\t'+wpatch[:-1]+'\t'+wrsa[:-1]+'\t'+wdpx[:-1]+'\t'+wdist[:-1]+'\n')
                    pnum += 1
                del(working_patch)
        else: print '/n*****/n*****/n*****residue '+str(i)+' is HOH but asa_data2[i]='#+str(asa_data2[i])+'/n*****/n*****'
    print "length of final_patches ",str(len(final_patches))
    print "End of program\n"
    




#################################################################
# Main                                                          #
# Get cmd line arguments and call the main processing function  #
#################################################################
def main(argv):
    pdbfile =''
    datafile=''
    size = 20

    # Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"p:d:s:h",["pdb_filename,data_filename,patch_size,help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-p","--pdb"): pdbfile = arg
        elif opt in ("-d","--data"): datafile = arg
        elif opt in ("-s","--size"): size = Decimal(arg)
   
    if pdbfile !='' and datafile !='':
        processFiles(pdbfile,datafile,size)
    else:
        usage()
        sys.exit()


if __name__ == '__main__':
    main(sys.argv[1:])