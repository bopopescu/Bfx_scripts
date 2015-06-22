#!python
import sys,os,getopt,re
from operator import itemgetter
from decimal import *


def process_files(cluster_input,individual_input):
    getcontext().prec = 6
    raw_patches = {}
    residue_result = {}

    f1 = open(cluster_input,'r')
    f2 = open(individual_input,'r')

    print "\n\nReading PEASE data...."
    c = 1
    for line in f1.readlines():
        a = line.split()
        if a[0] != "chain":
            p = a[2]
            b = re.sub("[A-Z]",'',a[1]).split(',')
            raw_patches[c] = {}
            for idx in b:
                raw_patches[c][int(idx)] = Decimal(p)
            c += 1
    f1.close()

    for line in f2.readlines():
        a =line.split()
        if a[0] != 'chain':
            residue_result[int(a[1])] = [Decimal(a[3]),a[2]]
    f2.close()

    patch_calc = {}
    for idx in raw_patches.keys():
        for pos in raw_patches[idx].keys():
            if patch_calc.has_key(pos): patch_calc[pos] = Decimal(patch_calc[pos])+(Decimal(raw_patches[idx][pos])*Decimal(residue_result[pos][0]))
            else: patch_calc[pos] = Decimal(raw_patches[idx][pos])*Decimal(residue_result[pos][0])

    maxval = 0
    for pos in residue_result.keys():
        if patch_calc.has_key(pos):
            if patch_calc[pos] > maxval: maxval = patch_calc[pos]
        else: patch_calc[pos] = 0

    print patch_calc
    for pos in patch_calc.keys():
        patch_calc[pos] = (Decimal(patch_calc[pos])/Decimal(maxval))*100

    f3 = open("epitope_map.txt","w")
    positions = patch_calc.keys()
    positions.sort()
    title_line =""
    val_line = ""
    
    for pos in positions:
        title_line += str(pos)+"\t"
        if str(patch_calc[pos]).find('+') != -1: val_line += "0\t"
        else:  val_line += str(patch_calc[pos])+"\t"
    title_line += "\n"
    val_line += "\n"
    f3.write(title_line+val_line)
    f3.close()


def main(argv):
    cluster_input =''
    individual_input = ''
    
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"c:i:",["cluster:indiv:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-c":
            cluster_input = arg    
        elif opt == "-i":
            individual_input = arg
        else:
            print "There was an error, not right args\n"
            print opts
            sys.exit()
    if cluster_input != '' and individual_input != '': 
        process_files(cluster_input,individual_input)
       



if __name__ == '__main__':
    main(sys.argv[1:])