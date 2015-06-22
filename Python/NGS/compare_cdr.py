import sys,os,getopt
from operator import itemgetter
from decimal import *


def process_files(xoma_cdr_input, ngs_cdr_input,N_ngs):
    xoma_cdrs = {}
    ngs_cdrs ={}
    xcdr_freq = {}
    not_found = []
    heavy_family = {}
    studies = {}
    unseen_cdrs = {}
    getcontext().prec = 3

    f1 = open(xoma_cdr_input,'r')
    f2 = open(ngs_cdr_input,'r')


    print "\n\nReading xoma cdrs...."
    for line in f1.readlines():
        a = line.split()
        cdr = a[4]
        if xcdr_freq.has_key(cdr): xcdr_freq[cdr] = xcdr_freq[cdr]+1
        else: xcdr_freq[cdr] = 1
        if heavy_family.has_key(cdr):
            if heavy_family[cdr].find(a[2]) == -1: 
                heavy_family[cdr] = heavy_family[cdr]+","+a[2]
        else:
            heavy_family[cdr] = a[2]
        if studies.has_key(cdr):
            if studies[cdr].find(a[0]) == -1:
                studies[cdr] = studies[cdr]+", "+str(a[0])
        else: studies[cdr] = str(a[0])
        xoma_cdrs[a[1]] = [cdr,a[0],a[2],a[3]]
    f1.close()



    print "Reading ngs cdrs..."
    for line in f2.readlines():
        a = line.split()
        ngs_cdrs[a[0]] = a[1]
    f2.close()

    f3 = open("not_seen_by_xoma.txt","w")
    for cdr in ngs_cdrs.keys():
    #    print cdr, "is present? ",str(xcdr_freq.has_key(cdr))
        if not xcdr_freq.has_key(cdr): 
            f3.write(cdr+"\t"+str(ngs_cdrs[cdr])+"\t"+str((Decimal(ngs_cdrs[cdr])/Decimal(N_ngs))*100)+"\n")
            unseen_cdrs[cdr] = [str(ngs_cdrs[cdr]),str((Decimal(ngs_cdrs[cdr])/Decimal(N_ngs))*100)]
    f3.close()

    N_xoma = len(xoma_cdrs)
    print "N_ngs = ",N_ngs," N_xoma = ",N_xoma

    
    print "There are ",str(len(xcdr_freq))," unique CDRs seen by XOMA out of a total of ",len(xoma_cdrs)," sequences"

    print "\nXoma_freq\t\tNGS_freq\t\tHeavy family\tStudies\tCDR\n"
    #for cdr in xcdr_freq.keys():
    for cdr, index in sorted(xcdr_freq.items(), key=itemgetter(1), reverse=True): #order the CDR3s by abundance
        if ngs_cdrs.has_key(cdr):
            print "\t".join([str((Decimal(xcdr_freq[cdr])/Decimal(N_xoma))*100),"\t",str((Decimal(ngs_cdrs[cdr])/Decimal(N_ngs))*100),"\t\t",heavy_family[cdr],studies[cdr],cdr])
        else:
            not_found.append(cdr)

    print "\nThe following CDR's were found by XOMA but not seen by NGS\n"
    for x in not_found: 
        print x,"\t",studies[x]

    f4 = open("not_seen_by_xoma_DNA.txt","w")
    f5 = open("bins.dna.cdr3.csv","r")
    for line in f5.readlines():
        a = line.split(',')
        if not xcdr_freq.has_key(a[0]):
            f4.write(a[0]+"\t"+str(((Decimal(a[2])/Decimal(N_ngs))*100))+"\t"+a[1]+"\n")
    f4.close()
    f5.close()



def main(argv):
    xoma_cdr_input =''
    ngs_cdr_input = ''
    N_ngs = ''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"x:s:n:",["xoma:seq:number"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-x":
            xoma_cdr_input = arg    
        elif opt == "-s":
            ngs_cdr_input = arg
        elif opt == "-n":
            N_ngs = arg
        else:
            print "There was an error, not right args\n"
            print opts
            sys.exit()
    if xoma_cdr_input != '' and ngs_cdr_input != '': 
        process_files(xoma_cdr_input, ngs_cdr_input, N_ngs)
       
   

if __name__ == '__main__':
    main(sys.argv[1:])