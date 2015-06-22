import sys,os,getopt

xoma_cdrs = {}
xcdr_freq = {}

def extractCDR_freq (f):
    print "\n\nReading xoma cdrs...."
    for line in f.readlines():
        a = line.split()
        cdr = a[4]
        if xcdr_freq.has_key(cdr): xcdr_freq[cdr] = xcdr_freq[cdr]+1
        else: xcdr_freq[cdr] = 1
        xoma_cdrs[a[1]] = [cdr,a[0],a[2],a[3]]
    f.close()

def process_files(xoma_cdr_input):



    f1 = open(xoma_cdr_input,'r')
    f2 = open(xoma_cdr_input[:-4]+'.uniq.txt','w')

    extractCDR_freq(f1)

    print "There are ",str(len(xcdr_freq))," unique CDRs seen by XOMA out of a total of ",len(xoma_cdrs)," sequences"
    for cdr in xcdr_freq.keys():
        f2.write('C'+cdr+'W'+'\n')  # Adding C and W at beginning and end to match ab_toolbox CDR calls
    f2.close()
    print "process complete\n\n"

def main(argv):
    xoma_cdr_input =''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:",["filename:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-f":
            xoma_cdr_input = arg    
        else:
            print "There was an error, not right args\n"
            print opts
            sys.exit()
    if xoma_cdr_input != '': 
        process_files(xoma_cdr_input)

if __name__ == '__main__':
    main(sys.argv[1:])