import sys,os,glob,getopt,subprocess

def processFiles(patchlist,seqs):
    datapath='patchfiles/'

    f = open(patchlist,'r')
    for line in f:
        a = line.rstrip().split('\t')
        process = subprocess.Popen(["python","pepa_seq_prep7.5.py","-s",seqs,"-p",datapath+a[0],"-q",datapath+a[1],"-d",a[2]])
        process.wait()
         
            

def main(argv):
    patchlist=''
    seqs=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"p:s:",["patchlist_idx,allseqs_filename"])
    except getopt.GetoptError:
        print "There was an error\n"
        sys.exit(2)
    print "batch opts = ",opts
    for opt, arg in opts:
        if opt == "-p": patchlist = arg
        elif opt =="-s": seqs = arg
        
    if patchlist !='' and seqs !='':
        processFiles(patchlist,seqs)
    else:
        print "There was an error\n"
        sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])