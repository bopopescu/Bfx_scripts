import sys,os,glob,getopt,subprocess

def processFiles(gp41,gp120,seqs):
    datapath='input/'
    path=""
    binpath = "~/Documents/Monogram/Terri/patch_analysis/multiclade_project/final/data/"
    for datafile in glob.glob (os.path.join(datapath,"*")):            # Read all filenames in directory
        a = datafile.split('/')[1].split('_')
        mab = a[0]
    #    print "In "+os.getcwd()
        process = subprocess.Popen(["python","pepa_seq_prep6.py","-i",datafile,"-s",seqs,"-p",gp41,"-q",gp120,"-m",mab])
        process.wait()
         
            

def main(argv):
    gp41=''
    gp120=''
    seqs=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"p:q:s:",["gp41patchfilename,gp120patchfilename,seq_filename"])
    except getopt.GetoptError:
        print "There was an error\n"
        sys.exit(2)
    print "batch opts = ",opts
    for opt, arg in opts:
        if opt == "-p": gp41 = arg
        elif opt =="-q": gp120 = arg
        elif opt =="-s": seqs = arg
        
    if gp41 !='' and gp120 !='' and seqs !='':
        processFiles(gp41,gp120,seqs)
    else:
        print "There was an error\n"
        sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])