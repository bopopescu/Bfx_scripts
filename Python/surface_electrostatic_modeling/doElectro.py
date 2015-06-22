import os,sys,getopt, subprocess
from chimera import runCommand as rc 
from chimera import replyobj

def processFiles(abname,pdb,ab_dx,imagedir,rootdir):
    # base_dir = os.getcwd()
    # files = os.listdir(base_dir)
    # for x in files: 
        # --ff=AMBER --apbs-input --nodebump --verbose
        #subprocess.call(["pdb2pqr","--ff","AMBER","--apbs-input","--nodebump",x,x[:-4]+'.pqr'])
        #subprocess.call(["apbs",x[:-4]+'.in'])

    # Chimera commands to color by electrostatic potential
    # surface
    # scolor #0 volume #2 perPixel true offset 1.4 cmap -5,red:0,white:5,blue
    # colorkey .8,.04 .98,.05 -5 red 0 white 5 blue
        #if x[-4:] == '.pdb':
            #ab = base_dir+'/'+x
            #ab_dx = base_dir+'/'+x[:-4]+'.pqr.dx'
    replyobj.status("Beginning to process ab")
    rc("open " + pdb)
    rc("open " + ab_dx)
    rc("surface")
    rc("color darkgrey")
    rc("set dcColor black")
    rc("background solid white")
    rc("scolor #0 volume #1 perPixel true offset 1.4 cmap -5,red:0,white:5,blue")
    #rc("colorkey 0.8,0.4 0.98,0.05 -5 red 0 white 5 blue")
    rc("colorkey 0.8,0.08 0.98,0.05 -5 red 0 white 5 blue")
    rc("turn 0,0,1 -75")
    rc("windowsize 1050 700")
    rc("focus")
    png_name = abname+"_left.png"
    replyobj.status("Writing PNG file")
    rc("cd "+imagedir)
    rc("copy file "+png_name+" supersample 3")
    rc("turn 0,1,0 190")
    rc("focus")
    png_name = abname+"_right.png"
    replyobj.status("Writing PNG file")
    rc("copy file "+png_name+" supersample 3")
    rc("turn 1,0,0 90")
    rc("focus")
    png_name = abname+"_top.png"
    replyobj.status("Writing PNG file")
    rc("copy file "+png_name+" supersample 3")
    rc("reset")
    rc("close session")
    replyobj.status("Closing session")
    rc("close all")
    replyobj.status("Finished processing...exiting")
    sys.exit()

def main(argv):
    abname = pdb = ab_dx = imagedir = rootdir = ''
    #Read command line arguments
    print "\\\\\\\\\\\ WELCOME to chimeraImageGenerator \\\\\\\\\\\\n"
    try:
        opts, args = getopt.getopt(argv,"f:e:d:r:",["pdb,dx,imagedir,rootdir"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-f": 
            abname = arg[:-4]
            pdb = arg
        elif opt == "-e":
            ab_dx = arg
        elif opt == "-d":
            imagedir = arg
        elif opt == "-r":
            rootdir = arg
    if abname !="" and pdb!="" and ab_dx!="" and imagedir !="" and rootdir!="":
        print ">>>>>>>>>>>> Got chimera args\nabname: ",abname,"\npdb: ",pdb,"\nab_dx: ",ab_dx,"\nimagedir: ",imagedir,"\nrootdir: ",rootdir
        processFiles(abname,pdb,ab_dx,imagedir,rootdir)
    else:
            print "There was an error, not right args\n"
            sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])