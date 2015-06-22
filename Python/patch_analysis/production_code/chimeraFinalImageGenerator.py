
"""
chimeraFinalImageGenerator.py
Script to be run via Chimera to generate 3D images of HIV gp120 and gp41

Created by Mark Evans on 2012-09-06.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,getopt,string,subprocess, random
from decimal import *
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

def callChimeraGP120(epitope,fn,dirname,rootdir):
    base_dir = os.getcwd()
    trimer = rootdir+'/gp120CPHmodel_trimer.pdb'
    monomer = rootdir+'/gp120_CPHModel.pdb'
    print "trimer path= ",trimer
    print "monomer path= ",monomer
    print "epitope = ",epitope

    # Create 1st trimer image
    replyobj.status("Processing trimer") # show what file we're working on
    rc("open " + trimer)
    rc("preset apply publication 3")    # make everything look nice
    rc("turn 0,0,1 15")                  # tweak the rotation a bit
    rc("focus")                          # center/zoom 
    rc("split #0")                       # split model into individually addressable sub models
    rc("surf")                           # generate surface
    rc("color dark gray")
    rc("surftransp 70")                  # make the surface a little bit see-through
    rc("focus")                          # readjust zoom
    # Color epitopes
    if epitope != 'No Result' and epitope != '':
        print "@@@@@@@@@@ cmd = color red #0.1 :",epitope
        rc("color red #0.1 :"+epitope)
        rc("color red #0.2 :"+epitope)
        rc("color red #0.3 :"+epitope)
    
    # save image to a file that ends in .png rather than .pdb
    rc("cd "+dirname)
    rc("windowsize 1050 700") # 2100 1400
    png_name = fn + "-final.trimer-top.png"
    rc("copy file " + png_name + " supersample 3")
    
    # Rotate trimer for second image
    rc("turn 0,1,0 90")
    rc("turn 0,0,1 90")
    rc("turn 0,1,0 30")
    rc("~modeldisp #0.3")               # hide the third molec
    rc("focus")
    rc("color light blue #0.1")         # color one chain light blue for contrast
    if epitope != 'No Result' and epitope != '':
        rc("color red #0.1 :"+epitope)
    
    png_name = fn+"-final.trimer-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("reset")
    rc("close session")

    # Open monomer model
    rc("open "+monomer)
    rc("turn 1,0,0 125")
    rc("preset apply publication 3")    # make everything look nice
    rc('color tan')
    rc("surface")
    rc("surftransp 70")
    rc("focus")
    if epitope != 'No Result' and epitope != '':
        rc("color red #0 :"+epitope)
    png_name = fn+"-final.front-side.png"
    rc("windowsize 1050 700 ")
    rc("copy file " + png_name + " supersample 3")
    rc("turn 0,1,0 180")
    png_name = fn+"-final.back-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("reset")
    rc("close all")

    return


def callChimeraGP41(epitope,fn,dirname,rootdir):
    base_dir = os.getcwd()
    trimer = rootdir+'/gp41_cphtrimer.pdb'
    monomer = rootdir+'/gp41_cphmodel.pdb'
    print "trimer path= ",trimer
    print "monomer path= ",monomer
    print "epitope = ",epitope

    # Create 1st trimer image
    replyobj.status("Processing trimer") # show what file we're working on
    rc("open " + trimer)
    rc("windowsize 1050 700")
    rc("preset apply publication 3")    # make everything look nice
    rc("turn 0,1,0 90")                  # tweak the rotation a bit
    rc("surf")                           # generate surface
    rc("focus")
    rc("color tan")
    rc("surftransp 60")                  # make the surface a little bit see-through
    rc("scale 0.70")                    # Fit into window
    rc("ribcolor rosy brown #0.1")  # rosey brown
    rc("ribcolor dark khaki #0.2")  # dark khaki
    # Color epitopes
    combo=''
    if epitope != 'No Result' and epitope != '':
        print "@@@@@@@@@@ cmd = color red #0.1 :",epitope
        rc("color red #0.1 :"+epitope)
        rc("color red #0.2 :"+epitope)
        rc("color red #0.3 :"+epitope)
        combo += epitope+','
    
    combo = combo[:-1]

    # Show ball & stick for epitopes
    rc("show #0.1 :"+combo)
    rc("show #0.2 :"+combo)
    rc("show #0.3 :"+combo)
    rc("represent bs #0.1 :"+combo)
    rc("represent bs #0.2 :"+combo)
    rc("represent bs #0.3 :"+combo)

    # save image to a file that ends in .png rather than .pdb
    rc("cd "+dirname)
    
    png_name = fn + "-final.trimer-side.png"
    rc("copy file " + png_name + " supersample 3")
    
    # Rotate trimer for second image
    rc("turn 0,1,0 90")
    png_name = fn+"-final.trimer-top.png"
    
    rc("copy file " + png_name + " supersample 3")
    rc("reset")
    rc("close session")

    # Open monomer model
    rc("open "+monomer)
    rc("surface")
    rc("turn 0,0,1 -115")
    rc("turn 0,1,0 90")
    rc("scale 0.70")
    rc("preset apply publication 3")    # make everything look nice
    rc('color tan')
    rc("surftransp 60")
    combo = ''
    if epitope != 'No Result' and epitope != '':
        rc("color red #0 :"+epitope)
        combo += epitope+','
    combo = combo[:-1]

    # Show ball & stick for epitopes
    rc("show #0 :"+combo)
    rc("represent bs #0 :"+combo)

    png_name = fn+"-final.front-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("turn 0,1,0 180")
    png_name = fn+"-final.back-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("reset")
    rc("close all")

    return



def processFiles(molec,epitope,fn,dirname,rootdir):
    print "\\\\\\\\ Writing images to ",dirname,'\n'
    if molec == 'gp120': callChimeraGP120(epitope.replace('_',' '),fn,dirname,rootdir)
    elif molec == 'gp41': callChimeraGP41(epitope.replace('_',' '),fn,dirname,rootdir)
    sys.exit()

def main(argv):
    molec=epitope=fn=dirname=rootdir = ''
    #Read command line arguments
    print "\\\\\\\\\\\ WELCOME to chimeraFinalImageGenerator \\\\\\\\\\\\n"
    try:
        opts, args = getopt.getopt(argv,"m:e:f:d:r:",["molec,epitope,fn,dirname,rootdir"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-m": 
            molec = arg
        elif opt == "-e":
            epitope = arg
        elif opt == "-f":
            fn = arg
        elif opt == "-d":
            dirname = arg
        elif opt == "-r":
            rootdir = arg
    if molec !="" and epitope!="" and fn!="" and dirname !="" and rootdir!="":
        processFiles(molec,epitope,fn,dirname,rootdir)
    else:
            print "There was an error, not right args\n"
            sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])