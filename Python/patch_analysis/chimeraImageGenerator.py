
"""
chimeraImageGenerator.py
Script to be run via Chimera to generate 3D images of HIV gp120 and gp41

Created by Mark Evans on 2012-08-06.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,getopt,string,subprocess, random
from decimal import *
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages


def callChimeraGP120(pri_epitope,sec_epitope,fn,dirname,rootdir):
    base_dir = os.getcwd()
    trimer = rootdir+'/gp120CPHmodel_trimer.pdb'
    monomer = rootdir+'/gp120_CPHModel.pdb'
    print "trimer path= ",trimer
    print "monomer path= ",monomer
    print "pri_epitope = ",pri_epitope
    print "sec_epitope = ",sec_epitope
    print "destination directory = ",dirname

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
    if pri_epitope != 'No Result' and pri_epitope != '':
        print "@@@@@@@@@@ cmd = color red #0.1 :",pri_epitope
        rc("color red #0.1 :"+pri_epitope)
        rc("color red #0.2 :"+pri_epitope)
        rc("color red #0.3 :"+pri_epitope)
    if sec_epitope != 'No Result' and sec_epitope != '':
        print "@@@@@@@@@@ cmd = color blue #0.1 :",sec_epitope
        rc("color blue #0.1 :"+sec_epitope)
        rc("color blue #0.2 :"+sec_epitope)
        rc("color blue #0.3 :"+sec_epitope)
    # save image to a file that ends in .png rather than .pdb
    rc("cd "+dirname)
    rc("windowsize 1050 700") # 2100 1400
    png_name = fn + ".trimer-top.png"
    rc("copy file " + png_name + " supersample 3")
    
    # Rotate trimer for second image
    rc("turn 0,1,0 90")
    rc("turn 0,0,1 90")
    rc("turn 0,1,0 30")
    rc("~modeldisp #0.3")               # hide the third molec
    rc("focus")
    rc("color light blue #0.1")         # color one chain light blue for contrast
    if pri_epitope != 'No Result' and pri_epitope != '':
        rc("color red #0.1 :"+pri_epitope)
    if sec_epitope != 'No Result' and sec_epitope != '':
        rc("color blue #0.1 :"+sec_epitope)
    png_name = fn+".trimer-side.png"
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
    if pri_epitope != 'No Result' and pri_epitope != '':
        rc("color red #0 :"+pri_epitope)
    if sec_epitope != 'No Result' and sec_epitope != '':
        rc("color blue #0 :"+sec_epitope)
    png_name = fn+".front-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("turn 0,1,0 180")
    png_name = fn+".back-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("reset")
    rc("close all")

    return


def callChimeraGP41(pri_epitope,sec_epitope,fn,dirname,rootdir):
    base_dir = os.getcwd()
    trimer = rootdir+'/gp41_cphtrimer.pdb'
    monomer = rootdir+'/gp41_cphmodel.pdb'
    print "trimer path= ",trimer
    print "monomer path= ",monomer
    print "pri_epitope = ",pri_epitope
    print "sec_epitope = ",sec_epitope

    # Create 1st trimer image
    replyobj.status("Processing trimer") # show what file we're working on
    rc("open " + trimer)
    rc("windowsize 1050 700") # 2100 1400
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
    if pri_epitope != 'No Result' and pri_epitope != '':
        print "@@@@@@@@@@ cmd = color red #0.1 :",pri_epitope
        rc("color red #0.1 :"+pri_epitope)
        rc("color red #0.2 :"+pri_epitope)
        rc("color red #0.3 :"+pri_epitope)
        combo += pri_epitope+','
    if sec_epitope != 'No Result' and sec_epitope != '':
        print "@@@@@@@@@@ cmd = color blue #0.1 :",sec_epitope
        rc("color blue #0.1 :"+sec_epitope)
        rc("color blue #0.2 :"+sec_epitope)
        rc("color blue #0.3 :"+sec_epitope)
        combo += sec_epitope+","
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
    
    png_name = fn + ".trimer-side.png"
    rc("copy file " + png_name + " supersample 3")
    
    # Rotate trimer for second image
    rc("turn 0,1,0 90")
    png_name = fn+".trimer-top.png"
    
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
    if pri_epitope != 'No Result' and pri_epitope != '':
        rc("color red #0 :"+pri_epitope)
        combo += pri_epitope+','
    if sec_epitope != 'No Result' and sec_epitope != '':
        rc("color blue #0 :"+sec_epitope)
        combo += sec_epitope+','
    combo = combo[:-1]

    # Show ball & stick for epitopes
    rc("show #0 :"+combo)
    rc("represent bs #0 :"+combo)

    png_name = fn+".front-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("turn 0,1,0 180")
    png_name = fn+".back-side.png"
    rc("windowsize 1050 700")
    rc("copy file " + png_name + " supersample 3")
    rc("reset")
    rc("close all")

    return


def processFiles(molec,pri_res,sec_res,fn,dn,rootdir):
    print "\\\\\\\\ Writing images to ",dn,'\n'
    if molec == 'gp120': callChimeraGP120(pri_res.replace('_',' '),sec_res.replace('_',' '),fn,dn,rootdir)
    elif molec == 'gp41': callChimeraGP41(pri_res.replace('_',' '),sec_res.replace('_',' '),fn,dn,rootdir)
    sys.exit()

def main(argv):
    molec=pri_res=sec_res=fn=dn=rootdir = ''
    #Read command line arguments
    print "\\\\\\\\\\\ WELCOME to chimeraImageGenerator \\\\\\\\\\\\n"
    try:
        opts, args = getopt.getopt(argv,"m:p:s:f:d:r:",["molec,pri_res,sec_res,fn,dn,rootdir"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-m": 
            molec = arg
        elif opt == "-p":
            pri_res = arg
        elif opt == "-s":
            sec_res = arg
        elif opt == "-f":
            fn = arg
        elif opt == "-d":
            dn = arg
        elif opt == "-r":
            rootdir = arg
    if molec !="" and pri_res!="" and sec_res !="" and fn!="" and dn !="" and rootdir!="":
        processFiles(molec,pri_res,sec_res,fn,dn,rootdir)
    else:
            print "There was an error, not right args\n"
            sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])