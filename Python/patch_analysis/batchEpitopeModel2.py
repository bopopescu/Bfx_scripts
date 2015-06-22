
"""
batchEpitopeModel.py
Loop through all epitope definitions and models them in 3D

Created by Mark Evans on 2012-06-22.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,getopt,string
from decimal import *
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages


def callChimera(pri_epitope,sec_epitope,fn,dirname):
    #trimer = '/Users/talon34/Dropbox/Bali/pgt121-123/models/gp120CPHmodel_trimer.pdb'
    trimer = '/Users/Bali/Documents/Monogram/Neut/Patch_analysis/chimera_integration/models/gp120CPHmodel_trimer.pdb'
    #monomer = '/Users/talon34/Dropbox/Bali/pgt121-123/models/gp120_CPHModel.pdb'
    monomer = '/Users/Bali/Documents/Monogram/Neut/Patch_analysis/chimera_integration/models/gp120_CPHModel.pdb'

    print "pri_epitope = ",pri_epitope
    print "sec_epitope = ",sec_epitope

    # Create 1st trimer image
    replyobj.status("Processing trimer") # show what file we're working on
    rc("open " + trimer)
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
    rc("copy file " + png_name + " supersample 3")
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
    rc("copy file " + png_name + " supersample 3")
    rc("turn 0,1,0 180")
    png_name = fn+".back-side.png"
    rc("copy file " + png_name + " supersample 3")

    rc("close all")

    return


def resultHTMLFile(ab,fn,dirname,pri,sec):
    f = open(dirname+fn+".idx.html",'w')
    f.write("<html><body><h3>Epitope Images for "+ab+" on gp120</h3>")
    code = """<p>
                <b>Primary Epitope: """+pri+"""<br>
                 Secondary Epitope: """+sec+"""</b></p>
                    <table>
                    <tr><th>Trimer top</th><th>Trimer Side</th><th>Monomer Front</th><th>Monomer Back</th></tr>
                    <tr>
                        <td><a href='"""+fn+""".trimer-top.png'><img src='"""+fn+""".trimer-top.png' width='200'></a></td>
                        <td><a href='"""+fn+""".trimer-side.png'><img src='"""+fn+""".trimer-side.png' width='200'></a></td>
                        <td><a href='"""+fn+""".front-side.png'><img src='"""+fn+""".front-side.png' width='200'></a></td>
                        <td><a href='"""+fn+""".back-side.png'><img src='"""+fn+""".back-side.png' width='200'></a></td>
                    </tr></table>
                </body></html>"""
    f.write(code)
    return


##############
# Main
###################
datadir='/Users/Bali/Documents/Monogram/Neut/Patch_analysis/chimera_integration/images/'
replyobj.status("Beginning script...\n ")

#files = os.listdir(datadir)     # get all of the result file names
EPITOPES = open("epitopes.txt","r")
HTMLIDX = open(datadir+"index.html",'w')

HTMLIDX.write("<html><body>")
HTMLIDX.write("<h3>Epitope Analysis Results</h3><P>\n")
HTMLIDX.write("<table width='100%' border=1>\n<col width='20'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<tr><th>Antibody</th><th>Result File</th><th>Trimer top view</th><th>Primary Epitope</th><th>Secondary Epitope</th></tr>\n")
# Epitope file should be in the following format, tab-delimited
# antibody_name  result_file    pri_epitope   sec_epitope   pri_epitope_num   sec_epitope_num

data = {}

for line in EPITOPES:
    vals = line.rstrip().split('\t')
    data[vals[0]] = vals

mabs = data.keys()
mabs.sort()

for ab in mabs:
    fn = ab +"_"+data[ab][1][:-4]
    print "######################\n# Processing ",ab,"\n########################"
    code =''
    if data[ab][2] != 'No Result' and data[ab][3] != 'No Result':
        os.mkdir(datadir+fn)
        dn = datadir+fn+"/"
        code ="<tr><td nowrap='nowrap'><a href='"+fn+"/"+fn+".idx.html'>"+ab+"</a></td><td>"+data[ab][1]+"</td><td><img src='"+fn+"/"+fn+".trimer-top.png' width='200'></td><td>"+data[ab][2]+"</td><td>"+data[ab][3]+"</td></tr>\n"
        callChimera(data[ab][4],data[ab][5],fn,dn)
        resultHTMLFile(ab,fn,dn,data[ab][2],data[ab][3])
    else:
        code ="<tr><td nowrap='nowrap'>"+ab+"</td><td>"+data[ab][1]+"</td><td></td><td>"+data[ab][2]+"</td><td>"+data[ab][3]+"</td></tr>\n"
    print "+++++++++ Code for ",ab," = ",code
    HTMLIDX.write(code)
HTMLIDX.write("</table></body></html>")
HTMLIDX.close()
        
print "\n\n******** Analysis complete\n"
