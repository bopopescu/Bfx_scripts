
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


def callChimeraGP120(pri_epitope,sec_epitope,fn,dirname,rootdir):
    base_dir = os.getcwd()
    trimer = rootdir+'/gp120CPHmodel_trimer.pdb'
    monomer = rootdir+'/gp120_CPHModel.pdb'
    print "trimer path= ",trimer
    print "monomer path= ",monomer
    print "pri_epitope = ",pri_epitope
    print "sec_epitope = ",sec_epitope

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
replyobj.status("Beginning script...\n ")
molec = ''
data = {}

# Epitope file should be in the following format, tab-delimited
# patch_width  dpx  antibody_name  result_file    pri_epitope   sec_epitope   pri_epitope_num   sec_epitope_num
EPITOPES = open("chimera_epitope_input.txt","r")
for line in EPITOPES:
    vals = line.rstrip().split('\t')
    a = vals[0].split('_')[0]
    d = vals[1].replace('dpx','').replace('.','')
    key = a+'A'+d+'dpx'
    if data.has_key(key):
        data[key][vals[2]] = vals
    else:
        data[key] = {}
        data[key][vals[2]] = vals
    molec = vals[2].split('_')[0]

rootdir=os.getcwd()
datadir = rootdir+'/images/'





dirs = data.keys()
dirs.sort()
for d in dirs:
    if not os.path.exists(datadir+d): os.mkdir(datadir+d)
    path = datadir+d+'/'
    
    HTMLIDX = open(path+"index.html",'w')
    HTMLIDX.write("<html><body>")
    HTMLIDX.write("<h3>"+molec+" Epitope Analysis Results</h3><P>\n")
    HTMLIDX.write("<table width='100%' border=1>\n<col width='20'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<tr><th>Antibody</th><th>Result File</th><th>Trimer top view</th><th>Primary Epitope</th><th>Secondary Epitope</th></tr>\n")

    mabs = data[d].keys()
    mabs.sort()
    for ab in mabs:
        fn = ab +"_"+data[d][ab][3][:-4]
        print "######################\n# Processing ",ab,"\n########################"
        code =''
        if not (data[d][ab][4] == 'No Result' and data[d][ab][5] == 'No Result'):
            os.mkdir(path+fn)
            dn = path+fn+"/"
            pri_res = sec_res = ''
            if data[d][ab][4] == 'No Result': pri_res = 'No Result'
            else: pri_res = data[d][ab][6]
            if data[d][ab][5] == 'No Result': sec_res = 'No Result'
            else: sec_res = data[d][ab][7]
            code ="<tr><td nowrap='nowrap'><a href='"+fn+"/"+fn+".idx.html'>"+ab+"</a></td><td>"+data[d][ab][3]+"</td><td><img src='"+fn+"/"+fn+".trimer-top.png' width='200'></td><td>"+data[d][ab][4]+"</td><td>"+data[d][ab][5]+"</td></tr>\n"
            if molec == 'gp120': callChimeraGP120(pri_res,sec_res,fn,dn,rootdir)
            elif molec == 'gp41': callChimeraGP41(pri_res,sec_res,fn,dn,rootdir)
            resultHTMLFile(ab,fn,dn,data[d][ab][4],data[d][ab][5])
        else:
            code ="<tr><td nowrap='nowrap'>"+ab+"</td><td>"+data[d][ab][3]+"</td><td></td><td>"+data[d][ab][4]+"</td><td>"+data[d][ab][5]+"</td></tr>\n"
        print "+++++++++ Code for ",ab," = ",code
        HTMLIDX.write(code)
    HTMLIDX.write("</table></body></html>")
    HTMLIDX.close()
        
print "\n\n******** Analysis complete\n"
