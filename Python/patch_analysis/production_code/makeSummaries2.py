"""
makeSummaries.py
Creates alternative summary view, linking previously generated images into table view for easy review

Created by Mark Evans on 2012-08-02
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

2012-08-09  Changed input file to use chimera_epitope_input, added creation of main index HTML file
"""


import sys,os,getopt,string
from decimal import *


def resultHTMLFile(abset,molec,ab,view):
    vws = {'Trimer Top View':'trimer-top','Trimer Side View':'trimer-side','Monomer Front View':'front-side','Monomer Back View':'back-side'}
    side = vws[view]
    my_abset = {}
    dpx = abset.keys()
    dpx.sort()
    area = abset[dpx[0]].keys()
    area.sort()

    for d in dpx:
        my_abset[d] = {}
        for a in area:
            if abset[d][a] != 'No result':
                r = abset[d][a]                                                                         # REPLACE env <> cph dending on what you are making
                my_abset[d][a]= "<img width='314 height='217' src='http://gamera.virologic.com/neut/results/cph/"+molec+"/"+a+"A"+d.replace('.','')+"dpx/"+ab+"_"+r[:-4]+"/"+ab+"_"+r[:-4]+"."+side+".png' />"
            else: my_abset[d][a] = 'No result'
    code = """<html>
    <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <title>"""+ab+"""</title>
    </head>
    <body>
    <p></p>
    <table width="100%" cellspacing="1" border="0" align="center">
      <tbody>
        <tr>
          <th rowspan="2" colspan="2" style="background-color: #b4e1e1; text-align: center; vertical-align: middle;">"""+ab+""" Result<br />"""+view+"""</th>
          <th rowspan="1" colspan="6" style="text-align: center; vertical-align: middle; background-color: #ededbd;">Patch Diameter</th>
        </tr>
        <tr>
          <th style="text-align: center; vertical-align: middle; background-color: #ededbd;">6 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; background-color: #ededbd;">8 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; background-color: #ededbd;">10 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; background-color: #ededbd;">12 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; background-color: #ededbd;">15 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; background-color: #ededbd;">20 Angstroms</th>
          <td></td>
          <td></td>
        </tr>
        <tr>
          <th rowspan="4" colspan="1" style="text-align: center; vertical-align: middle; background-color: #90f090;">Atomic Depth</th>
          <th style="text-align: center; vertical-align: middle; margin-left: -29px; background-color: #90f090; ">0 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['0']['6']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['0']['8']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['0']['10']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['0']['12']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['0']['15']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['0']['20']+"""</th>
        </tr>
        <tr>
          <th style="text-align: center; vertical-align: middle; background-color: #90f090; ">1.5 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['1.5']['6']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['1.5']['8']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['1.5']['10']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['1.5']['12']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['1.5']['15']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227px; width: 320px;">"""+my_abset['1.5']['20']+"""</th>
          <td style=""></td>
        </tr>
        <tr>
          <th style="text-align: center; vertical-align: middle; background-color: #90f090; ">2.5 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['2.5']['6']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['2.5']['8']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['2.5']['10']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['2.5']['12']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['2.5']['15']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['2.5']['20']+"""</th>
          <td style=""></td>
        </tr>
        <tr>
          <th style="text-align: center; vertical-align: middle; background-color: #90f090; ">3.0 Angstroms</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['3.0']['6']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['3.0']['8']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['3.0']['10']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['3.0']['12']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['3.0']['15']+"""</th>
          <th style="text-align: center; vertical-align: middle; height: 227;">"""+my_abset['3.0']['20']+"""</th>
          <td style=""></td>
        </tr>
      </tbody>
    </table>
    <p><br />
    </p>
    </body>
    </html>"""
    del(my_abset)
    return code


def createIndexHtml(mabs,molec):
    f2 = open(molec+'_index.html','w')
    mabs.sort()
    code = """<html><head>
                  <meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">
                  <title>Antibody Summaries</title></head><body>
                  <h2>Antibody """+molec+""" Summary Displays</h2>

                  <table border="1" cellpadding="2" cellspacing="2" width="100%">
                  <tbody>
                    <tr>
                        <th valign="top">Antibody<br></th>
                        <th valign="top">Trimer Top View<br></th>
                        <th valign="top">Trimer Side View<br></th>
                        <th valign="top">Monomer Front View<br></th>
                        <th valign="top">Monomer Back View<br></th>
                    </tr>"""
    for m in mabs:
        ab = m.split('_')[1]
        code += """ <tr>
                        <td align="center" valign="top">"""+ab+"""<br></td>
                        <td align="center" valign="top"><a href='"""+m+'_'+molec+"""-1.html'>View 1</a></td>
                        <td align="center" valign="top"><a href='"""+m+'_'+molec+"""-2.html'>View 2</a></td>
                        <td align="center" valign="top"><a href='"""+m+'_'+molec+"""-3.html'>View 3</a></td>
                        <td align="center" valign="top"><a href='"""+m+'_'+molec+"""-4.html'>View 4</a></td>
                    </tr>"""
    code += """</tbody></table><br><br></body></html>"""
    f2.write(code)
    f2.close()
    return



def processFiles(chimera_epitope_input):
    f1 = open(chimera_epitope_input,'r')
    
    rootdir=os.getcwd()

    epitope_data = {}
    ekeys = []
    mabs = []

    # Load chimera epitope result file into dictionary
    for line in f1:
        vals = line.rstrip().split('\t')
        epitope_data[vals[3]] = vals  # Use result file name as index, e.g. result853.txt
    f1.close()

    # get dictionary keys
    ekeys = epitope_data.keys()

    # get all the mabs that are in the file
    for k in ekeys:
        if epitope_data[k][2] not in mabs: mabs.append(epitope_data[k][2])

    # Begin processing epitope data per antibody
    mabs.sort()
    views = ['Trimer Top View','Trimer Side View','Monomer Front View','Monomer Back View']
    for ab in mabs:
        # Hold all the residue information for this ab regardless of analysis combo
        # Need to change this matrix if analysis combos change
        abset = {'0':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'},
                 '1.5':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'},
                 '2.5':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'},
                 '3.0':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'}} 
        molec = code = ""
        for k in ekeys:
            if epitope_data[k][2] == ab:
                molec = epitope_data[k][2].split('_')[0]
                area = epitope_data[k][0].split('_')[0]
                dpx = epitope_data[k][1][3:]
                if epitope_data[k][4] != 'No Result' or epitope_data[k][5] != 'No Result':
                    abset[dpx][area] = epitope_data[k][3]			
        c = 1

        for v in views:
            code = resultHTMLFile(abset,molec,ab,v)
            f3 = open(ab+'_'+molec+'-'+str(c)+'.html','w')
            c += 1
            f3.write(code)
            f3.close()

    createIndexHtml(mabs,molec)





def main(argv):
    chimera_epitope_input=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:",["file:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-f":
            chimera_epitope_input = arg    
            processFiles(chimera_epitope_input)
        else:
            print "There was an error, not right args\n"
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])