"""
calcFinalEpitope.py
Process all epitopes indentified and find highest ranking common residues for a "final" epitope call

Created by Mark Evans on 2012-07-27.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

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


def processFiles(main_epitope_file):
    f1 = open(main_epitope_file,'r')
    
    rootdir=os.getcwd()

    epitope_data = {}
    ekeys = []
    mabs = []

    # Load main epitope result file into dictionary
    for line in f1:
        vals = line.rstrip().split('\t')
        if vals[0] != 'Dir1': epitope_data[vals[7]] = vals
    f1.close()

    # get dictionary keys
    ekeys = epitope_data.keys()

    # get all the mabs that are in the file
    for k in ekeys:
        if epitope_data[k][6] not in mabs: mabs.append(epitope_data[k][6])

    # Begin processing epitope data per antibody
    mabs.sort()
    views = ['Trimer Top View','Trimer Side View','Monomer Front View','Monomer Back View']
    for ab in mabs:
        code = ""
        abset = {'0':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'},
                 '1.5':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'},
                 '2.5':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'},
                 '3.0':{'6':'No result','8':'No result','10':'No result','12':'No result','15':'No result','20':'No result'}} # hold all the residue information for this ab regardless of analysis combo
        molec = ''
        for k in ekeys:
            if epitope_data[k][6] == ab:
                if ab == 'gp41_pg16': print epitope_data[k]
                molec = epitope_data[k][5]
                area = epitope_data[k][3].split('_')[0]
                dpx = epitope_data[k][4][3:]
                if epitope_data[k][8] != 'No Result' or epitope_data[k][9] != 'No Result':
                    if ab =='gp41_pg16': print "************* Found one ",epitope_data[k][7]
                    abset[dpx][area] = epitope_data[k][7]			
        c = 1

        for v in views:
            code = resultHTMLFile(abset,molec,ab,v)
            f3 = open(ab+'_'+molec+'-'+str(c)+'.html','w')
            c += 1
            f3.write(code)
            f3.close()





def main(argv):
    main_epitope_file=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:",["file:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-f":
            main_epitope_file = arg    
            processFiles(main_epitope_file)
        else:
            print "There was an error, not right args\n"
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])