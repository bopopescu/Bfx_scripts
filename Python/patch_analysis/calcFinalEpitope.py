"""
calcFinalEpitope.py
Process all epitopes indentified and find highest ranking common residues for a "final" epitope call

Created by Mark Evans on 2012-07-27.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""


import sys,os,getopt,string
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

def resultHTMLFile(ab,epitope,path,molec):
	code = """<tr>
				<th bgcolor="#ffffcc" valign="top">"""+ab+"""<br></th>
				<td valign="top">"""+epitope+"""<br></td>
				<td align="center" valign="middle"><a href='"""+path+'/'+ab+'-'+molec+"""-final.trimer-top.png'><img src='"""+path+'/'+ab+'-'+molec+"""-final.trimer-top.png' width='317'></a><br></td>
				<td align="center" valign="middle"><a href='"""+path+'/'+ab+'-'+molec+"""-final.trimer-side.png'><img src='"""+path+'/'+ab+'-'+molec+"""-final.trimer-top.png' width='317'></a></td>
				<td align="center" valign="middle"><a href='"""+path+'/'+ab+'-'+molec+"""-final.monomer-front.png'><img src='"""+path+'/'+ab+'-'+molec+"""-final.monomer-back.png' width='317'></a></td>
				<td align="center" valign="middle"><a href='"""+path+'/'+ab+'-'+molec+"""-final.monomer-front.png'><img src='"""+path+'/'+ab+'-'+molec+"""-final.monomer-back.png' width='317'></a></td>
			</tr>"""
	return code


def processFiles(master_epitope_file):
	f1 = open(master_epitope_file,'r')
	f2 = open('final_epitope_results.txt','w')
	rootdir=os.getcwd()
	datadir = rootdir+'/images/'
	epitope_data = {}
	ekeys = []
	mabs = []
	getcontext().prec=2
	code = """<html><head>
			  <meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">
			  <title>Antibody Summary Results</title>
			  </head><body>
			  <table border="0" cellpadding="0" cellspacing="0" width="100%">
			  <tbody>
				<tr>
					<th bgcolor="#ccffff" valign="top">Antibody<br> </th>
					<th bgcolor="#ccffff" valign="top">Predicted Epitope<br></th>
					<th bgcolor="#ccffff" valign="top">Trimer Top View<br></th>
					<th bgcolor="#ccffff" valign="top">Trimer Side View<br></th>
					<th bgcolor="#ccffff" valign="top">Monomer Front View<br></th>
					<th bgcolor="#ccffff" valign="top">Monomer Side View<br></th>
				</tr>"""
	# Load master epitope result file into dictionary
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
	for ab in mabs:
		abset = {'residues':{},'count':{}}	# hold all the residue information for this ab regardless of analysis combo
		c = 0
		molec = ''
		for k in ekeys:
			if epitope_data[k][6] == ab:
				c += 1		# keep count how many analyses were done for this ab
				molec = epitope_data[k][5]
				pe = epitope_data[k][8]		# primary epitope residues
				se = epitope_data[k][9]		# secondary epitope residues
				p = []						# merge primary and secondary into this
				
				if pe != 'No Result':	p.extend(pe.replace('"','').split(','))
				if se != 'No Result':	p.extend(se.replace('"','').split(','))
				
				if len(p) != 0:
					for x in p:
						xd = int(x[1:])	# residue position integer
						if not abset['residues'].has_key(xd): abset['residues'][xd] = x
						if not abset['count'].has_key(xd): abset['count'][xd] = 1
						else: abset['count'][xd] += 1
		freq = {}
		for pos in abset['count']:
			fq = Decimal(abset['count'][pos])/Decimal(c)*Decimal(100)
			if freq.has_key(fq): freq[fq].append(pos)
			else: freq[fq] = [pos]
		fq_vals = freq.keys()
		if len(fq_vals) != 0: N = max(fq_vals)
		else: print "WARNING: fq_vals is empty. abset[count]=",abset['count'], 'k=',k,' ab=',ab
		norm_freq = {}
		for fq in freq:
			fn = Decimal(fq)/Decimal(N)*Decimal(100)
			if fn >= 50: norm_freq[fn] = fq
		
		final_epitope = chimera_epitope = ''
		fe = []
		for x in norm_freq:
			fe.extend(freq[norm_freq[x]])
		fe.sort()
		for x in fe:
	#		print "x=",x," abset[residues]=",abset['residues']
			final_epitope = final_epitope + abset['residues'][x]+','
			chimera_epitope = chimera_epitope + str(x)+','
		final_epitope = final_epitope[:-1]
		chimera_epitope = chimera_epitope[:-1]

		f2.write(ab+'\t'+molec+'\t'+final_epitope+'\n')
		if molec =='gp120' and len(fq_vals) != 0:
			path = datadir+ab+'/'
			os.mkdir(path)
			callChimeraGP120(chimera_epitope,ab+'-'+molec,path,rootdir)
			code += resultHTMLFile(ab,final_epitope,path,molec)
		elif molec == 'gp41' and len(fq_vals) != 0:
			path = datadir+ab+'/'
			os.mkdir(path)
			callChimeraGP41(chimera_epitope,ab+'-'+molec,path,rootdir)
			code += resultHTMLFile(ab,final_epitope,path,molec)

	f2.close()
	code += """</tbody></table></body></html>"""
	f3 = open(datadir+'/index.html','w')
	f3.write(code)
	f3.close()





def main(argv):
    master_epitope_file=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:",["file:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-f":
            master_epitope_file = arg    
            processFiles(master_epitope_file)
        else:
            print "There was an error, not right args\n"
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])