"""
calcFinalEpitope.py
Process all epitopes indentified and find highest ranking common residues for a "final" epitope call

Created by Mark Evans on 2012-07-27.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""


import sys,os,getopt,string,subprocess
from decimal import *
from multiprocessing import Pool
import multiprocessing
from operator import itemgetter
import logging, logging.handlers


###########################################################################################
# Set up global logging object
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

# This log handler writes to console
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# This log handler writes to file
fh = logging.FileHandler(filename='final_image_generator.log', mode='a')
fh.setLevel(logging.DEBUG)

# this log handler emails anything that is error or worse
#sh = logging.handlers.SMTPHandler('autofax.virologic.com','logger@gamera.virologic.com',['Evansm8@labcorp.com','talon.sensei@gmail.com'],'HIVWH Update Error Log')
#sh.setLevel(logging.ERROR)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)
fh.setFormatter(formatter)
#sh.setFormatter(formatter)

#add ch to logger
LOG.addHandler(ch)
LOG.addHandler(fh)
#LOG.addHandler(sh)

###########################################################################################


def resultHTMLFile(ab,epitope,path,molec):
    code = """<tr>
                <th bgcolor="#ffffcc" valign="middle">"""+ab+"""<br></th>
                <td align="center" valign="middle"><a href='"""+path+ab+'-'+molec+"""-final.trimer-top.png'><img src='"""+path+ab+'-'+molec+"""-final.trimer-top.png' width='317'></a><br></td>
                <td align="center" valign="middle"><a href='"""+path+ab+'-'+molec+"""-final.trimer-side.png'><img src='"""+path+ab+'-'+molec+"""-final.trimer-side.png' width='317'></a></td>
                <td align="center" valign="middle"><a href='"""+path+ab+'-'+molec+"""-final.front-side.png'><img src='"""+path+ab+'-'+molec+"""-final.front-side.png' width='317'></a></td>
                <td align="center" valign="middle"><a href='"""+path+ab+'-'+molec+"""-final.back-side.png'><img src='"""+path+ab+'-'+molec+"""-final.back-side.png' width='317'></a></td>
                <td align="left" valign="middle" width="10%">"""+epitope+"""<br></td>
            </tr>"""
    return code


###########################################################
# ProcessImages                                           #
# Calls Chimera and executes chimeraImageGenerator script #
###########################################################
def processImages(chimera_args):
    molecule = chimera_args[0].replace(' ','_')
    epitope = chimera_args[1].replace(' ','_')
    fn = chimera_args[2].replace(' ','_')
    dirname = chimera_args[3]
    rootdir = chimera_args[4]
    script_path = os.getcwd()
    script_path += "/chimeraFinalImageGenerator.py"
    LOG.info("##################")
    LOG.info("SCRIPT PATH: "+script_path)
    LOG.info("VARS: "+molecule+", "+str(chimera_args))
    try:
        # This line works in terminal
        #process = subprocess.Popen(['chimera','--nogui','--nosilent','--script=chimeraImageGenerator.py -m gp120 -p 152,155,157,160,162 -s 132,154,158,159,165,166,167,168,169 -f gp120_pgt-143_result77 -d /home/bfx_analysis/model_test/images/10A15dpx/gp120_pgt-143_result77/ -r /home/bfx_analysis/model_test'])
        LOG.info(multiprocessing.current_process().name+" calling Chimera for "+fn)
        process = subprocess.Popen(['chimera','--nogui','--nosilent','--script='+script_path+' -m '+molecule+' -e '+epitope+' -f '+fn+' -d '+dirname+' -r '+rootdir])
        process.wait()
    except:
        LOG.info("!!!!!!!!! Chimera subprocess failed for "+fn+" in "+multiprocessing.current_process().name+" !!!!!!!!!")
    LOG.info(multiprocessing.current_process().name+" completed Chimera subprocess")
    return


############################
# Pool process initializer #
############################
def start_process():
    LOG.info('Starting '+multiprocessing.current_process().name)
    return


def processFiles(chimera_epitope_input):
    f1 = open(chimera_epitope_input,'r')
    
    rootdir=os.getcwd()
    datadir = rootdir+'/images/final_results/'
    if not os.path.exists(datadir): os.mkdir(datadir)
    epitope_data = {}
    chimera_args = []
    ekeys = []
    mabs = []
    final_results=""
    getcontext().prec=2
    code = """<html><head>
                <meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">
                <title>Antibody Summary Results</title>
                </head><body>
                <table border="0" cellpadding="0" cellspacing="0" width="100%">
                <tbody>
                <tr>
                    <th bgcolor="#ccffff" valign="top">Antibody<br> </th>
                    <th bgcolor="#ccffff" valign="top">Trimer Top View<br></th>
                    <th bgcolor="#ccffff" valign="top">Trimer Side View<br></th>
                    <th bgcolor="#ccffff" valign="top">Monomer Front View<br></th>
                    <th bgcolor="#ccffff" valign="top">Monomer Side View<br></th>
                    <th align="left" bgcolor="#ccffff" valign="middle">Predicted Epitope<br></th>
                </tr>"""
    # Load master epitope result file into dictionary
    for line in f1:
        vals = line.rstrip().split('\t')
        epitope_data[vals[3]] = vals  # was 7
    f1.close()

    # get dictionary keys
    ekeys = epitope_data.keys()

    # get all the mabs that are in the file
    for k in ekeys:
        if epitope_data[k][2] not in mabs: mabs.append(epitope_data[k][2])  # was 6

    # Begin processing epitope data per antibody
    for ab in mabs:
        abset = {'residues':{},'count':{}}	# hold all the residue information for this ab regardless of analysis combo
        c = 0
        molec = ''
        for k in ekeys:
            if epitope_data[k][2] == ab: # was 6
                c += 1		                # keep count how many analyses were done for this ab
                molec = epitope_data[k][2].split('_')[0]    # was 5
                pe = epitope_data[k][4]		# primary epitope residues   was 8
                se = epitope_data[k][5]		# secondary epitope residues   was 9
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

        final_results += ab+'\t'+molec+'\t'+final_epitope+'\n'
        if len(fq_vals) != 0:
            path = datadir+ab+'/'
            os.mkdir(path)
            chimera_args.append((molec,chimera_epitope,ab+'-'+molec,path,rootdir))
            code += resultHTMLFile(ab,final_epitope,ab+'/',molec)

    
    code += """</tbody></table></body></html>"""
    f2 = open(molec+'_final_epitope_results.txt','w')
    f2.write(final_results)
    f2.close()
    f3 = open(datadir+'/'+molec+'_index.html','w')
    f3.write(code)
    f3.close()

    # Multiprocessing Section
    #########################################
    multiprocessing.log_to_stderr(logging.INFO)  # set logging to info level rather than DEBUG
    pool_size = multiprocessing.cpu_count()*2  # number 4 is arbitrary, example used 2. Trying to maximize throughput without using all RAM and crashing
    pool = Pool(processes = pool_size, initializer = start_process, maxtasksperchild = 4,)
    pool.map(processImages,chimera_args)
    pool.close()    # no more tasks
    pool.join()     # wrap up current tasks

    LOG.info("\n\n******** Image generation complete\n")
    sys.exit()






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