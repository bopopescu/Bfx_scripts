
"""
batchEpitopeModel.py
Loop through all epitope definitions and models them in 3D

Created by Mark Evans on 2012-06-22.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

2012-08-06  This version functions as a wrapper for chimeraImageGenerator.py so that image generation will be multiprocessed
"""
#
#
# This version crashed gamera because it sucked up all RAM, still not fixed. Going to try switching to Pool in v.6 instead
#
#
import sys,os,getopt,string,subprocess, random
from decimal import *
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter
import logging, logging.handlers

#####################################################################################
# glbl                                                                              #
# Creates a class containing global variables that can be accessed across processes # 
#####################################################################################
class glbl:
   jobs =[]
   model_input={}
   molec=''

###########################################################################################
# Set up global logging object
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

# This log handler writes to console
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# This log handler writes to file
fh = logging.FileHandler(filename='chimera_generator.log', mode='a')
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
#glbl.LOG.addHandler(sh)

###########################################################################################





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


def processImages(nid):
    pri_res = sec_res = fn = dn = rootdir =''
    while 1:
        try: mid = nid.get(False)                                                             # Get next query id from Queue generator, 
        except: 
            print " CPU not in use"
            break
        pri_res = glbl.model_input[glbl.molec][mid][0]
        sec_res = glbl.model_input[glbl.molec][mid][1]
        fn = glbl.model_input[glbl.molec][mid][2]
        dn = glbl.model_input[glbl.molec][mid][3]
        rootdir = glbl.model_input[glbl.molec][mid][4]
        script_path = os.getcwd()
        script_path += "/chimeraImageGenerator.py"
        print "##################\n"
        print "SCRIPT PATH: ",script_path,"\n"
        print "VARS: ",glbl.model_input[glbl.molec][mid],"\n"
        try:
            # This line works in terminal
            #process = subprocess.Popen(['chimera','--nogui','--nosilent','--script=chimeraImageGenerator.py -m gp120 -p 152,155,157,160,162 -s 132,154,158,159,165,166,167,168,169 -f gp120_pgt-143_result77 -d /home/bfx_analysis/model_test/images/10A15dpx/gp120_pgt-143_result77/ -r /home/bfx_analysis/model_test'])

            process = subprocess.Popen(['chimera','--nogui','--nosilent','--script='+script_path+' -m '+glbl.molec+' -p '+pri_res.replace(' ','_')+' -s '+sec_res.replace(' ','_')+' -f '+fn.replace(' ','_')+' -d '+dn+' -r '+rootdir]) #' -m '+glbl.molec+' -p '+pri_res+' -s '+sec_res+' -f '+fn+' -d '+dn+' -r '+rootdir+
            process.wait()
        except:
            print "Chimera subprocess failed"
        nid.task_done()
        return


##############
# Main
###################
def processFiles(input_file):
    # input file should be chimera_epitope_input.txt for which ever model you are building
    print "Beginning script...\n "
    glbl.molec = ''
    data = {}
    glbl.model_input={'gp120':[],'gp41':[]}

    # Epitope file should be in the following format, tab-delimited
    # patch_width  dpx  antibody_name  result_file    pri_epitope   sec_epitope   pri_epitope_num   sec_epitope_num
    EPITOPES = open(input_file,"r")
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
        glbl.molec = vals[2].split('_')[0]

    rootdir=os.getcwd()
    datadir = rootdir+'/images/'

    dirs = data.keys()
    dirs.sort()
    for d in dirs:
        if not os.path.exists(datadir+d): os.mkdir(datadir+d)
        path = datadir+d+'/'
            
        HTMLIDX = open(path+"index.html",'w')
        HTMLIDX.write("<html><body>")
        HTMLIDX.write("<h3>"+glbl.molec+" Epitope Analysis Results</h3><P>\n")
        HTMLIDX.write("<table width='100%' border=1>\n<col width='20'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<tr><th>Antibody</th><th>Result File</th><th>Trimer top view</th><th>Primary Epitope</th><th>Secondary Epitope</th></tr>\n")

        mabs = data[d].keys()
        mabs.sort()
        for ab in mabs:
            fn = ab +"_"+data[d][ab][3][:-4]
 #           print "######################\n# Processing ",ab,"\n########################"
            code =''
            if not (data[d][ab][4] == 'No Result' and data[d][ab][5] == 'No Result'):
                if not os.path.exists(path+fn):os.mkdir(path+fn)
                dn = path+fn+"/"
                pri_res = sec_res = ''
                if data[d][ab][4] == 'No Result': pri_res = 'No Result'
                else: pri_res = data[d][ab][6]
                if data[d][ab][5] == 'No Result': sec_res = 'No Result'
                else: sec_res = data[d][ab][7]
                code ="<tr><td nowrap='nowrap'><a href='"+fn+"/"+fn+".idx.html'>"+ab+"</a></td><td>"+data[d][ab][3]+"</td><td><img src='"+fn+"/"+fn+".trimer-top.png' width='200'></td><td>"+data[d][ab][4]+"</td><td>"+data[d][ab][5]+"</td></tr>\n"
                if glbl.molec == 'gp120': glbl.model_input['gp120'].append((pri_res,sec_res,fn,dn,rootdir))
                elif glbl.molec == 'gp41': glbl.model_input['gp41'].append((pri_res,sec_res,fn,dn,rootdir))
                resultHTMLFile(ab,fn,dn,data[d][ab][4],data[d][ab][5])
            else:
                code ="<tr><td nowrap='nowrap'>"+ab+"</td><td>"+data[d][ab][3]+"</td><td></td><td>"+data[d][ab][4]+"</td><td>"+data[d][ab][5]+"</td></tr>\n"
 #           print "+++++++++ Code for ",ab," = ",code
            HTMLIDX.write(code)
        HTMLIDX.write("</table></body></html>")
        HTMLIDX.close()
    print "HTML generation finished....beginning image generation...\n"

    for x in glbl.model_input['gp120']:
        print x[3]
    print "length of model_input = ",len(glbl.model_input[glbl.molec])

    # Multiprocessing Section
    #########################################
    multiprocessing.log_to_stderr(logging.INFO)  # set logging to info level rather than DEBUG
    manager = Manager()                                          # creates shared memory manager object
    nextid = Queue()                                             # Create Queue object to serve as shared id generator across processes
    for qid in range(len(glbl.model_input[glbl.molec])): nextid.put(qid)                             # Load the ids to be tested into the Queue
    LOG.info('There are '+str(len(glbl.model_input[glbl.molec]))+' items in queue')
    #for x in range(0,multiprocessing.cpu_count()):               # Create one process per logical CPU
    for x in range(0,len(glbl.model_input[glbl.molec])):
        LOG.info('Starting process '+str(x))
        p = Process(target=processImages, args=(nextid,)) # Assign process to processCBR function, passing in the Queue and shared dictionary
        glbl.jobs.append(p)                                      # Add the process to a list of running processes
        p.start()                                                # Start process running
    for j in glbl.jobs:
        j.join()                                              # For each process, join them back to main, blocking on each one until finished
    

    print "\n\n******** Image generation complete\n"


def main(argv):
    input_file=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"f:",["input_file"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-f": 
            input_file = arg
            processFiles(input_file)
        else:
            print "There was an error, not right args\n"
            sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])