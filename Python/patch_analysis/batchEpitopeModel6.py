
"""
batchEpitopeModel.py
Loop through all epitope definitions and models them in 3D

Created by Mark Evans on 2012-06-22.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

2012-08-06  This version functions as a wrapper for chimeraImageGenerator.py so that image generation will be multiprocessed.
"""

import sys,os,getopt,string,subprocess, random
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
#LOG.addHandler(sh)

###########################################################################################


########################################################
# resultHTMLFile                                       #
# Creates HTML index file for each antibody result set #
########################################################
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


###########################################################
# ProcessImages                                           #
# Calls Chimera and executes chimeraImageGenerator script #
###########################################################
def processImages(chimera_args):
    pri_res = chimera_args[0].replace(' ','_')
    sec_res = chimera_args[1].replace(' ','_')
    fn = chimera_args[2].replace(' ','_')
    dn = chimera_args[3]
    rootdir = chimera_args[4]
    molecule = chimera_args[5]
    script_path = os.getcwd()
    script_path += "/chimeraImageGenerator.py"
    LOG.info("##################")
    LOG.info("SCRIPT PATH: "+script_path)
    LOG.info("VARS: "+molecule+", "+str(chimera_args))
    try:
        # This line works in terminal
        #process = subprocess.Popen(['chimera','--nogui','--nosilent','--script=chimeraImageGenerator.py -m gp120 -p 152,155,157,160,162 -s 132,154,158,159,165,166,167,168,169 -f gp120_pgt-143_result77 -d /home/bfx_analysis/model_test/images/10A15dpx/gp120_pgt-143_result77/ -r /home/bfx_analysis/model_test'])
        LOG.info(multiprocessing.current_process().name+" calling Chimera for "+fn)
        process = subprocess.Popen(['chimera','--nogui','--nosilent','--script='+script_path+' -m '+molecule+' -p '+pri_res+' -s '+sec_res+' -f '+fn+' -d '+dn+' -r '+rootdir])
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


########################################################################################
# Main                                                                                 #
# input file should be chimera_epitope_input.txt for which ever model you are building #
########################################################################################
def processFiles(input_file):
    LOG.info("Beginning script...")
    data = {}
    rootdir=os.getcwd()
    datadir = rootdir+'/images/'
    model_input={'gp120':[],'gp41':[]}
    molec = ''

    # Load epitope data.  Epitope file should be in the following format, tab-delimited
    # patch_width  dpx  antibody_name  result_file    pri_epitope   sec_epitope   pri_epitope_num   sec_epitope_num
    ################################################################################################################
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
        molec = vals[2].split('_')[0]


    # Create HTML index files and identify combinations to generate images for
    ###########################################################################
    dirs = data.keys()
    dirs.sort()
    for d in dirs:
        if not os.path.exists(datadir+d): os.mkdir(datadir+d)   # make directory if doesn't exist
        path = datadir+d+'/'                                    # define path to write index file to
        
        # This HTML page is the main index.html for all the abs in this diam / dpx combo
        HTMLIDX = open(path+"index.html",'w')
        HTMLIDX.write("<html><body>")
        HTMLIDX.write("<h3>"+molec+" Epitope Analysis Results</h3><P>\n")
        HTMLIDX.write("<table width='100%' border=1>\n<col width='20'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<col width='100'>\n<tr><th>Antibody</th><th>Result File</th><th>Trimer top view</th><th>Primary Epitope</th><th>Secondary Epitope</th></tr>\n")

        # Get list of antibodies to process
        mabs = data[d].keys()
        mabs.sort()
        for ab in mabs:
            fn = ab +"_"+data[d][ab][3][:-4]                                            # create filename which contains antibody name in it
            code =''
            if not (data[d][ab][4] == 'No Result' and data[d][ab][5] == 'No Result'):   # Make sure there will be an image generated for this diam / dpx combo
                if not os.path.exists(path+fn):os.mkdir(path+fn)
                dn = path+fn+"/"                                                        # create directory path
                pri_res = sec_res = ''
                if data[d][ab][4] == 'No Result': pri_res = 'No Result'                 # store primary epitope residues
                else: pri_res = data[d][ab][6]
                if data[d][ab][5] == 'No Result': sec_res = 'No Result'                 # store secondary epitope residues
                else: sec_res = data[d][ab][7]
                code ="<tr><td nowrap='nowrap'><a href='"+fn+"/"+fn+".idx.html'>"+ab+"</a></td><td>"+data[d][ab][3]+"</td><td><img src='"+fn+"/"+fn+".trimer-top.png' width='200'></td><td>"+data[d][ab][4]+"</td><td>"+data[d][ab][5]+"</td></tr>\n"
                
                # store ChimeraImageGenerator script arguments in a dictionary array for the multiprocessing step
                ################################################################################################## 
                if molec == 'gp120': model_input['gp120'].append((pri_res,sec_res,fn,dn,rootdir,molec))
                elif molec == 'gp41': model_input['gp41'].append((pri_res,sec_res,fn,dn,rootdir,molec))
                
                # create individual index.html files for each antibody folder
                ##############################################################
                resultHTMLFile(ab,fn,dn,data[d][ab][4],data[d][ab][5])
            else:
                code ="<tr><td nowrap='nowrap'>"+ab+"</td><td>"+data[d][ab][3]+"</td><td></td><td>"+data[d][ab][4]+"</td><td>"+data[d][ab][5]+"</td></tr>\n"
 
            HTMLIDX.write(code)
        HTMLIDX.write("</table></body></html>")
        HTMLIDX.close()

    LOG.info("HTML generation finished....beginning image generation...")
    LOG.info("Chimera should be called "+str(len(model_input[molec]))+" times")

    # Multiprocessing Section
    #########################################
    multiprocessing.log_to_stderr(logging.INFO)  # set logging to info level rather than DEBUG
    pool_size = multiprocessing.cpu_count() * 4  # number 4 is arbitrary, example used 2. Trying to maximize throughput without using all RAM and crashing
    pool = Pool(processes = pool_size, initializer = start_process, maxtasksperchild = 4,)
    pool.map(processImages,model_input[molec])
    pool.close()    # no more tasks
    pool.join()     # wrap up current tasks

    LOG.info("\n\n******** Image generation complete\n")
    sys.exit()


################
# Main process #
################
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