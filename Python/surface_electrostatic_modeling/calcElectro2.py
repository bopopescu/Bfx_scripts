"""
calcElectro.py
Process all epitopes indentified and find highest ranking common residues for a "final" epitope call

Revised 2014-08-06 
Created by Mark Evans on 2014-07-27.
Copyright (c) 2014 __XOMA_(US)_LLC__. All rights reserved.

"""


import sys,os,getopt,string,subprocess
from decimal import *
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter
from time import localtime,strftime
import logging, logging.handlers


###########################################################################################
# Set up global logging object
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

# This log handler writes to console
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# This log handler writes to file
fh = logging.FileHandler(filename='electro_calculations.log', mode='a')
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

#############
# Usage     #
#############
def usage():
    code =  "\n\n#############\nLoop through all pdb files in specified directory and calculate electrostatics and generate model images colored by electrostatic property\n"
    code += "Usage: $ python calcElectro.py -d data_directory_path \n\n"
    code += "-d [data_dir] path to data directory containing pdb files to process\n"
    code += "-h help\n\n#############\n\n"
    print code
    return


########################
# Returns current time #
########################
def gt():
    return strftime("%d %b %Y %H:%M:%S",localtime())


#####################################################################################
# glbl                                                                              #
# Creates a class containing global variables that can be accessed across processes # 
#####################################################################################
class glbl:
   jobs =[]
   model_input = {}


#################################################
# doPDB2PQR                                     #
# Calls PDB2PQR and generates APBS input file   #
#################################################
def doPDB2PQR(pdb_id):
    pid =''
    while 1:
        try: pid = pdb_id.get(False)
        except:
            break

        basename = glbl.model_input[pid][0]
        LOG.info("################## Doing PDB2PQR")
        try:
            # Run PQR2PDB first
            LOG.info(multiprocessing.current_process().name+" calling PDB2PQR for "+basename)
            process = subprocess.call(["pdb2pqr","--ff","AMBER","--apbs-input","--nodebump", pid, basename+'.pqr'])
        except:
            LOG.info(multiprocessing.current_process().name+" !!!!!!!!! PDB2PQR unexpected error for "++str(sys.exc_info()[0])+"\n"+str(sys.exc_info()[1])+"\n"+str(sys.exc_info()[2]))
        LOG.info(multiprocessing.current_process().name+" completed PDB2PQR subprocess")
    return


######################################################
# doAPBS                                             #
# Calls APBS and generates electrostatics .dx file   #
######################################################
def doAPBS(pdb_id):
    pid =''
    while 1:
        try: pid = pdb_id.get(False)
        except:
            break

        basename = glbl.model_input[pid][0]
        LOG.info("################## Doing APBS")
        try:
            # Run APBS
            LOG.info(multiprocessing.current_process().name+" calling APBS for " + basename)
            process = subprocess.call(["apbs", basename + '.in'])
        except:
            LOG.info(multiprocessing.current_process().name+" !!!!!!!!! APBS unexpected error "+str(sys.exc_info()[0])+"\n"+str(sys.exc_info()[1])+"\n"+str(sys.exc_info()[2]))
        LOG.info(multiprocessing.current_process().name+" completed APBS subprocess")
    return


#########################################################################
# doChimera                                                             #
# Calls Chimera and loads .pdb and .dx files then creates 3 .png images #
#########################################################################
def doChimera(pdb_id):
    pid =''
    while 1:
        try: pid = pdb_id.get(False)
        except:
            break
      
        basename = glbl.model_input[pid][0]
        rootdir = glbl.model_input[pid][1]
        imagedir = glbl.model_input[pid][2]
        script_path = os.getcwd()
        script_path += "/doElectro.py"
        LOG.info("################## Doing Chimera")
        LOG.info("SCRIPT PATH: "+script_path)
        LOG.info("Imagedir PATH: "+imagedir)
        LOG.info("Rootdir: "+rootdir)
        try:
            # Run Chimera and generate images
            LOG.info(multiprocessing.current_process().name+" calling Chimera for "+basename)
            path =  script_path+" -f "+ pid +" -e "+basename+".pqr.dx"+ " -d "+imagedir+" -r "+rootdir
            LOG.info("Chimera args path: "+path)
            process = subprocess.call(['chimera','--bgopacity', '--multisample', '--script='+path])
        except:
            LOG.info(multiprocessing.current_process().name+" !!!!!!!!! Chimera unexpected error "+str(sys.exc_info()[0])+"\n"+str(sys.exc_info()[1])+"\n"+str(sys.exc_info()[2]))
        LOG.info(multiprocessing.current_process().name+" completed Chimera subprocess")
    return


################################################
# doMultiprocessing                            #
# Starts multiprocess Queue to execute apps in #
################################################
def doMultiprocess(app):
    applist = {'pdb2pqr':doPDB2PQR,
               'apbs':doAPBS,
               'chimera':doChimera}
    multiprocessing.log_to_stderr(logging.INFO)               # set logging to info level rather than DEBUG
    pdb_ids = glbl.model_input.keys()
    manager = Manager()                                       # creates shared memory manager object
    nextPDBid = Queue()                                       # Create Queue object to serve as shared id generator across processes
    for pid in pdb_ids: nextPDBid.put(pid)                    # Load the ids to be tested into the Queue
    for x in range(0,multiprocessing.cpu_count()):            # Create one process per logical CPU
        p = Process(target=applist[app], args=(nextPDBid,))   # Assign process to app function, passing in the Queue
        glbl.jobs.append(p)                                   # Add the process to a list of running processes
        p.start()                                             # Start process running
    for j in glbl.jobs:
        j.join()  
    return     


##############################################################
# processFiles                                               #
# Initiate processing of all pdb files in the directory,     #
# to calculate electrostatics and generate 3D images of them #
##############################################################
def processFiles(pdb_files_dir):
    LOG.info("Beginning script...")
    start_time = gt()
    all_files = os.listdir(pdb_files_dir)
    rootdir=os.getcwd()
    datadir = rootdir+'/images/'
    if not os.path.exists(datadir): os.mkdir(datadir)


    # Begin processing pdb files
    for fn in all_files:
        if fn[-4:] == '.pdb':
            basename = fn[:-4]
            glbl.model_input[fn] = (basename,rootdir,datadir)

    doMultiprocess('pdb2pqr')
    doMultiprocess('apbs')
    doMultiprocess('chimera')
    finish_time = gt()
    LOG.info("\n\n******** Image generation complete\n")
    LOG.info("Start time: "+start_time)
    LOG.info("Finish Time: "+finish_time)
    sys.exit()


def main(argv):
    pdb_files_dir=''

    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"d:h:",["dir:help:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-d":
            pdb_files_dir = arg    
            processFiles(pdb_files_dir)
        elif opt == "-h":
            usage()
            sys.exit()
        else:
            print "There was an error, not right args\n"
            usage()
            sys.exit()
       

if __name__ == '__main__':
    main(sys.argv[1:])