#!/usr/local/bin/ python2.7
# encoding: utf-8
"""
hcv_comm_query.py

Created by Mark Evans on 2011-10-03.
Copyright (c) 2011 __Monogram Biosciences__. All rights reserved.

Revised 10.06.2011   Modified Oracle query to add accession and visit_no
Revised 10.13.2011   Added checking for duplicate accession and mutation counting
Revised 02.06.2012   Fixed bugs relating to hcv_mutations.total_count updates and deleting revised records
"""

import sys
import os
os.environ["LD_LIBRARY_PATH"]="/usr/lib/oracle/11.2/client64/lib"
os.environ["ORACLE_HOME"]="/usr/lib/oracle/11.2/client64"
os.environ["TNS_ADMIN"]="/usr/lib/oracle/11.2/client64"
import cx_Oracle as cxo
import psycopg2 as pg
from time import localtime, strftime

# HCV_RAMS FOR BOC, TVR
NS3_RAMS = {'V36A':'','V36I':'','V36M':'','V36G':'','V36C':'','V36L':'','Q41R':'','F43C':'','F43S':'','T54A':'','T54S':'','T54G':'','T54C':'','V55A':'','V55I':'','V107I':'','I132V':'',
            'R155G':'','R155I':'','R155K':'','R155M':'','R155Q':'','R155T':'','R155P':'','R155LL':'','R155S':'','A156S':'','A156T':'','A156V':'','A156G':'','A156F':'','A156N':'',
            'A156I':'','V158I':'','D168N':'','D168Y':'','I170A':'','I170F':'','I170T':'','I170L':'','M175L':''}

NS4A_RAMS = {}


def connectPG(LOG):
   # Define connection string
  # con = "host='gamera.virologic.com' dbname='mgrm' user='mark' password='mark'"
   con = "host='localhost' dbname='mgrm' user='mark' password='mark'"
   t = strftime("%d %b %Y %H:%M:%S",localtime())
   try:
      pg_con = pg.connect(con)       # get a connection or raise exception
      pg_cursor = pg_con.cursor()    # connection returns a cursor object which is used to perform queries
      LOG.write(t+"\tSuccessful connection made to postgres\n")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(t+"\tPostgres Databse connection failed!\n ->%s" % (exceptionValue))
      sys.exit(t+"Postgres Database connection failed!\n ->%s" % (exceptionValue))
   return pg_cursor, pg_con



def connectORA(LOG):
   t = strftime("%d %b %Y %H:%M:%S",localtime())  # Get current time
   try:
      ip = '10.196.1.140'
      port = 1521
      SID = 'vrl1'
      dsn_tns = cxo.makedsn(ip,port,SID)
      con = cxo.connect('medusa_ro','xfinity1',dsn_tns)
      #con = cxo.connect('medusa_ro/xfinity1@failover_db.virologic.com/vrl1') # get a connection or raise exception
      cxo_cursor = con.cursor()                # connection returns a cursor object which is used to perform queries
      LOG.write(t+"\tSuccessful connection made to Oracle\n")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(t+"\tOracle Databse connection failed!\n ->%s" % (exceptionValue))
      sys.exit(t+"Oracle Database connection failed!\n ->%s" % (exceptionValue))
   return cxo_cursor



def prepareMedusaSQL(last_update):
   sql = """SELECT distinct 
   ur.ALIQUOT_ID,
   ur.STATUS,
   ur.project,
   to_char(ur.COLLECTED_DATE,'yyyy-MM-dd') as collected_date,
   to_char(ur.REP_DATE,'yyyy-mm-dd') as rep_date,
   ur.PATIENT_DOB DOB,
   ur.PATIENT_SEX SEX,
   ur.BATCH_ID,
   ur.PATIENT_LNAME||', '||ur.PATIENT_FNAME PName, 
   ur.patient_id PID,
   ur.REF_PHYS_LNAME||', '||ur.REF_PHYS_FNAME DOC_Name,
   pc.state client_state, 
   pc.first_name client_name,
   gt.sub_type,
   NS3.NS3_MUTs, 
   NS4.NS4A_MUTs,
   (CASE 
      when TVR3.TVR_RESULT=0 then 'Sensitive' 
      when TVR3.TVR_RESULT=1 then 'Resistant' 
      when TVR3.TVR_RESULT=2 then 'RP' END) TVR_RES,
   TVR3.TVR_NS3_RAMS, 
   TVR4.TVR_NS4A_RAMS,
   (CASE 
      when BOC3.BOC_RESULT=0 then 'Sensitive' 
      when BOC3.BOC_RESULT=1 then 'Resistant' 
      when BOC3.BOC_RESULT=2 then 'RP' END) BOC_RES,
   BOC3.BOC_NS3_RAMS, 
   BOC4.BOC_NS4A_RAMS,
   ur.VRL_ACCESSION,
   ur.VISIT_NO
   
   FROM 
   MEDUSA.SAMPLES s, 
   MEDUSA.GT_HEADERS gt,
   MEDUSA.UR_REPORT_HEADERS_VIEW ur, 
   MEDUSA.PATIENT_CONTACTS pc, 
   MEDUSA.PATIENT_HEADERS pt, 
   (select sampleid NS3_id, REGIONNAME NS3_r, MUTATIONLIST NS3_MUTs  
      from MEDUSA.SAMPLE_REGIONS
      where REGIONNAME='NS3') NS3,
   (select sampleid NS4_id, REGIONNAME NS4_r, MUTATIONLIST NS4A_MUTs  
      from MEDUSA.SAMPLE_REGIONS
      where upper(REGIONNAME)='NS4A') NS4,
   (select gt_header_id boc3_id, RESULT BOC_RESULT, COMMENTS BOC_CMNT,RESISTANCE_MUTATIONS BOC_NS3_RAMs, region_name boc3_region 
      from MEDUSA.GT_DRUG_RESULTS
      where drug_abrv_name='BOC') BOC3,
   (select gt_header_id boc4_id, RESISTANCE_MUTATIONS BOC_NS4A_RAMs, region_name boc4_region 
      from MEDUSA.GT_DRUG_RESULTS
      where drug_abrv_name='BOC') BOC4,
   (select gt_header_id TVR3_id, RESULT TVR_RESULT, COMMENTS TVR_CMNT,RESISTANCE_MUTATIONS TVR_NS3_RAMs, region_name tvr3_region 
      from MEDUSA.GT_DRUG_RESULTS
      where drug_abrv_name='TVR') TVR3,
   (select gt_header_id TVR4_id, RESISTANCE_MUTATIONS TVR_NS4A_RAMs, region_name tvr4_region 
      from MEDUSA.GT_DRUG_RESULTS
      where drug_abrv_name='TVR') TVR4
      
   WHERE 
   --ur.project='00073' and
   to_char(ur.REP_DATE,'yyyy-mm-dd') > '%s' and
   s.id = NS3.NS3_id and 
   s.id = NS4.NS4_id and 
   s.aliquotid = gt.aliquot_id and 
   ur.aliquot_id = gt.aliquot_id and 
   gt.REPORT_LOCKS = 1 and 
   gt.id = BOC3.BOC3_ID and 
   gt.id = BOC4.BOC4_ID and 
   gt.id = TVR3.TVR3_ID and 
   gt.id = TVR4.TVR4_ID and 
   NS3.NS3_r = TVR3.tvr3_region and 
   NS3.NS3_r = BOC3.boc3_region and 
   NS4.NS4_r = TVR4.tvr4_region and 
   NS4.NS4_r = BOC4.boc4_region and 
   pc.patient_header_id = pt.id and 
   pc.contact_type = 'CLIENT' and 
   pt.report_id = gt.report_id""" % (last_update)
   
   return sql



def checkSingleMut(mut,gene,pgc,LOG,pg_con):
   ram_list = {'ns3':NS3_RAMS,'ns4a':NS4A_RAMS}
   ram=''
   pos=''
   pos = mut[1:-1]
   # Search hcv_mutations table to see if this mutation has been seen before
   try:
      sql = "select mut_id,total_count from wh.hcv_mutations where mutation='%s' and gene='%s'" % (mut,gene.upper())
      pgc.execute(sql)
#      print pgc.statusmessage
#      print pgc.query
      result = pgc.fetchone()
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not execute postgres query!\n"+sql+"\n ->%s" % (exceptionValue))
      sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute postgres query!\n ->%s" % (exceptionValue))
   
   # Mutation is new, so insert it into the table
   if str(result) == 'None':
      if ram_list[gene].has_key(mut): ram='yes'
      else: ram='no'
      try:
         sql2 = "insert into wh.hcv_mutations(mut_key,mutation,gene,total_count,ram,pos) values('%s','%s','%s',%d,'%s',%d)" % (gene.upper()+':'+mut,mut,gene.upper(),1,ram,int(pos))
         pgc.execute(sql2)
         pg_con.commit()

      except:
         exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
         LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not execute postgres query!\n"+sql2+"\n ->%s" % (exceptionValue))
         sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute postgres query!\n ->%s" % (exceptionValue))
   else:
      # Mutation already exists, so increment the counter and update the table
      count = int(result[1]) + 1
      try:
         sql2 = "update wh.hcv_mutations set total_count=%d where mut_id=%d" % (count,result[0])
         pgc.execute(sql2)
         pg_con.commit()

      except:
         exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
         LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not execute postgres query!\n"+sql2+"\n ->%s" % (exceptionValue))
         sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute postgres query!\n ->%s" % (exceptionValue))
   return
   
   
def countMuts(mut_sum,gene,pgc,LOG,pg_con):
   if str(mut_sum) !='None':
      muts = mut_sum.split(', ')
      for m in muts:
         # Check for single mutants first
         if m.find('/') == -1: checkSingleMut(m,gene,pgc,LOG,pg_con)
         
         # Have multiple mutants
         else:
            muts = m.split('/')   # ['I50A','C','H']
            mut_count = len(muts)
            
            # check first mutation in list
            if muts[0][0] != muts[0][-1]: checkSingleMut(muts[0],gene,pgc,LOG,pg_con)
            
            # Generate the rest of mutations in list and check them.
            for i in range(1,len(muts)):
                  m2 = muts[0][:-1] + muts[i] # automatically create I50C, I50H etc
                  checkSingleMut(m2,gene,pgc,LOG,pg_con)
   return




def main():
   # Create or append to update log file
   ######################################
  # LOG = open('/home/mevans/logs/hcv_comm_query_update.log','a')
   LOG = open('hcv_update.log','w')
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" beginning postgres hcv_sum update\n")
   
   
   # Create database connections
   ##############################
   pgc, pg_con = connectPG(LOG)
   cxc = connectORA(LOG)
   
   
   # Get the last update date from Postgres
   #########################################
   try:
      pgc.execute("select max(rep_date) as rep_date from wh.hcv_sum")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not execute postgres query!\n ->%s" % (exceptionValue))
      sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute postgres query!\n ->%s" % (exceptionValue))
      
   last_update = pgc.fetchone()[0]
   if last_update == None: last_update = '2010-01-01'
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" Last update to hcv_sum was "+str(last_update)+"\n")
   
   
   # Prepare Medusa SQL to retrieve incremental update
   ####################################################
   medusa_sql = prepareMedusaSQL(last_update)
   try:
      cxc.execute(medusa_sql)
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not execute Oracle query!\n ->%s" % (exceptionValue))
      sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute Oracle query!\n ->%s" % (exceptionValue))
      
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " retrieved update from medusa\n")
   
   
   # Loop through each new record and insert it into Postgres table
   #################################################################
   try:
      result = cxc.fetchall()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " retrieved "+str(len(result))+" records from medusa\n")
      sql =''
      if len(result) != 0:
         for aliquot_id,status,project,collected_date,rep_date,dob,sex,batch_id,pt_name,pid,doc_name,client_state,client_name,subtype,ns3_muts,ns4a_muts,tvr_res,tvr_ns3_rams,tvr_ns4a_rams,boc_res,boc_ns3_rams,boc_ns4a_rams,accession,visit in result:
            # Check if this accession has been entered before. 
            sql = "select uid,accession_id from wh.hcv_sum where accession_id='%s'" % (accession)
            try:
               pgc.execute(sql)
               r = pgc.fetchone()
            except:
               exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
               LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\tCould not execute Postgres sql!\n"+sql+"\n ->%s" % (exceptionValue))
               sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute Postgres sql\n ->%s" % (exceptionValue))
            # check the results. If yes, delete previous record and proceed to insert this one
            if str(r) != 'None':
               sql2 = "delete from wh.hcv_sum where uid=%s" % (r[0])
               try:
                  pgc.execute(sql2)
               except:
                  exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
                  LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\tCould not execute Postgres deletion!\n"+sql+"\n ->%s" % (exceptionValue))
                  sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute Postgres deletion\n ->%s" % (exceptionValue))
            else:
               # Since this record did not exist before, it is okay to proceed to process and count the mutations that are present
               countMuts(ns3_muts,'ns3',pgc,LOG,pg_con)
               countMuts(ns4a_muts,'ns4a',pgc,LOG,pg_con)
               
            # Now we can proceed to insert the actual record itself, since the preprocessing validation is finished
            sql = "insert into wh.hcv_sum (aliquot_id,status,project,collected_date,rep_date,dob,sex,batch_id,pt_name,pid,doc_name,client_state,client_name,subtype,ns3_sum,ns4a_sum,tvr_res,tvr_ns3_rams,tvr_ns4a_rams,boc_res,boc_ns3_rams,boc_ns4a_rams,accession_id,visit_no) "
            sql += "values('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" % (aliquot_id,status,project,collected_date,rep_date,dob,sex,batch_id,pt_name.replace("'",""),pid,doc_name.replace("'",""),client_state,client_name.replace("'",""),subtype,ns3_muts,ns4a_muts,tvr_res,tvr_ns3_rams,tvr_ns4a_rams,boc_res,boc_ns3_rams,boc_ns4a_rams,accession,visit)
            LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " Attempting insert: "+sql+"\n")
            try:
               pgc.execute(sql)
            except:
               exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
               LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\tCould not execute Postgres insert!\n"+sql+"\n ->%s" % (exceptionValue))
               sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute Postgres insert!\n ->%s" % (exceptionValue))
      else:
         LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " No records have been added to Oracle, therefore update cancelled\n")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not read Oracle results!\n ->%s" % (exceptionValue))
      sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not read Oracle results!\n ->%s\n->%s\n->%s\n" % (exceptionType,exceptionValue,exceptionTraceback))

   pg_con.commit()
   
   # Clean up and exit
   ####################
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" completed updating postgres\n")
   cxc.close()
   pgc.close()
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" Update finished successfully\n++++++++++++++++++++++++++++++\n")
   LOG.close()

if __name__ == '__main__':
   main()

