#!/usr/bin/env python
# encoding: utf-8
"""
calcmut.py

Created by Mark Evans on 2011-04-29.
Revised 2011-05-06
Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""

import sys
import os
from decimal import *
from time import localtime, strftime
os.environ["LD_LIBRARY_PATH"]="/usr/lib/oracle/11.2/client64/lib"
os.environ["ORACLE_HOME"]="/usr/lib/oracle/11.2/client64"
os.environ["TNS_ADMIN"]="/usr/lib/oracle/11.2/client64"
import cx_Oracle as cxo
import psycopg2 as pg

#aa = {A':'',C':'',D':'',E':'',F':'',G':'',H':'',I':'',K':'',L':'',M':'',N':'',P':'',Q':'',R':'',S':'',T':'',V':'',W':'',Y':''}

# HCV_RAMS FOR BOC, TVR
NS3_RAMS = {'V36A':'','V36I':'','V36M':'','V36G':'','V36C':'','V36L':'','Q41R':'','F43C':'','F43S':'','T54A':'','T54S':'','T54G':'','T54C':'','V55A':'','V55I':'','V107I':'','I132V':'',
            'R155G':'','R155I':'','R155K':'','R155M':'','R155Q':'','R155T':'','R155P':'','R155LL':'','R155S':'','A156S':'','A156T':'','A156V':'','A156G':'','A156F':'','A156N':'',
            'A156I':'','V158I':'','D168N':'','D168Y':'','I170A':'','I170F':'','I170T':'','I170L':'','M175L':''}

NS4A_RAMS = {}

def connectPG(LOG):
   # Define connection string
   con = "host='gamera.virologic.com' dbname='mgrm' user='mark' password='mark'"
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




def main():
   # Create or append to update log file
   ######################################
   LOG = open('hcv_count_mut.log','a')
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" beginning postgres retrieve data\n")
   
   
   # Create database connections
   ##############################
   pgc, pg_con = connectPG(LOG)
   
   
   # Get the last update date from Postgres
   #########################################
   try:
      pgc.execute("select ns3_sum,ns4a_sum from wh.hcv_sum")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not execute postgres query!\n ->%s" % (exceptionValue))
      sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute postgres query!\n ->%s" % (exceptionValue))
      
   
   # Loop through each new record and insert it into Postgres table
   #################################################################
   # seen = {'A123T':['NS3:A123T','A123T','NS3',1,'yes',123]}
   seen = {}
   ram=''
   pos=''
   try:
      result = pgc.fetchall()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " retrieved "+str(len(result))+" records from gamera\n")
      if len(result) != 0:
         for ns3_sum,ns4a_sum in result:
            ns3 = ns3_sum.split(', ')
            ns4a = ns4a_sum.split(', ')
            for val in ns3:
               # Check for single mutants first
               if val.find('/') == -1: 

                  pos = val[1:-1]

                  if NS3_RAMS.has_key(val): ram='yes'
                  else: ram='no'
                  if seen.has_key(val):
                     seen[val][3] = seen[val][3] + 1
                  else: seen[val] = ['NS3:'+val,val,'NS3',1,ram,int(pos)]

               # Have multiple mutants
               else:
                  val2 = val.split('/')   # ['I50A','C','H']
                  mut_count = len(val2)
                  if val2[0][0] != val2[0][-1]:
                     # check first mutation in list
                     
                     pos = val2[0][1:-1]
                     
                     if NS3_RAMS.has_key(val2[0]): ram='yes'
                     else: ram='no'
                     if seen.has_key(val2[0]):
                        seen[val2[0]][3] = seen[val2[0]][3] + 1
                     else: seen[val2[0]] = ['NS3:'+val2[0],val2[0],'NS3',1,ram,int(pos)]
                  # Generate the rest of mutations in list and check them.
                  for i in range(1,len(val2)):
                     m = val2[0][:-1] + val2[i] # automatically create I50C, I50H etc
                     
                     pos = m[1:-1]
                     
                     if NS3_RAMS.has_key(m): ram='yes'
                     else: ram='no'
                     if seen.has_key(m):
                        seen[m][3] = seen[m][3] + 1
                     else: seen[m] = ['NS3:'+m,m,'NS3',1,ram,int(pos)]

            for val in ns4a:
               # Check for single mutants first
               if val != 'None':
                  if val.find('/') == -1: 
                  
                     pos = val[1:-1]
                  
                     if NS4A_RAMS.has_key(val): ram='yes'
                     else: ram='no'
                     if seen.has_key(val):
                        seen[val][3] = seen[val][3] + 1
                     else: 
                        try:
                           seen[val] = ['NS4A:'+val,val,'NS4A',1,ram,int(pos)]
                        except:
                           print val,"FAILED at ",pos,"\n"

               # Have multiple mutants
                  else:
                     val2 = val.split('/')   # ['I50A','C','H']
                     mut_count = len(val2)
                     if val2[0][0] != val2[0][-1]:
                        # check first mutation in list
                     
                        pos = val2[0][1:-1]
                     
                        if NS4A_RAMS.has_key(val2[0]): ram='yes'
                        else: ram='no'
                        if seen.has_key(val2[0]):
                           seen[val2[0]][3] = seen[val2[0]][3] + 1
                        else: seen[val2[0]] = ['NS4A:'+val2[0],val2[0],'NS4A',1,ram,int(pos)]
                        # Generate the rest of mutations in list and check them.
                     for i in range(1,len(val2)):
                        m = val2[0][:-1] + val2[i] # automatically create I50C, I50H etc
                     
                        pos = m[1:-1]
                     
                        if NS4A_RAMS.has_key(m): ram='yes'
                        else: ram='no'
                        if seen.has_key(m):
                           seen[m][3] = seen[m][3] + 1
                        else: seen[m] = ['NS4A:'+m,m,'NS4A',1,ram,int(pos)]
               
      else:
         LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " No records have been added to Oracle, therefore update cancelled\n")
      
      for r in seen:
         sql=''
         try:
            sql = """insert into hcv_mutations(mut_key,mutation,gene,total_count,ram,pos) values ('%s','%s','%s',%d,'%s',%d)""" % (seen[r][0],seen[r][1],seen[r][2],int(seen[r][3]),seen[r][4],int(seen[r][5]))
         except:
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            print "FAILED",r,":",seen[r],"\n",sys.exc_info()
         LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+ " Attempting insert: "+sql+"\n")
         try:
            pgc.execute(sql)
         except:
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\tCould not execute Postgres insert!\n ->%s" % (exceptionValue))
            sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not execute Postgres insert!\n ->%s" % (exceptionValue))
            
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+"\Could not read Postgres results!\n ->%s\n->%s\n->%s" % (exceptionType,exceptionValue,exceptionTraceback))
      sys.exit(strftime("%d %b %Y %H:%M:%S",localtime())+"Could not read Postgres results!\n ->%s\n->%s\n->%s" % (exceptionType,exceptionValue,exceptionTraceback))

   pg_con.commit()
   
   # Clean up and exit
   ####################
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" completed updating postgres\n")
   pgc.close()
   LOG.write(strftime("%d %b %Y %H:%M:%S",localtime())+" Update finished successfully\n++++++++++++++++++++++++++++++\n")
   LOG.close()
   
   
if __name__ == '__main__':
   main()

