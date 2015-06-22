#!/usr/bin/env python
# encoding: utf-8
"""
assignTestNum2.py

Created by Mark Evans on 2011-10-21.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from time import localtime, strftime
import psycopg2 as pg


# Connect to Postgres
#########################
def connectPG(LOG):
   # Define connection string
   #   con = "host='gamera.virologic.com' dbname='mgrm' user='mark' password='mark'"
   con = "host='localhost' dbname='mgrm' user='mark' password='mark'"
   
   try:
      pg_con = pg.connect(con)       # get a connection or raise exception
      pg_cursor = pg_con.cursor()    # connection returns a cursor object which is used to perform queries
      LOG.write(gt()+"\tSuccessful connection made to postgres\n")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(gt()+"\tPostgres Databse connection failed!\n ->%s" % (exceptionValue))
      sys.exit(gt()+"\tPostgres Database connection failed!\n ->%s" % (exceptionValue))
   return pg_cursor, pg_con


# execute Postgres SQL, do not return result
#############################################
def execPG(pgc,pgc_sql,vals,LOG):
   try:
      pgc.execute(pgc_sql,vals)
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(gt()+"\tCould not execute Postgres sql!\n"+pgc_sql+"\n ->%s" % (exceptionValue))
      return "Error"
   return "ok"



# Returns current time
#########################
def gt():
   return strftime("%d %b %Y %H:%M:%S",localtime())



def main():
   filename = raw_input("What is name of patient_test_export file? ")
   f1 = open(filename,'r')
   
   buffer_size = 0
   LOG = open("test_count_log.txt","w", buffer_size)
   LOG.write("***********************\n"+gt()+" Beginning update process...\n")

   # Create database connections
   ##############################
   pgc, pg_con = connectPG(LOG)


   LOG.write(gt()+"\tstarted counting tests....\n")
   
   c = 1
   cid = ''
   sql = "update wh.patient set test_num=%s where patient_id=%s"

   for line in f1:
      a = line.split('\t')
      a[-1].rstrip()  # if did not export with trailing null column, use this
 #     del a[-1]  # if exported with last column as null, use this
      LOG.write(gt()+"\t\tupdating "+a[1]+"\n")
      if a[1] != cid:
         cid = a[1]
         c = 1
      elif a[1]==cid:
         c = c + 1
      execPG(pgc,sql,[c,a[0]],LOG)

   pg_con.commit()
   f1.close()
   LOG.write(gt()+"\tfinished \n")



if __name__ == '__main__':
   main()

