#!/usr/bin/env python
# encoding: utf-8
"""
load_mysql_data.py

Created by Mark Evans on 2011-06-30.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""
import os,sys
import psycopg2 as pg
import pymysql as ms
from time import localtime, strftime

# Connect to Postgres
#########################
def connectPG(LOG):
   # Define connection string
	con = "host='gamera.virologic.com' dbname='mgrm' user='mark' password='mark'"
	#con = "host='localhost' dbname='mgrm' user='mark' password='mark'"
   
	try:
		pg_con = pg.connect(con)       # get a connection or raise exception
		pg_cursor = pg_con.cursor()    # connection returns a cursor object which is used to perform queries
		LOG.write(gt()+"\tSuccessful connection made to postgres\n")
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tPostgres Database connection failed!\n ->%s" % (exceptionValue))
		sys.exit(gt()+"\tPostgres Database connection failed!\n ->%s" % (exceptionValue))
	return pg_cursor, pg_con

# Connect to MySQL Databases
################################
def connectMySQL(LOG):
	try:
		con = ms.connect(host='kong.virologic.com',user='seq_pipe',passwd='seqpipe',db='GEN_PIPE_DB')  # mgrm samples
		kong = con.cursor()
		LOG.write(gt()+"\tSuccessful connection made to Kong\n")
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tKong Database connection failed!\n ->%s" % (exceptionValue))
		sys.exit(gt()+"\tKong Database connection failed!\n ->%s" % (exceptionValue))
	try:
		con2 = ms.connect(host='kongdb.virologic.com',user='root',passwd='olympics',db='GENORULES')	# Labcorp samples
		kong_db = con2.cursor()
		LOG.write(gt()+"\tSuccessful connection made to KongDB\n")
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tKongDB Database connection failed!\n ->%s" % (exceptionValue))
		sys.exit(gt()+"\tKongDB Database connection failed!\n ->%s" % (exceptionValue))
	return kong,kong_db


# Returns current time
#########################
def gt():
	return strftime("%d %b %Y %H:%M:%S",localtime())


# execute Postgres SQL, do not return result
#############################################
def execPG(pgc,pgc_sql,vals,LOG):
	try:
		pgc.execute(pgc_sql,vals)
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tCould not execute Postgres sql!\nsql="+pgc_sql+"\nvals="+str(vals)+"\n ->%s" % (exceptionValue))
		return "Error"
	return "ok"


# Execute any MySQL query and return a result
###############################################
def getMySQL(msc,msc_sql,vals,LOG):
	try:
		msc.execute(msc_sql,vals)
		msc_result = msc.fetchall()
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tCould not execute MySQL query!\n"+msc_sql+"\n ->%s\n" % (exceptionValue))
		return "Error"
	LOG.write(gt()+ "\t\tRetrieved update from MySQL\n")
	return msc_result


def loadTables(pgc,kong,kong_db,LOG):
	LOG.write(gt()+"\tUpdating wh.gt_rule_result_rows....\n")
	# Update MGRM data

	# Get Kong records > last_update
	sql = """select SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,'mgrm' from GT_RULE_RESULT_ROWS where CREATED_DATE < '2011-10-23 01:00:00'"""
	mgrm = kong.execute(sql) #getMySQL(kong,sql,'',LOG)
	if mgrm != "Error":
		# Verify there are records to insert
		if len(mgrm) != 0:
			LOG.write(gt()+"\t\tWriting MGRM update to wh.gt_rule_result_rows table...\n")
			# Insert new records in wh.gt_rule_result_rows
			for SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src in mgrm:
				sql = "insert into wh.gt_rule_result_rows(SAMPLE_ID,DRUG_ABRV,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
				vals = (str(SAMPLE_ID),str(DRUG_ABRV_NAME),str(DRUG_CLASS),str(COMMENTS),str(RAMs),RESULT,RESISTANT_FLG,MUTATION_SCORE,str(CREATED_DATE),RULE_VERSION,ENGINE_VERSION,str(data_src))
				msr = execPG(pgc,sql,vals,LOG)
			LOG.write(gt()+"\twh.gt_rule_result_rows part 1 update was successful..."+str(len(mgrm))+" records added\n")
			
		else:
			LOG.write(gt()+"\twh.gt_rule_result_rows update part 1 was skipped, there are no records to add today\n")
	else:
		LOG.write(gt()+"\tFailed to retrieve update from KONG, wh.gt_rule_result_rows update part 1 was not done\n")
	
	# Update LABC data
	# Get KongDB records > last_update
	sql = """select SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,'labc' from GT_RULE_RESULT_ROWS where CREATED_DATE < '2011-10-23 01:00:00'"""
	labc = kong_db.execute(sql) #getMySQL(kong_db,sql,'',LOG)
	if labc != "Error":
		# Verify there are records to insert
		if len(labc) != 0:
			LOG.write(gt()+"\t\tWriting LABC update to wh.gt_rule_result_rows table...\n")
			# Insert new records in wh.gt_rule_result_rows
			for SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src in labc:
				sql = "insert into wh.gt_rule_result_rows(SAMPLE_ID,DRUG_ABRV,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
				vals = (str(SAMPLE_ID),str(DRUG_ABRV_NAME),str(DRUG_CLASS),str(COMMENTS),str(RAMs),RESULT,RESISTANT_FLG,MUTATION_SCORE,str(CREATED_DATE),RULE_VERSION,ENGINE_VERSION,str(data_src))
				msr = execPG(pgc,sql,vals,LOG)
			LOG.write(gt()+"\twh.gt_rule_result_rows part 2 update was successful..."+str(len(labc))+" records added\n\t\t\twh.gt_rule_result_rows update was successful...\n")
			return
		else:
			LOG.write(gt()+"\twh.gt_rule_result_rows update part 2 was skipped, there are no records to add today\n")
			return
	else:
		LOG.write(gt()+"\tFailed to retrieve update from KONG_DB, wh.gt_rule_result_rows update part 2 was not done\n")
		return
	
	return



def main():
	buffer_size = 0
	LOG = open("table_load_log.txt","w", buffer_size)  # /home/mevans/logs/hivdb_update_log.txt
	LOG.write("***********************\n"+gt()+" Beginning update process...\n")

	# Create database connections
	##############################
	pgc, pg_con = connectPG(LOG)
	kong, kong_db = connectMySQL(LOG)

   # Begin processing tables
	loadTables(pgc,kong,kong_db,LOG)	

  	# Commit all changes to the database
 	pg_con.commit()
	print "Database update complete\n"
	LOG.write(gt()+"\tDatabase update completed\n##############################\n")
	LOG.close()

if __name__ == '__main__':
	main()