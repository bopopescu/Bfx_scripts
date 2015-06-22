#!/usr/bin/env python
# encoding: utf-8
"""
updateHIVdb.py

Created by Mark Evans on 2011-06-30.
Copyright (c) 2011 __Monogram Biosciences__. All rights reserved.

Updated 01.04.2012  Added function to update etag_sum table. Modified connectORA to add etag login, 
					updated hivdb_update_data.py with new sql
Updated 11.02.2011	Added sendMail function to email log file upon failure, other minor fixes
"""

import sys, os, re
os.environ["LD_LIBRARY_PATH"]="/usr/lib/oracle/11.2/client64/lib"
os.environ["ORACLE_HOME"]="/usr/lib/oracle/11.2/client64"
os.environ["TNS_ADMIN"]="/usr/lib/oracle/11.2/client64"
import cx_Oracle as cxo
import psycopg2 as pg
import pymysql as ms
from time import localtime, strftime
from decimal import *
from hivdb_update_data import *

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
		#sys.exit(gt()+"\tPostgres Database connection failed!\n ->%s" % (exceptionValue))
		return 'failed','',''
	return 'ok', pg_cursor, pg_con


# Connect to Oracle
#########################
def connectORA(LOG):
   
   try:
      ip = '10.196.1.140'
      port = 1521
      SID = 'vrl1'
      dsn_tns = cxo.makedsn(ip,port,SID)
      #con = cxo.connect('vlps/vlps@failover_db.virologic.com/vrl1') # get a connection or raise exception
      vlps = cxo.connect('vlps','vlps',dsn_tns)
      vlps_cursor = vlps.cursor()                # connection returns a cursor object which is used to perform queries
      LOG.write(gt()+"\tSuccessful connection made to Oracle as vlps\n")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(gt()+"\tvlps Oracle Database connection failed!\n ->%s" % (exceptionValue))
      #sys.exit(gt()+"\tOracle Database connection failed!\n ->%s" % (exceptionValue))
      return 'failed','',''
   try:
   	  etag = cxo.connect('etag','etag',dsn_tns)
   	  etag_cursor = etag.cursor()
   	  LOG.write(gt()+"\tSuccessful connection made to Oracle as etag\n")
   except:
      exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
      LOG.write(gt()+"\tetag Oracle Database connection failed!\n ->%s" % (exceptionValue))
      return 'failed','',''
   return 'ok', vlps_cursor, etag_cursor


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
		#sys.exit(gt()+"\tKong Database connection failed!\n ->%s" % (exceptionValue))
		return 'failed', '',''
	try:
		con2 = ms.connect(host='kongdb.virologic.com',user='root',passwd='olympics',db='GENORULES')	# Labcorp samples
		kong_db = con2.cursor()
		LOG.write(gt()+"\tSuccessful connection made to KongDB\n")
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tKongDB Database connection failed!\n ->%s" % (exceptionValue))
		#sys.exit(gt()+"\tKongDB Database connection failed!\n ->%s" % (exceptionValue))
		return 'failed', '', ''
	return 'ok', kong, kong_db


# Returns current time
#########################
def gt():
	return strftime("%d %b %Y %H:%M:%S",localtime())


# Execute Oracle SQL, do not return result
############################################
def execORA(cxc,cxc_sql,vals,LOG):
	try:
		cxc.execute(cxc_sql,vals)
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tCould not execute Oracle query!\n"+cxc_sql+"\n ->%s\n" % (exceptionValue))
		return "Error"
	LOG.write(gt()+ "\tOracle sql executed successfully\n")
	return "ok"

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


# Execute any Oracle query and return result
##############################################
def getOracle(cxc,cxc_sql,vals,LOG):
	try:
		#cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		cxc.execute(cxc_sql,vals)
		cxc_result = cxc.fetchall()
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tCould not execute Oracle query!\n"+cxc_sql+"\n ->%s\n" % (exceptionValue))
		return "Error"
	LOG.write(gt()+ "\t\tRetrieved update from Oracle\n")
	return cxc_result


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


# Execute any Postgres query and return result
#################################################
def getPG(pgc,pgc_sql,vals,LOG):
	try:
		pgc.execute(pgc_sql,vals)
		pgc_result = pgc.fetchall()
	except:
		exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
		LOG.write(gt()+"\tCould not execute Postgres sql!\n"+pgc_sql+"\n ->%s" % (exceptionValue))
		return "Error"
	
	return pgc_result


# Update wh.Drugs table
############################
def updateDrugs(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.drugs....\n")
	cxr = getOracle(cxc,update_drugs_sql,'',LOG)
	if cxr != "Error":
		sql = "delete from wh.drugs"
		pgr = execPG(pgc,sql,'',LOG)
		if pgr != "Error":
			LOG.write(gt()+"\tWriting update to wh.drugs table...\n")
			for drug_id,drug_abrv,dname,fname,bname,drug_class,mfc,dec,hsc,lc,clc,uc,drug_flg in cxr:
				sql = "insert into wh.drugs values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
				vals = (drug_id,drug_abrv,dname,fname,bname,drug_class,mfc,dec,hsc,lc,clc,uc,drug_flg)
				stat = execPG(pgc,sql,vals,LOG)
				if stat == 'Error': 
					LOG.write(gt()+"\tFailed to insert row in Postgres, wh.Drugs update is incomplete\n")
					return 'bad'
			LOG.write(gt()+"\twh.drugs table was updated successfully..."+str(len(cxr))+" records added\n")
			return 'ok'
	else:
		LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.Drugs was not updated\n")
		return 'bad'


# Update wh.Projects table
#############################
def updateProjects(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.project....\n")
	cxr = getOracle(cxc,update_projects_sql,'',LOG)
	if cxr != "Error":
		sql = "delete from wh.project"
		pgr = execPG(pgc,sql,'',LOG)
		if pgr != "Error":
			LOG.write(gt()+"\tWriting update to wh.project table...\n")
			for project_id,project_code,project_name,project_type in cxr:
				sql = "insert into wh.project values(%s,%s,%s,%s)"
				vals = (project_id,project_code,project_name,project_type)
				stat = execPG(pgc,sql,vals,LOG)
				if stat == 'Error': 
					LOG.write(gt()+"\tFailed to insert row in Postgres, wh.Project update is incomplete\n")
					return
			LOG.write(gt()+"\twh.project table was updated successfully..."+str(len(cxr))+" records added\n")
			return 'ok'
	else:
		LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.Project was not updated\n")
		return 'bad'


# Update wh.Mutations table
#############################
def updateMutations(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.mutations....\n")
	# Get last update date from wh.mutations
	sql = "select max(created_date) from wh.mutations" #max(to_char(to_timestamp(created_date,'DD-MON-YYYY HH24:MI:SS'),'YYYY-MM-DD HH24:MI:SS')) from wh.mutations"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_mutations_sql%last_update 
		cxr = getOracle(cxc,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:
				LOG.write(gt()+"\t\tWriting update to wh.mutations table...\n")
				# Insert new records in wh.mutations
				for mutation_id,mutation_type,mutation,position,wt,mut,insertion_flg,created_date in cxr:
					sql = "insert into wh.mutations values (%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (mutation_id,mutation_type,mutation,position,wt,mut,insertion_flg,created_date)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\twh.mutations was updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'
			else:
				LOG.write(gt()+"\twh.mutations update was skipped, there are no records to add today\n")
				return 'ok'

		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.Mutations was not updated\n")
			return 'bad'
		
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.Mutations was not updated\n")
		return 'bad'
	return 'ok'


# Check if mutant is a RAM
############################
def checkRAMS(mutant, mut_type):
	global summary
	pattern = re.compile(r'(\w(\d+)([\w+\^*]))')	# Regular expression pattern to identify letter|number(s)|letter(s) or ^ or * pattern
	
	match = re.findall(pattern,mutant)
	p = match[0][1]		# position
	r = match[0][2]		# mutation or insertion amino acids

	score = 0
	if p=='69' and (mut_type == 'tams' or mut_type == 'nams' or mut_type == 'nnrti'):
		score = 1
		if summary[mut_type] != '': summary[mut_type] = summary[mut_type] + ", "+ mutant
		else: summary[mut_type] = mutant
		return score
		
	if muts[mut_type].has_key(mutant):
		score = 1
		if summary[mut_type] != '': summary[mut_type] = summary[mut_type] + ", "+ mutant
		else: summary[mut_type] = mutant
	return score


# Record all mutants present in sequence
##########################################
def seenMut(mutant,mut_type):
   global seen
   if mut_type == 'tams' or mut_type == 'nams' or mut_type == 'nnrti': mut_type = 'rt'
   if seen[mut_type].has_key(mutant): seen[mut_type][mutant] = seen[mut_type][mutant] + 1
   else: seen[mut_type][mutant] = 1
   return


# Evaluate each mutant
############################
def processMuts(mutation_summary,mut_type,check_freq):
	mut_score = 0
      
	# Process mutation list for sequences
	for mut in mutation_summary:
		if check_freq == 'yes': seenMut(mut,mut_type)

		# Check for single mutants first
		if mut.find('/') == -1: mut_score = mut_score + checkRAMS(mut,mut_type)
      
		# Have multiple mutants
		else:
			mut2 = mut.split('/')   # ['I50A','C','H']
			mut_count = len(mut2)
			freq = Decimal(1)/Decimal(mut_count)

			if mut2[0][0] != mut2[0][-1]:
				# check first mutation in list
				if check_freq == 'yes': seenMut(mut2[0], mut_type)
				s = checkRAMS(mut2[0], mut_type)
				if s == 1: mut_score = mut_score + freq
            
			# Generate the rest of mutations in list and check them.
			for i in range(1,len(mut2)):
				m = mut2[0][:-1] + mut2[i] # automatically create I50C, I50H etc
				if check_freq == 'yes': seenMut(m, mut_type)       # Add the individual mutations of a mixture to the seen total
				s = checkRAMS(m, mut_type)
				if s == 1: mut_score = mut_score + freq
	return mut_score



# Update wh.seq table
########################
def updateSeq(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.seq....\n")
	sql = "select max(rep_date) from wh.seq" # max(to_char(to_timestamp(rep_date,'DD-MON-YYYY HH24:MI:SS'),'YYYY-MM-DD HH24:MI:SS')) from wh.seq"
	pgr = getPG(pgc,sql,'',LOG)

	if pgr != "Error":
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_seq_sql % (pgr[0][0],pgr[0][0]) #last_update 
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		cxr = getOracle(cxc,sql,'',LOG)

		if cxr != "Error":
			
			# Verify there are records to insert
			if len(cxr) != 0:
				data_ids = []
				for r in cxr:
					data_ids.append(r[0])
				sql = "select distinct seq_id from wh.seq where seq_id in %s"
				r_ids = tuple(x for x in data_ids)
				pgresult = getPG(pgc,sql,(r_ids,),LOG)

				dup_ids = []  # use this to capture ids that are dups and will be dropped, so that we can not include them when we do mut_freq update
				if pgresult !="Error":
					if len(pgresult) !=0:
						for x in pgresult:
							sql = "delete from wh.seq where seq_id=%s" % x[0]
							bad_id = x[0]
							execPG(pgc,sql,'',LOG)						
							dup_ids.append(x[0])
				else:
					LOG.write(gt()+"\twh.seq and wh.mut_freq updates was skipped, there was a problem checking for existing wh.seq records\n")
					return	'bad'


				# Insert new records in wh.seq
				getcontext().prec=2     # Sets precision for decimal to 2 digits
				global seen, summary, scores, muts
				muts = {'pri':PI_RAMS,'tams':TAMS,'nams':NAMS,'nnrti':NNRTI_RAMS,'int':INT_RAMS,'rt2':RT2_RAMS}
				
				LOG.write(gt()+"\t\tWriting update to wh.seq and wh.mut_freq table...\n")
				for row in cxr:
					# Begin processing record (Oracle export for building HIVWH - contains pri_summary and rt_summary fields)
					vals =[]
					for x in row: vals.extend([x]) # Convert tuple from db into array
					seen = {'pri':{},'rt':{},'int':{},'rt2':{}}     # Record of all protease mutations that have been seen, whether they are RAMs or not
					summary = {'tams':'','nams':'','nnrti':'','pri':'','int':'','rt2':''}
					scores = {'tams':Decimal(0),'nams':Decimal(0),'nnrti':Decimal(0),'pri':Decimal(0),'int':Decimal(0),'rt2':Decimal(0)}
					wt=''
					
					check_freq='yes'
					if vals[0] in dup_ids: check_freq = 'no'
					
					if vals[5] is not None and len(vals[5]) > 0:
						pri_sum = vals[5].replace(' ','').replace('"','').split(',')
						scores['pri'] = processMuts(pri_sum,'pri',check_freq)

					if vals[6] is not None and len(vals[6]) > 0:
						rt_sum = vals[6].replace(' ','').replace('"','').split(',')
						scores['tams'] = processMuts(rt_sum,'tams',check_freq)
						scores['nams'] = processMuts(rt_sum,'nams',check_freq)
						scores['nnrti'] = processMuts(rt_sum,'nnrti',check_freq)

					if vals[7] is not None and len(vals[7]) > 0:
						int_sum = vals[7].replace(' ','').replace('"','').split(',')
						scores['int'] = processMuts(int_sum,'int',check_freq)

					if vals[8] is not None and len(vals[8]) > 0:
						rt2_summary = vals[8].replace(' ','').replace('"','').split(',')
						scores['rt2'] = processMuts(rt2_summary,'rt2',check_freq)

					if scores['tams'] + scores['nams'] + scores['nnrti'] + scores['pri'] + scores['int'] == 0: wt='y'
					else: wt = 'n'

					if vals[5] is not None: vals[5] = vals[5].replace('"','')
					if vals[6] is not None: vals[6] = vals[6].replace('"','')
					if vals[7] is not None: vals[7] = vals[7].replace('"','')
					if vals[8] is not None: vals[8] = vals[8].replace('"','')
      			
					# a = data_file_id,id,vrl_Accession,aliquot_id,nucleotide_seq,pri_summary,rt_summary,in_summary,rt2_summary,assay_type,rep_date,project
					# joined to: wt,tam_sum,nam_sum,nnrti_sum,pri_sum,int_sum,rt2_sum,tam_score,nam_Score,nnrti_Score,pri_score,int_score,rt2_score
					vals.extend([wt,summary['tams'],summary['nams'],summary['nnrti'],summary['pri'],summary['int'],summary['rt2'],scores['tams'],scores['nams'],scores['nnrti'],scores['pri'],scores['int'],scores['rt2']])
					sql = "insert into wh.seq values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
					execPG(pgc,sql,vals,LOG)

					# Also update wh.mut_freq table concurrently, since numbers come from these calculations
					for x in seen['pri']:	# Update protease mutant_freq
						m = x.replace('/','')
						sql = "select mut_freq, id from wh.mut_freq where mutation=%s and mut_type='Protease'"
						pgr2 = getPG(pgc,sql,[m,],LOG)

						if len(pgr2) == 0:
							sql = "insert into wh.mut_freq(mutation,mutation_sep,mut_freq,mut_type) values(%s,%s,%s,%s)"
							vals2 = (m,x,seen['pri'][x],'Protease')
							execPG(pgc,sql,vals2,LOG)
						else:
							sql = "update wh.mut_freq set mut_freq=%s where id=%s"
							c = int(pgr2[0][0])+1
							sid = pgr2[0][1]
							execPG(pgc,sql,(c,sid),LOG)

					for y in seen['rt']:		# Update RT mutant_freq
						m = y.replace('/','')
						sql = "select mut_freq, id from wh.mut_freq where mutation=%s and mut_type='RT'"
						pgr2 = getPG(pgc,sql,[m,],LOG)

						if len(pgr2) == 0:
							sql = "insert into wh.mut_freq(mutation,mutation_sep,mut_freq,mut_type) values(%s,%s,%s,%s)"
							vals2 = (m,y,seen['rt'][y],'RT')
							execPG(pgc,sql,vals2,LOG)
						else:
							sql = "update wh.mut_freq set mut_freq=%s where id = %s"
							c = int(pgr2[0][0])+1
							sid = pgr2[0][1]
							execPG(pgc,sql,(c,sid),LOG)

					for z in seen['int']:	# Update INT mutatnt_freq
						m = z.replace('/','')
						sql = "select mut_freq, id from wh.mut_freq where mutation=%s and mut_type='Int'"
						pgr2 = getPG(pgc,sql,[m,],LOG)
						
						if len(pgr2) == 0:
							sql = "insert into wh.mut_freq(mutation,mutation_sep,mut_freq,mut_type) values(%s,%s,%s,%s)"
							vals2 = (m,x,seen['int'][z],'Int')
							execPG(pgc,sql,vals2,LOG)
						else:
							sql = "update wh.mut_freq set mut_freq=%s where id= %s"
							c = int(pgr2[0][0])+1
							sid = pgr2[0][1]
							execPG(pgc,sql,(c,sid),LOG)

					for w in seen['rt2']:	# Update RT2 mutant_freq
						m = w.replace('/','')
						sql = "select mut_freq, id from wh.mut_freq where mutation=%s and mut_type='RT2'"
						pgr2 = getPG(pgc,sql,[m,],LOG)

						if len(pgr2) == 0:
							sql = "insert into wh.mut_freq(mutation,mutation_sep,mut_freq,mut_type) values(%s,%s,%s,%s)"
							vals2 = (m,x,seen['rt2'][w],'RT2')
							execPG(pgc,sql,vals2,LOG)
						else:
							sql = "update wh.mut_freq set mut_freq=%s where id= %s"
							c = int(pgr2[0][0])+1
							sid = pgr2[0][1]
							execPG(pgc,sql,(c,sid),LOG)

				LOG.write(gt()+"\twh.seq and wh.mut_freq were updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'
			else:
				LOG.write(gt()+"\twh.seq and wh.mut_freq updates was skipped, there are no records to add today\n")
				return 'ok'
		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.Seq was not updated\n")
			return 'bad'
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.Seq was not updated\n")
		return 'bad'
	return 'ok'


# Update wh.Seq_Mutations table
#############################
def updateSeqMutations(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.seq_mutations....\n")
	# Get last update date from wh.seq_mutations
	sql = "select max(rep_date) from wh.seq_mutations"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_seq_mutations_sql%last_update 
		cxr = getOracle(cxc,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:

				# check for duplicate gt_report_id
				data_ids = []
				for r in cxr:
					data_ids.append(r[0])

				sql = "select distinct gt_report_id from wh.seq_mutations where gt_report_id in %s"
				r_ids = tuple(x for x in data_ids)
				pgresult = getPG(pgc,sql,(r_ids,),LOG)

				dup_ids = []  # use this to capture ids that are dups and will be dropped, so that we can not include them when we do mut_freq update
				if pgresult !="Error":
					if len(pgresult) !=0:
						for x in pgresult:
							sql = "delete from wh.seq_mutations where gt_report_id=%s" % x[0]
							execPG(pgc,sql,'',LOG)						
							dup_ids.append(x[0])
				else:
					LOG.write(gt()+"\twh.seq_mutations update was skipped, there was a problem checking for existing wh.seq_mutations records\n")
					return 'bad'
				
				LOG.write(gt()+"\t\tWriting update to wh.seq_mutations table...\n")
				# Insert new records in wh.seq_mutations
				for gt_report_id,seq_id,mutation_id,mutation_type,mutation,rep_date in cxr:
					sql = "insert into wh.seq_mutations(gt_report_id,seq_id,mutation_id,mutation_type,mutation,rep_date) values (%s,%s,%s,%s,%s,%s)"
					vals = (gt_report_id,seq_id,mutation_id,mutation_type,mutation,rep_date)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\twh.seq_mutations was updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'
			else:
				LOG.write(gt()+"\twh.seq_mutations update was skipped, there are no records to add today\n")
				return 'ok'

		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.seq_mutations was not updated\n")
			return 'bad'
		
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.seq_mutations was not updated\n")
		return 'bad'
	return 'ok'



# Update wh.PSGT_Result table
#############################
def updatePSGTResults(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.psgt_result....\n")
	# Get last update date from wh.psgt_result
	sql = "select max(rep_date) from wh.psgt_result" # max(to_char(to_timestamp(rep_date,'DD-MON-YYYY HH24:MI:SS'),'YYYY-MM-DD HH24:MI:SS')) from wh.psgt_result"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_psgtresult_sql%last_update 
		cxr = getOracle(cxc,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:

				LOG.write(gt()+"\t\tchecking for duplicates...\n")
				# check for duplicate gt_report_ids
				data_ids = []
				for r in cxr:
					data_ids.append(r[0])

				sql = "select distinct pg_report_id from wh.psgt_result where pg_report_id in %s"
				r_ids = tuple(x for x in data_ids)
				pgresult = getPG(pgc,sql,(r_ids,),LOG)

				dup_ids = []  # use this to capture ids that are dups and will be dropped, so that we can not include them when we do mut_freq update
				if pgresult !="Error":
					if len(pgresult) !=0:
						for x in pgresult:
							sql = "delete from wh.psgt_result where pg_report_id=%s" % x[0]
							execPG(pgc,sql,'',LOG)						
							dup_ids.append(x[0])
				else:
					LOG.write(gt()+"\twh.psgt_result update was skipped, there was a problem checking for existing wh.psgt_result records\n")
					return 'bad'

				# Insert new records in wh.psgt_result
				LOG.write(gt()+"\t\tWriting update to wh.psgt_result table...\n")
				for pg_report_id,vrl_accession,pt_report_id,gt_report_id,test_code,project_code,drug_abrv,resistant_flg,result_flg,rep_date in cxr:
					sql = "insert into wh.psgt_result(pg_report_id,vrl_accession,pt_report_id,gt_report_id,test_code,project_code,drug_abrv,resistant_flg,result_flg,rep_date) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (pg_report_id,vrl_accession,pt_report_id,gt_report_id,test_code,project_code,drug_abrv,resistant_flg,result_flg,rep_date)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\twh.psgt_result was updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'
			else:
				LOG.write(gt()+"\twh.psgt_result update was skipped, there are no records to add today\n")
				return 'ok'

		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.psgt_result was not updated\n")
			return 'bad'
		
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.psgt_result was not updated\n")
		return 'bad'
	return 'ok'



# Update wh.Geno_Report and wh.Geno_Rep_Rows table
####################################################
def updateGenoReport(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.geno_report and wh.geno_rep_row....\n")
	# Get last update date from wh.geno_report
	sql = "select max(rep_date) from wh.geno_report" # max(to_char(to_timestamp(rep_date,'DD-MON-YYYY HH24:MI:SS'),'YYYY-MM-DD HH24:MI:SS')) from wh.geno_report"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_geno_report_sql%last_update 
		cxr = getOracle(cxc,sql,'',LOG)
		dup_ids = []  # use this to capture ids that are dups and will be dropped, so that we can not include them when we do mut_freq update

		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:
				# check for duplicate gt_report_ids
				data_ids = []
				for r in cxr:
					data_ids.append(r[0])
				sql = "select distinct gt_report_id from wh.geno_report where gt_report_id in %s"
				r_ids = tuple(x for x in data_ids)
				pgresult = getPG(pgc,sql,(r_ids,),LOG)
		
				if pgresult !="Error":
					if len(pgresult) !=0:
						for x in pgresult:
							sql = "delete from wh.geno_report where gt_report_id=%s" % x[0]
							bad_id = x[0]
							execPG(pgc,sql,'',LOG)						
							dup_ids.append(x[0])
				else:
					LOG.write(gt()+"\twh.geno_report and wh.geno_rep_row updates was skipped, there was a problem checking for existing wh.geno_report records\n")
					return 'bad'
				
				LOG.write(gt()+"\t\tWriting update to wh.geno_report table...\n")
				# Insert new records in wh.geno_report
				for gt_report_id,vrl_accession,aliquot_id,subtype,rule_ver,rep_date,test_code,assay_type in cxr:
					sql = "insert into wh.geno_report values (%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (gt_report_id,vrl_accession,aliquot_id,subtype,rule_ver,rep_date,test_code,assay_type)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\t\twh.geno_report was updated successfully..."+str(len(cxr))+" records added\n")
				
			else:
				LOG.write(gt()+"\twh.geno_report update was skipped, there are no records to add today\n")
				return 'ok'

		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.geno_report was not updated\n")
			return 'bad'
		sql = update_geno_rep_row_sql%last_update
		cxr = getOracle(cxc,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:

				# remove rep_rows that have the same dup gt_report_id's
				if len(dup_ids) !=0:
					for x in dup_ids:
						sql = "delete from wh.geno_rep_row where gt_report_id=%s" % x
						execPG(pgc,sql,'',LOG)
				
				LOG.write(gt()+"\t\tWriting update to wh.geno_rep_row table...\n")
				# Insert new records in wh.geno_rep_row
				for report_id,gt_report_id,drug_id,result,reported_drug_flg,major_mut_sum,resistant_flg,drug_abrv in cxr:
					sql = "insert into wh.geno_rep_row values (%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (report_id,gt_report_id,drug_id,result,reported_drug_flg,major_mut_sum,resistant_flg,drug_abrv)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\twh.geno_rep_row was updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'

			else:
				LOG.write(gt()+"\twh.geno_rep_row update was skipped, there are no records to add today\n")
				return 'ok'
		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.geno_rep_row was not updated\n")
			return 'bad'		
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.geno_report and wh.geno_rep_row were not updated\n")
		return 'bad'
	return 'ok'



# Update wh.Pheno_Report and wh.Pheno_Rep_Rows table
######################################################
def updatePhenoReport(pgc,cxc,LOG):
	LOG.write(gt()+"\tUpdating wh.pheno_report and wh.pheno_rep_row....\n")
	# Get last update date from wh.pheno_report
	sql = "select max(rep_date) from wh.pheno_report" # max(to_char(to_timestamp(rep_date,'DD-MON-YYYY HH24:MI:SS'),'YYYY-MM-DD HH24:MI:SS')) from wh.pheno_report"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_pheno_report_sql%last_update 
		cxr = getOracle(cxc,sql,'',LOG)
		dup_ids = []  # use this to capture ids that are dups and will be dropped, so that we can not include them when we do mut_freq update
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:
				# check for duplicate pt_report_ids
				data_ids = []
				for r in cxr:
					data_ids.append(r[0])
				sql = "select distinct pt_report_id from wh.pheno_report where pt_report_id in %s"
				r_ids = tuple(x for x in data_ids)
				pgresult = getPG(pgc,sql,(r_ids,),LOG)

				
				if pgresult !="Error":
					if len(pgresult) !=0:
						for x in pgresult:
							sql = "delete from wh.pheno_report where pt_report_id=%s" % x[0]
							bad_id = x[0]
							execPG(pgc,sql,'',LOG)						
							dup_ids.append(x[0])
				else:
					LOG.write(gt()+"\twh.pheno_report and wh.pheno_rep_row updates was skipped, there was a problem checking for existing wh.pheno_report records\n")
					return 'bad'
				
				LOG.write(gt()+"\t\tWriting update to wh.pheno_report table...\n")
				# Insert new records in wh.pheno_report
				for pt_report_id,vrl_accession,aliquot_id,test_code,assay_type,rc_status,reported_rc,rep_date,project_code in cxr:
					sql = "insert into wh.pheno_report values (%s,%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (pt_report_id,vrl_accession,aliquot_id,test_code,assay_type,rc_status,reported_rc,rep_date,project_code)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\t\twh.pheno_report was updated successfully..."+str(len(cxr))+" records added\n")
				
			else:
				LOG.write(gt()+"\twh.pheno_report update was skipped, there are no records to add today\n")
				return 'ok'

		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.pheno_report was not updated\n")
			return 'bad'
		sql = update_pheno_rep_row_sql%last_update
		cxr = getOracle(cxc,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:
				
				# remove rep_rows that have the same dup pt_report_id's
				if len(dup_ids) !=0:
					for x in dup_ids:
						sql = "delete from wh.pheno_rep_row where pt_report_id=%s" % x
						execPG(pgc,sql,'',LOG)						
				
				LOG.write(gt()+"\tWriting update to wh.pheno_rep_row table...\n")
				# Insert new records in wh.pheno_rep_row
				for report_id,pt_report_id,ic50,ic95,fold_res,fold_change,resistant_flg,result_flg,reported_drug_flg,drug_class,drug_abrv in cxr:
					sql = "insert into wh.pheno_rep_row values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (report_id,pt_report_id,ic50,ic95,fold_res,fold_change,resistant_flg,result_flg,reported_drug_flg,drug_class,drug_abrv)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\twh.pheno_rep_row was updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'
			else:
				LOG.write(gt()+"\twh.pheno_rep_row update was skipped, there are no records to add today\n")
				return 'ok'
		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.pheno_rep_row was not updated\n")
			return 'bad'		
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.pheno_report and wh.pheno_rep_row were not updated\n")
		return 'bad'
	return 'ok'



# Remove newline. p= pattern, x= string to search
def rn(p,x):
	m = p.search(x)
	if m: x = x.replace('\\','')
	return x



# Update wh.Patient table
#############################
def updatePatient(pgc,cxc,LOG):
	global temp_pth
	LOG.write(gt()+"\tUpdating wh.patient....\n")
	# Get last update date from wh.patient
	sql = "select max(received_date) from wh.patient"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_patient_sql%last_update 
		cxr = getOracle(cxc,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:
				# Drop and recreate patient_load table
				# Using patient_load is not strictly necessary, but protects patient table in case there is a problem with load
				pxr = execPG(pgc,create_pt_load,'',LOG)
				LOG.write(gt()+"\t\t created patient_load table\n")
				# Need to search fields for hidden carriage returns that would break the load, so make patterns here
				p1 = re.compile(r'\\')

				LOG.write(gt()+"\t\tWriting update to wh.patient_load table...\n")
				# Insert new records in wh.patient_load
				for vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs in cxr:
					sql = "insert into wh.patient_load(vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
					client_name = rn(p1,client_name)
					vals = (vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs)
					pxr2 = execPG(pgc,sql,vals,LOG)

				LOG.write(gt()+"\t\twh.patient_load was updated successfully..."+str(len(cxr))+" records added\n")
			
				LOG.write(gt()+"\t\tWorking with patient_load table...\n")
				# Begin prepping patient_load table for transfer to patient table
				sql = "update wh.patient_load set pt_hash = md5(lastname||firstname||sex||dob) where pt_hash is null"
				r1 = execPG(pgc,sql,'',LOG)

				if r1 != 'Error':
					sql = "update wh.patient_load set pt_hash = md5(lastname||sex||dob) where pt_hash is null"
					r2 = execPG(pgc,sql,'',LOG)
				
					if r2 != 'Error':
						sql = "update wh.patient_load set pt_hash = md5(lastname||dob) where pt_hash is null"
						r3 = execPG(pgc,sql,'',LOG)
						
						# Nothing else is uniquely identifiable, so delete it. Should be very small number
						if r3 !='Error':
							sql = "delete from wh.patient_load where pt_hash is null"
							r4 = execPG(pgc,sql,'',LOG)
							LOG.write(gt()+"\t\tpt_hash keys have been generated for patient_load, assigning test numbers...\n")

						else:
							LOG.write(gt()+"\twh.patient update was skipped, there was an error generating pt_hash\n")
							return 'bad'
					else:
						LOG.write(gt()+"\twh.patient update was skipped, there was an error generating pt_hash\n")
						return 'bad'
			
				# Generate test_numbers
				LOG.write(gt()+"\t\tUpdating patient test numbers...\n")
				f = open("temp_patient_data","w")
				sql = "copy (select pt_temp_id,pt_hash from wh.patient_load order by pt_hash,collected_date) to STDOUT"
				try:
					pgc.copy_expert(sql,f)
					f.close()
				except:
					LOG.write(gt()+"\t\tCould not execute COPY from patient_load. Update aborted for wh.patient\n")
					return 'bad'

				f2 = open(temp_pth,"w")
				sql = "copy (select distinct pt_hash,max(test_num) from patient group by pt_hash) to STDOUT"
				try:
					pgc.copy_expert(sql,f2)
					f2.close()
				except:
					LOG.write(gt()+"\t\tCould not execute COPY from patient. Update aborted for wh.patient\n")
					return 'bad'
				
				f3 = open(temp_pth,"r")
				existing = {}
				for line in f3:
					a = line.split('\t')
					existing[a[0]] = int(a[1].rstrip())
				f3.close()

				f4 = open("temp_patient_data","r")
				c = 1
				cid = ''
				sql = "update wh.patient_load set test_num=%s where pt_temp_id=%s"

				for line in f4:
					a = line.split('\t')
					a[-1].rstrip()  # if did not export with trailing null column, use this
	 				
					if existing.has_key(a[1]) and cid != a[1]:
						c = existing[a[1]] + 1
						cid = a[1]
					elif existing.has_key(a[1]) and cid == a[1]: c = c + 1
					else:
						if a[1] != cid:
							cid = a[1]
							c = 1
						elif a[1]==cid: c = c + 1
					execPG(pgc,sql,[c,a[0]],LOG)
				LOG.write(gt()+"\t\tFinished assigning test_numbers and loading wh.patient_load\n")

				pt_rec = 0
				# Make the transfer to wh.patient
				sql = "select count(pt_hash) from wh.patient_load"
				r5 = getPG(pgc,sql,'',LOG)
				pt_rec = r5[0][0]
				try:
					r6 = execPG(pgc,load_patient_sql,'',LOG)
				except:
					LOG.write(gt()+"\t\tCould not transfer data from patient_load to wh.patient. Update of wh.patient aborted\n")
					return 'bad'
				LOG.write(gt()+"\twh.patient was loaded successfully..."+str(pt_rec)+" of "+str(len(cxr))+" records added.\n\t\t\tAny that remainded could not make pt_hash and were dropped\n")
				execPG(pgc,"drop table patient_load",'',LOG)
				f4.close()
				os.remove(temp_pth)
				os.remove('temp_patient_data')
				return 'ok'

			else:
				LOG.write(gt()+"\twh.patient update was skipped, there are no records to add today\n")
				return 'ok'
		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.patient was not updated\n")
			return 'bad'
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.patient was not updated\n")
		return 'bad'
	return 'ok'



# Update wh.GT_Rule_Result_Rows table from MySQL
##################################################
def updateGTRuleResultRows(pgc,kong,kong_db,LOG):
	LOG.write(gt()+"\tUpdating wh.gt_rule_result_rows....\n")
	# Update MGRM data
	# Get last update date from wh.gt_rule_result_rows and data_src='mgrm'
	sql = "select max(created_date) from wh.gt_rule_result_rows where data_src='mgrm'"
	pgr = getPG(pgc,sql,'',LOG)
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update from MGRM table was "+pgr[0][0]+"\n")
		
		# Get Kong records > last_update
		mgrm = getMySQL(kong,kong_update_sql,pgr[0][0],LOG)
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
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.gt_rule_result_rows update part 1 was not done\nUpdate is aborted\n")
		return 'bad'
	
	# Update LABC data
	# Get last update date from wh.gt_rule_result_rows and data_src='labc'
	sql = "select max(created_date) from wh.gt_rule_result_rows where data_src='labc'"
	pgr2 = getPG(pgc,sql,'',LOG)
	if pgr2 != "Error":
		last_update = (pgr2[0][0],)
		LOG.write(gt()+"\t\tLast update from LABC table was "+pgr2[0][0]+"\n")
		
		# Get KongDB records > last_update
		labc = getMySQL(kong_db,kongdb_update_sql,pgr2[0][0],LOG)
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
				return 'ok'
			else:
				LOG.write(gt()+"\twh.gt_rule_result_rows update part 2 was skipped, there are no records to add today\n")
				return 'ok'
		else:
			LOG.write(gt()+"\tFailed to retrieve update from KONG_DB, wh.gt_rule_result_rows update part 2 was not done\n")
			return 'bad'
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.gt_rule_result_rows update part 2 was not done\n")
		return 'bad'
	return 'ok'


# update wh.etag_sum from Oracle etag tables
###############################################
def updateEtag(pgc,etg,LOG):
	LOG.write(gt()+"\tUpdating wh.etag_sum....\n")
	sql = "select max(rep_date) from wh.etag_sum"
	pgr = getPG(pgc,sql,'',LOG)
	getcontext().prec=2     # Sets precision for decimal to 2 digits
	if pgr != "Error":
		last_update = (pgr[0][0],)
		LOG.write(gt()+"\t\tLast update was "+pgr[0][0]+"\n")
		# Get Oracle records > last_update
		# date format causes problems here, so I use normal string formatting instead of relying on
		# the cx_Oracle cursor to format.  Also avoids using this ==>	cxc.execute('ALTER SESSION SET NLS_DATE_FORMAT = \'YYYY-MM-DD HH24:MI:SS\'')
		sql = update_etag_sql%last_update 
		cxr = getOracle(etg,sql,'',LOG)
		if cxr != "Error":
			# Verify there are records to insert
			if len(cxr) != 0:
				LOG.write(gt()+"\t\tWriting update to wh.etag_sum table...\n")
				# Insert new records in wh.etag_sum
				for vrl_acc,dob,col_date,rcv_date,rep_date,h2t,h2d,h2t_r,h2d_r,h2t_ta,h2d_ta,h2t_ss,h2d_ss,h2t_pt,h2d_pt,h2t_id,h2d_id,state,country in cxr:
					sql = "insert into wh.etag_sum (vrl_accession,patient_dob,collected_date,rcv_date,rep_date,h2t,h2d,h2t_result,h2d_result,h2t_tumor_area,h2d_tumor_area,h2t_section_size,h2d_section_size,h2t_pct_tumor,h2d_pct_tumor,h2t_batch_id,h2d_batch_id,client_state,client_country) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
					vals = (vrl_acc,dob,col_date,rcv_date,rep_date,Decimal(str(h2t)),Decimal(str(h2d)),h2t_r,h2d_r,Decimal(str(h2t_ta)),Decimal(str(h2d_ta)),Decimal(str(h2t_ss)),Decimal(str(h2d_ss)),Decimal(str(h2t_pt)),Decimal(str(h2d_pt)),h2t_id,h2d_id,state,country)
					pxr = execPG(pgc,sql,vals,LOG)
				LOG.write(gt()+"\twh.etag_sum was updated successfully..."+str(len(cxr))+" records added\n")
				return 'ok'
			else:
				LOG.write(gt()+"\twh.etag_sum update was skipped, there are no records to add today\n")
				return 'ok'

		else:
			LOG.write(gt()+"\tFailed to retrieve update from Oracle, wh.etag_sum was not updated\n")
			return 'bad'
		
	else:
		LOG.write(gt()+"\tFailed to retrieve last_update from Postgres, wh.etag_sum was not updated\n")
		return 'bad'

	return 'ok'


# Send mail upon failure to update
######################################	
def sendMail(log_pth):
	import smtplib
	from email.mime.text import MIMEText
	fp = open(log_pth, 'rb')
	msg = MIMEText(fp.read())	# read logfile into message body
	fp.close()
	me = 'gamera_server@monogrambio.com'
	you = 'bioinformatics@monogrambio.com'
	msg['Subject'] = 'HIVdb Warehouse update failed last night'
	msg['From'] = me
	msg['To'] = you
	# Send the message via our own SMTP server, but don't include the envelope header.
	s = smtplib.SMTP('autofax.virologic.com')
	s.sendmail(me, [you], msg.as_string())
	s.quit()
	return


# Check status of proceedure
#################################
def checkStatus(status):
	global notify
	if status =='bad': notify = 'yes'
	return


# End program if can't make connection
########################################
def endProgram(LOG,log_pth):
	LOG.write(gt()+"\tDatabase update failed because could not make database connections\nProcess terminated\n")
	LOG.close()
	sendMail(log_pth)
	sys.exit()



# Main Program
###############################
def main():
	global temp_pth, notify
	#temp_pth="/Users/Bali/Documents/Monogram/Genemachine/update_dev/existing.txt"		# Use for testing on Bali 
	temp_pth="/opt/bin/existing.txt"													# Use for production on Gamera
	log_pth = "/home/mevans/logs/hivdb_update_log.txt"									# Use for production on Gamera
	#log_pth = "dev_log.txt"															# Use for testing on Bali
	buffer_size = 0
	LOG = open(log_pth,"w",buffer_size)					
	LOG.write("***********************\n"+gt()+" Beginning update process...\n")

	notify = 'no'	# flag to email logfile if there are problems
	status = 'ok'	# flag to be set by each update process
	cxc = ''
	kong = ''
	kong_db = ''

	# Create database connections
	##############################
	status, pgc, pg_con = connectPG(LOG)	#** Don't forget to switch Postgres connection string when going from production to testing & vice versa
	if status == 'ok':	status, cxc, etg = connectORA(LOG)
	else: endProgram(LOG,log_pth)
	if status == 'ok':	status, kong, kong_db = connectMySQL(LOG)
	else: endProgram(LOG,log_pth)

	if status == 'ok':
		# Begin processing tables
  		checkStatus(updateDrugs(pgc,cxc,LOG))
  		checkStatus(updateProjects(pgc,cxc,LOG))
  		checkStatus(updateMutations(pgc,cxc,LOG))
  		checkStatus(updateSeq(pgc,cxc,LOG)) # updates both Seq and Mut_freq tables
  		checkStatus(updateSeqMutations(pgc,cxc,LOG))
		checkStatus(updatePSGTResults(pgc,cxc,LOG))
		checkStatus(updateGenoReport(pgc,cxc,LOG))  # updates both Geno_Report and Geno_Rep_Row tables
		checkStatus(updatePhenoReport(pgc,cxc,LOG)) # updates both Pheno_Report and Pheno_Rep_Row tables
		checkStatus(updatePatient(pgc,cxc,LOG))
		checkStatus(updateGTRuleResultRows(pgc,kong,kong_db,LOG))
		checkStatus(updateEtag(pgc,etg,LOG))

  		# Commit all changes to the database
 		pg_con.commit()
		print "Database update complete\n"
		LOG.write(gt()+"\tDatabase update completed\n##############################\n")
		LOG.close()

		# Check for errors and email logfile if ther are any
		if notify == 'yes':
			sendMail(log_pth)
			sys.exit()
			
	else: endProgram(LOG,log_pth)

if __name__ == '__main__':
	main()

