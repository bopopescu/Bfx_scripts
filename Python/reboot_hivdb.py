#!/usr/bin/env python
# encoding: utf-8
"""
updateHIVdb.py

Created by Mark Evans on 2012-03-16.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""


import sys, os, re, csv
from types import *
os.environ["LD_LIBRARY_PATH"]="/usr/lib/oracle/11.2/client64/lib"
os.environ["ORACLE_HOME"]="/usr/lib/oracle/11.2/client64"
os.environ["TNS_ADMIN"]="/usr/lib/oracle/11.2/client64"
import cx_Oracle as cxo
import psycopg2 as pg
import pymysql as ms
from time import localtime, strftime
from cStringIO import StringIO
from decimal import *

# NRTI TAMS: M41L,D67N,K70R,L210W,T215F,T215Y,K219
TAMS = {'M41L':'',
        'D67N':'',
        'K70R':'',
        'L210W':'',
        'T215F':'','T215Y':'',
        'K219A':'','K219C':'','K219D':'','K219E':'','K219F':'','K219G':'','K219H':'','K219I':'','K219L':'','K219M':'','K219N':'','K219P':'','K219Q':'','K219R':'','K219S':'','K219T':'','K219V':'','K219W':'','K219Y':'',
        }
        
# NRTI NAMS: K65R,T69,K70E,L74,V75A,V75M,V75S,V75T,Y115F,Q151M,M184
NAMS = {'K65R':'',
        'T69A':'','T69C':'','T69D':'','T69E':'','T69F':'','T69G':'','T69H':'','T69I':'','T69K':'','T69L':'','T69M':'','T69N':'','T69P':'','T69Q':'','T69R':'','T69S':'','T69V':'','T69W':'','T69Y':'','T69^':'',
        'K70E':'',
        'L74A':'','L74C':'','L74D':'','L74E':'','L74F':'','L74G':'','L74H':'','L74I':'','L74K':'','L74M':'','L74N':'','L74P':'','L74Q':'','L74R':'','L74S':'','L74T':'','L74V':'','L74W':'','L74Y':'',
        'V75A':'','V75M':'','V75S':'','V75T':'',
        'Y115F':'',
        'Q151M':'',
        'M184A':'','M184C':'','M184D':'','M184E':'','M184F':'','M184G':'','M184H':'','M184I':'','M184K':'','M184L':'','M184N':'','M184P':'','M184Q':'','M184R':'','M184S':'','M184T':'','M184V':'','M184W':'','M184Y':''
        }
        
# NRTI RAMS (TAMS + NAMS): M41L,K65R,D67N,T69,K70R,K70E,L74,V75A,V75M,V75S,V75T,Y115F,Q151M,M184,L210W,T215F,T215Y,K219
NRTI_RAMS = {'M41L':'',
             'K65R':'',
             'D67N':'',
             'T69':'',
             'K70R':'','K70E':'',
             'L74A':'','L74C':'','L74D':'','L74E':'','L74F':'','L74G':'','L74H':'','L74I':'','L74K':'','L74M':'','L74N':'','L74P':'','L74Q':'','L74R':'','L74S':'','L74T':'','L74V':'','L74W':'','L74Y':'',
             'V75A':'','V75M':'','V75S':'','V75T':'',
             'L210W':'',
             'Y115F':'',
             'Q151M':'',
             'M184A':'','M184C':'','M184D':'','M184E':'','M184F':'','M184G':'','M184H':'','M184I':'','M184K':'','M184L':'','M184N':'','M184P':'','M184Q':'','M184R':'','M184S':'','M184T':'','M184V':'','M184W':'','M184Y':'',
             'T215F':'','T215Y':'',
             'K219A':'','K219C':'','K219D':'','K219E':'','K219F':'','K219G':'','K219H':'','K219I':'','K219L':'','K219M':'','K219N':'','K219P':'','K219Q':'','K219R':'','K219S':'','K219T':'','K219V':'','K219W':'','K219Y':'',
             }
             
# NNRTI RAMS: A98G, L100I, K101E, K101P, K103N, K103S, V106A, V106M, Y181, Y188, G190, P225, F227, M230L, P236L
NNRTI_RAMS = {'A98G':'', 
              'L100I':'', 
              'K101E':'','K101P':'', 
              'K103N':'','K103S':'', 
              'V106A':'','V106M':'', 
              'Y181A':'','Y181C':'','Y181D':'','Y181E':'','Y181F':'','Y181G':'','Y181H':'','Y181I':'','Y181K':'','Y181L':'','Y181M':'','Y181N':'','Y181P':'','Y181Q':'','Y181R':'','Y181S':'','Y181T':'','Y181V':'','Y181W':'',
              'Y188A':'','Y188C':'','Y188D':'','Y188E':'','Y188F':'','Y188G':'','Y188H':'','Y188I':'','Y188K':'','Y188L':'','Y188M':'','Y188N':'','Y188P':'','Y188Q':'','Y188R':'','Y188S':'','Y188T':'','Y188V':'','Y188W':'',
              'G190A':'','G190C':'','G190D':'','G190E':'','G190F':'','G190H':'','G190I':'','G190K':'','G190L':'','G190M':'','G190N':'','G190P':'','G190Q':'','G190R':'','G190S':'','G190T':'','G190V':'','G190W':'','G190Y':'',
              'P225A':'','P225C':'','P225D':'','P225E':'','P225F':'','P225G':'','P225H':'','P225I':'','P225K':'','P225L':'','P225M':'','P225N':'','P225Q':'','P225R':'','P225S':'','P225T':'','P225V':'','P225W':'','P225Y':'',
              'F227A':'','F227C':'','F227D':'','F227E':'','F227G':'','F227H':'','F227I':'','F227K':'','F227L':'','F227M':'','F227N':'','F227P':'','F227Q':'','F227R':'','F227S':'','F227T':'','F227V':'','F227W':'','F227Y':'',
              'M230L':'', 
              'P236L':''
              }
              
# PI RAMS: D30,V32,M46,I47,G48,I50,I54,V82A,V82F,V82S,V82T,V82C,V82G,V82L,V82M,I84,N88,L90
PI_RAMS = {'D30A':'','D30C':'','D30E':'','D30F':'','D30G':'','D30H':'','D30I':'','D30K':'','D30L':'','D30M':'','D30N':'','D30P':'','D30Q':'','D30R':'','D30S':'','D30T':'','D30V':'','D30W':'','D30Y':'',
           'V32A':'','V32C':'','V32D':'','V32E':'','V32F':'','V32G':'','V32H':'','V32I':'','V32K':'','V32L':'','V32M':'','V32N':'','V32P':'','V32Q':'','V32R':'','V32S':'','V32T':'','V32W':'','V32Y':'',
           'M46A':'','M46C':'','M46D':'','M46E':'','M46F':'','M46G':'','M46H':'','M46I':'','M46K':'','M46L':'','M46N':'','M46P':'','M46Q':'','M46R':'','M46S':'','M46T':'','M46V':'','M46W':'','M46Y':'',
           'I47A':'','I47C':'','I47D':'','I47E':'','I47F':'','I47G':'','I47H':'','I47K':'','I47L':'','I47M':'','I47N':'','I47P':'','I47Q':'','I47R':'','I47S':'','I47T':'','I47V':'','I47W':'','I47Y':'',
           'G48A':'','G48C':'','G48D':'','G48E':'','G48F':'','G48H':'','G48I':'','G48K':'','G48L':'','G48M':'','G48N':'','G48P':'','G48Q':'','G48R':'','G48S':'','G48T':'','G48V':'','G48W':'','G48Y':'',
           'I50A':'','I50C':'','I50D':'','I50E':'','I50F':'','I50G':'','I50H':'','I50K':'','I50L':'','I50M':'','I50N':'','I50P':'','I50Q':'','I50R':'','I50S':'','I50T':'','I50V':'','I50W':'','I50Y':'',
           'I54A':'','I54C':'','I54D':'','I54E':'','I54F':'','I54G':'','I54H':'','I54K':'','I54L':'','I54M':'','I54N':'','I54P':'','I54Q':'','I54R':'','I54S':'','I54T':'','I54V':'','I54W':'','I54Y':'',
           'V82A':'','V82F':'','V82S':'','V82T':'','V82C':'','V82G':'','V82L':'','V82M':'',
           'I84A':'','I84C':'','I84D':'','I84E':'','I84F':'','I84G':'','I84H':'','I84K':'','I84L':'','I84M':'','I84N':'','I84P':'','I84Q':'','I84R':'','I84S':'','I84T':'','I84V':'','I84W':'','I84Y':'',
           'N88A':'','N88C':'','N88D':'','N88E':'','N88F':'','N88G':'','N88H':'','N88I':'','N88K':'','N88L':'','N88M':'','N88P':'','N88Q':'','N88R':'','N88S':'','N88T':'','N88V':'','N88W':'','N88Y':'',
           'L90A':'','L90C':'','L90D':'','L90E':'','L90F':'','L90G':'','L90H':'','L90I':'','L90K':'','L90M':'','L90N':'','L90P':'','L90Q':'','L90R':'','L90S':'','L90T':'','L90V':'','L90W':'','L90Y':'',
           }

# INT RAMS: Y143R,Y143H,Y143C,Q148H,Q148K,Q148R,N155H
INT_RAMS = {'E92Q':'',
            'F121C':'','F121N':'','F121Y':'',
            'Y143R':'','Y143H':'','Y143C':'',
            'Q148H':'','Q148K':'','Q148R':'',
            'N155H':'','N155S':''
            }

# RT2-specific RAMS only, not including overlap with RT
RT2_RAMS = {'N348I':'',
            'T369I':'','T369V':''
            }


#######################
# Connect to Postgres #
#######################
def connectPG(LOG):
   # Define connection string
    #con = "host='gamera.virologic.com' dbname='mgrm' user='mark' password='mark'"
    con = "host='localhost' dbname='mgrm' user='mark' password='mark'"
   
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


#####################
# Connect to Oracle #
#####################
def connectORA(LOG):
   
   try:
      ip = '10.196.1.140'
      port = 1521
      SID='vrl1'
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


##############################
# Connect to MySQL Databases #
##############################
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
        con2 = ms.connect(host='kongdb.virologic.com',user='root',passwd='olympics',db='GENORULES') # Labcorp samples
        kong_db = con2.cursor()
        LOG.write(gt()+"\tSuccessful connection made to KongDB\n")
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        LOG.write(gt()+"\tKongDB Database connection failed!\n ->%s" % (exceptionValue))
        #sys.exit(gt()+"\tKongDB Database connection failed!\n ->%s" % (exceptionValue))
        return 'failed', '', ''
    return 'ok', kong, kong_db


##############################################
# Execute any Oracle query and return result #
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


##############################################
# execute Postgres SQL, do not return result #
##############################################
def execPG(pgc,pgc_sql,vals,LOG):
    try:
        pgc.execute(pgc_sql,vals)
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        LOG.write(gt()+"\tCould not execute Postgres sql!\nsql="+pgc_sql+"\nvals="+str(vals)+"\n ->%s" % (exceptionValue))
        return "Error"
    return "ok"


################################################
# Execute any Postgres query and return result #
################################################
def getPG(pgc,pgc_sql,vals,LOG):
    try:
        pgc.execute(pgc_sql,vals)
        pgc_result = pgc.fetchall()
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        LOG.write(gt()+"\tCould not execute Postgres sql!\n"+pgc_sql+"\n ->%s" % (exceptionValue))
        return "Error"
    
    return pgc_result


##############################################
# execute Postgres SQL, do not return result #
##############################################
def execPGCopyExpert(pgc,fileObj,sql,tablename,LOG):
    if sql =='': sql = "copy "+tablename+" from STDIN with null ''"
    try:
        fileObj.seek(0)
        pgc.copy_expert(sql,fileObj)
        LOG.write(gt()+"\t\tWrote data to "+tablename+"\n")
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        LOG.write(gt()+"\tCould not execute Postgres copy_expert for tablename="+tablename+"!\n ->%s" % (exceptionValue))
        return "Error"
    return "ok"


###############################################
# Execute any MySQL query and return a result #
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


####################################
# Send mail upon failure to update #
####################################
def sendMail(log_pth):
    import smtplib
    from email.mime.text import MIMEText
    fp = open(log_pth, 'rb')
    msg = MIMEText(fp.read())   # read logfile into message body
    fp.close()
    #me = 'gamera_server@monogrambio.com'
    me = 'bali_server@monogrambio.com'
    #you = 'bioinformatics@monogrambio.com'
    you = 'mevans@monogrambio.com'
    #msg['Subject'] = 'HIVdb Warehouse update failed last night'
    msg['Subject'] = 'Bali HIVdb Warehouse update failed last night'
    msg['From'] = me
    msg['To'] = you
    # Send the message via our own SMTP server, but don't include the envelope header.
    #s = smtplib.SMTP('localhost')
    s = smtplib.SMTP('autofax.virologic.com')
    s.sendmail(me, [you], msg.as_string())
    s.quit()
    return


##############################
# Check status of proceedure #
##############################
def checkStatus(status):
    global notify
    if status =='bad': notify = 'yes'
    return


########################
# Returns current time #
########################
def gt():
    return strftime("%d %b %Y %H:%M:%S",localtime())


###################################################
# Remove newline. p= pattern, x= string to search #
###################################################
def rn(p,x):
    m = p.search(x)
    if m: x = x.replace('\\','')
    return x


#########################################################################
# Catch None type, convert to empty string to prevent str concat errors #
#########################################################################
def catchNone(x):
    if x is None: return ""
    else: return x


###########################################################################################
#######                           Oracle Export SQL                               #########
###########################################################################################

#-- Export data for drugs table
def exportDrugsData():
    sql = """SELECT
    d.id AS drug_id,
    d.abrv_name AS drug_abrv,
    d.generic_name,
    d.full_name,
    d.def_brand_name,
    dc.name AS drug_class,
    d.max_fold_change,
    d.dec,
    d.ref_range_div AS hyper_susc_cutoff,
    d.ref_range_mult AS lower_cutoff,
    d.clinical_data_flg AS clin_lower_cut,
    d.upper_cutoff,
    d.boosted_drug_flg
    FROM 
    vlps.drugs d, 
    vlps.drug_classes dc
    WHERE d.drug_class_id=dc.id """ # and rownum < 51

    return sql


#-- Export data for projects table
def exportProjectData():
    sql ="""SELECT
    p.id AS project_id,
    p.code AS project_code,
    p.project_name,
    p.project_type

    FROM 
    reportmanager.projects p"""

    return sql


#-- Export data for mutation table
#-- 07.18.11 added creation date to allow updates
def exportMutationsData():
    sql ="""SELECT
    m.id AS mutation_id,
    (CASE m.type 
          WHEN 1 THEN 'RT'
          WHEN 2 THEN 'Prt'
          WHEN 3 THEN 'Int'
          WHEN 4 THEN 'RT2'
          ELSE TO_CHAR(m.type)
    END) AS "mutation_type",
    m.abrv_name AS mutation,
    m.position,
    m.orig_protein AS wt,
    m.new_protein AS mut,
    m.insertion_flg,
    TO_CHAR(m.CREATED_DATE,'DD-MON-YYYY HH24:MI:SS')

    FROM 
    genotype.mutations m"""
    
    return sql


#-- Export data for geno_rep_rows table
def exportGenoRepRowData():
    sql ="""SELECT
    r.id AS report_id,
    g.id AS gt_report_id,
    d.id AS drug_id,
    r.result, 
    r.reported_drug_flg,
    r.major_mutation_sum,
    r.resistant_flg,
    d.abrv_name AS drug_abrv

    FROM  
    genotype.gt_reports g, 
    genotype.gt_rep_rows r, 
    vlps.drugs d

    WHERE r.gt_report_id  = g.id
    AND   g.status_code = 1
    AND   g.rejected_report_flg =0
    AND   r.drug_abrv_name = d.abrv_name"""

    return sql


#-- Export data for geno_report table
#-- 07.18.11 added timestamp to date to make update easier
def exportGenoReportData():
    sql ="""SELECT
    id AS gt_report_id,
    vrl_accession,
    aliquot_id,
    subtype,
    rule_version_used AS rule_ver,
    TO_CHAR(rep_date,'DD-MON-YYYY HH24:MI:SS') AS rep_date,
    test_code,
    assay_type

    FROM  gt_reports g
    WHERE g.status_code = 1
    AND   g.rejected_report_flg =0"""

    return sql


#-- Export data for pheno_rep_row table
def exportPhenoRepRowData():
    sql ="""SELECT
    r.id AS rpt_row_id,
    r.pt_report_id AS report_id,
    r.ic50,
    r.ic95,
    r.fold_res,
    CASE 
        WHEN r.fold_change IS NULL THEN
            CASE 
                WHEN r.ic50 LIKE '>%' THEN r.max_fold_change 
                WHEN r.ic50 LIKE '<%' THEN '0.3'
                ELSE r.ic50
            END
        ELSE r.fold_change 
    END AS fold_change,
    r.resistant_flg,   --binary call: 0, 1, 2=sens; 3, 4=reduced susc
    r.result_flg,      --3 category call: 0=sens; 1=resistant; 2=partially sens
    r.reported_drug_flg,
    r.drug_class,
    r.abrv_name as drug_abrv

    FROM  
    vlps.pt_reports p, 
    vlps.pt_rep_rows r

    WHERE r.pt_report_id = p.id
    AND   p.rejected_report_flg = 0
    AND   p.status_code = 1 and rownum < 1001
    --AND   g.project in ('00073','00094','00203','00218') --commercial PJs only
    """

    return sql


#-- Export data for pheno_report table
#-- 07.18.11 added timestamp to date to make update easier
def exportPhenoReportData():
    sql ="""SELECT
    p.id as pt_report_id,
    p.vrl_accession,
    p.aliquot_id,
    p.test_code,
    p.assay_type,
    c.rc_status,
    c.reported_rc,
    TO_CHAR(p.rep_date,'YYYY-MM-DD HH24:MI:SS'),
    p.project as project_code

    FROM 
    vlps.pt_reports p, vlps.rc c

    WHERE   
    p.status_code = 1
    AND p.rejected_report_flg =0
    --commercial PJs only
    --AND   p.project in ('00073','00094','00203','00218')
    AND   c.ID (+) = p.RC_ID and rownum < 1001"""

    return sql


#-- Export data for PSGT_Result
#-- Must be queried from PSGT tablespace
#-- 07.18.11 formatted rep_date
def exportPSGTResultData():
    sql ="""SELECT
    p.pg_report_id,
    p.accession_id AS vrl_accession,
    pt.id AS pt_report_id,
    g.id AS gt_report_id,
    p.test_code,
    p.project_code,
    pr.drug_abrv_name AS drug_abrv,
    pr.resistant_flg,
    pr.result_flg,
    TO_CHAR(p.rep_date,'DD-MON-YYYY HH24:MI:SS')

    FROM 
    psgt.pg_reports p,
    psgt.pg_rep_rows pr,
    vlps.pt_reports pt,
    genotype.gt_reports g

    WHERE
    pr.pg_report_id=p.pg_report_id
    AND current_report_flg=1 
    AND rep_status IN ('FINAL','CORRECTED','AMENDED')
    AND p.accession_id = pt.vrl_accession
    AND p.accession_id = g.vrl_accession and rownum <1001"""

    return sql


#-- Export data for seq_mutations table
def exportSeqMutationsData():
    sql ="""SELECT
    g.id AS gt_report_id,
    g.data_file_id AS seq_id,
    u.id AS mutation_id,
    (CASE u.type 
        WHEN 1 THEN 'RT'
        WHEN 2 THEN 'Prt'
        WHEN 3 THEN 'Int'
        WHEN 4 THEN 'RT2'
        ELSE TO_CHAR(u.type)
    END) AS "mutation_type",
    u.abrv_name AS mutation,
    g.rep_date

    FROM 
    genotype.gt_reports g, 
    genotype.data_file_mutation m, 
    genotype.mutations u

    WHERE 
    g.data_file_id = m.data_file_id
    AND m.mutation_id = u.id
    AND g.status_code = 1
    AND g.rejected_report_flg =0 and rownum < 1001"""

    return sql


#- Update for etag_sum table. Connect via etag schema to battle
def exportEtagData():
    sql ="""SELECT distinct
    h2t.accession_id,
    to_char(h2t.patient_dob,'YYYY-MM-DD') patient_dob,
    to_char(h2t.collected_date,'YYYY-MM-DD') collected_date,
    to_char(h2t.rcv_date,'YYYY-MM-DD') rcv_date,
    to_char(h2t.rep_date,'YYYY-MM-DD') rep_date,
    (case when h2t.result_strvalue = '<Min' then '-9999' when h2t.result_strvalue = '>Max' then '99999' else h2t.result_strvalue end) H2T,
    (case when h2d.result_strvalue = '<Min' then '-9999' when h2d.result_strvalue = '>Max' then '99999' else h2d.result_strvalue end) H2D,
    h2t.H2T_Res,
    h2d.H2D_Res,
    h2t.area H2T_TA,
    h2d.area H2D_TA,
    h2t.sectionsize H2T_SS,
    h2d.sectionsize H2D_SS,
    h2t.percent_tumor H2T_P_Tumor,
    h2d.percent_tumor H2D_P_Tumor,
    h2t.batch_id H2T_batch,
    h2d.batch_id H2D_batch,
    h2t.client_state,
    h2t.client_country
    FROM
    (select 
        r.accession_id,
        r.patient_lname,
        r.patient_id,
        r.patient_dob,
        r.collected_date,
        r.rcv_date,
        r.rep_date,
        r.client_state,
        r.client_country,
        t.area,
        t.sectionsize,
        t.percent_tumor,
        rw.*, 
        rwt.result_category H2T_Res
        from 
        etag.tumor_measurement t,
        etag.worklist_rows w,
        etag.etag_rep_rows rw,
        etag.etag_reports r, 
        etag.etag_result_categories rwt
        where 
        rw.target_name='H2T' 
        and t.aliquot_id=w.tm_best_aliquot 
        and rw.aliquot_id=w.aliquot_id
        and r.id=rw.etag_report_id 
        and project_code='02233' 
        and status_code=1 
        and rw.result_category_id=rwt.id) h2t,
        
    (select 
        r.accession_id,
        t.area,
        t.sectionsize,
        t.percent_tumor,
        rw.*, 
        rwd.result_category H2D_Res
        from 
        etag.tumor_measurement t,
        etag.worklist_rows w,
        etag.etag_rep_rows rw,
        etag.etag_reports r, 
        etag.etag_result_categories rwd
        where 
        rw.target_name='H2D' 
        and t.aliquot_id=w.tm_best_aliquot 
        and rw.aliquot_id=w.aliquot_id
        and r.id=rw.etag_report_id 
        and project_code='02233' 
        and status_code=1 
        and rw.result_category_id=rwd.id) h2d
        
    WHERE
    h2t.accession_id = h2d.accession_id
    and to_char(h2t.rep_date,'YYYY-MM-DD') > '1969-07-20'
    order by collected_date, rep_date""" 

    return sql


#-- Export data for hivwh.Patient, but still need to create the pt_hash
#-- Export from Genotype schema
#-- 07.18.11 added order by clause, formatted dates
def exportPatientData():
    sql ="""SELECT
    r.vl_accession_id as vrl_accession,
    p.lastname as pt_lastname,
    p.firstname as pt_firstname,
    p.middlename as pt_middlename,
    p.sex as sex,
    p.dob as DOB,
    p.postal_code as pt_ZIP,
    p.province as pt_state,
    O.CODE as client_code,
    O.NAME as client_name,
    A.COUNTRY_NAME as client_country,
    A.CITY as client_city,
    A.PROVINCE_NAME as client_state,
    A.POST_CODE as client_zip,
    r.project_code,
    r.doctor_code,
    r.client_notes,
    r.visit_number,
    r.site_number,
    r.patient_viroload as viroload,
    r.account_class,
    r.collected_date,
    r.received_date,
    pn.vl_test_code,
    pn.ul_panel_code,
    td.product_name,
    td.product_type,
    pn.vl_reported_drugs
    FROM REPORTMANAGER.ORGANIZATION_VIEW O, 
         REPORTMANAGER.CUSTOMER_ADDRESSES_VIEW A,
        BRIDGEHEAD.patients p, 
        BRIDGEHEAD.requests r,
        BRIDGEHEAD.panels pn,
        REPORTMANAGER.test_definitions td
    where A.ADDRESS_TYPE = 1001 
    and A.organization_id = o.id
    and o.code=r.client_code
    and p.ul_patient_id = r.ul_patient_id
    and r.vl_accession_id = pn.vl_accession_id
    and pn.vl_test_code = td.code
    and r.received_date > to_date('1969-07-20','YYYY-MM-DD HH24:MI:SS') order by r.received_date,r.vl_accession_id,pn.vl_test_code""" # and rownum < 1001
    
    return sql


#-- Export data for hivwh.sequence
#-- 07.18.11 updated query to add rep_date and project
def exportSequenceData():
    sql ="""SELECT
    d.data_file_id,
    g.id,
    g.vrl_Accession,
    g.aliquot_id,
    d.nucleotide_seq,
    g.pri_summary,
    g.rt_summary,
    g.in_summary,
    g.rt2_summary,
    g.assay_type,
    TO_CHAR(g.rep_date,'DD-MON-YYYY HH24:MI:SS'),
    g.project

    FROM genotype.gt_reports g, genotype.data_files_nucleotide_seq d
    
    WHERE d.data_file_id = g.data_file_id
    AND g.status_code =1
    AND g.rejected_report_flg = 0
    UNION
    SELECT
    g.data_file_id,
    g.id,
    g.vrl_accession,
    g.aliquot_id,
    null as nucleotide_seq,
    g.pri_summary,
    g.rt_summary,
    g.in_summary,
    g.rt2_summary,
    g.assay_type,
    TO_CHAR(g.rep_date,'DD-MON-YYYY HH24:MI:SS'),
    g.project
    
    FROM genotype.gt_reports g
    
    WHERE g.data_file_id < 25972
    AND g.status_code =1
    AND g.rejected_report_flg=0 and rownum < 1001"""
    
    return sql


def exportGTRuleResultRowsData():
    sql="""select 
    SAMPLE_ID,
    DRUG_ABRV_NAME,
    DRUG_CLASS,
    COMMENTS,
    RAMs,
    RESULT,
    RESISTANT_FLG,
    MUTATION_SCORE,
    CREATED_DATE,
    RULE_VERSION,
    ENGINE_VERSION,
    "%s" 
    
    FROM 
    GT_RULE_RESULT_ROWS 

    WHERE 
    SAMPLE_ID !='dummy id' 
    AND CREATED_DATE > '1969-07-20' """

    return sql



###########################################################################################
########                              Build SQL                                  ##########
###########################################################################################

def createDrugsTable():
    sql ="""DROP TABLE wh.drugs CASCADE;

    CREATE TABLE wh.drugs
    (
        drug_id           integer primary key,
        drug_abrv         text not null,
        generic_name      text,
        full_name         text,
        def_brand_name    text,
        drug_class        text,
        max_fold_change   numeric,
        dec               numeric,
        hyper_susc_cutoff numeric,
        lower_cutoff      numeric,
        clinical_data_flg numeric,
        upper_cutoff      numeric,
        boosted_drug_flg  numeric
    );
    
    CREATE INDEX drug_abrv_idx ON wh.drugs(drug_abrv);"""
    return sql


def createProjectTable():
    sql="""DROP TABLE wh.project CASCADE;

    CREATE TABLE wh.project
    (
        project_id integer,
        project_code text not null,
        project_name text,
        project_type numeric
    );
    
    CREATE INDEX project_proj_code_idx ON wh.project(project_code);"""
    return sql


def createMutationsTable():
    #-- 07.18.11  Added created_date to mutations table
    sql="""DROP TABLE wh.mutations CASCADE;

    CREATE TABLE wh.mutations
    (
        mutation_id    integer primary key,
        mutation_type  text,
        mutation       text,
        position       integer,
        wt             text,
        mut            text,
        insertion_flg  integer,
        created_date   text
    );
    
    CREATE INDEX mutation_mut_idx ON wh.mutations(mutation);"""    
    return sql


def createGenoRepRowTable():
    sql="""DROP TABLE wh.geno_rep_row CASCADE;

    CREATE TABLE wh.geno_rep_row
    (
        report_id         integer  primary key,
        gt_report_id      integer not null,
        drug_id           integer not null,
        result            numeric,
        reported_drug_flg numeric,
        major_mut_sum     text,
        resistant_flg     numeric,
        drug_abrv         text
    );
    
    CREATE INDEX geno_rep_row_idx ON wh.geno_rep_row(gt_report_id,report_id);"""
    return sql


def createGenoReportTable():
    sql="""DROP TABLE wh.geno_report CASCADE;

    CREATE TABLE wh.geno_report
    (
        gt_report_id   integer primary key,
        vrl_accession  text not null,
        aliquot_id     text not null,
        subtype        text,
        rule_ver       numeric,
        rep_date       text,
        test_code      text,
        assay_type     integer
    );
    
    CREATE INDEX geno_rpt_idx ON wh.geno_report(vrl_accession,aliquot_id,gt_report_id);"""
    return sql


def createPhenoRepRowTable():
    sql="""DROP TABLE wh.pheno_rep_row CASCADE;

    CREATE TABLE wh.pheno_rep_row
    (
        report_id         integer primary key,
        pt_report_id      integer not null,
        ic50              text,
        ic95              text,
        fold_res          text,
        fold_change       text,
        resistant_flg     integer,
        result_flg        integer,
        reported_drug_flg numeric,
        drug_class        text,
        drug_abrv         text
    );
    
    CREATE INDEX pheno_rep_row_idx ON wh.pheno_rep_row(pt_report_id,report_id,drug_abrv);"""
    return sql


def createPhenoReportTable():
    sql="""DROP TABLE wh.pheno_report CASCADE;

    CREATE TABLE wh.pheno_report
    (
        pt_report_id   integer primary key,
        vrl_accession  text not null,
        aliquot_id     text not null,
        test_code      text,
        assay_type     integer,
        rc_status      integer,
        reported_rc    numeric,
        rep_date       text,
        project_code   text
    );
    
    CREATE INDEX pheno_rpt_idx ON wh.pheno_report(vrl_accession,aliquot_id,pt_report_id);"""
    return sql


def createPSGTResultTable():
    sql="""DROP TABLE wh.psgt_result CASCADE;

    CREATE TABLE wh.psgt_result
    (
        id             serial primary key,
        pg_report_id   integer,
        vrl_accession  text,
        pt_report_id   numeric not null,
        gt_report_id   numeric not null,
        test_code      text,
        project_code   text,
        drug_abrv      text,
        resistant_flg  integer,
        result_flg     integer,
        rep_date       text
    );
    
    CREATE INDEX psgt_result_idx on wh.psgt_result(vrl_accession,pt_report_id,gt_report_id);"""
    return sql


def createSeqMutationsTable():
    sql="""DROP TABLE wh.seq_mutations CASCADE;

    CREATE TABLE wh.seq_mutations
    (
        id             serial primary key,
        gt_report_id   integer,
        seq_id         integer,
        mutation_id    integer,
        mutation_type  text,
        mutation       text,
        rep_date       text
    );
    
    CREATE INDEX seq_mutation_idx ON wh.seq_mutations(gt_report_id,seq_id,mutation_id,mutation);"""
    return sql


def createETAGSumTable():
    sql ="""DROP TABLE wh.etag_sum CASCADE;

    CREATE TABLE wh.etag_sum
    (
        eid                 serial primary key,
        vrl_accession       text,
        patient_dob         text,
        collected_date      text,
        rcv_date            text,
        rep_date            text,
        h2t                 numeric,
        h2d                 numeric,
        h2t_result          text,
        h2d_result          text,
        h2t_tumor_area      numeric,
        h2d_tumor_area      numeric,
        h2t_section_size    numeric,
        h2d_section_size    numeric,
        h2t_pct_tumor       numeric,
        h2d_pct_tumor       numeric,
        h2t_batch_id        numeric,
        h2d_batch_id        numeric,
        client_state        text,
        client_country      text
    );
    
    CREATE INDEX etag_sum_idx ON wh.etag_sum(vrl_accession,rcv_date,rep_date);"""
    return sql


def createPatientTable():
    sql="""DROP TABLE wh.patient CASCADE;

    CREATE TABLE wh.patient
    (
        
        vrl_accession  text not null,
        lastname       text,
        firstname      text,
        middlename     text,
        sex            text,
        dob            text,
        pt_zip         text,
        pt_state       text,
        client_code    text,
        client_name    text,
        client_country text,
        client_city    text,
        client_state   text,
        client_zip     text,
        project_code   text,
        doctor_code    text,
        client_notes   text,
        visit_num      integer,
        site_num       text,
        viroload       text,
        account_class  text,
        collected_date text,
        received_date  text,
        test_code      text,
        panel_code     text,
        product_name   text,
        product_type   text,
        reported_drugs text,
        pt_hash        text,
        patient_id     serial primary key,
        test_num       int4
    );
    
    CREATE INDEX patient_idx on wh.patient(vrl_accession,pt_hash,project_code);"""
    return sql


def createPatientLoadTable():
    sql="""CREATE TABLE wh.patient_load
    (
        vrl_accession text,
        lastname text,
        firstname text,
        middlename text,
        sex text,
        dob text,
        pt_ZIP text,
        pt_state text,
        client_code text,
        client_name text,
        client_country text,
        client_city text,
        client_state text,
        client_zip text,
        project_code text,
        doctor_code text,
        client_notes text,
        visit_num int4,
        site_num text,
        viroload text,
        account_class text,
        collected_date text,
        received_date text,
        test_code text,
        panel_code text,
        product_name text,
        product_type text,
        reported_drugs text,
        pt_hash text,
        pt_temp_id serial,
        test_num int4 
    );"""
    return sql


def createSeqTable():
    #-- 07.18.11  Added rep_date, project fields to seq table 
    sql="""DROP TABLE wh.seq CASCADE;

    CREATE TABLE wh.seq
    (
        seq_id         integer primary key,
        gt_report_id   integer,
        aliquot_id     text,
        vrl_accession  text,
        nucleotide_seq text,
        pri_sum        text,
        rt_sum         text,
        int_sum        text,
        rt2_sum        text,
        assay_type     integer,
        rep_date       text,
        project        text,
        wt_flg         text,
        tams_sum       text,
        nams_sum       text,
        nnrti_sum      text,
        pri_mut_sum    text,
        int_mut_sum    text,
        rt2_mut_sum    text,
        tam_score      numeric,
        nam_score      numeric,
        nnrti_score    numeric,
        pri_score      numeric,
        int_score      numeric,
        rt2_score      numeric
    );
    
    CREATE INDEX seq_idx ON wh.seq(seq_id,gt_report_id,vrl_accession,aliquot_id,pri_sum,rt_sum,int_sum,pri_mut_sum,tams_sum,nams_sum,nnrti_sum,int_mut_sum);"""
    return sql


def createMutFreqTable():
    sql="""DROP TABLE wh.mut_freq CASCADE;

    CREATE TABLE wh.mut_freq
    (
        id             serial primary key,
        mutation       text,
        mutation_sep   text,
        mut_freq       integer,
        mut_type       text
    );
    
    CREATE INDEX mut_freq_idx ON wh.mut_freq(mutation,mutation_sep);"""
    return sql


def createGTRuleResultRowsTable():
    sql=""" DROP TABLE wh.gt_rule_result_rows CASCADE;

    CREATE TABLE wh.gt_rule_result_rows
    (
        result_id      serial primary key,
        sample_id      text,
        drug_abrv      text,
        drug_class     text,
        comments       text,
        rams           text,
        result         numeric,
        resistant_flg  numeric,
        mutation_score numeric,
        created_date   text,
        rule_version   numeric,
        engine_version numeric,
        data_src       text
    );
    
    CREATE INDEX gt_rule_result_rows_idx ON (sample_id,drug_abrv,drug_class);"""
    return sql


def createGenoResultSumView():
    sql="""CREATE VIEW wh.geno_result_sum AS
    SELECT DISTINCT
    ab.project_code,
    a.*,
    seq_id,
    wt_flg,
    pri_sum,
    rt_sum,
    int_sum,
    rt2_sum,
    tams_sum,
    nams_sum,
    nnrt_sum,
    pri_mut_sum,
    int_mut_sum,
    rt2_mut_sum,
    tam_score,
    nam_score,
    nnrti_score,
    pri_score,
    int_score,
    rt2_score,

    geno_result_3TC,
    geno_result_ABC,
    geno_result_d4T,
    geno_result_ddI,
    geno_result_FTC,
    geno_result_TFV,
    geno_result_ZDV,
    geno_result_DLV,
    geno_result_EFV,
    geno_result_ETR,
    geno_result_NVP,
    geno_result_AMP,
    geno_result_AMPr,
    geno_result_ATV,
    geno_result_ATVr,
    geno_result_DRVr,
    geno_result_IDV,
    geno_result_IDVr,
    geno_result_LPVr,
    geno_result_NFV,
    geno_result_RTV,
    geno_result_SQV,
    geno_result_SQVr,
    geno_result_TPVr,
    geno_result_RPV,

    mmut_sum_3TC,
    mmut_sum_ABC,
    mmut_sum_d4T,
    mmut_sum_ddI,
    mmut_sum_FTC,
    mmut_sum_TFV,
    mmut_sum_ZDV,
    mmut_sum_DLV,
    mmut_sum_EFV,
    mmut_sum_ETR,
    mmut_sum_NVP,
    mmut_sum_AMP,
    mmut_sum_AMPr,
    mmut_sum_ATV,
    mmut_sum_ATVr,
    mmut_sum_DRVr,
    mmut_sum_IDV,
    mmut_sum_IDVr,
    mmut_sum_LPVr,
    mmut_sum_NFV,
    mmut_sum_RTV,
    mmut_sum_SQV,
    mmut_sum_SQVr,
    mmut_sum_TPVr,
    mmut_sum_RPV,

    rflg_3TC,
    rflg_ABC,
    rflg_d4T,
    rflg_ddI,
    rflg_FTC,
    rflg_TFV,
    rflg_ZDV,
    rflg_DLV,
    rflg_EFV,
    rflg_ETR,
    rflg_NVP,
    rflg_AMP,
    rflg_AMPr,
    rflg_ATV,
    rflg_ATVr,
    rflg_DRVr,
    rflg_IDV,
    rflg_IDVr,
    rflg_LPVr,
    rflg_NFV,
    rflg_RTV,
    rflg_SQV,
    rflg_SQVr,
    rflg_TPVr,
    rflg_RPV

    from 
    wh.geno_report a
    left outer join (select distinct vrl_accession, project_code from wh.patient)   as ab on (a.vrl_accession = ab.vrl_accession)
    left outer join (select distinct gt_report_id, seq_id, wt_flg, pri_sum,rt_sum,int_sum,rt2_sum,coalesce(tams_sum,'') as tams_sum, coalesce(nams_sum,'') as nams_sum, coalesce(nnrti_sum,'') as nnrt_sum, coalesce(pri_mut_sum,'') as pri_mut_sum, coalesce(int_mut_sum,'') as int_mut_sum,coalesce(rt2_mut_sum,'') as rt2_mut_sum, tam_score, nam_score, nnrti_score, pri_score,int_score,rt2_score from wh.seq ) as b on (a.gt_report_id = b.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_3TC , coalesce(major_mut_sum,'') as mmut_sum_3TC , resistant_flg as rflg_3TC  from wh.geno_rep_row where drug_abrv='3TC')     as c on (a.gt_report_id = c.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_ABC , coalesce(major_mut_sum,'') as mmut_sum_ABC , resistant_flg as rflg_ABC  from wh.geno_rep_row where drug_abrv='ABC')     as d on (a.gt_report_id = d.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_d4T , coalesce(major_mut_sum,'') as mmut_sum_d4T , resistant_flg as rflg_d4T  from wh.geno_rep_row where drug_abrv='d4T')     as e on (a.gt_report_id = e.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_ddI , coalesce(major_mut_sum,'') as mmut_sum_ddI , resistant_flg as rflg_ddI  from wh.geno_rep_row where drug_abrv='ddI')     as f on (a.gt_report_id = f.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_FTC , coalesce(major_mut_sum,'') as mmut_sum_FTC , resistant_flg as rflg_FTC  from wh.geno_rep_row where drug_abrv='FTC')     as g on (a.gt_report_id = g.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_TFV , coalesce(major_mut_sum,'') as mmut_sum_TFV , resistant_flg as rflg_TFV  from wh.geno_rep_row where drug_abrv='TFV')     as h on (a.gt_report_id = h.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_ZDV , coalesce(major_mut_sum,'') as mmut_sum_ZDV , resistant_flg as rflg_ZDV  from wh.geno_rep_row where drug_abrv='ZDV')     as i on (a.gt_report_id = i.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_DLV , coalesce(major_mut_sum,'') as mmut_sum_DLV , resistant_flg as rflg_DLV  from wh.geno_rep_row where drug_abrv='DLV')     as j on (a.gt_report_id = j.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_EFV , coalesce(major_mut_sum,'') as mmut_sum_EFV , resistant_flg as rflg_EFV  from wh.geno_rep_row where drug_abrv='EFV')     as k on (a.gt_report_id = k.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_ETR , coalesce(major_mut_sum,'') as mmut_sum_ETR , resistant_flg as rflg_ETR  from wh.geno_rep_row where drug_abrv='ETR')     as l on (a.gt_report_id = l.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_NVP , coalesce(major_mut_sum,'') as mmut_sum_NVP , resistant_flg as rflg_NVP  from wh.geno_rep_row where drug_abrv='NVP')     as m on (a.gt_report_id = m.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_AMP , coalesce(major_mut_sum,'') as mmut_sum_AMP , resistant_flg as rflg_AMP  from wh.geno_rep_row where drug_abrv='AMP')     as n on (a.gt_report_id = n.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_AMPr, coalesce(major_mut_sum,'') as mmut_sum_AMPr, resistant_flg as rflg_AMPr from wh.geno_rep_row where drug_abrv='AMP/r')   as o on (a.gt_report_id = o.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_ATV , coalesce(major_mut_sum,'') as mmut_sum_ATV , resistant_flg as rflg_ATV  from wh.geno_rep_row where drug_abrv='ATV')     as p on (a.gt_report_id = p.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_ATVr, coalesce(major_mut_sum,'') as mmut_sum_ATVr, resistant_flg as rflg_ATVr from wh.geno_rep_row where drug_abrv='ATV/r')   as q on (a.gt_report_id = q.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_DRVr, coalesce(major_mut_sum,'') as mmut_sum_DRVr, resistant_flg as rflg_DRVr from wh.geno_rep_row where drug_abrv='DRV/r')   as r on (a.gt_report_id = r.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_IDV , coalesce(major_mut_sum,'') as mmut_sum_IDV , resistant_flg as rflg_IDV  from wh.geno_rep_row where drug_abrv='IDV')     as s on (a.gt_report_id = s.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_IDVr, coalesce(major_mut_sum,'') as mmut_sum_IDVr, resistant_flg as rflg_IDVr from wh.geno_rep_row where drug_abrv='IDV/r')   as t on (a.gt_report_id = t.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_LPVr, coalesce(major_mut_sum,'') as mmut_sum_LPVr, resistant_flg as rflg_LPVr from wh.geno_rep_row where drug_abrv='LPV/r')   as u on (a.gt_report_id = u.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_NFV , coalesce(major_mut_sum,'') as mmut_sum_NFV , resistant_flg as rflg_NFV  from wh.geno_rep_row where drug_abrv='NFV')     as v on (a.gt_report_id = v.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_RTV , coalesce(major_mut_sum,'') as mmut_sum_RTV , resistant_flg as rflg_RTV  from wh.geno_rep_row where drug_abrv='RTV')     as w on (a.gt_report_id = w.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_SQV , coalesce(major_mut_sum,'') as mmut_sum_SQV , resistant_flg as rflg_SQV  from wh.geno_rep_row where drug_abrv='SQV')     as x on (a.gt_report_id = x.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_SQVr, coalesce(major_mut_sum,'') as mmut_sum_SQVr, resistant_flg as rflg_SQVr from wh.geno_rep_row where drug_abrv='SQV/r')   as y on (a.gt_report_id = y.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_TPVr, coalesce(major_mut_sum,'') as mmut_sum_TPVr, resistant_flg as rflg_TPVr from wh.geno_rep_row where drug_abrv='TPV/r')   as z on (a.gt_report_id = z.gt_report_id)
    left outer join (select distinct gt_report_id, result as geno_result_RPV , coalesce(major_mut_sum,'') as mmut_sum_RPV , resistant_flg as rflg_RPV  from wh.geno_rep_row where drug_abrv='RPV')     as aa on (a.gt_report_id = aa.gt_report_id)"""

    return sql


def createPhenoResultSumView():
    sql="""CREATE VIEW wh.pheno_result_sum2 AS
    SELECT DISTINCT
    a.pt_report_id,
    a.vrl_Accession,
    a.aliquot_id,
    a.test_code,
    a.assay_type,
    a.rc_status,
    a.reported_rc,
    a.rep_date,
    a.project_code,

    FC_3TC,
    FC_ABC,
    FC_d4T,
    FC_ddI,
    FC_FTC,
    FC_TFV,
    FC_ZDV,
    FC_DLV,
    FC_EFV,
    FC_ETR,
    FC_NVP,
    FC_AMP,
    FC_AMPr,
    FC_ATV,
    FC_ATVr,
    FC_DRVr,
    FC_IDV,
    FC_IDVr,
    FC_LPVr,
    FC_NFV,
    FC_RTV,
    FC_SQV,
    FC_SQVr,
    FC_TPVr,
    FC_RPV,

    pheno_result_3TC,
    pheno_result_ABC,
    pheno_result_d4T,
    pheno_result_ddI,
    pheno_result_FTC,
    pheno_result_TFV,
    pheno_result_ZDV,
    pheno_result_DLV,
    pheno_result_EFV,
    pheno_result_ETR,
    pheno_result_NVP,
    pheno_result_AMP,
    pheno_result_AMPr,
    pheno_result_ATV,
    pheno_result_ATVr,
    pheno_result_DRVr,
    pheno_result_IDV,
    pheno_result_IDVr,
    pheno_result_LPVr,
    pheno_result_NFV,
    pheno_result_RTV,
    pheno_result_SQV,
    pheno_result_SQVr,
    pheno_result_TPVr,
    pheno_result_RPV,

    IC50_3TC,
    IC50_ABC,
    IC50_d4T,
    IC50_ddI,
    IC50_FTC,
    IC50_TFV,
    IC50_ZDV,
    IC50_DLV,
    IC50_EFV,
    IC50_ETR,
    IC50_NVP,
    IC50_AMP,
    IC50_AMPr,
    IC50_ATV,
    IC50_ATVr,
    IC50_DRVr,
    IC50_IDV,
    IC50_IDVr,
    IC50_LPVr,
    IC50_NFV,
    IC50_RTV,
    IC50_SQV,
    IC50_SQVr,
    IC50_TPVr,
    IC50_RPV,

    rflg_3TC,
    rflg_ABC,
    rflg_d4T,
    rflg_ddI,
    rflg_FTC,
    rflg_TFV,
    rflg_ZDV,
    rflg_DLV,
    rflg_EFV,
    rflg_ETR,
    rflg_NVP,
    rflg_AMP,
    rflg_AMPr,
    rflg_ATV,
    rflg_ATVr,
    rflg_DRVr,
    rflg_IDV,
    rflg_IDVr,
    rflg_LPVr,
    rflg_NFV,
    rflg_RTV,
    rflg_SQV,
    rflg_SQVr,
    rflg_TPVr,
    rflg_RPV

    from 
    wh.pheno_report a 
    left outer join (select distinct pt_report_id, fold_change as FC_3TC,  result_flg as pheno_result_3TC , coalesce(ic50,'') as IC50_3TC , resistant_flg as rflg_3TC  from wh.pheno_rep_row where drug_abrv='3TC')   as b on (a.pt_report_id=b.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_ABC,  result_flg as pheno_result_ABC , coalesce(ic50,'') as IC50_ABC , resistant_flg as rflg_ABC  from wh.pheno_rep_row where drug_abrv='ABC')   as c on (a.pt_report_id=c.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_d4T,  result_flg as pheno_result_d4T , coalesce(ic50,'') as IC50_d4T , resistant_flg as rflg_d4T  from wh.pheno_rep_row where drug_abrv='d4T')   as d on (a.pt_report_id=d.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_ddI,  result_flg as pheno_result_ddI , coalesce(ic50,'') as IC50_ddI , resistant_flg as rflg_ddI  from wh.pheno_rep_row where drug_abrv='ddI')   as e on (a.pt_report_id=e.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_FTC,  result_flg as pheno_result_FTC , coalesce(ic50,'') as IC50_FTC , resistant_flg as rflg_FTC  from wh.pheno_rep_row where drug_abrv='FTC')   as f on (a.pt_report_id=f.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_TFV,  result_flg as pheno_result_TFV , coalesce(ic50,'') as IC50_TFV , resistant_flg as rflg_TFV  from wh.pheno_rep_row where drug_abrv='TFV')   as g on (a.pt_report_id=g.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_ZDV,  result_flg as pheno_result_ZDV , coalesce(ic50,'') as IC50_ZDV , resistant_flg as rflg_ZDV  from wh.pheno_rep_row where drug_abrv='ZDV')   as h on (a.pt_report_id=h.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_DLV,  result_flg as pheno_result_DLV , coalesce(ic50,'') as IC50_DLV , resistant_flg as rflg_DLV  from wh.pheno_rep_row where drug_abrv='DLV')   as i on (a.pt_report_id=i.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_EFV,  result_flg as pheno_result_EFV , coalesce(ic50,'') as IC50_EFV , resistant_flg as rflg_EFV  from wh.pheno_rep_row where drug_abrv='EFV')   as j on (a.pt_report_id=j.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_ETR,  result_flg as pheno_result_ETR , coalesce(ic50,'') as IC50_ETR , resistant_flg as rflg_ETR  from wh.pheno_rep_row where drug_abrv='ETR')   as k on (a.pt_report_id=k.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_NVP,  result_flg as pheno_result_NVP , coalesce(ic50,'') as IC50_NVP , resistant_flg as rflg_NVP  from wh.pheno_rep_row where drug_abrv='NVP')   as l on (a.pt_report_id=l.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_AMP,  result_flg as pheno_result_AMP , coalesce(ic50,'') as IC50_AMP , resistant_flg as rflg_AMP  from wh.pheno_rep_row where drug_abrv='AMP')   as m on (a.pt_report_id=m.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_AMPr, result_flg as pheno_result_AMPr, coalesce(ic50,'') as IC50_AMPr, resistant_flg as rflg_AMPr from wh.pheno_rep_row where drug_abrv='AMP/r') as n on (a.pt_report_id=n.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_ATV,  result_flg as pheno_result_ATV , coalesce(ic50,'') as IC50_ATV , resistant_flg as rflg_ATV  from wh.pheno_rep_row where drug_abrv='ATV')   as o on (a.pt_report_id=o.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_ATVr, result_flg as pheno_result_ATVr, coalesce(ic50,'') as IC50_ATVr, resistant_flg as rflg_ATVr from wh.pheno_rep_row where drug_abrv='ATV/r') as p on (a.pt_report_id=p.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_DRVr, result_flg as pheno_result_DRVr, coalesce(ic50,'') as IC50_DRVr, resistant_flg as rflg_DRVr from wh.pheno_rep_row where drug_abrv='DRV/r') as q on (a.pt_report_id=q.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_IDV,  result_flg as pheno_result_IDV , coalesce(ic50,'') as IC50_IDV , resistant_flg as rflg_IDV  from wh.pheno_rep_row where drug_abrv='IDV')   as r on (a.pt_report_id=r.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_IDVr, result_flg as pheno_result_IDVr, coalesce(ic50,'') as IC50_IDVr, resistant_flg as rflg_IDVr from wh.pheno_rep_row where drug_abrv='IDV/r') as s on (a.pt_report_id=s.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_LPVr, result_flg as pheno_result_LPVr, coalesce(ic50,'') as IC50_LPVr, resistant_flg as rflg_LPVr from wh.pheno_rep_row where drug_abrv='LPV/r') as t on (a.pt_report_id=t.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_NFV,  result_flg as pheno_result_NFV , coalesce(ic50,'') as IC50_NFV , resistant_flg as rflg_NFV  from wh.pheno_rep_row where drug_abrv='NFV')   as u on (a.pt_report_id=u.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_RTV,  result_flg as pheno_result_RTV , coalesce(ic50,'') as IC50_RTV , resistant_flg as rflg_RTV  from wh.pheno_rep_row where drug_abrv='RTV')   as v on (a.pt_report_id=v.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_SQV,  result_flg as pheno_result_SQV , coalesce(ic50,'') as IC50_SQV , resistant_flg as rflg_SQV  from wh.pheno_rep_row where drug_abrv='SQV')   as w on (a.pt_report_id=w.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_SQVr, result_flg as pheno_result_SQVr, coalesce(ic50,'') as IC50_SQVr, resistant_flg as rflg_SQVr from wh.pheno_rep_row where drug_abrv='SQV/r') as x on (a.pt_report_id=x.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_TPVr, result_flg as pheno_result_TPVr, coalesce(ic50,'') as IC50_TPVr, resistant_flg as rflg_TPVr from wh.pheno_rep_row where drug_abrv='TPV/r') as y on (a.pt_report_id=y.pt_report_id)
    left outer join (select distinct pt_report_id, fold_change as FC_RPV,  result_flg as pheno_result_RPV , coalesce(ic50,'') as IC50_RPV , resistant_flg as rflg_RPV  from wh.pheno_rep_row where drug_abrv='RPV')   as z on (a.pt_report_id=z.pt_report_id)"""

    return sql


def createIntegrasePhenoResultSumView():
    sql="""CREATE VIEW wh.integrase_pheno_result_sum AS 
    SELECT DISTINCT 
    a.pt_report_id, 
    a.vrl_accession, 
    a.aliquot_id, 
    a.test_code, 
    a.assay_type, 
    a.rc_status, 
    a.reported_rc, 
    a.rep_date, 
    a.project_code, 
    b.fc_ral, 
    c.fc_evg, 
    d.fc_572, 
    e.fc_744, 
    b.pheno_result_ral, 
    c.pheno_result_evg, 
    d.pheno_result_572, 
    e.pheno_result_744, 
    b.ic50_ral, 
    c.ic50_evg, 
    d.ic50_572, 
    e.ic50_744, 
    b.rflg_ral, 
    c.rflg_evg, 
    d.rflg_572, 
    e.rflg_744 

    FROM 
    (
        (
            (
                (pheno_report a LEFT JOIN 
                    (SELECT DISTINCT 
                        pheno_rep_row.pt_report_id, 
                        pheno_rep_row.fold_change AS fc_ral, 
                        pheno_rep_row.result_flg AS pheno_result_ral, 
                        COALESCE(pheno_rep_row.ic50, ''::text) AS ic50_ral, 
                        pheno_rep_row.resistant_flg AS rflg_ral 
                        FROM 
                        pheno_rep_row 
                        WHERE 
                        (pheno_rep_row.drug_abrv = 'RAL'::text)
                    ) b ON 
                    (
                        (a.pt_report_id = b.pt_report_id)
                    )
                ) LEFT JOIN 
                    (SELECT DISTINCT 
                        pheno_rep_row.pt_report_id, 
                        pheno_rep_row.fold_change AS fc_evg, 
                        pheno_rep_row.result_flg AS pheno_result_evg, 
                        COALESCE(pheno_rep_row.ic50, ''::text) AS ic50_evg, 
                        pheno_rep_row.resistant_flg AS rflg_evg 
                        FROM 
                        pheno_rep_row 
                        WHERE 
                        (pheno_rep_row.drug_abrv = 'EVG'::text)
                    ) c ON 
                    (
                        (a.pt_report_id = c.pt_report_id)
                    )
                ) LEFT JOIN 
                    (SELECT DISTINCT 
                        pheno_rep_row.pt_report_id, 
                        pheno_rep_row.fold_change AS fc_572, 
                        pheno_rep_row.result_flg AS pheno_result_572, 
                        COALESCE(pheno_rep_row.ic50, ''::text) AS ic50_572, 
                        pheno_rep_row.resistant_flg AS rflg_572 
                        FROM 
                        pheno_rep_row 
                        WHERE 
                        (pheno_rep_row.drug_abrv = '572'::text)
                    ) d ON 
                    (
                        (a.pt_report_id = d.pt_report_id)
                    )
                ) LEFT JOIN 
                    (SELECT DISTINCT 
                        pheno_rep_row.pt_report_id, 
                        pheno_rep_row.fold_change AS fc_744, 
                        pheno_rep_row.result_flg AS pheno_result_744, 
                        COALESCE(pheno_rep_row.ic50, ''::text) AS ic50_744, 
                        pheno_rep_row.resistant_flg AS rflg_744 
                        FROM 
                        pheno_rep_row 
                        WHERE 
                        (pheno_rep_row.drug_abrv = '744'::text)
                    ) e ON 
                    (
                        (a.pt_report_id = e.pt_report_id)
                    )
                ) WHERE (a.assay_type = 2)"""
    return sql




# General function to add data from Oracle query to Postgres table
###################################################################
def addData(cxc,pgc,pg_con,tablename,sql,insert_sql,LOG):
    LOG.write(gt()+"Beginning to fill "+tablename+"\n")
    f = StringIO()
    csv_writer = csv.writer(f, delimiter='\t')
    data = getOracle(cxc,sql,'',LOG)
    if data != 'Error':
        LOG.write(gt()+"\t\tRetrieved data from Oracle...\n")
        csv_writer.writerows(data)
        execPGCopyExpert(pgc,f,insert_sql,tablename,LOG)
    else: 
        LOG.write(gt()+"\tFailed to get query results from Oracle, "+tablename+" not updated\n")
        return 'bad'
    f.close()
    pg_con.commit()
    LOG.write(gt()+"\t\t"+str(len(data))+" records were added to "+tablename+"\n")
    return 'ok'


# Load patient data into table, calculate unique patient ids and assign test numbers
#####################################################################################
def loadPatientData(pgc,pg_con,cxc,LOG):
    LOG.write(gt()+"\tBeginning to process data for wh.patient....\n")

    # Get all patient data from Oracle
    ###################################
    cxr = getOracle(cxc,exportPatientData(),'',LOG)

    if cxr != "Error":

        # Verify there are records to insert
        if len(cxr) != 0:

            # Drop and recreate patient_load table
            # Using patient_load is not strictly necessary, but protects patient table in case there is a problem with load
            #################################################################################################################
            execPG(pgc,createPatientLoadTable(),'',LOG)
            pg_con.commit()
            LOG.write(gt()+"\t\tCreated patient_load table\n")
            
            # Insert new records in wh.patient_load
            #####################################
            LOG.write(gt()+"\t\tWriting update to wh.patient_load table...\n")
            
            p1 = re.compile(r'\\') # Need to search fields for hidden carriage returns that would break the load, so make patterns here
            counter = 0
            for vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs in cxr:
                sql = "insert into wh.patient_load(vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                client_name = rn(p1,client_name)
                vals = (vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs)
                execPG(pgc,sql,vals,LOG)
                counter += 1
                if counter == 100000:
                    pg_con.commit()
                    LOG.write(gt()+"\t\tLoaded 100,000 rows into wh.patient_load\n")
                    counter = 0

            pg_con.commit()
            LOG.write(gt()+"\t\twh.patient_load was updated successfully..."+str(len(cxr))+" records added\n")
            
            # Begin prepping patient_load table for transfer to patient table
            ###################################################################
            LOG.write(gt()+"\t\tWorking with patient_load table...\n")
            # Create unique patient id's for every patient. ID is a hash
            # This is the best id possible
            sql = "update wh.patient_load set pt_hash = md5(lastname||firstname||sex||dob) where pt_hash is null"
            r1 = execPG(pgc,sql,'',LOG)
            LOG.write(gt()+"\t\tCreated level 1 pt_hash\n")

            if r1 != 'Error':
                # This is next best id possible
                sql = "update wh.patient_load set pt_hash = md5(lastname||sex||dob) where pt_hash is null"
                r2 = execPG(pgc,sql,'',LOG)
                pg_con.commit()
                LOG.write(gt()+"\t\tCreated level 2 pt_hash\n")

                if r2 != 'Error':
                    # This is the last best id possible
                    sql = "update wh.patient_load set pt_hash = md5(lastname||dob) where pt_hash is null"
                    r3 = execPG(pgc,sql,'',LOG)
                    pg_con.commit()
                    LOG.write(gt()+"\t\tCreated level 3 pt_hash\n")

                    # Nothing else is uniquely identifiable, so delete it. Should be very small number, ~150
                    if r3 !='Error':
                        sql = "delete from wh.patient_load where pt_hash is null"
                        r4 = execPG(pgc,sql,'',LOG)
                        pg_con.commit()
                        LOG.write(gt()+"\t\tRemoved entires where not possible to create a hash\n")
                        LOG.write(gt()+"\t\tpt_hash keys have been generated for patient_load, assigning test numbers...\n")

                    else:
                        LOG.write(gt()+"\twh.patient update was skipped, there was an error generating pt_hash\n")
                        return 'bad'
                else:
                    LOG.write(gt()+"\twh.patient update was skipped, there was an error generating pt_hash\n")
                    return 'bad'
                pg_con.commit()

                # Generate initial test_numbers
                ################################
                LOG.write(gt()+"\t\tUpdating patient test numbers...\n")
                f = StringIO()
                sql = "copy (select pt_temp_id,pt_hash from wh.patient_load order by pt_hash,collected_date,pt_temp_id) to STDOUT"
                try:
                    pgc.copy_expert(sql,f)
                    LOG.write(gt()+"\t\tGot ids from patient_load\n")
                except:
                    LOG.write(gt()+"\t\tCould not execute COPY from patient_load. Update aborted for wh.patient\n")
                    return 'bad'

                # Begin assigning test numbers, writeing to a temp file obj
                f.seek(0)   # reset the IO buffer to the beginning for reading
                c = 0
                cid = ''
                sql = "update wh.patient_load set test_num=%s where pt_temp_id=%s"

                tmp = StringIO()    # This object will hold all the new numbers
                count = 0
                for line in f:
                    a = line.replace("\\N","").rstrip().split('\t') # Remove null values, trailing whitespace
                    if a[1] != cid:
                        cid = a[1]
                        c = 1
                        tmp.write(a[0]+"\t"+str(c)+"\n")
                    elif a[1]==cid:
                        c = c + 1
                        tmp.write(a[0]+"\t"+str(c)+"\n")
                    count += 1
                    
                    if count == 100000:
                        LOG.write(gt()+"\t\tAssigned test numbers to 100,000 records\n")
                        count = 0
                
                # reset tmp buffer to beginning and create a temp table to transfer data to
                tmp.seek(0)
                execPG(pgc,"create table wh.tmp_test_count (pt_id integer,test integer);",'',LOG)
                pg_con.commit()
                execPGCopyExpert(pgc,tmp,'wh.tmp_test_count',LOG)
                pg_con.commit()

                # Update patient_load with data from the temp table
                execPG(pgc,"update wh.patient_load set test_num = a.test from tmp_test_count a where a.pt_id=pt_temp_id",'',LOG)
                LOG.write(gt()+"\t\tFinished putting test num back in patient_load\n")
                pg_con.commit()

                # Drop the temp table
                execPG(pgc,"drop table wh.tmp_test_count",'',LOG)
                LOG.write(gt()+"\t\tFinished assigning test_numbers and moving to load wh.patient_load\n")

                # Recreate Patient table
                #########################
                execPG(pgc,createPatientTable(),'',LOG)
                pg_con.commit()

                # Make the transfer to wh.patient and drop load table
                ######################################################
                pt_rec = 0
                sql = "select count(pt_hash) from wh.patient_load"
                r5 = getPG(pgc,sql,'',LOG)
                pt_rec = r5[0][0]
                load_patient_sql="""insert into wh.patient(vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs,pt_hash,test_num)
                                    select vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs,pt_hash,test_num from wh.patient_load"""
                try:
                    r6 = execPG(pgc,load_patient_sql,'',LOG)
                except:
                    LOG.write(gt()+"\t\tCould not transfer data from patient_load to wh.patient. Update of wh.patient aborted\n")
                    return 'bad'
                LOG.write(gt()+"\twh.patient was loaded successfully..."+str(pt_rec)+" of "+str(len(cxr))+" records added.\n\t\t\tAny that remainded could not make pt_hash and were dropped\n")
                execPG(pgc,"drop table patient_load",'',LOG)
                f.close()
                tmp.close()
                pg_con.commit()
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



############################
# Check if mutant is a RAM #
############################
def checkRAMS(mutant, mut_type):
    global summary
    pattern = re.compile(r'(\w(\d+)([\w+\^*]))')    # Regular expression pattern to identify letter|number(s)|letter(s) or ^ or * pattern
    
    match = re.findall(pattern,mutant)
    p = match[0][1]     # position
    r = match[0][2]     # mutation or insertion amino acids

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


##########################################
# Record all mutants present in sequence #
##########################################
def seenMut(mutant,mut_type):
    global seen
    if mut_type == 'tams' or mut_type == 'nams' or mut_type == 'nnrti': mut_type = 'rt'
    if seen[mut_type].has_key(mutant): seen[mut_type][mutant] = seen[mut_type][mutant] + 1
    else: seen[mut_type][mutant] = 1
    return


########################
# Evaluate each mutant #
########################
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


################################################################
# Load seq data into table, clculate mut_freq and load as well #
################################################################
def loadSeqData(pgc,pg_con,cxc,LOG):
    LOG.write(gt()+"\tBeginning to process data for wh.seq....\n")
    
    # Get all seq data from Oracle
    ###################################
    cxr = getOracle(cxc,exportSequenceData(),'',LOG)

    if cxr != "Error":
        
        # Drop and recreate seq and mut_freq tables
        ############################################
        execPG(pgc,createSeqTable(),'',LOG)
        pg_con.commit()
        execPG(pgc,createMutFreqTable(),'',LOG)
        pg_con.commit()
        LOG.write(gt()+"\t\tCreated wh.seq and wh.mut_freq tables\n")

        # Insert new records in wh.seq
        ###############################
        getcontext().prec=2     # Sets precision for decimal to 2 digits
        global seen, summary, scores, muts
        muts = {'pri':PI_RAMS,'tams':TAMS,'nams':NAMS,'nnrti':NNRTI_RAMS,'int':INT_RAMS,'rt2':RT2_RAMS}
                
        LOG.write(gt()+"\t\tLoading data to wh.seq and wh.mut_freq table...\n")
            
        # Begin processing record (Oracle export for building HIVWH - contains pri_summary and rt_summary fields)
        for row in cxr:        
            vals =[]
            pri_sum = []
            rt_sum = []
            int_sum = []
            rt2_summary = []
            for x in row: vals.extend([x]) # Convert tuple from db into array
            seen = {'pri':{},'rt':{},'int':{},'rt2':{}}     # Record of all protease mutations that have been seen, whether they are RAMs or not
            summary = {'tams':'','nams':'','nnrti':'','pri':'','int':'','rt2':''}
            scores = {'tams':Decimal(0),'nams':Decimal(0),'nnrti':Decimal(0),'pri':Decimal(0),'int':Decimal(0),'rt2':Decimal(0)}
            wt=''
                    
            check_freq='yes'
                    
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

            # Begin processing individual mutations frequencies
            ####################################################
            if check_freq == 'yes':
                # Create list of all mutation summary fields so we can do one big query and try to speed up things
                indivmuts = ""
                indiv_muts = {}
                for i in (5,6,7,8):
                    if catchNone(vals[i]) != '':
                        indivmuts = indivmuts+vals[i]

                #indivmuts = catchNone(vals[5]) + catchNone(vals[6]) + catchNone(vals[7]) +catchNone(vals[8])
                indivmuts = "'"+"','".join(indivmuts.replace(' ','').replace('"','').replace('/','').split(','))+"'"
                               
                sql = "select mutation, mutation_sep, mut_freq,mut_type,id from wh.mut_freq where mutation in (%s)"
                mutfreq_dbr = getPG(pgc,sql%indivmuts,'',LOG)
                if mutfreq_dbr != "Error":
                    mutids = []  # pull out ids, will want to delete these later so we can use copy_from instead of having to update them
                    if len(mutfreq_dbr) != 0:
                        # Load dict with query results
                        for record in mutfreq_dbr: 
                            indiv_muts[record[0]] = {'mutation_sep':record[1],'mut_freq':record[2],'mut_type':record[3],'id':str(record[4])}
                            mutids.append(str(record[4]))

                    # Evaluate mutations wrt freq. If already present, increment value. If not, add them to dict.
                    for x in pri_sum:
                        y = x.replace('/','')
                        if indiv_muts.has_key(y) and indiv_muts[y]['mut_type'] == 'Protease':
                            indiv_muts[y]['mut_freq'] = str(int(indiv_muts[y]['mut_freq']) + 1)
                        else:
                            indiv_muts[y]= {'mutation_sep':x,'mut_freq':'1','mut_type':'Protease'}

                    for x in rt_sum:
                        y = x.replace('/','')
                        if indiv_muts.has_key(y) and indiv_muts[y]['mut_type'] == 'RT':
                            indiv_muts[y]['mut_freq'] = str(int(indiv_muts[y]['mut_freq']) + 1)
                        else:
                            indiv_muts[y]= {'mutation_sep':x,'mut_freq':'1','mut_type':'RT'}

                    for x in int_sum:
                        y = x.replace('/','')
                        if indiv_muts.has_key(y) and indiv_muts[y]['mut_type'] == 'Int':
                            indiv_muts[y]['mut_freq'] = str(int(indiv_muts[y]['mut_freq']) + 1)
                        else:
                            indiv_muts[y]= {'mutation_sep':x,'mut_freq':'1','mut_type':'Int'}

                    for x in rt2_summary:
                        y = x.replace('/','')
                        if indiv_muts.has_key(y) and indiv_muts[y]['mut_type'] == 'RT2':
                            indiv_muts[y]['mut_freq'] = str(int(indiv_muts[y]['mut_freq']) + 1)
                        else:
                            indiv_muts[y]= {'mutation_sep':x,'mut_freq':'1','mut_type':'RT2'}

                    # Get ready to write output to WH. First clear out records for ones that were already there.
                    if len(mutids) != 0:
                        mutids_str = ','.join(mutids)
                        sql = "delete from wh.mut_freq where id in (%s)"
                        execPG(pgc,sql%mutids_str,'',LOG)

                    # Build result string to copy back to WH
                    data = ""
                    for m in indiv_muts:
                        data = data + "\t".join([m, indiv_muts[m]['mutation_sep'], indiv_muts[m]['mut_freq'], indiv_muts[m]['mut_type']])+"\n"

                    # Write mut_freq data back to table
                    f = StringIO(data)  # create a file-like object
                    pgc.copy_from(f,'wh.mut_freq',columns=('mutation','mutation_sep','mut_freq','mut_type'))

                    del(mutids)
                else:
                    LOG.write(gt()+"\tThere was an error updating wh.mut_freq...update cancelled")
                    return 'bad'

        pg_con.commit()
        LOG.write(gt()+"\twh.seq and wh.mut_freq were loaded successfully..."+str(len(cxr))+" records added\n")
        return 'ok'  
    else:
        if cxr == 'Error': 
            LOG.write(gt()+"\tFailed to retrieve records from Oracle, wh.Seq was not loaded\n")
            return 'bad'
    return 'ok'


###########################################################
# Create and load wh.GT_Rule_Result_Rows table from MySQL #
###########################################################
def loadGTRuleResultRows(pgc,pg_con,LOG):
    status, kong, kong_db = connectMySQL(LOG)
    LOG.write(gt()+"\tCreating wh.gt_rule_result_rows....\n")

    # Drop and recreate seq and mut_freq tables
    ############################################
    execPG(pgc,createGTRuleResultRowsTable(),'',LOG)
    pg_con.commit()

    # Update MGRM data
    # Get Kong records 
    ######################
    mgrm = getMySQL(kong,exportGTRuleResultRowsData(),'mgrm',LOG)
    if mgrm != "Error":
        # Verify there are records to insert
        if len(mgrm) != 0:
            LOG.write(gt()+"\t\tWriting MGRM data to wh.gt_rule_result_rows table...\n")
            # Insert new records in wh.gt_rule_result_rows
            for SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src in mgrm:
                sql = "insert into wh.gt_rule_result_rows(SAMPLE_ID,DRUG_ABRV,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                vals = (str(SAMPLE_ID),str(DRUG_ABRV_NAME),str(DRUG_CLASS),str(COMMENTS),str(RAMs),RESULT,RESISTANT_FLG,MUTATION_SCORE,str(CREATED_DATE),RULE_VERSION,ENGINE_VERSION,str(data_src))
                msr = execPG(pgc,sql,vals,LOG)
            LOG.write(gt()+"\twh.gt_rule_result_rows part 1 creation was successful..."+str(len(mgrm))+" records added\n")
            
        else:
            LOG.write(gt()+"\twh.gt_rule_result_rows population part 1 failed, no records were retrieved from Kong\n")
    else:
        LOG.write(gt()+"\Error retrieving records from KONG, wh.gt_rule_result_rows update part 1 was not done\n")
    
    # Update LABC data
    # Get KongDB records 
    #######################
    labc = getMySQL(kong_db,exportGTRuleResultRowsData(),'labc',LOG)
    if labc != "Error":
        # Verify there are records to insert
        if len(labc) != 0:
            LOG.write(gt()+"\t\tWriting LABC data to wh.gt_rule_result_rows table...\n")
            # Insert new records in wh.gt_rule_result_rows
            for SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src in labc:
                sql = "insert into wh.gt_rule_result_rows(SAMPLE_ID,DRUG_ABRV,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,data_src) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                vals = (str(SAMPLE_ID),str(DRUG_ABRV_NAME),str(DRUG_CLASS),str(COMMENTS),str(RAMs),RESULT,RESISTANT_FLG,MUTATION_SCORE,str(CREATED_DATE),RULE_VERSION,ENGINE_VERSION,str(data_src))
                msr = execPG(pgc,sql,vals,LOG)
            LOG.write(gt()+"\twh.gt_rule_result_rows part 2 creation was successful..."+str(len(labc))+" records added\n\t\t\twh.gt_rule_result_rows loading was successful...\n")
            return 'ok'
        else:
            LOG.write(gt()+"\twh.gt_rule_result_rows update part 2 failed, no records were retrieved from KongDB\n")
            return 'ok'
    else:
        LOG.write(gt()+"\Error retrieving records from KONG_DB, wh.gt_rule_result_rows update part 2 was not done\n")
        return 'bad'

    pg_con.commit()
    return 'ok'







################
# Main Program #
################
def main():
    print "hello"
    #log_pth = "/home/mevans/logs/hivdb_update_log.txt"               # Use for production on Gamera
    log_pth = "/opt/bin/logs/reboot_hivdb_log.txt"                    # Use for testing on Bali
    buffer_size = 0
    LOG = open(log_pth,"w",buffer_size)   
    print  "\n***********************\n***********************\n"+gt()+" Beginning reboot process...\n"             
    LOG.write("\n***********************\n***********************\n"+gt()+" Beginning reboot process...\n")

    notify = 'no'   # flag to email logfile if there are problems
    status = 'ok'   # flag to be set by each update process

    # Create database connections
    ##############################
    status, pgc, pg_con = connectPG(LOG)    #** Don't forget to switch Postgres connection string when going from production to testing & vice versa
    status, cxc, etg = connectORA(LOG)

    if status == 'ok':
        # Rebuild structure
        # Add the data that does not require any additional manipulation
#        checkStatus(execPG(pgc,createDrugsTable(),'',LOG))
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.drugs',exportDrugsData(),'',LOG))

#        checkStatus(execPG(pgc,createProjectTable(),'',LOG))
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.project',exportProjectData(),'',LOG))

#        checkStatus(execPG(pgc,createMutationsTable(),'',LOG))
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.mutations',exportMutationsData(),'',LOG))

#        checkStatus(execPG(pgc,createGenoRepRowTable(),'',LOG)) 
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.geno_rep_row',exportGenoRepRowData(),'',LOG))

#        checkStatus(execPG(pgc,createGenoReportTable(),'',LOG)) 
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.geno_report',exportGenoReportData(),'',LOG))

#        checkStatus(execPG(pgc,createPhenoRepRowTable(),'',LOG)) 
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.pheno_rep_row',exportPhenoRepRowData(),'',LOG))

#        checkStatus(execPG(pgc,createPhenoReportTable(),'',LOG)) 
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.pheno_report',exportPhenoReportData(),'',LOG))

#        isql = "copy psgt_result (pg_report_id,vrl_accession,pt_report_id,gt_report_id,test_code,project_code,drug_abrv,resistant_flg,result_flg,rep_date) from STDIN with null ''"
#        checkStatus(execPG(pgc,createPSGTResultTable(),'',LOG))
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.psgt_result',exportPSGTResultData(),isql,LOG))
#        checkStatus(execPG(pgc,'update wh.psgt_result set test_code = wh.geno_report.test_code from wh.geno_report where wh.psgt_result.gt_report_id=wh.geno_report.gt_report_id','',LOG))
#        pg_con.commit()
#        checkStatus(execPG(pgc,'update wh.psgt_result set project_code = wh.pheno_report.project_code from wh.pheno_report where wh.psgt_result.project_code is null and wh.psgt_result.pt_report_id=wh.pheno_report.pt_report_id','',LOG))
#        pg_con.commit()

#        isql = "copy seq_mutations (gt_report_id,seq_id,mutation_id,mutation_type,mutation,rep_date) from STDIN with null ''"
#        checkStatus(execPG(pgc,createSeqMutationsTable(),'',LOG))
#        pg_con.commit()
#        checkStatus(addData(cxc,pgc,pg_con,'wh.seq_mutations',exportSeqMutationsData(),isql,LOG))

#        isql = "copy etag_sum (vrl_accession,patient_dob,collected_date,rcv_date,rep_date,h2t,h2d,h2t_result,h2d_result,h2t_tumor_area,h2d_tumor_area,h2t_section_size,h2d_section_size,h2t_pct_tumor,h2d_pct_tumor,h2t_batch_id,h2d_batch_id,client_state,client_country) from STDIN with null ''"
#        checkStatus(execPG(pgc,createETAGSumTable(),'',LOG))
#        pg_con.commit()
#        checkStatus(addData(etg,pgc,pg_con,'wh.etag_sum',exportEtagData(),isql,LOG))

#        checkStatus(execPG(pgc,createGTRuleResultRowsTable(),'',LOG))


    
        
        # Add the more complicated stuff
#        checkStatus(loadPatientData(pgc,pg_con,cxc,LOG))
#        checkStatus(loadSeqData(pgc,pg_con,cxc,LOG))
#        checkStatus(loadGTRuleResultRows(pgc,pg_con,LOG))

        # Create views on the new tables
        checkStatus(execPG(pgc,createGenoResultSumView(),'',LOG))
        pg_con.commit()
        checkStatus(execPG(pgc,createPhenoResultSumView(),'',LOG))
        pg_con.commit()
        checkStatus(execPG(pgc,createIntegrasePhenoResultSumView(),'',LOG))
        pg_con.commit()


    LOG.close()
    print "\nReboot complete\n\n"




if __name__ == '__main__':
    main()