#!/usr/bin/env python
# encoding: utf-8
"""
hivdb_update_data.py

Created by Mark Evans on 2011-10-14.

Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""


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


#-----------------------------------------------------------------
#-- Update wh.drugs table
#-- Just re-export entire thing and reload wh table
#-----------------------------------------------------------------
update_drugs_sql = """SELECT
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

WHERE d.drug_class_id=dc.id"""

#-----------------------------------------------------------------
#-- Update wh.projects table
#-- Just re-export entire thing and reload wh table
#-----------------------------------------------------------------
update_projects_sql="""SELECT
p.id AS project_id,
p.code AS project_code,
p.project_name,
p.project_type

FROM 
reportmanager.projects p"""

#-----------------------------------------------------------------
#-- Update wh.mutations table
#-- created_date is already in YYYY-MM-DD format
#-----------------------------------------------------------------
update_mutations_sql="""SELECT 
m.id AS mutation_id, 
(CASE m.type  
          WHEN 1 THEN 'RT' 
          WHEN 2 THEN 'Prt' 
          WHEN 3 THEN 'Int' 
          WHEN 4 THEN 'RT2' 
          ELSE TO_CHAR(m.type) 
END) AS mutation_type, 
m.abrv_name AS mutation, 
m.position, 
m.orig_protein AS wt, 
m.new_protein AS mut, 
m.insertion_flg, 
m.CREATED_DATE   

FROM 
genotype.mutations m 

WHERE 
m.created_date > TO_DATE('%s','YYYY-MM-DD HH24:MI:SS')"""  #  --  <=== replace with max(created_date) from wh.mutations

#-----------------------------------------------------------------
#-- Update wh.seq table
#-- rep_date is already in YYYY-MM-DD format
#-----------------------------------------------------------------
update_seq_sql="""select distinct
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
g.rep_date,      
g.project

from genotype.gt_reports g, genotype.data_files_nucleotide_seq d
where d.data_file_id = g.data_file_id
and g.status_code =1
and g.rejected_report_flg = 0
and g.rep_date > to_date('%s','YYYY-MM-DD HH24:MI:SS')

union
select distinct      
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
g.rep_date,
g.project
from genotype.gt_reports g
where g.data_file_id < 25972
and g.status_code =1
and g.rejected_report_flg=0
and g.rep_date > to_date('%s','YYYY-MM-DD HH24:MI:SS')"""  #--   <=== replace with max(rep_date) from wh.seq


#-----------------------------------------------------------------
#-- Update seq_mutations table
#-- rep_date is already in YYYY-MM-DD format
#-----------------------------------------------------------------
update_seq_mutations_sql = """SELECT
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
AND g.rejected_report_flg =0
and g.rep_date > to_date('%s','YYYY-MM-DD HH24:MI:SS')"""  #--  <==== replace with max(rep_date) from wh.geno_report


#-----------------------------------------------------------------
#-- Update wh.PSGT_Result table
#-- Must be queried from PSGT tablespace
#-- rep_Date is already in YYYY-MM-DD format
#-----------------------------------------------------------------
update_psgtresult_sql = """SELECT
p.pg_report_id,
p.accession_id AS vrl_accession,
pt.id AS pt_report_id,
g.id AS gt_report_id,
p.test_code,
p.project_code,
pr.drug_abrv_name AS drug_abrv,
pr.resistant_flg,
pr.result_flg,
p.rep_date

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
AND p.accession_id = g.vrl_accession
AND p.rep_date > TO_DATE('%s','YYYY-MM-DD HH24:MI:SS')"""  #--  <====  replace with max(rep_date) from wh.psgt_result


#-----------------------------------------------------------------
#-- Update geno_rep_rows table
#-----------------------------------------------------------------
update_geno_rep_row_sql="""SELECT
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
AND   r.drug_abrv_name = d.abrv_name
AND g.rep_date > to_date('%s','YYYY-MM-DD HH24:MI:SS')"""  #--  <=== replace with max(rep_date) from wh.geno_report


#-----------------------------------------------------------------
#-- Update geno_report table
#-- rep_date is already in YYYY-MM-DD format
#-----------------------------------------------------------------
update_geno_report_sql="""SELECT distinct
id AS gt_report_id,
vrl_accession,
aliquot_id,
subtype,
rule_version_used AS rule_ver,
rep_date,
test_code,
assay_type

FROM  gt_reports g
WHERE g.status_code = 1
AND   g.rejected_report_flg =0
AND g.rep_date > to_date('%s','YYYY-MM-DD HH24:MI:SS')"""  #--  <=== replace with max(rep_date) from wh.geno_report


#-----------------------------------------------------------------
#-- Update for pheno_report table
#-- rep_date is already in YYYY-MM-DD format
#-----------------------------------------------------------------
update_pheno_report_sql="""SELECT
p.id as pt_report_id,
p.vrl_accession,
p.aliquot_id,
p.test_code,
p.assay_type,
c.rc_status,
c.reported_rc,
p.rep_date,
p.project as project_code

FROM 
vlps.pt_reports p, vlps.rc c

WHERE   
p.status_code = 1
AND p.rejected_report_flg =0
AND   c.ID (+) = p.RC_ID
AND p.REP_DATE > to_date('%s','YYYY-MM-DD HH24:MI:SS')"""  #--  <==== replace with max(rep_date) from wh.pheno_report


#-----------------------------------------------------------------
#-- Update for pheno_rep_rows table
#-----------------------------------------------------------------
update_pheno_rep_row_sql="""SELECT
r.id AS rpt_row_id,
r.pt_report_id AS report_id,
r.ic50,
r.ic95,
r.fold_res,
CASE 
    WHEN r.fold_change IS NULL THEN
         CASE 
             WHEN r.ic50 LIKE '>%%' THEN r.max_fold_change 
           WHEN r.ic50 LIKE '<%%' THEN '0.3'
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
AND   p.status_code = 1
AND   p.rep_date > to_date('%s','YYYY-MM-DD HH24:MI:SS')"""  #--  <=== replace with max(rep_date) from wh.pheno_report


#-----------------------------------------------------------------
#-- Create patient_load table and sequence
#-----------------------------------------------------------------
create_pt_load="""CREATE TABLE wh.patient_load
   (vrl_accession text,
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
    test_num int4 )"""


#-----------------------------------------------------------------
#-- Update for patient table
#----------------------------------------------------------------- 
#-- Get last patient received date
#-- select max(to_char(to_timestamp(received_date,'DD-MON-YYYY HH24:MI:SS'),'YYYY-MM-DD HH24:MI:SS')) from wh.patient  e.g. 2011-07-19 10:40:00
#-- should convert all dates to be in this format already, sorts the best
update_patient_sql="""SELECT
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
and r.received_date > to_date('%s','YYYY-MM-DD HH24:MI:SS') order by r.received_date,r.vl_accession_id,pn.vl_test_code"""  #-- <== replace date with last date in wh.patients


#--------------------------------------------------------------------
#-- Transfer data from patient_load to patient table
#--------------------------------------------------------------------
load_patient_sql="""insert into wh.patient(vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs,pt_hash,test_num)
 select vrl_accession,lastname,firstname,middlename,sex,dob,pt_zip,pt_state,client_code,client_name,client_country,client_city,client_state,client_zip,project_code,doctor_code,client_notes,visit_num,site_num,viroload,account_class,collected_date,received_date,test_code,panel_code,product_name,product_type,reported_drugs,pt_hash,test_num from wh.patient_load"""


#-------------------------------------------------------------------
#- Update SQL from Kong for GT_RULE_RESULT_ROWS (MGRM)
#-------------------------------------------------------------------
kong_update_sql="""select SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,'mgrm' from GT_RULE_RESULT_ROWS where SAMPLE_ID !='dummy id' and CREATED_DATE > %s """


#--------------------------------------------------------------------
#- Update SQL from KongDB for GT_RULE_RESULT_ROWS (LABC)
#--------------------------------------------------------------------
kongdb_update_sql="""select SAMPLE_ID,DRUG_ABRV_NAME,DRUG_CLASS,COMMENTS,RAMs,RESULT,RESISTANT_FLG,MUTATION_SCORE,CREATED_DATE,RULE_VERSION,ENGINE_VERSION,'labc' from GT_RULE_RESULT_ROWS where SAMPLE_ID !='dummy id' and CREATED_DATE > %s """


#--------------------------------------------------------------------
#- Update for etag_sum table. Connect via etag schema to battle
#--------------------------------------------------------------------
update_etag_sql = """SELECT distinct
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
and to_char(h2t.rep_date,'YYYY-MM-DD') > '%s'
order by collected_date, rep_date"""   #-- <== replace date with last rep_date in wh.etag_sum


