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
TO_CHAR(m.CREATED_DATE,'DD-MON-YYYY HH24:MI:SS') as created_date 

FROM 
genotype.mutations m 

WHERE 
m.created_date > TO_DATE('%s','YYYY-MM-DD HH24:MI:SS')"""  #  --  <=== replace with max(created_date) from wh.mutations

#-----------------------------------------------------------------
#-- Update wh.seq table
#-----------------------------------------------------------------
update_seq_sql="""select
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
select
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
















