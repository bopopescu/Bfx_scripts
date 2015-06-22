#!/usr/bin/env python
# encoding: utf-8
"""
getpepstats.py

Created by Mark Evans on 2011-01-05.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os, subprocess
import psycopg2 as pg
import random


# Connect to database. Should be in 'try' block but complained about not having been declared
# so we are going at risk ;-)
def pgConnect():
   conn = pg.connect("dbname='mgrm' user='mark' host='localhost' password='mark'")
   return conn.cursor()


# Calls pepstat, returns output filename
def runPepstats(tmpid):
   fileloc = tmpid+".seq"
   outfile = tmpid+".pepstat"
   result = os.popen("pepstats -sequence "+fileloc+" -outfile "+outfile+" -sformat1 fasta -sprotein1 -aadata Eamino.dat -mwdata Emolwt.dat -termini -nomono -auto") # os.popen is deprecated but easier
   return outfile


# Parses pepstat output into two dictionaries, one for the amino acid DayhoffStat and one for everything else
def parsePepstats(tmp):
   p = open(tmp,"r")
   result = {'mw':0,'charge':0,'iep':0,'Tiny':0,'Small':0,'Aliphatic':0,'Aromatic':0,'Non-polar':0,'Polar':0,'Charged':0,'Basic':0,'Acidic':0}
   DayhoffStat = {'A':0,'B':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'J':0,'K':0,'L':0,'M':0,'N':0,'O':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'U':0,'V':0,'W':0,'X':0,'Y':0,'Z':0}
   lines = p.readlines()
   p.close()
   
   result['mw']     = lines[2].split()[3]
   result['charge'] = lines[3].split()[7]
   result['iep']    = lines[4].split()[3]
   for x in range(38,47):
      result[lines[x].split()[0]] = lines[x].split()[3]
   for y in range(10,36):
      DayhoffStat[lines[y].split()[0]] = lines[y].split()[5]
   del DayhoffStat['J']   #
   del DayhoffStat['O']   # These are never used, but inserted just to make the code simple
   del DayhoffStat['U']   #
   
   return (result,DayhoffStat)
   
   
def main():
   # Establish database connection and fetch sequence record set
   #############################################################
   cur = pgConnect()
   sql = """select v3id,mgrmseq from master_v3_commercial_set_152"""
   cur.execute(sql)
   rows = cur.fetchall()
   
   # Create output files and begin generating and parsing pepstat output
   #####################################################################
   RESULTS = open("pepstats_result.txt","w")
   RESULTS.write("v3id\tmw\tcharge\tiep\tTiny\tSmall\tAliphatic\tAromatic\tNon-polar\tPolar\tCharged\tBasic\tAcidic\n")
   DAY = open("pepstats_dayhoff.txt","w")
   DAY.write("v3id\tA\tB\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tX\tY\tZ\n")
   for row in rows:
      # Create temp file for pepstat to read
      t = str(random.random())[2:9] # generate temp id
      f = open(t+".seq","w")
      f.write(">"+row[0]+"\n"+row[1])
      f.close()
      
      # Run pepstats, parse results, catch returned dictionaries with results
      r,d = parsePepstats(runPepstats(t))
      
      # Parse results for output
      dkeys = d.keys()
      dkeys.sort()
      dvals = []
      for i in dkeys: dvals.append(d[i])
      
      # Write output
      RESULTS.write(row[0]+"\t"+"\t".join([r['mw'],r['charge'],r['iep'],r['Tiny'],r['Small'],r['Aliphatic'],r['Aromatic'],r['Non-polar'],r['Polar'],r['Charged'],r['Basic'],r['Acidic'],])+"\n")
      DAY.write(row[0]+"\t"+"\t".join(dvals)+"\n")
      
      # Cleanup temp files
      os.remove(t+".seq")
      os.remove(t+".pepstat")
   RESULTS.close()
   DAY.close()

if __name__ == '__main__':
   main()

