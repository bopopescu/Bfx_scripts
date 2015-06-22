#!/usr/bin/env python
# encoding: utf-8
"""
createSampleXML.py

Created by Mark Evans on 2011-07-26.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys






def main(argv=None):
   print"Input file format is tab-delimited: sample_id \t Subtype \t NS3_muts \t NS4A_muts\n"
   file1 = raw_input("Enter file with HCV mutation data: ")
   file2 = file1+".xml"
   
   f1 = open(file1,'r')
   f2 = open(file2,'w')
   
   for line in f1:
      line = line.replace('"','')
      line = line.rstrip("\n")
      cols = line.split("\t")
      xmlstring = "   <sample name='"+cols[0].upper()+"'>\n"\
                  "      <NS3Mutations> "+cols[2].upper()+"</NS3Mutations>\n"\
                  "      <NS4AMutations> "+cols[3].upper()+"</NS4AMutations>\n"\
                  "      <Subtype>"+cols[1].upper()+"</Subtype>\n"\
                  "   </sample>\n"
      f2.write(xmlstring)
   f1.close()
   f2.close()
   print "Finished formatting.  Please add drug rules next\n"

if __name__ == "__main__":
   sys.exit(main())
