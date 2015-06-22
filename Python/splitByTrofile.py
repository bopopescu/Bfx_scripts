from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

filename = raw_input("Enter name of seq file to parse by Trop call ")
f1 = open('trofileR5.seq', "w")
f2 = open('trofileDX.seq', "w")
f3 = open('esR5.seq', 'w')
f4 = open('esDX.seq','w')
f5 = open('bothR5.seq','w')
f6 = open('bothDX.seq','w')
f7 = open('mixedR5DX.seq','w')
for record in SeqIO.parse(filename,"fasta",generic_dna):
   ac = record.id.split("|")
   if ac[2] =='R5': f1.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[3] =='R5': f3.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[2] !='R5' and ac[2] !='': f2.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[3] !='R5' and ac[3] !='': f4.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[2] =='R5' and ac[3] =='R5': f5.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[2] !='R5' and ac[2] !='' and ac[3]!='R5' and ac[3]!='': f6.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if (ac[2] =='R5' and ac[3]!='R5' and ac[3]!='') or (ac[2]!='R5' and ac[2]!='' and ac[3]=='R5'): f7.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')

print "Process complete\n\n"
