from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

filename = raw_input("Enter filename of master seq file to split: ")
f1 = open('actg.seq', "w")
f2 = open('toro2.seq', "w")
f3 = open('ltm.seq', 'w')
f4 = open('zabrina1.seq','w')
f5 = open('zabrina2.seq','w')
f6 = open('zepto.seq','w')

for record in SeqIO.parse(filename,"fasta",generic_dna):
   ac = record.id.split("|")
   if ac[1] =='ACTG':     f1.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[1] =='TOROII':   f2.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[1] =='LTM':      f3.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[1] =='Zabrina1': f4.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[1] =='Zabrina2': f5.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')
   if ac[1] =='Zepto':    f6.write(">"+"|".join(ac)+'\n'+record.seq.tostring()+'\n')

print "Proces complete\n\n"
   

