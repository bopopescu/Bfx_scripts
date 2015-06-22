from Bio import SeqIO
from Bio.Seq import Seq

pr = SeqIO.to_dict(SeqIO.parse('pr.txt','fasta'))
rt1 = SeqIO.to_dict(SeqIO.parse('rt1.txt','fasta'))
rt2 = SeqIO.to_dict(SeqIO.parse('rt2.txt','fasta'))

f = open('merged.txt','w')

for acc in pr.keys():
	f.write('>'+acc+'\n'+pr[acc].seq.tostring()+rt1[acc].seq.tostring()+rt2[acc].seq.tostring()+'\n')

f.close()