from Bio import SeqIO



ns3 = SeqIO.to_dict(SeqIO.parse('mgrm_seq.txt.ns3.dna.seq','fasta'))
ns4a = SeqIO.to_dict(SeqIO.parse('mgrm_seq.txt.ns4a.dna.seq','fasta'))

set1 = {}
set2 = {}

for y in ns3:
	set1[y.split('|')[0]] = ns3[y]

for z in ns4a:
	set2[z.split('|')[0]] = ns4a[z]

f = open('mgrm_seq_joined.seq','w')

ids = set1.keys()

for x in ids:
	if set2.has_key(x):
		f.write('>'+set1[x].id+'\n'+set1[x].seq.tostring()+set2[x].seq.tostring()+'\n')