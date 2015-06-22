from Bio import SeqIO
from Bio.Seq import Seq

fn = raw_input("What is the name of the dna file to QC?")

seqs={}

f = open(fn+".clean.txt",'w')

def findAmbig(seq):
    for v in ['N','V','B','H','D','M','K','W','S','Y','R','-','.']:
        if v in seq: return 1
    return 0

c = 1

for record in SeqIO.parse(fn,'fasta'):
    s = record.seq.tostring().upper()
    if findAmbig(s) == 0 and seqs.has_key(s)==0:
        a = record.id.replace('_','.').split('.')
        rid = a[0]+'.'+str(c)
        f.write('>'+rid+'\n'+s+'\n')
        seqs[s] = record.id
        c += 1

f.close()

