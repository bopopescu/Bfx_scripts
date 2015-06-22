import sys
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq


def renameFile(filename,idx):
	writepath = 'renamed/'

	if filename.find('.seq') != -1:
		print "SEQ file in ...",os.getcwd(),filename
		for record in SeqIO.parse(filename,'fasta'):
			a = record.id.split('-')
			b = a[1].split('_')
			print a
			print b
			ordernum = a[0][3:]	
			print "ordernum ",ordernum
			newfilename=""
			if idx.has_key(ordernum):
				if idx[ordernum].has_key(b[0]):
					newfilename = a[0]+'-'+idx[ordernum][b[0]][1]+"_"+b[2]+"_"+idx[ordernum][b[0]][0]+".seq"
					newdefline = a[0]+'-'+idx[ordernum][b[0]][1]+"_"+b[2]+"_"+idx[ordernum][b[0]][0]+".ab1"
					print newfilename
					f2 = open(writepath+newfilename,"w")
					f2.write(">"+newdefline+"\n"+record.seq.tostring())
					f2.close()

	elif filename.find('.ab1') != -1:
		a = filename.split('-')
		b = a[1].split('_')
		ordernum = a[0][3:]
		newfilename = ""
		if idx.has_key(ordernum):
			if idx[ordernum].has_key(b[0]):
				newfilename = a[0]+'-'+idx[ordernum][b[0]][1]+"_"+b[2]+"_"+idx[ordernum][b[0]][0]+"_"+b[3][3:]
				print newfilename
		os.rename(filename,newfilename)
		os.system("mv "+newfilename+" "+writepath+newfilename)
	return


def main():
	indexfile = open('index.txt','r')
	idx = {}
   	for line in indexfile:
   		x = line.rstrip().split("\t")
   		if idx.has_key(x[0]):
   			idx[x[0]][x[1]] = (x[2],x[3])
   		else:
   			idx[x[0]] = {x[1]:(x[2],x[3])}
   	print idx
	path = ""
	os.mkdir('renamed')
	for infile in glob.glob (os.path.join(path,"*")):
		if infile.find(".py") == -1:
			if infile.find(".txt") == -1:
				renameFile(infile,idx)
					
	print "finished"

   






if __name__ == '__main__':
   main()