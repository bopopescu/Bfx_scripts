# purpose is to parse tab-delimited BLAST output
# and evaluate hits for genotype consistency, i.e. try to identify potential mixed
# genotypes

import sys
from decimal import Decimal

filename = raw_input("What is filename of tab-delimited BLAST file to parse? ")
f1 = open(filename,'r')
f2 = open(filename+".result",'w')  # Put everything in here, both possible and probable
f3 = open(filename+".urgent",'w')  # Put probable in here 

current_seq =''
current_geno=''
mgrm_geno = ''
geno_class=''
count = 0
best_geno =''
best_pct=''
mix_count=0
total_mix=0
code = ''
flg = 'no'

# Loop through each Blast record
for record in f1:
	vals = record.rstrip().split('\t')
	a = vals[0].split('|')
	count += 1
	
	# get the query id for the first set of results
	if current_seq == '' or a[0] != current_seq:
		count = 1
		flg = 'no'
		code = "#####################################\n"
		code = code + "count="+str(count)+str(vals)+"\n"
		
		# Since this is the top hit for this query, set everything based on this record
		current_seq = a[0]
		mgrm_geno = a[2]
		current_geno = vals[1].split('.')[0]
		geno_class = current_geno[:-1]
		best_geno = current_geno
		best_pct = str(vals[2])

		# This script was intended to check previous subtype calls made by MGRM, so this section can be eliminated
		if mgrm_geno != current_geno: 
			print '**** WRONG CALL: '+'\t'.join([current_seq,mgrm_geno,current_geno,vals[2],vals[3]])
			f2.write('**** WRONG CALL: '+code+'\n')
	else:
		if count < 5:
			code = code + "count="+str(count)+str(vals)+"\n"
#			print "count="+str(count),vals
			current_geno = vals[1].split('.')[0]

			# if the subtype changes
			if mgrm_geno != current_geno:
				if Decimal(best_pct)-Decimal(vals[2]) <=10:
					if Decimal(best_pct)-Decimal(vals[2]) <=5:
						code = code+ '!!!!!> MIXTURE: '+'\t'.join([current_seq,mgrm_geno,best_geno,best_pct,' hit# '+str(count),current_geno,vals[2],vals[3]])+'\n\n'
						print '!!!!!> MIXTURE: '+'\t'.join([current_seq,mgrm_geno,best_geno,best_pct,' hit# '+str(count),current_geno,vals[2],vals[3]])
						mix_count += 1
						f3.write(code)
					else:
						code = code+ '===> Possible mixture: '+'\t'.join([current_seq,mgrm_geno,best_geno,best_pct,' hit# '+str(count),current_geno,vals[2],vals[3]])+'\n\n'
						print '===> Possible mixture: '+'\t'.join([current_seq,mgrm_geno,best_geno,best_pct,' hit# '+str(count),current_geno,vals[2],vals[3]])
					total_mix += 1
					f2.write(code)
					flg = 'yes'
					
					count = 5
	if count >= 5: code = ''
	if flg == 'yes':
		if vals[1].split('.')[0][:-1] != geno_class:
			f2.write(">>>>>>>> Genotype class switch occured at \n"+"       "+str(vals)+"\n\n")
			flg = 'no'
#	record_count += 1
#	if record_count == 100: sys.exit()

# These last two lines could be eliminated.
print "*********** Total potential mixtures: "+str(total_mix)+" and total probable mixtures: "+str(mix_count)
f2.write("\n*********** Total potential mixtures: "+str(total_mix)+" and total probable mixtures: "+str(mix_count))
