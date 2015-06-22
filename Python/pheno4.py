from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import os,math

#*******************************************************
#* doHMMPfam                                           *
#* Takes sequence input and searches V3 HMM's          *
#* using hmmpfam. Parses result and returns array with *
#* phenotype scores and E-values                       *
#* @param sequence string                              *
#* @return array of phenotype scores and E-values      *
#*******************************************************/
def doHMMPfam(seq):
   infile = open("tmp_input.seq",'w')
   infile.write(">sequence\n"+seq)
   infile.close()	
   bin_path="/usr/local/bin/hmmpfam"  # change depending on system being run on
   
   # Run sequence against proper chain HMMs
 #  result = os.popen(bin_path+" -A0 pheno.hmm tmp_input.seq",'r')
   result = os.popen(bin_path+" -A0 all_pheno.hmm tmp_input.seq",'r')
   hmm ={}
   read ='no'
   count = 1
   
   for line in result:
      if read=='no':
         if (count == 17):
            words = line.split()
            read ='yes'
            hmm[words[0]] = [float(words[1]),float(words[2])]
            count +=1
         else: count += 1
      elif read=='yes':
         if len(line) != 1:
            words = line.split()
            hmm[words[0]] = [float(words[1]),float(words[2])]
            count +=1
         else: read='no'

   os.system('rm tmp_input.seq')
   
   return hmm

#############################################################################################################
# assocCoeff                                                                                                #
# Function to calculate association coefficient, true pos percent, false pos percent, and accuracy          #
# from a dictionary of data                                                                                 #
# @input: dictionary, index of actual call, index of data to be evaluated, pheno to be predicted (DM or R5) #
# @return dictionary of cutoff:[tn,fn,tp,fp,phi,tpp,fpp,acc]                                                #
#############################################################################################################
def assocCoeff(data,ac,sc,predictor):
   # ac = actual call index, sc = stat calculated call index, predictor = pheno a positive score predicts
   # set variables
   acstats={}
   pheno = {'DM':['DM','R5'],
            'R5':['R5','DM']}
   for x in range(50,100):
      fp=0.0
      fn=0.0
      tp=0.0
      tn=0.0
      pid=''
      for v3id in data:
         # if predictor matches mgrm call, this is the positive direction (tp,fp)
         if predictor == 'R5':
            # check if score is below current cutoff x and call pheno for this round
            if float(data[v3id][sc]) <= x: pid=pheno[predictor][1]
            else: pid=pheno[predictor][0]
         elif predictor == 'DM':
            if float(data[v3id][sc]) >= x: pid=pheno[predictor][0]
            else: pid=pheno[predictor][1]

         if   pid == 'DM' and pid == data[v3id][ac]: tp = tp + 1     # true positive call
         elif pid == 'DM' and pid != data[v3id][ac]: fp = fp + 1     # false positive
         elif pid == 'R5' and pid == data[v3id][ac]: tn = tn + 1     # true negative
         elif pid == 'R5' and pid != data[v3id][ac]: fn = fn + 1     # false negative

         pid=''
      try:
         # phi = sqrt(((tp*tn - fp*fn)^2/(tp+fp)(fn+tn)(tp+fn)(fp+tn)))  association coefficient 
         float(phi)
     #    phi = math.sqrt((math.pow(tp*tn - fp*fn,2)/(tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
      except:
         phi = 0.0
      try: tpp = tp/(tp+fn)  # true positive percentage
      except: tpp = 0.0
      try: fpp = fp/(fp+tn)  # false positive percentage
      except: fpp = 0.0
      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
      except: acc = 0.0
      try: spec = tn/(tn+fp)  # specificity for X4
      except: spec = 0.0
      try: sens = tp/(tp+fn)  # sensitivity for X4
      except: sens = 0.0
      acstats[str(x)] = {'tn':tn,'fn':fn,'tp':tp,'fp':fp,'phi':phi,'tpp':tpp,'fpp':fpp,'acc':acc,'spec':spec,'sens':sens}
   return (acstats)

######################################### Main
filename = raw_input("Enter the name of the sequence file (.seq)? ")
outfile = open(filename+".out",'w')

# use with pheno.hmm
#outfile.write("v3id\tTF\tES\tHMM\tHMM_R5_score\tHMM_R5_Eval\tHMM_DX_score\tHMM_DX_Eval\n")
# Use with all_pheno.hmm
outfile.write("v3id\tTF\tES\tHMM\tBestHit\tbothR5_score\tbothR5_Eval\tbothDX_score\tbothDX_Eval\tesR5_score\tesR5_Eval\tesDX_score\tesDX_Eval\tmixedR5DX_score\tmixedR5X4_Eval\ttrofileR5_score\ttrofileR5_Eval\ttrofileDX_score\ttrofileDX_Eval\n")
print "Beginning HMM scans...\n"
for record in SeqIO.parse(filename,"fasta",generic_dna):
   ac = record.id.split("|")
   stats = doHMMPfam(record.seq.tostring())
   
   # Choose between all R5 and all DM to make HMM call for now
   H=''
   if stats['bothR5'][1] < stats['bothDX'][1] and stats['bothR5'][1] != stats['bothDX'][1]: H='R5'
   else: H='DM'

   # Find the hmm with the lowest eval and report it as the best hit for diagnostics
   foo =100
   foo_id=''
   for var in stats:
      if stats[var][1] < foo: 
         foo = stats[var][1]
         foo_id = var

   # Use with pheno.hmm
 #  code = "%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\n" % (ac[0],ac[2],ac[3],H,stats['bothR5'][0],stats['bothR5'][1],stats['bothDX'][0],stats['bothDX'][1])
   # Use with all_pheno.hmm
   code = "%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" % (ac[0],ac[2],ac[3],H,foo_id,stats['bothR5'][0],stats['bothR5'][1],stats['bothDX'][0],stats['bothDX'][1],stats['esR5'][0],stats['esR5'][1],stats['esDX'][0],stats['esDX'][1],stats['mixedR5DX'][0],stats['mixedR5DX'][1],stats['trofileR5'][0],stats['trofileR5'][1],stats['trofileDX'][0],stats['trofileDX'][1])

   outfile.write(code)
outfile.close()
print "Finished writing HMM scan data.....\n"

# Statistically Analyze the HMM output 
###########################################
infile = open(filename+".out","r")
v = [('bothR5',1,'R5'),('bothDX',2,'DM'),('esR5',3,'R5'),('esDX',4,'DM'),('mixedR5DX',5,'DM'),('trofileR5',6,'R5'),('trofileDX',7,'DM')]

print "Reading in hmm results....\n"
hmmresults = {}
for line in infile.readlines():
   a = line.split('\t')
   if a[0] !='v3id':
      mgrm = ''
      # Consolidate Mgrm calls into a single R5/DM call
      if a[1] == 'DM' or a[1]=='X4' or a[2]=='DM' or a[2]=='X4': mgrm = 'DM'
      else: mgrm = 'R5'
      # key = v3id, array = [mgrm_pheno,bothR5_score, bothDX_score, esR5_score, esDX_score, mixedR5DX_Score,trofileR5,trofileDX] 
      hmmresults[a[0]]= [mgrm,a[5],a[7],a[9],a[11],a[13],a[15],a[17]]

for a,b,c in v:
   f = open(filename+"."+str(a)+".stats","w")
   stats = assocCoeff(hmmresults,0,int(b),str(c))
   f.write("Cutoff\tTn\tFn\tTp\tFp\tPhi\tTPP\tFPP\tAccuracy\tSpecificity\tSensitivity\n")
   cutoffs = stats.keys()
   cutoffs.sort()
   for c in cutoffs:
      f.write(str(c)+"\t"+"\t".join([str(stats[c]['tn']),str(stats[c]['fn']),str(stats[c]['tp']),str(stats[c]['fp']),str(stats[c]['phi']),str(stats[c]['tpp']),str(stats[c]['fpp']),str(stats[c]['acc']),str(stats[c]['spec']),str(stats[c]['sens'])])+"\n")
   f.close()
print "Finished analyzing HMM data\n"