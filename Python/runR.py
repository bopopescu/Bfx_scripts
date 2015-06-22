from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects.lib import grid
from rpy2.robjects.vectors import DataFrame
import os, os.path, sys

def scatterPlot(stat_data,x,y,title):
   plot(stat_data.rx2(x+"_score"),stat_data.rx2(y+"_score"),pch=16,xlim=base.c(30,100),ylim=base.c(30,100),xlab=x+' score',ylab=y+' score',main=title,col=stat_data.rx2('TF'))
   robjects.r.abline(0,1)
   robjects.r.legend("bottomright",legend=base.c("R5","DM","X4"),col=base.c(2,1,3),pch=16)
   return

def createGraphSeries(cohort):
   hm = ['bothR5','bothDX','esR5','esDX','trofileR5','trofileDX','mixedR5DX']
   stat_data = DataFrame.from_csvfile(cohort+'.seq.out', sep = "\t")
   
   # Begin witing graphs to pdf
   grdevices.pdf(file=cohort+".stats.pdf",width=7,height=7)
   
   # graph 1 scatter plot
   scatterPlot(stat_data,'bothR5','bothDX',cohort.upper()+' cohort\n bothR5/bothDX correlation')
   scatterPlot(stat_data,'esR5','esDX',cohort.upper()+' cohort\n esR5/esDX correlation')
   scatterPlot(stat_data,'trofileR5','trofileDX',cohort.upper()+' cohort\n trofileR5/trofileDX correlation')
   scatterPlot(stat_data,'bothR5','esR5',cohort.upper()+' cohort\n bothR5/esR5 correlation')
   scatterPlot(stat_data,'trofileR5','esR5',cohort.upper()+' cohort\n trofileR5/esR5 correlation')
   scatterPlot(stat_data,'bothDX','esDX',cohort.upper()+' cohort\n bothDX/esDX correlation')
   scatterPlot(stat_data,'trofileDX','esDX',cohort.upper()+' cohort\n trofileDX/esDX correlation')
   scatterPlot(stat_data,'bothDX','mixedR5DX',cohort.upper()+' cohort\n bothDX/mixedR5DX correlation')

   for hmm in hm:
      hmm_scores = DataFrame.from_csvfile(cohort+'.seq.'+hmm+'.stats',sep="\t")
      # graph 3
      plot(hmm_scores.rx2("Cutoff"),hmm_scores.rx2("Accuracy"),type='o',main='Accuracy vs. Cutoff '+cohort.upper()+' Cohort\n'+hmm+'.hmm',xlab='Cutoff',ylab='Accuracy')
      # graph 4
      plot(hmm_scores.rx2("FPP"),hmm_scores.rx2("TPP"),type='o',xlim=base.c(0,1),ylim=base.c(0,1),main=cohort.upper()+' Cohort ROC\n'+hmm+'.hmm',xlab='% False Pos',ylab='% False Neg')
      # graph 5
      plot(hmm_scores.rx2("Cutoff"),hmm_scores.rx2("Phi"),type='o',xlim=base.c(50,100),ylim=base.c(0,0.75),main=cohort.upper()+' Cohort Association Coeff\n'+hmm+'.hmm',xlab='Cutoff',ylab='Association Coefficient')
   # close pdf file
   grdevices.dev_off()
   return


################## Main ############

# create basic R functions
rprint = robjects.globalenv.get("print")
stats = importr('stats')
grdevices = importr('grDevices')
base = importr('base')
palette = grdevices.palette()
r = robjects.r
plot = robjects.r.plot
ch = []  # names of cohorts
hm = [] # HMM names
statsfiles = []
outfiles = []

# Read in data from external files...seems to attach dataframes automatically
print "Checking directory for .stats and .out files....\n"
for top, middle, bottom in os.walk('.'):
   for filename in bottom:
      if filename.find(".stats") != -1:  statsfiles.append(filename)
      if filename.find(".out") != -1: outfiles.append(filename)
check = 0
if len(statsfiles) > 0: check += 1
if len(outfiles) > 0: check += 2
if check == 0:
   print "There were no .stats or .out files in this directory. Ending program...\n"
   sys.exit()
if check == 1:
   print "There were no .out files in this directory.  Ending program...\n"
   sys.exit()
if check == 2:
   print "THere were no .stats files in this directory.  Ending program...\n"
   sys.exit()

# Load cohort names and hmm names from this directory into arrays
#for x in statsfiles:
#   y = x.split('.')
#   if y[2] not in hm: hm.append(y[2])
for x in outfiles:
   y = x.split('.')
   if y[0] not in ch: ch.append(y[0])
print "Beginning to process data files...\n"
# Begin reading files and doing statistics
for cohort in ch:
   print "Working with cohort "+cohort+" data...\n"
   createGraphSeries(cohort)

print "Completed data analysis...\n"
