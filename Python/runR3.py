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

def createGraphSeries(cohort,t,sfiles,ofile):
   hm = ['bothR5','bothDX','esR5','esDX','trofileR5','trofileDX','mixedR5DX']
   
   # Begin witing graphs to pdf
   grdevices.pdf(file=cohort+".stats.pdf",width=7,height=7)
   
   if t == 'both' or t =='out':
      stat_data = DataFrame.from_csvfile(ofile, sep = "\t")
      # graph 1 scatter plot
      scatterPlot(stat_data,'bothR5','bothDX',cohort.upper()+' cohort\n bothR5/bothDX correlation')
      scatterPlot(stat_data,'esR5','esDX',cohort.upper()+' cohort\n esR5/esDX correlation')
      scatterPlot(stat_data,'trofileR5','trofileDX',cohort.upper()+' cohort\n trofileR5/trofileDX correlation')
      scatterPlot(stat_data,'bothR5','esR5',cohort.upper()+' cohort\n bothR5/esR5 correlation')
      scatterPlot(stat_data,'trofileR5','esR5',cohort.upper()+' cohort\n trofileR5/esR5 correlation')
      scatterPlot(stat_data,'bothDX','esDX',cohort.upper()+' cohort\n bothDX/esDX correlation')
      scatterPlot(stat_data,'trofileDX','esDX',cohort.upper()+' cohort\n trofileDX/esDX correlation')
      scatterPlot(stat_data,'bothDX','mixedR5DX',cohort.upper()+' cohort\n bothDX/mixedR5DX correlation')

   if t == 'both' or t =='stats':
      for hmm in hm:
         if sfiles.has_key(hmm):
            f = sfiles[hmm]
            hmm_scores = DataFrame.from_csvfile(f,sep="\t")
            # graph 3
            plot(hmm_scores.rx2("Cutoff"),hmm_scores.rx2("Accuracy"),type='o',main='Accuracy vs. Cutoff '+cohort.upper()+' Cohort\n'+hmm+'.hmm',xlab='Cutoff',ylab='Accuracy')
            # graph 4
            plot(hmm_scores.rx2("FPP"),hmm_scores.rx2("TPP"),type='o',xlim=base.c(0,1),ylim=base.c(0,1),main=cohort.upper()+' Cohort ROC\n'+hmm+'.hmm',xlab='% False Pos',ylab='% False Neg')
            # graph 5
            plot(hmm_scores.rx2("Cutoff"),hmm_scores.rx2("Phi"),type='o',xlim=base.c(50,100),ylim=base.c(0,0.75),main=cohort.upper()+' Cohort Association Coeff\n'+hmm+'.hmm',xlab='Cutoff',ylab='Association Coefficient')
            # graph 6
            plot(hmm_scores.rx2("Cutoff"),hmm_scores.rx2("Specificity"),type='o',main='Specificity vs. Cutoff '+cohort.upper()+' Cohort\n'+hmm+'.hmm',xlab='Cutoff',ylab='Specificity')
            # graph 7
            plot(hmm_scores.rx2("Cutoff"),hmm_scores.rx2("Sensitivity"),type='o',main='Sensitivity vs. Cutoff '+cohort.upper()+' Cohort\n'+hmm+'.hmm',xlab='Cutoff',ylab='Sensitivity')
            # graph 8
            plot(hmm_scores.rx2("Sensitivity"),hmm_scores.rx2("Specificity"),type='o',main='Sensitivity vs. Specificity '+cohort.upper()+' Cohort\n'+hmm+'.hmm',xlab='Sensitivity',ylab='Specificity')
   
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
sf = {} # statfile names with cohort keys
ouf = {} # outfile names with cohort keys

# Read in data from external files...seems to attach dataframes automatically
print "Checking directory for .stats and .out files....\n"
filelist = os.listdir(os.getcwd())
for filename in filelist:
   if filename.find(".stats") != -1:  
      if len(filename)-6 == filename.find(".stats"): statsfiles.append(filename)
   if filename.find(".out") != -1: 
      if len(filename)-4 == filename.find(".out"): outfiles.append(filename)
check = 0
if len(statsfiles) > 0: check += 1
if len(outfiles) > 0: check += 2
if check == 0:
   print "There were no .stats or .out files in this directory. Ending program...\n"
   sys.exit()
if check == 1:
   print "There were no .out files in this directory.  Some graphs may be missing. Continuing program...\n"
   for x in statsfiles:
      y = x.split('.')
      if y[0] not in ch: 
         ch.append(y[0])
         sf[y[0]] = {y[len(y)-2]:x}
      else: sf[y[0]][y[len(y)-2]] = x
   print "Beginning to process data files...\n"
   print str(sf)
   # Begin reading files and doing statistics
   for cohort in ch:
      print "Working with cohort "+cohort+" data...\n"
      createGraphSeries(cohort,'stats',sf[cohort],'')
if check == 2:
   print "THere were no .stats files in this directory.  Some graphs may be missing.  Continuing program...\n"
   for x in outfiles:
      y = x.split('.')
      if y[0] not in ch: 
         ch.append(y[0])
         ouf[y[0]]=x
   print "Beginning to process data files...\n"
   # Begin reading files and doing statistics
   for cohort in ch:
      print "Working with cohort "+cohort+" data...\n"
      createGraphSeries(cohort,'out','',ouf[cohort])
if check > 2:
   for x in statsfiles:
        y = x.split('.')
        if y[0] not in ch: 
           ch.append(y[0])
           sf[y[0]] = {y[len(y)-2]:x}
        else: sf[y[0]][y[len(y)-2]] = x
   for x in outfiles:
      y = x.split('.')
      if y[0] not in ch: 
         ch.append(y[0])
         ouf[y[0]]=x
   print "Beginning to process data files...\n"
   # Begin reading files and doing statistics
   for cohort in ch:
      print "Working with cohort "+cohort+" data...\n"
      createGraphSeries(cohort,'both',sf[cohort],ouf[cohort])

print "Completed data analysis...\n"
