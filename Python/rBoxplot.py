from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects.lib import grid
from rpy2.robjects.vectors import DataFrame
import os, os.path, sys
import math, datetime
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.rlike.container as rlc
import rpy2


################## Main ############

# create basic R functions
rprint = robjects.globalenv.get("print")
stats = importr('stats')
grdevices = importr('grDevices')
base = importr('base')
palette = grdevices.palette()
r = robjects.r
plot = robjects.r.plot

f = open('drugfc_sum.txt','r')
dsumFC ={}
dsumY ={}

for line in f.xreadlines():
   a = line.split()
   drug = a[1]
   val = a[2]
   yr = a[4]
  # print str(a)
   try:
      if dsumFC.has_key(drug):
         dsumFC[drug]['Fold_Change'].append(math.log10(float(val)))
         dsumY[drug]['Year'].append(yr)
      else:
         dsumFC[drug]= {'Fold_Change': [math.log10(float(val)),]}
         dsumY[drug]= {'Year': [yr,]}
   except:
      print "FAILURE: dsumFC="+str(dsumFC)+"\n\ndsumY="+str(dsumY)
      sys.exit()
drugs = dsumFC.keys()

for x in drugs:
   if '/' in x:
      grdevices.pdf(file=x.upper().replace('/','_')+".pdf",width=7,height=7)
   else:
      grdevices.pdf(file=x+".pdf",width=7,height=7)
   od = rlc.OrdDict([('Fold_Change',robjects.FloatVector(dsumFC[x]['Fold_Change'])),('Year',robjects.FactorVector(dsumY[x]['Year']))])
   dataf = robjects.DataFrame(od)
   gp3 = ggplot2.ggplot(dataf)
   pp3 = gp3 + ggplot2.scale_fill_brewer(palette='BrBG',name="Year")+ ggplot2.aes_string(x='Year',y='Fold_Change',fill='factor(Year)') +  ggplot2.geom_boxplot() + ggplot2.opts(title =  x+" Yearly Trend")
  # pp3 = gp3 + ggplot2.scale_colour_hue(h=base.c(180,270),name="Year")+ ggplot2.aes_string(x='Year',y='Fold_Change',fill='factor(Year)') +  ggplot2.geom_boxplot() + ggplot2.opts(title =  x+" Yearly Trend")
   #+ ggplot2.scale_y_log10()
   pp3.plot()
   grdevices.dev_off()
   
f.close()
print "\nfinished\n"


