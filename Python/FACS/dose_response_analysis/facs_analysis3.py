import os, os.path, sys
import math, datetime
from time import localtime, strftime

import rpy2
from rpy2.robjects.packages import importr
from rpy2 import robjects as ro
from rpy2.robjects import Formula

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# Example of log concentration values used for their dose-response curves.
# Conc = 25,8.33,2.77,0.92,0.30,0.10,0.03,0.01

WELLS = {'HORIZ':{1:'A01',2:'A02',3:'A03',4:'A04',5:'A05',6:'A06',7:'A07',8:'A08',9:'A09',10:'A10',11:'A11',12:'A12',
                  13:'B01',14:'B02',15:'B03',16:'B04',17:'B05',18:'B06',19:'B07',20:'B08',21:'B09',22:'B10',23:'B11',24:'B12',
                  25:'C01',26:'C02',27:'C03',28:'C04',29:'C05',30:'C06',31:'C07',32:'C08',33:'C09',34:'C10',35:'C11',36:'C12',
                  37:'D01',38:'D02',39:'D03',40:'D04',41:'D05',42:'D06',43:'D07',44:'D08',45:'D09',46:'D10',47:'D11',48:'D12',
                  49:'E01',50:'E02',51:'E03',52:'E04',53:'E05',54:'E06',55:'E07',56:'E08',57:'E09',58:'E10',59:'E11',60:'E12',
                  61:'F01',62:'F02',63:'F03',64:'F04',65:'F05',66:'F06',67:'F07',68:'F08',69:'F09',70:'F10',71:'F11',72:'F12',
                  73:'G01',74:'G02',75:'G03',76:'G04',77:'G05',78:'G06',79:'G07',80:'G08',81:'G09',82:'G10',83:'G11',84:'G12',
                  85:'H01',86:'H02',87:'H03',88:'H04',89:'H05',90:'H06',91:'H07',92:'H08',93:'H09',94:'H10',95:'H11',96:'H12'},
         'VERT':{   1:1,2:13,3:25,4:37,5:49,6:61,7:73,8:85,
                    9:2,10:14,11:26,12:38,13:50,14:62,15:74,16:86,
                    17:3,18:15,19:27,20:39,21:51,22:63,23:75,24:87,
                    25:4,26:16,27:28,28:40,29:52,30:64,31:76,32:88,
                    33:5,34:17,35:29,36:41,37:53,38:65,39:77,40:89,
                    41:6,42:18,43:30,44:42,45:54,46:66,47:78,48:90,
                    49:7,50:19,51:31,52:43,53:55,54:67,55:79,56:91,
                    57:8,58:20,59:32,60:44,61:56,62:68,63:80,64:92,
                    65:9,66:21,67:33,68:45,69:57,70:69,71:81,72:93,
                    73:10,74:22,75:34,76:46,77:58,78:70,79:82,80:94,
                    81:11,82:23,83:35,84:47,85:59,86:71,87:83,88:95,
                    89:12,90:24,91:36,92:48,93:60,94:72,95:84,96:96}}



def main():
    
    # Set up constants for this script. These variables will eventually be derived from a user input form
    # but for now are fixed to help development
    ######################################################################################################
    layout = 8
    maxvalue = 0
    multiple = 12  # if layout = 8, multiple = 12 and if layout = 12, multiple=8
    #conc = '25,8.33,2.77,0.92,0.30,0.10,0.03,0.01'  # Values are in ug/ml
    conc = '100,25,6.25,1.56,0.39,0.10,0.02,0.01'  # Values are in nM
    #conditions = {'tranc_130529_40IgGs_Lep_P1.csv':"Lep+",'tranc_130529_40IgGs_NoLep_P1.csv':"Lep-"}
    conditions = {'cleaned_data.csv':"Cond1",}
    clonelist = ['Clone 1','Clone 2','Clone 3','Clone 4','Clone 5','Clone 6','Clone 7','Clone 8','Clone 9','Clone 10','Clone 11','Control']
    #comparisons = [('Median__ss_293F','Median__ss_293F_H8'),('Median__ss_CHOK1','Median__ss_CHO_M8')]
    comparisons = [('cells_CHOk1_Median_FL4_H_FL2_Width','cells_hInsR_Median_FL4_H_FL2_Width')]
    datasets = {}
    antigens = {}
    multiple = 0
    conc2 =[]

    # Determine sample layout pattern, 8-point (vertical) titration or 12-point (horizontal) titration
    if layout == 8: 
        layout = 'VERT'
        multiple = 12
    elif layout == 12: 
        layout = 'HORIZ'
        multiple = 8

    
    # <=== Most likely delete in final version
    #ismore = 'start'
    #while ismore =='start' or ismore=='yes':
    #    filename = raw_input('Enter name of file to analyse: ')
    #    cond = raw_input('What condition is this file (e.g. +Lep) ')
    #    conditions[filename] = cond
    #    check = raw_input('Is there another file? (y/n) ')
    #    if check =='y' or check == 'Y' or check =='yes': ismore = 'yes'
    #    else: ismore = 'no'
    # ======>
    
    
    # Begin data processing by reading and reformating original input file(s) 
    # into a single temp file that is more compatible with R
    ##########################################################################
    tmp = open('tmp.txt','w')
    files = conditions.keys()
    
    # Create list with conc distributed properly to match layout to be used later in temp file
    # If we do this, we don't have to worry so much about getting right index later because values are 'matched' with Conc
    for x in conc.split(','): conc2.extend(((x+',')*multiple).rstrip(',').split(','))
    
    # Get antigens from the column headers of the original data files
    for fn in files:
        f = open(fn,'r')
        for line in f:
            vals = line.rstrip().replace(' ','_').replace(':','_').replace('/','_').replace('-','_').replace('"','').replace('.','_').split(',') # replace funky chars that R might not like later
            if vals[0]=='Sample' and vals[-2:len(vals)-1][0].find('Ratio') == 0 and vals[-1:len(vals)][0].find('Ratio')==0:
                antigens[conditions[fn]] = vals[1:-2]  # Just grab Ags, not sample or Ratio columns
                f.close()
                break
            elif  vals[0]=='Sample':
                antigens[conditions[fn]] = vals[1:]  # failsafe, just grab everyting except first Sample column
                f.close()
                break

    # Parse data and write to a tmp file in modified format
    tmp.write('Conditions,Clones,Conc,'+','.join(antigens[antigens.keys()[0]])+'\n')
    for fn in files:
        f = open(fn,'r')
        c=1  # Had to add this 'fulllength' check because some datasets might not have the two 'Ratio' columns at the end
        for line in f:
            vals = line.rstrip().split(',')
            if vals[0]=='Sample' and vals[-1:len(vals)][0].find('Ratio') != 0: c ='fulllength'
            if vals[0]!='Sample' and vals[0] !='Mean' and vals[0] != 'StdDev':
                idx = int(vals[0].split(': ')[0])
                if layout=='VERT':
                    if datasets.has_key(fn):
                        if c != 'fulllength': datasets[fn][idx] = vals[1:-2]
                        else: datasets[fn][idx] = vals[1:]
                        tmp.write(conditions[fn]+',')
                        tmp.write(clonelist[idx%multiple-1]+',')
                        tmp.write(conc2[idx-1]+',')
                        if c != 'fulllength': tmp.write(','.join(vals[1:-2])+'\n')
                        else: tmp.write(','.join(vals[1:])+'\n')
                    else:
                        if c != 'fulllength': datasets[fn] = {idx:vals[1:-2]}
                        else: datasets[fn] = {idx:vals[1:]}
                        tmp.write(conditions[fn]+',')
                        tmp.write(clonelist[idx%multiple-1]+',')
                        tmp.write(conc2[idx-1]+',')
                        if c != 'fulllength': tmp.write(','.join(vals[1:-2])+'\n')
                        else: tmp.write(','.join(vals[1:])+'\n')
        f.close()
    tmp.close()


    # R functions
    ############################################
    grdevices = importr('grDevices')
    r = ro.r
    plot = r.plot
    drc = importr('drc')

    # Import data from tmp file ito R dataframe object
    clone_data = r['read.csv'](os.path.join(os.getcwd(),'tmp.txt'))

    # Create subsets of clone_data by condition and then by clones
    conds = conditions.values()
    subsets = {'Conditions':{},}
    for c in conds:
        subsets['Conditions'][c] = {'data':'','clones':{}}
        sub = clone_data.rx(clone_data.rx2("Conditions").ro == c, True)
        subsets['Conditions'][c]['data'] = sub
        for clone in clonelist:
            subsets['Conditions'][c]['clones'][clone] = sub.rx(sub.rx2("Clones").ro == clone, True)
    
    # Create aliases for analysis functions in th edrc package because rPY2 doesn't like the () in LL.4()
    L4 = SignatureTranslatedAnonymousPackage("""L4eval<-LL.4(names=c("Slope","Lower Limit","UpperLimit","EC50"))""","L4eval")
    L5 = SignatureTranslatedAnonymousPackage("""L5eval<-LL.5(names=c("Skew","Lower limit","Upper limit","Slope","EC50"))""","L5eval")
    W2 = SignatureTranslatedAnonymousPackage("""W2eval<-W2.4()""","W2eval") # This function doesn't always converge in data, need to trap error if used

    print subsets
    # Perform curve-fit analysis and generate graphs for each parental/transfected antigen pair by clone and condition
    path = "/opt/local/var/media/facstool/"
    graphfile = path+"Rplots.pdf"
    jpg_names =[]
    for d in ('no','yes'):
        if d == 'no': grdevices.pdf(onefile=True,file=graphfile,width=10,height=10)  # writes multiple graphs per single pdf file
        for c in subsets['Conditions']:
            for clone in subsets['Conditions'][c]['clones']:
                for ag1,ag2 in comparisons:
                    if d == 'yes': 
                        grdevices.jpeg(file=path+c+"_"+clone+'_'+ag1+'_'+ag2+'_graph.jpg',width=10,height=10,units='in',res=300)
                        jpg_names.append(path+c+"_"+clone+'_'+ag1+'_'+ag2+'_graph.jpg')
                    # set max value of yaxis from the data, add 500 so it graphs nicely
                    maxvalue = 0  
                    print "c =",c," clone = ",clone," ag1=",ag1," ag2=",ag2  
                    if max(subsets['Conditions'][c]['clones'][clone].rx2(ag1)) >= max(subsets['Conditions'][c]['clones'][clone].rx2(ag2)): 
                        maxvalue = max(subsets['Conditions'][c]['clones'][clone].rx2(ag1))+100
                    else: maxvalue = max(subsets['Conditions'][c]['clones'][clone].rx2(ag2))+100

                    # Begin creating graph
                    #grdevices.pdf(onefile=True,file=c+"_"+ag2+"_"+clone+".pdf",width=10,height=10) # Use this instead to write one file per graph
                                    
                    # Define Y~X, which columns to graph for X and Y by name, since R uses column names
                    formula1 = Formula(ag1+' ~ Conc')
                    formula2 = Formula(ag2+' ~ Conc')
                    
                    # Calculate dose response model
                    p = drc.drm(formula = formula2,data=subsets['Conditions'][c]['clones'][clone], fct=L5.L5eval)
                    p2 = drc.drm(formula = formula1, data=subsets['Conditions'][c]['clones'][clone], fct=L5.L5eval)
                    r.summary(p) # dumb to need to do this line twice, but required to kick summary out of sci. note. mode
                    summary = r.summary(p)
                    summary2 = r.summary(p2)
                    r.ED(p,r.c(10,90))
                    ec = r.ED(p,r.c(10,90))
                    ec2 = r.ED(p2,r.c(10,90))
                
                    
                                    
                    # Set graph margin parameters
                    r.par(mar=r.c(5.1, 6, 4.1, 12), oma=r.c(14,0,0,0), xpd=True, bty='l')
                    
                    # Draw the graphs
                    plot(p,type='all',ylim=r.c(0,maxvalue),main=c+" ("+ag2+" \nvs.\n "+ag1+") "+clone,lwd=2,col="blue",ylab="MFI",log='x')
                    plot(p2,type='all',ylim=r.c(0,maxvalue),add=True,lwd=2,col="red",pch=2,log='x')

                    # Draw the legend
                    r.legend("topright",r.c(ag2,ag1),inset=r.c(-0.34,0.15),pch=r.c(1,2),col=r.c('blue','red'))
                    r.box("figure",col='black')
                    #r.box("outer", lty="solid", col="green")
                    # Add Summary statistic and timestamp
                    r.mtext("Blue Curve", side=1, line=1,at=0.01, adj=0, outer=True, cex=0.75)
                    r.mtext(str(summary), side=1, line=13,at=0.01,adj=0, outer=True, cex=0.75)
                    r.mtext("EC90: "+str(ec[1]), side=1, line=10, at=0.01,adj=0,outer=True, cex=0.75, col='brown')
                    r.mtext("EC10: "+str(ec[0]), side=1, line=10, at=0.2, adj=0,outer=True, cex=0.75,col='brown')
                    r.mtext("Red Curve",  side=1, line=1,at=0.5, adj=0, outer=True, cex=0.75)
                    r.mtext(str(summary2),side=1, line=13,at=0.5, adj=0, outer=True, cex=0.75)
                    r.mtext("EC90: "+str(ec2[1]), side=1, line=10, at=0.5,adj=0,outer=True, cex=0.75, col='brown')
                    r.mtext("EC10: "+str(ec2[0]), side=1, line=10, at=0.7, adj=0,outer=True, cex=0.75,col='brown')
                    r.mtext(r.c(strftime("%d %b %Y %H:%M",localtime())), cex=0.75, line=13, side=1, adj=1, outer=True)

                    # Close individual graph file
                    if d =='yes': grdevices.dev_off()
        if d == 'no': grdevices.dev_off()    # close multi-graph file

    print "\n\n>>> Analysis Complete\n\n"


if __name__ == '__main__':
    main()