import os, os.path, sys
import math, datetime

from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects.lib import grid
from rpy2.robjects.vectors import DataFrame
from rpy2.robjects import Formula
import rpy2.robjects.lib.ggplot2 as ggplot2
import rpy2.rlike.container as rlc
import rpy2


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
    
    #layout  = raw_input('Is this for 8 or 12 point titration? ')
    #conc = raw_input("Enter the concentrations as comma seperated list: ")        
    #conditions ={}
    layout = 8
    conc = '25,8.33,2.77,0.92,0.30,0.10,0.03,0.01'
    conditions = {'tranc_130529_40IgGs_Lep_P1.csv':"+Lep",'tranc_130529_40IgGs_NoLep_P1.csv':"-Lep"}
    clonelist = ['Clone 1','Clone 2','Clone 3','Clone 4','Clone 5','Clone 6','Clone 7','Clone 8','Clone 9','Clone 10','Clone 11','Control']

    #ismore = 'start'
    #while ismore =='start' or ismore=='yes':
    #    filename = raw_input('Enter name of file to analyse: ')
    #    cond = raw_input('What condition is this file (e.g. +Lep) ')
    #    conditions[filename] = cond
    #    check = raw_input('Is there another file? (y/n) ')
    #    if check =='y' or check == 'Y' or check =='yes': ismore = 'yes'
    #    else: ismore = 'no'
    
    files = conditions.keys()
    datasets = {}
    antigens = {}
    for fn in files:
        f = open(fn,'r')
        for line in f:
            vals = line.rstrip().split(',')
            if vals[0]=='Sample' and vals[-2:len(vals)-1][0].find('Ratio') == 0 and vals[-1:len(vals)][0].find('Ratio')==0:
                antigens[conditions[fn]] = vals[1:-2]  # Just grab Ags, not sample or Ratio columns
            elif  vals[0]=='Sample':
                antigens[conditions[fn]] = vals[1:]  # failsafe, just grab everyting except first Sample column
            elif vals[0] !='Mean' and vals[0] != 'StdDev':
                idx = vals[0].split(': ')
                if datasets.has_key(fn):
                    datasets[fn][int(idx[0])] = vals[1:-2]
                else:
                    datasets[fn] = {int(idx[0]):vals[1:-2]}

    print datasets
    if layout == 8: layout = 'VERT'
    elif layout == 12: layout = 'HORIZ'

    # { condition: { clone_num: { ag: [vals] }}}
    clone_data = {}

    for fn in files:
        clone_data[conditions[fn]] = {}
        for clone in clonelist:
            clone_data[conditions[fn]][clone] = {}
            for ag in antigens[conditions[fn]]:
                clone_data[conditions[fn]][clone][ag] = []

        clone_count = 1
        val_keys = datasets[fn].keys()
        val_keys.sort()
        for i in val_keys:
            if layout=='VERT' and i%12 != 0:
                idx = i%12
                data = datasets[fn][i]
                for ag in antigens[conditions[fn]]:
                    clone_data[conditions[fn]][clonelist[idx-1]][ag].append(data[antigens[fn].index(ag)])
                clone_count += 1
            elif layout=='VERT' and i%12 == 0:
                idx = i%12
                data = datasets[fn][i]
                for ag in antigens[conditions[fn]]:
                    clone_data[conditions[fn]][clonelist[idx-1]][ag].append(data[antigens[fn].index(ag)])
                clone_count = 1
            
    tmp = open('tmp.txt','w')
    tmp.write('Condition,Clone,'+','.join(antigens[clone_data.keys[0]])+'\n')
    for cdn in clone_data.keys():
        for cln in clone_analysis[cdn].keys():
            ags = clone_analysis[cdn][cln].keys()
            tmp.write(cdn+','+cln+',')
            ags.sort()
  #          for ag in ags:


    # R functions
   #  rprint = robjects.globalenv.get("print")
   #  grdevices = importr('grDevices')
   #  base = importr('base')
   #  palette = grdevices.palette()
   #  r = robjects.r
   #  plot = robjects.r.plot
   #  drc = importr('drc')

   #  Yvector = []
   #  for x in  clone_data['+Lep']['Control']['Median: ss CHOK1']: Yvector.append(float(x))
   #  Xvector = []
   #  for x in conc.split(','): Xvector.append(float(x))

   #  robjects.FloatVector(Xvector)
   #  robjects.FloatVector(Yvector)
   #  od = rlc.OrdDict([('Concentration',robjects.FloatVector(Xvector)),('MFI',robjects.FloatVector(Yvector))])
   #  grdevices.pdf(file="test.pdf",width=7,height=7)
   #  dataf = robjects.DataFrame(od)

   #  formula = Formula('MFI ~ Concentration')
   #  p = drc.drm(formula = formula, data=od, fct='l5')

   #  gp3 = ggplot2.ggplot(dataf)
   #  pp3 = gp3 + ggplot2.scale_fill_brewer(palette='BrBG',name="Year")+ ggplot2.aes_string(x='MFI',y='Concentration') + ggplot2.opts(title =  " Control MFI")
    
   # # pp3.plot()
   #  #plot(dataf)
   #  plot(p,broken=TRUE)
   #  grdevices.dev_off()

    
 #   for x in clone_data.iterkeys():
 #       print x," : ",clone_data[x]





if __name__ == '__main__':
    main()