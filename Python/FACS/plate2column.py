# 384 Well Conversion Utility
# This script takes a 96-well Softmax format file exported from Softmax or Platespitter for 384 well experiment
# and transforms it into a column file format where Plate names are column headers
#
# To Do: Need to combine like plate names/antigens into single column instead of having seperate columns
# To Do: Change the Well column to be Sample columm with format 1: 01A01 ... 96L: 01H12 ... 01: 02A01
#
# by Mark Evans
# 
# Created 08.27.2013
#

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


filename = raw_input("what is the file to parse: ")

f1 = open(filename,'r')
f2 = open(filename[:-4]+".list.txt",'w')

plates = {}
plate = ""
for line in f1.readlines():
    a = line.rstrip().split('\t')
    
    if a[0] == 'Plate:': plate = a[1]

    if len(a) > 1 and a[1] == '0':
        plates[plate] = a[2:14]
    elif len(a) > 1 and a[0] != '~End' and a[0]==a[1]=='':
        plates[plate].extend(a[2:14])

pkeys = plates.keys()
pkeys.sort()
column_values = zip(*map(plates.get, sorted(plates)))
f2.write('Well \t'+'\t'.join(pkeys)+'\n')
c = 1
for vals in column_values:
    f2.write(WELLS['HORIZ'][c]+'\t'+'\t'.join(list(vals))+'\n')
    #f2.write(str(c)+'\t'+'\t'.join(list(vals))+'\n')
    c += 1

f1.close()
f2.close()