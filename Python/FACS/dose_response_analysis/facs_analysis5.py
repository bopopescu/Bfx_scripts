# facs_analysis5.py
# By Mark Evans
# Created 11.01.2013
# Revised 11.13.2013



import os, os.path, sys, subprocess
import math, datetime
from time import localtime, strftime
from decimal import *

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
# For Quad-well layout 384-well plates
# {HORIZ ABS pos:(394 WELL-H, 96-plates num, 96-horiz-pos, 96-horiz-well)}
QUAD_WELLS_384 = {        1:('A01','P1',1,'A01'),2:('A02','P2',1,'A01'),3:('A03','P1',2,'A02'),4:('A04','P2',2,'A02'),5:('A05','P1',3,'A03'),6:('A06','P2',3,'A03'),7:('A07','P1',4,'A04'),8:('A08','P2',4,'A04'),9:('A09','P1',5,'A05'),10:('A10','P2',5,'A05'),
                    11:('A11','P1',6,'A06'),12:('A12','P2',6,'A06'),13:('A13','P1',7,'A07'),14:('A14','P2',7,'A07'),15:('A15','P1',8,'A08'),16:('A16','P2',8,'A08'),17:('A17','P1',9,'A09'),18:('A18','P2',9,'A09'),19:('A19','P1',10,'A10'),20:('A20','P2',10,'A10'),
                    21:('A21','P1',11,'A11'),22:('A22','P2',11,'A11'),23:('A23','P1',12,'A12'),24:('A24','P2',12,'A12'),25:('B01','P3',1,'A01'),26:('B02','P4',1,'A01'),27:('B03','P3',2,'A02'),28:('B04','P4',2,'A02'),29:('B05','P3',3,'A03'),30:('B06','P4',3,'A03'),
                    31:('B07','P3',4,'A04'),32:('B08','P4',4,'A04'),33:('B09','P3',5,'A05'),34:('B10','P4',5,'A05'),35:('B11','P3',6,'A06'),36:('B12','P4',6,'A06'),37:('B13','P3',7,'A07'),38:('B14','P4',7,'A07'),39:('B15','P3',8,'A08'),40:('B16','P4',8,'A08'),
                    41:('B17','P3',9,'A09'),42:('B18','P4',9,'A09'),43:('B19','P3',10,'A10'),44:('B20','P4',10,'A10'),45:('B21','P3',11,'A11'),46:('B22','P4',11,'A11'),47:('B23','P3',12,'A12'),48:('B24','P4',12,'A12'),49:('C01','P1',13,'B01'),50:('C02','P2',13,'B01'),
                    51:('C03','P1',14,'B02'),52:('C04','P2',14,'B02'),53:('C05','P1',15,'B03'),54:('C06','P2',15,'B03'),55:('C07','P1',16,'B04'),56:('C08','P2',16,'B04'),57:('C09','P1',17,'B05'),58:('C10','P2',17,'B05'),59:('C11','P1',18,'B06'),60:('C12','P2',18,'B06'),
                    61:('C13','P1',19,'B07'),62:('C14','P2',19,'B07'),63:('C15','P1',20,'B08'),64:('C16','P2',20,'B08'),65:('C17','P1',21,'B09'),66:('C18','P2',21,'B09'),67:('C19','P1',22,'B10'),68:('C20','P2',22,'B10'),69:('C21','P1',23,'B11'),70:('C22','P2',23,'B11'),
                    71:('C23','P1',24,'B12'),72:('C24','P2',24,'B12'),73:('D01','P3',13,'B01'),74:('D02','P4',13,'B01'),75:('D03','P3',14,'B02'),76:('D04','P4',14,'B02'),77:('D05','P3',15,'B03'),78:('D06','P4',15,'B03'),79:('D07','P3',16,'B04'),80:('D08','P4',16,'B04'),
                    81:('D09','P3',17,'B05'),82:('D10','P4',17,'B05'),83:('D11','P3',18,'B06'),84:('D12','P4',18,'B06'),85:('D13','P3',19,'B07'),86:('D14','P4',19,'B07'),87:('D15','P3',20,'B08'),88:('D16','P4',20,'B08'),89:('D17','P3',21,'B09'),90:('D18','P4',21,'B09'),
                    91:('D19','P3',22,'B10'),92:('D20','P4',22,'B10'),93:('D21','P3',23,'B11'),94:('D22','P4',23,'B11'),95:('D23','P3',24,'B12'),96:('D24','P4',24,'B12'),97:('E01','P1',25,'C01'),98:('E02','P2',25,'C01'),99:('E03','P1',26,'C02'),100:('E04','P2',26,'C02'),
                    101:('E05','P1',27,'C03'),102:('E06','P2',27,'C03'),103:('E07','P1',28,'C04'),104:('E08','P2',28,'C04'),105:('E09','P1',29,'C05'),106:('E10','P2',29,'C05'),107:('E11','P1',30,'C06'),108:('E12','P2',30,'C06'),109:('E13','P1',31,'C07'),110:('E14','P2',31,'C07'),
                    111:('E15','P1',32,'C08'),112:('E16','P2',32,'C08'),113:('E17','P1',33,'C09'),114:('E18','P2',33,'C09'),115:('E19','P1',34,'C10'),116:('E20','P2',34,'C10'),117:('E21','P1',35,'C11'),118:('E22','P2',35,'C11'),119:('E23','P1',36,'C12'),120:('E24','P2',36,'C12'),
                    121:('F01','P3',25,'C01'),122:('F02','P4',25,'C01'),123:('F03','P3',26,'C02'),124:('F04','P4',26,'C02'),125:('F05','P3',27,'C03'),126:('F06','P4',27,'C03'),127:('F07','P3',28,'C04'),128:('F08','P4',28,'C04'),129:('F09','P3',29,'C05'),130:('F10','P4',29,'C05'),
                    131:('F11','P3',30,'C06'),132:('F12','P4',30,'C06'),133:('F13','P3',31,'C07'),134:('F14','P4',31,'C07'),135:('F15','P3',32,'C08'),136:('F16','P4',32,'C08'),137:('F17','P3',33,'C09'),138:('F18','P4',33,'C09'),139:('F19','P3',34,'C10'),140:('F20','P4',34,'C10'),
                    141:('F21','P3',35,'C11'),142:('F22','P4',35,'C11'),143:('F23','P3',36,'C12'),144:('F24','P4',36,'C12'),145:('G01','P1',37,'D01'),146:('G02','P2',37,'D01'),147:('G03','P1',38,'D02'),148:('G04','P2',38,'D02'),149:('G05','P1',39,'D03'),150:('G06','P2',39,'D03'),
                    151:('G07','P1',40,'D04'),152:('G08','P2',40,'D04'),153:('G09','P1',41,'D05'),154:('G10','P2',41,'D05'),155:('G11','P1',42,'D06'),156:('G12','P2',42,'D06'),157:('G13','P1',43,'D07'),158:('G14','P2',43,'D07'),159:('G15','P1',44,'D08'),160:('G16','P2',44,'D08'),
                    161:('G17','P1',45,'D09'),162:('G18','P2',45,'D09'),163:('G19','P1',46,'D10'),164:('G20','P2',46,'D10'),165:('G21','P1',47,'D11'),166:('G22','P2',47,'D11'),167:('G23','P1',48,'D12'),168:('G24','P2',48,'D12'),169:('H01','P3',37,'D01'),170:('H02','P4',37,'D01'),
                    171:('H03','P3',38,'D02'),172:('H04','P4',38,'D02'),173:('H05','P3',39,'D03'),174:('H06','P4',39,'D03'),175:('H07','P3',40,'D04'),176:('H08','P4',40,'D04'),177:('H09','P3',41,'D05'),178:('H10','P4',41,'D05'),179:('H11','P3',42,'D06'),180:('H12','P4',42,'D06'),
                    181:('H13','P3',43,'D07'),182:('H14','P4',43,'D07'),183:('H15','P3',44,'D08'),184:('H16','P4',44,'D08'),185:('H17','P3',45,'D09'),186:('H18','P4',45,'D09'),187:('H19','P3',46,'D10'),188:('H20','P4',46,'D10'),189:('H21','P3',47,'D11'),190:('H22','P4',47,'D11'),
                    191:('H23','P3',48,'D12'),192:('H24','P4',48,'D12'),193:('I01','P1',49,'E01'),194:('I02','P2',49,'E01'),195:('I03','P1',50,'E02'),196:('I04','P2',50,'E02'),197:('I05','P1',51,'E03'),198:('I06','P2',51,'E03'),199:('I07','P1',52,'E04'),200:('I08','P2',52,'E04'),
                    201:('I09','P1',53,'E05'),202:('I10','P2',53,'E05'),203:('I11','P1',54,'E06'),204:('I12','P2',54,'E06'),205:('I13','P1',55,'E07'),206:('I14','P2',55,'E07'),207:('I15','P1',56,'E08'),208:('I16','P2',56,'E08'),209:('I17','P1',57,'E09'),210:('I18','P2',57,'E09'),
                    211:('I19','P1',58,'E10'),212:('I20','P2',58,'E10'),213:('I21','P1',59,'E11'),214:('I22','P2',59,'E11'),215:('I23','P1',60,'E12'),216:('I24','P2',60,'E12'),217:('J01','P3',49,'E01'),218:('J02','P4',49,'E01'),219:('J03','P3',50,'E02'),220:('J04','P4',50,'E02'),
                    221:('J05','P3',51,'E03'),222:('J06','P4',51,'E03'),223:('J07','P3',52,'E04'),224:('J08','P4',52,'E04'),225:('J09','P3',53,'E05'),226:('J10','P4',53,'E05'),227:('J11','P3',54,'E06'),228:('J12','P4',54,'E06'),229:('J13','P3',55,'E07'),230:('J14','P4',55,'E07'),
                    231:('J15','P3',56,'E08'),232:('J16','P4',56,'E08'),233:('J17','P3',57,'E09'),234:('J18','P4',57,'E09'),235:('J19','P3',58,'E10'),236:('J20','P4',58,'E10'),237:('J21','P3',59,'E11'),238:('J22','P4',59,'E11'),239:('J23','P3',60,'E12'),240:('J24','P4',60,'E12'),
                    241:('K01','P1',61,'F01'),242:('K02','P2',61,'F01'),243:('K03','P1',62,'F02'),244:('K04','P2',62,'F02'),245:('K05','P1',63,'F03'),246:('K06','P2',63,'F03'),247:('K07','P1',64,'F04'),248:('K08','P2',64,'F04'),249:('K09','P1',65,'F05'),250:('K10','P2',65,'F05'),
                    251:('K11','P1',66,'F06'),252:('K12','P2',66,'F06'),253:('K13','P1',67,'F07'),254:('K14','P2',67,'F07'),255:('K15','P1',68,'F08'),256:('K16','P2',68,'F08'),257:('K17','P1',69,'F09'),258:('K18','P2',69,'F09'),259:('K19','P1',70,'F10'),260:('K20','P2',70,'F10'),
                    261:('K21','P1',71,'F11'),262:('K22','P2',71,'F11'),263:('K23','P1',72,'F12'),264:('K24','P2',72,'F12'),265:('L01','P3',61,'F01'),266:('L02','P4',61,'F01'),267:('L03','P3',62,'F02'),268:('L04','P4',62,'F02'),269:('L05','P3',63,'F03'),270:('L06','P4',63,'F03'),
                    271:('L07','P3',64,'F04'),272:('L08','P4',64,'F04'),273:('L09','P3',65,'F05'),274:('L10','P4',65,'F05'),275:('L11','P3',66,'F06'),276:('L12','P4',66,'F06'),277:('L13','P3',67,'F07'),278:('L14','P4',67,'F07'),279:('L15','P3',68,'F08'),280:('L16','P4',68,'F08'),
                    281:('L17','P3',69,'F09'),282:('L18','P4',69,'F09'),283:('L19','P3',70,'F10'),284:('L20','P4',70,'F10'),285:('L21','P3',71,'F11'),286:('L22','P4',71,'F11'),287:('L23','P3',72,'F12'),288:('L24','P4',72,'F12'),289:('M01','P1',73,'G01'),290:('M02','P2',73,'G01'),
                    291:('M03','P1',74,'G02'),292:('M04','P2',74,'G02'),293:('M05','P1',75,'G03'),294:('M06','P2',75,'G03'),295:('M07','P1',76,'G04'),296:('M08','P2',76,'G04'),297:('M09','P1',77,'G05'),298:('M10','P2',77,'G05'),299:('M11','P1',78,'G06'),300:('M12','P2',78,'G06'),
                    301:('M13','P1',79,'G07'),302:('M14','P2',79,'G07'),303:('M15','P1',80,'G08'),304:('M16','P2',80,'G08'),305:('M17','P1',81,'G09'),306:('M18','P2',81,'G09'),307:('M19','P1',82,'G10'),308:('M20','P2',82,'G10'),309:('M21','P1',83,'G11'),310:('M22','P2',83,'G11'),
                    311:('M23','P1',84,'G12'),312:('M24','P2',84,'G12'),313:('N01','P3',73,'G01'),314:('N02','P4',73,'G01'),315:('N03','P3',74,'G02'),316:('N04','P4',74,'G02'),317:('N05','P3',75,'G03'),318:('N06','P4',75,'G03'),319:('N07','P3',76,'G04'),320:('N08','P4',76,'G04'),
                    321:('N09','P3',77,'G05'),322:('N10','P4',77,'G05'),323:('N11','P3',78,'G06'),324:('N12','P4',78,'G06'),325:('N13','P3',79,'G07'),326:('N14','P4',79,'G07'),327:('N15','P3',80,'G08'),328:('N16','P4',80,'G08'),329:('N17','P3',81,'G09'),330:('N18','P4',81,'G09'),
                    331:('N19','P3',82,'G10'),332:('N20','P4',82,'G10'),333:('N21','P3',83,'G11'),334:('N22','P4',83,'G11'),335:('N23','P3',84,'G12'),336:('N24','P4',84,'G12'),337:('O01','P1',85,'H01'),338:('O02','P2',85,'H01'),339:('O03','P1',86,'H02'),340:('O04','P2',86,'H02'),
                    341:('O05','P1',87,'H03'),342:('O06','P2',87,'H03'),343:('O07','P1',88,'H04'),344:('O08','P2',88,'H04'),345:('O09','P1',89,'H05'),346:('O10','P2',89,'H05'),347:('O11','P1',90,'H06'),348:('O12','P2',90,'H06'),349:('O13','P1',91,'H07'),350:('O14','P2',91,'H07'),
                    351:('O15','P1',92,'H08'),352:('O16','P2',92,'H08'),353:('O17','P1',93,'H09'),354:('O18','P2',93,'H09'),355:('O19','P1',94,'H10'),356:('O20','P2',94,'H10'),357:('O21','P1',95,'H11'),358:('O22','P2',95,'H11'),359:('O23','P1',96,'H12'),360:('O24','P2',96,'H12'),
                    361:('P01','P3',85,'H01'),362:('P02','P4',85,'H01'),363:('P03','P3',86,'H02'),364:('P04','P4',86,'H02'),365:('P05','P3',87,'H03'),366:('P06','P4',87,'H03'),367:('P07','P3',88,'H04'),368:('P08','P4',88,'H04'),369:('P09','P3',89,'H05'),370:('P10','P4',89,'H05'),
                    371:('P11','P3',90,'H06'),372:('P12','P4',90,'H06'),373:('P13','P3',91,'H07'),374:('P14','P4',91,'H07'),375:('P15','P3',92,'H08'),376:('P16','P4',92,'H08'),377:('P17','P3',93,'H09'),378:('P18','P4',93,'H09'),379:('P19','P3',94,'H10'),380:('P20','P4',94,'H10'),
                    381:('P21','P3',95,'H11'),382:('P22','P4',95,'H11'),383:('P23','P3',96,'H12'),
                    384:('P24','P4',96,'H12')}

# For linear 96-well plate layout 384-well plates
LINEAR_WELLS_384 = { 1:('A01','P1',1,'A01'),2:('A02','P1',2,'A02'),3:('A03','P1',3,'A03'),4:('A04','P1',4,'A04'),5:('A05','P1',5,'A05'),6:('A06','P1',6,'A06'),7:('A07','P1',7,'A07'),8:('A08','P1',8,'A08'),9:('A09','P1',9,'A09'),
                    10:('A10','P1',10,'A10'),11:('A11','P1',11,'A11'),12:('A12','P1',12,'A12'),13:('A13','P2',1,'A01'),14:('A14','P2',2,'A02'),15:('A15','P2',3,'A03'),16:('A16','P2',4,'A04'),17:('A17','P2',5,'A05'),18:('A18','P2',6,'A06'),19:('A19','P2',7,'A07'),
                    20:('A20','P2',8,'A08'),21:('A21','P2',9,'A09'),22:('A22','P2',10,'A10'),23:('A23','P2',11,'A11'),24:('A24','P2',12,'A12'),25:('B01','P1',1,'B01'),26:('B02','P1',2,'B02'),27:('B03','P1',3,'B03'),28:('B04','P1',4,'B04'),29:('B05','P1',5,'B05'),
                    30:('B06','P1',6,'B06'),31:('B07','P1',7,'B07'),32:('B08','P1',8,'B08'),33:('B09','P1',9,'B09'),34:('B10','P1',10,'B10'),35:('B11','P1',11,'B11'),36:('B12','P1',12,'B12'),37:('B13','P2',1,'B01'),38:('B14','P2',2,'B02'),39:('B15','P2',3,'B03'),
                    40:('B16','P2',4,'B04'),41:('B17','P2',5,'B05'),42:('B18','P2',6,'B06'),43:('B19','P2',7,'B07'),44:('B20','P2',8,'B08'),45:('B21','P2',9,'B09'),46:('B22','P2',10,'B10'),47:('B23','P2',11,'B11'),48:('B24','P2',12,'B12'),49:('C01','P1',1,'C01'),
                    50:('C02','P1',2,'C02'),51:('C03','P1',3,'C03'),52:('C04','P1',4,'C04'),53:('C05','P1',5,'C05'),54:('C06','P1',6,'C06'),55:('C07','P1',7,'C07'),56:('C08','P1',8,'C08'),57:('C09','P1',9,'C09'),58:('C10','P1',10,'C10'),59:('C11','P1',11,'C11'),
                    60:('C12','P1',12,'C12'),61:('C13','P2',1,'C01'),62:('C14','P2',2,'C02'),63:('C15','P2',3,'C03'),64:('C16','P2',4,'C04'),65:('C17','P2',5,'C05'),66:('C18','P2',6,'C06'),67:('C19','P2',7,'C07'),68:('C20','P2',8,'C08'),69:('C21','P2',9,'C09'),
                    70:('C22','P2',10,'C10'),71:('C23','P2',11,'C11'),72:('C24','P2',12,'C12'),73:('D01','P1',1,'D01'),74:('D02','P1',2,'D02'),75:('D03','P1',3,'D03'),76:('D04','P1',4,'D04'),77:('D05','P1',5,'D05'),78:('D06','P1',6,'D06'),79:('D07','P1',7,'D07'),
                    80:('D08','P1',8,'D08'),81:('D09','P1',9,'D09'),82:('D10','P1',10,'D10'),83:('D11','P1',11,'D11'),84:('D12','P1',12,'D12'),85:('D13','P2',1,'D01'),86:('D14','P2',2,'D02'),87:('D15','P2',3,'D03'),88:('D16','P2',4,'D04'),89:('D17','P2',5,'D05'),
                    90:('D18','P2',6,'D06'),91:('D19','P2',7,'D07'),92:('D20','P2',8,'D08'),93:('D21','P2',9,'D09'),94:('D22','P2',10,'D10'),95:('D23','P2',11,'D11'),96:('D24','P2',12,'D12'),97:('E01','P1',1,'E01'),98:('E02','P1',2,'E02'),99:('E03','P1',3,'E03'),
                    100:('E04','P1',4,'E04'),101:('E05','P1',5,'E05'),102:('E06','P1',6,'E06'),103:('E07','P1',7,'E07'),104:('E08','P1',8,'E08'),105:('E09','P1',9,'E09'),106:('E10','P1',10,'E10'),107:('E11','P1',11,'E11'),108:('E12','P1',12,'E12'),109:('E13','P2',1,'E01'),
                    110:('E14','P2',2,'E02'),111:('E15','P2',3,'E03'),112:('E16','P2',4,'E04'),113:('E17','P2',5,'E05'),114:('E18','P2',6,'E06'),115:('E19','P2',7,'E07'),116:('E20','P2',8,'E08'),117:('E21','P2',9,'E09'),118:('E22','P2',10,'E10'),119:('E23','P2',11,'E11'),
                    120:('E24','P2',12,'E12'),121:('F01','P1',1,'F01'),122:('F02','P1',2,'F02'),123:('F03','P1',3,'F03'),124:('F04','P1',4,'F04'),125:('F05','P1',5,'F05'),126:('F06','P1',6,'F06'),127:('F07','P1',7,'F07'),128:('F08','P1',8,'F08'),129:('F09','P1',9,'F09'),
                    130:('F10','P1',10,'F10'),131:('F11','P1',11,'F11'),132:('F12','P1',12,'F12'),133:('F13','P2',1,'F01'),134:('F14','P2',2,'F02'),135:('F15','P2',3,'F03'),136:('F16','P2',4,'F04'),137:('F17','P2',5,'F05'),138:('F18','P2',6,'F06'),139:('F19','P2',7,'F07'),
                    140:('F20','P2',8,'F08'),141:('F21','P2',9,'F09'),142:('F22','P2',10,'F10'),143:('F23','P2',11,'F11'),144:('F24','P2',12,'F12'),145:('G01','P1',1,'G01'),146:('G02','P1',2,'G02'),147:('G03','P1',3,'G03'),148:('G04','P1',4,'G04'),149:('G05','P1',5,'G05'),
                    150:('G06','P1',6,'G06'),151:('G07','P1',7,'G07'),152:('G08','P1',8,'G08'),153:('G09','P1',9,'G09'),154:('G10','P1',10,'G10'),155:('G11','P1',11,'G11'),156:('G12','P1',12,'G12'),157:('G13','P2',1,'G01'),158:('G14','P2',2,'G02'),159:('G15','P2',3,'G03'),
                    160:('G16','P2',4,'G04'),161:('G17','P2',5,'G05'),162:('G18','P2',6,'G06'),163:('G19','P2',7,'G07'),164:('G20','P2',8,'G08'),165:('G21','P2',9,'G09'),166:('G22','P2',10,'G10'),167:('G23','P2',11,'G11'),168:('G24','P2',12,'G12'),169:('H01','P1',1,'H01'),
                    170:('H02','P1',2,'H02'),171:('H03','P1',3,'H03'),172:('H04','P1',4,'H04'),173:('H05','P1',5,'H05'),174:('H06','P1',6,'H06'),175:('H07','P1',7,'H07'),176:('H08','P1',8,'H08'),177:('H09','P1',9,'H09'),178:('H10','P1',10,'H10'),179:('H11','P1',11,'H11'),
                    180:('H12','P1',12,'H12'),181:('H13','P2',1,'H01'),182:('H14','P2',2,'H02'),183:('H15','P2',3,'H03'),184:('H16','P2',4,'H04'),185:('H17','P2',5,'H05'),186:('H18','P2',6,'H06'),187:('H19','P2',7,'H07'),188:('H20','P2',8,'H08'),189:('H21','P2',9,'H09'),
                    190:('H22','P2',10,'H10'),191:('H23','P2',11,'H11'),192:('H24','P2',12,'H12'),193:('I01','P3',1,'A01'),194:('I02','P3',2,'A02'),195:('I03','P3',3,'A03'),196:('I04','P3',4,'A04'),197:('I05','P3',5,'A05'),198:('I06','P3',6,'A06'),199:('I07','P3',7,'A07'),
                    200:('I08','P3',8,'A08'),201:('I09','P3',9,'A09'),202:('I10','P3',10,'A10'),203:('I11','P3',11,'A11'),204:('I12','P3',12,'A12'),205:('I13','P4',1,'A01'),206:('I14','P4',2,'A02'),207:('I15','P4',3,'A03'),208:('I16','P4',4,'A04'),209:('I17','P4',5,'A05'),
                    210:('I18','P4',6,'A06'),211:('I19','P4',7,'A07'),212:('I20','P4',8,'A08'),213:('I21','P4',9,'A09'),214:('I22','P4',10,'A10'),215:('I23','P4',11,'A11'),216:('I24','P4',12,'A12'),217:('J01','P3',1,'B01'),218:('J02','P3',2,'B02'),219:('J03','P3',3,'B03'),
                    220:('J04','P3',4,'B04'),221:('J05','P3',5,'B05'),222:('J06','P3',6,'B06'),223:('J07','P3',7,'B07'),224:('J08','P3',8,'B08'),225:('J09','P3',9,'B09'),226:('J10','P3',10,'B10'),227:('J11','P3',11,'B11'),228:('J12','P3',12,'B12'),229:('J13','P4',1,'B01'),
                    230:('J14','P4',2,'B02'),231:('J15','P4',3,'B03'),232:('J16','P4',4,'B04'),233:('J17','P4',5,'B05'),234:('J18','P4',6,'B06'),235:('J19','P4',7,'B07'),236:('J20','P4',8,'B08'),237:('J21','P4',9,'B09'),238:('J22','P4',10,'B10'),239:('J23','P4',11,'B11'),
                    240:('J24','P4',12,'B12'),241:('K01','P3',1,'C01'),242:('K02','P3',2,'C02'),243:('K03','P3',3,'C03'),244:('K04','P3',4,'C04'),245:('K05','P3',5,'C05'),246:('K06','P3',6,'C06'),247:('K07','P3',7,'C07'),248:('K08','P3',8,'C08'),249:('K09','P3',9,'C09'),
                    250:('K10','P3',10,'C10'),251:('K11','P3',11,'C11'),252:('K12','P3',12,'C12'),253:('K13','P4',1,'C01'),254:('K14','P4',2,'C02'),255:('K15','P4',3,'C03'),256:('K16','P4',4,'C04'),257:('K17','P4',5,'C05'),258:('K18','P4',6,'C06'),259:('K19','P4',7,'C07'),
                    260:('K20','P4',8,'C08'),261:('K21','P4',9,'C09'),262:('K22','P4',10,'C10'),263:('K23','P4',11,'C11'),264:('K24','P4',12,'C12'),265:('L01','P3',1,'D01'),266:('L02','P3',2,'D02'),267:('L03','P3',3,'D03'),268:('L04','P3',4,'D04'),269:('L05','P3',5,'D05'),
                    270:('L06','P3',6,'D06'),271:('L07','P3',7,'D07'),272:('L08','P3',8,'D08'),273:('L09','P3',9,'D09'),274:('L10','P3',10,'D10'),275:('L11','P3',11,'D12'),276:('L12','P3',12,'D12'),277:('L13','P4',1,'D01'),278:('L14','P4',2,'D02'),279:('L15','P4',3,'D03'),
                    280:('L16','P4',4,'D04'),281:('L17','P4',5,'D05'),282:('L18','P4',6,'D06'),283:('L19','P4',7,'D07'),284:('L20','P4',8,'D08'),285:('L21','P4',9,'D09'),286:('L22','P4',10,'D10'),287:('L23','P4',11,'D11'),288:('L24','P4',12,'D12'),289:('M01','P3',1,'E01'),
                    290:('M02','P3',2,'E02'),291:('M03','P3',3,'E03'),292:('M04','P3',4,'E04'),293:('M05','P3',5,'E05'),294:('M06','P3',6,'E06'),295:('M07','P3',7,'E07'),296:('M08','P3',8,'E08'),297:('M09','P3',9,'E09'),298:('M10','P3',10,'E10'),299:('M11','P3',11,'E11'),
                    300:('M12','P3',12,'E12'),301:('M13','P4',1,'E01'),302:('M14','P4',2,'E02'),303:('M15','P4',3,'E03'),304:('M16','P4',4,'E04'),305:('M17','P4',5,'E05'),306:('M18','P4',6,'E06'),307:('M19','P4',7,'E07'),308:('M20','P4',8,'E08'),309:('M21','P4',9,'E09'),
                    310:('M22','P4',10,'E10'),311:('M23','P4',11,'E11'),312:('M24','P4',12,'E12'),313:('N01','P3',1,'F01'),314:('N02','P3',2,'F02'),315:('N03','P3',3,'F03'),316:('N04','P3',4,'F04'),317:('N05','P3',5,'F05'),318:('N06','P3',6,'F06'),319:('N07','P3',7,'F07'),
                    320:('N08','P3',8,'F08'),321:('N09','P3',9,'F09'),322:('N10','P3',10,'F10'),323:('N11','P3',11,'F11'),324:('N12','P3',12,'F12'),325:('N13','P3',1,'F01'),326:('N14','P4',2,'F02'),327:('N15','P4',3,'F03'),328:('N16','P4',4,'F04'),329:('N17','P4',5,'F05'),
                    330:('N18','P4',6,'F06'),331:('N19','P4',7,'F07'),332:('N20','P4',8,'F08'),333:('N21','P4',9,'F09'),334:('N22','P4',10,'F10'),335:('N23','P4',11,'F11'),336:('N24','P4',12,'F12'),337:('O01','P3',1,'G01'),338:('O02','P3',2,'G02'),339:('O03','P3',3,'G03'),
                    340:('O04','P3',4,'G04'),341:('O05','P3',5,'G05'),342:('O06','P3',6,'G06'),343:('O07','P3',7,'G07'),344:('O08','P3',8,'G08'),345:('O09','P3',9,'G09'),346:('O10','P3',10,'G10'),347:('O11','P3',11,'G11'),348:('O12','P3',12,'G12'),349:('O13','P4',1,'G01'),
                    350:('O14','P4',2,'G02'),351:('O15','P4',3,'G03'),352:('O16','P4',4,'G04'),353:('O17','P4',5,'G05'),354:('O18','P4',6,'G06'),355:('O19','P4',7,'G07'),356:('O20','P4',8,'G08'),357:('O21','P4',9,'G09'),358:('O22','P4',10,'G10'),359:('O23','P4',11,'G11'),
                    360:('O24','P4',12,'G12'),361:('P01','P3',1,'H01'),362:('P02','P3',2,'H02'),363:('P03','P3',3,'H03'),364:('P04','P3',4,'H04'),365:('P05','P3',5,'H05'),366:('P06','P3',6,'H06'),367:('P07','P3',7,'H07'),368:('P08','P3',8,'H08'),369:('P09','P3',9,'H09'),
                    370:('P10','P3',10,'H10'),371:('P11','P3',11,'H11'),372:('P12','P3',12,'H12'),373:('P13','P4',1,'H01'),374:('P14','P4',2,'H02'),375:('P15','P4',3,'H03'),376:('P16','P4',4,'H04'),377:('P17','P4',5,'H05'),378:('P18','P4',6,'H06'),379:('P19','P4',7,'H07'),
                    380:('P20','P4',8,'H08'),381:('P21','P4',9,'H09'),382:('P22','P4',10,'H10'),383:('P23','P4',11,'H11'), 384:('P24','P4',12,'H12')}


#######################################################################
# split384To96                                                        #
# Splits well-quad format 384 plate into component 96-well plates 1-4 #
# input: 384 well FACS data in list format, 1 col per Ag              #
# return: dict plates[Ag][96-well plate num][data value]              #
#######################################################################
def split384To96(fn, plate_layout):
    wellmap = ""
    if plate_layout == 'linear': wellmap = LINEAR_WELLS_384
    elif plate_layout == 'quad': wellmap = QUAD_WELLS_384

    # Read in source file and parse
    plates = {}
    ratio_pos = 0
    ags =[]

    f1 = open(fn,'r')
    for line in f1.readlines():
        a = line.rstrip().split(',')
        if a[0] == 'Sample': 
            try: ratio_pos = a.index('Ratio')    # Want to exclude this column, so need to know where it is
            except ValueError: ratio_pos = 0
            if ratio_pos != 0: ags = a[1:ratio_pos]            # Get list of antigens from the column headers
            else: ags = a[1:]
            for x in ags: plates[x] = {'P1':[],'P2':[],'P3':[],'P4':[]}     # Preload dictionary with four 96-well plates
        
        # Handle data lines    
        if a[0] != 'Sample' and a[0] != 'StdDev' and a[0] != 'Mean':
            b = int(a[0].split(':')[0])     # Get out 384 well abs pos index
            
            # Use zip to pair Antigen and value on the fly
            if ratio_pos != 0:
                for val,ag in zip(a[1:ratio_pos],ags):
                    plates[ag][wellmap[b][1]].append(val)
            else:
                for val,ag in zip(a[1:],ags):
                    plates[ag][wellmap[b][1]].append(val)
    return (plates,ags)


##################################################################
# merge2List                                                     #
# Recombine split 384 column data back into list of 96-well data #
##################################################################
def merge2List(fn,data):
    # Combine individual 96-well plates into single list organized by antigen
    column_values =[]
    merged_plates = {}
    plates = data[0]
    ags = data[1]
    
    for ag in ags: 
        merged_plates[ag] = []
        for p in ('P1','P2','P3','P4'):
            merged_plates[ag].extend(plates[ag][p])
            
    # Create list of tuples containing all values across the columns for writing
    column_values = zip(*map(merged_plates.get, sorted(merged_plates)))
    
    # Create output
    output_fn = fn[:-4]+".list.txt"
    f2 = open(output_fn,'w')
    f2.write('Sample,'+','.join(sorted(ags))+'\n')
    c = 1
    well_count = 1
    plate_count = 1
    for vals in column_values:
        if well_count > 96: 
            well_count = 1
            plate_count +=1
        # Using sample ID format from FACScan so XabTracker will be happy
        # e.g. 01A01 ... 96L: 01H12 ... 01: 02A01
        sample = str(well_count)+': 0'+str(plate_count)+WELLS['HORIZ'][well_count]
        f2.write(sample+','+','.join(list(vals))+'\n')
        c += 1
        well_count += 1
    f2.close()
    return output_fn
    

############################################
# getAgList                                #
# Get list of Antigens from column headers #
############################################
def getAgList(files,conditions):
    antigens = {}
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
    return antigens


##########################################################################
# createDataset                                                          #
# Recombine all parsed data and antigens, conditions and clones          #
# into new data file to be used as input file for analysis               #
# This will work for duplicate samples as long as the user specifies all #
# of the duplicate names,e.g. Clone 1, Clone 1, Clone 2, Clone 2, etc    #
##########################################################################
def createDataset(files,conditions,antigens,clonelist,multiple,conc2):
    # Parse data and write to a tmp file in modified format
    tmp = open('tmp.txt','w')
    datasets = {}
    clone_set = zip(*[iter(clonelist)]*multiple)
    set_count = len(clone_set)  # counts how many plates in the files
    set_idx = 0  # This keeps track of # of 96-well plates, needed for 384-well plates
    idx_sum = 1
    tmp.write('Conditions,Clones,Conc,'+','.join(antigens[antigens.keys()[0]])+'\n')
    import pdb
    pdb.set_trace()

    for fn in files:
        f = open(fn,'r')
        c=1  # Had to add this 'fulllength' check because some datasets might not have the two 'Ratio' columns at the end
        for line in f:
            vals = line.rstrip().split(',')

            if vals[0]=='Sample' and vals[-1:len(vals)][0].find('Ratio') != 0: c ='fulllength'  # No ratio column, data cols only
            if vals[0]!='Sample' and vals[0] !='Mean' and vals[0] != 'StdDev':
                idx = int(vals[0].split(': ')[0])   # abs well pos
                if idx_sum == 96 and idx == 1: set_idx += 1               
                else: idx_sum = idx

                if datasets.has_key(fn):
                    if c != 'fulllength': datasets[fn][idx] = vals[1:-2]
                    else: datasets[fn][idx] = vals[1:]
                    
                    tmp.write(conditions[fn]+',')
                    tmp.write(clone_set[set_idx][idx%multiple-1]+',')
                    tmp.write(conc2[idx-1]+',')
                    
                    if c != 'fulllength': tmp.write(','.join(vals[1:-2])+'\n')
                    else: tmp.write(','.join(vals[1:])+'\n')
                    
                else:
                    if c != 'fulllength': datasets[fn] = {idx:vals[1:-2]}
                    else: datasets[fn] = {idx:vals[1:]}
                    
                    tmp.write(conditions[fn]+',')
                    tmp.write(clone_set[set_idx][idx%multiple-1]+',')
                    tmp.write(conc2[idx-1]+',')
                    
                    if c != 'fulllength': tmp.write(','.join(vals[1:-2])+'\n')
                    else: tmp.write(','.join(vals[1:])+'\n')
        f.close()
    tmp.close()
    return datasets, tmp.name


#########################################################################################################################
# createConcRange                                                                                                       #
# Create list with conc distributed properly to match layout to be used later in temp file                              #
# If we do this, we don't have to worry so much about getting right index later because values are 'matched' with Conc  #
# conc = string of comma-seperated concentration values                                                                 #
# multiple = 8 or 12, depending on layout value, indicating if vertical or horzontal samples                            #
# @rvalue = list of concentration values that have been properly repeated for the whole sample range on plate           #
#########################################################################################################################
def createConcRange(conc,multiple):
    conc2 = []
    for x in conc.replace(' ','').split(','): conc2.extend(((x+',')*multiple).rstrip(',').split(','))
    return conc2


########################################################################################################
# getPlatePattern                                                                                      #
# Return plate layout pattern and the multiple required to properly 'break' well positions by clone    #
# layout = 8 or 12 indicating 8-point titrattion (vertical layout) or 12-point titration (horz layout) #
# @rvalue = numeric multiplier and text indicating layout                                              #
########################################################################################################
def getPlatePattern(layout):
    multiple = 0
    if layout == '8': multiple = 12
    elif layout == '12': multiple = 8
    return multiple


####################################################################
# calculateModelFit                                                #
# Attempts to calculate DRC model using user-specified model first #
# If that fails to converge, it tries a 4-param fit. If that fails #
# it returns failure and the points will simply plotted            #
####################################################################
def calculateModelFit(formula,data, fct):
    R = ro.r
    drc = importr('drc')
    try:
        p = drc.drm(formula = formula,data=data, fct=fct)
        R.summary(p) # dumb to need to do this line twice, but required to kick summary out of sci. note. mode
        summary = R.summary(p)
        R.ED(p, R.c(10,50,90))
        ec = R.ED(p, R.c(10,50,90))
        return p,summary,ec
    except:
        L4 = SignatureTranslatedAnonymousPackage("""L4eval<-LL.4(names=c("Slope","Lower Limit","UpperLimit","EC50"))""","L4eval")
        try:
            p = drc.drm(formula = formula,data=data, fct=L4.L4eval)
            R.summary(p) # dumb to need to do this line twice, but required to kick summary out of sci. note. mode
            summary = R.summary(p)
            R.ED(p, R.c(10,50,90))
            ec = R.ED(p, R.c(10,50,90))
            return p,summary,ec
        except:
            #import pdb
            #pdb.set_trace()
            print "===================== Failed to converge ===============================\n"
            print formula
            print data
            print "========================================================================\n"
            return "Error",'Failed to converge',''


##########################################################
# parseSummary                                           #
# Parse the summary object that comes back from R        #
# Extract out upper/lower limits (or whatever they want) #
##########################################################
def parseSummary(summary):
    so = {}
    s = str(summary).split()
    var1 = var2 = var3 = var4 ="NA"
    
    if len(s) == 51:
        if s[23] != 'NA': 
            if s[23] != '0': var1 = str(Decimal(s[23]))[:str(Decimal(s[23])).index('.')+3]
            else: var1 = '0'
        if s[26] != 'NA': 
            if s[26] != '0': var2 = str(Decimal(s[26]))[:str(Decimal(s[26])).index('.')+5]
            else: var2 = '0'
        if s[29] != 'NA': 
            if s[29] != '0': var3 = str(Decimal(s[29]))[:str(Decimal(s[29])).index('.')+3]
            else: var3 = '0'
        if s[32] != 'NA': 
            if s[32] != '0': var4 = str(Decimal(s[32]))[:str(Decimal(s[32])).index('.')+5]
            else: var4 = '0'

        so['num_params'] = '5'
        so['Lower limit'] = (var1, var2)
        so['Upper limit'] = (var3, var4)

    elif len(s) == 44:
        if s[22] != 'NA': 
            if s[22] != '0': var1 = str(Decimal(s[22]))[:str(Decimal(s[22])).index('.')+3]
            else: var1 = '0'
        if s[25] != 'NA': 
            if s[25] != '0': var2 = str(Decimal(s[25]))[:str(Decimal(s[25])).index('.')+5]
            else: var2 = '0'
        if s[27] != 'NA': 
            if s[27] != '0': var3 = str(Decimal(s[27]))[:str(Decimal(s[27])).index('.')+3]
            else: var3 = '0'
        if s[30] != 'NA': 
            if s[30] != '0': var4 = str(Decimal(s[30]))[:str(Decimal(s[30])).index('.')+5]
            else: var4 = '0'

        so['num_params'] = '4'
        so['Lower limit'] = (var1, var2)
        so['Upper limit'] = (var3, var4)

    return so


#####################################################################################
# drawDRCCurves                                                                     #
# Plots the DRC model that was calculated as well as some of the summary statistics #
# Note: DRC.plot() params are different than the base R plot() function             #
#####################################################################################
def drawDRCCurves(ag1,ag2,p1,p2,summary1,summary2,ec1,ec2,maxvalue,units,clone,condition):
    R = ro.r
    plot = R.plot
    
    # Set graph margin parameters
    R.par(mar = R.c(5.1, 6, 4.1, 12), oma = R.c(11,0,0,0), xpd = True, bty = 'l')
                        
    # Draw the graphs
    plot(p1, type = 'all', ylim = R.c(0,maxvalue), main = condition+" "+ag2+" vs. "+ag1+"  "+clone, lwd = 2, col = "blue", ylab = "MFI", xlab = "Conc ("+units+")", log = 'x')
    plot(p2, type = 'all', ylim = R.c(0,maxvalue), add = True, lwd = 2, col = "red", pch = 2, log = 'x')
                        
    # Draw the legend
    R.legend("topright", R.c(ag2,ag1), inset = R.c(-0.34,0.15), pch = R.c(1,2), col = R.c('blue','red'))
    R.box("figure", col='black')

    # Add Summary statistic and timestamp
    # Blue Curve
    s = parseSummary(summary1)
    R.mtext(ag2+"  Blue Curve, "+s['num_params']+"-param fit", side=1, line=1, at=0.01, adj=0, outer=True, cex=0.75)
    R.mtext("Estimate", side=1, line=3, at=0.1, adj=0, outer=True, cex=0.75)
    R.mtext("p-value", side=1, line=3, at=0.18, adj=0, outer=True, cex=0.75)
    R.mtext("Lower limit", side=1, line=4, at=0.01, adj=0, outer=True, cex=0.75)
    R.mtext(s['Lower limit'][0], side=1, line=4, at=0.1, adj=0, outer=True, cex=0.75)
    R.mtext(s['Lower limit'][1], side=1, line=4, at=0.18, adj=0, outer=True, cex=0.75)
    R.mtext("Upper limit", side=1, line=5, at=0.01, adj=0, outer=True, cex=0.75)
    R.mtext(s['Upper limit'][0], side=1, line=5, at=0.1, adj=0, outer=True, cex=0.75)
    R.mtext(s['Upper limit'][1], side=1, line=5, at=0.18, adj=0, outer=True, cex=0.75)
    R.mtext("EC10:", side=1, line=6, at=0.01, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(str(ec1[0]), side=1, line=6, at=0.1, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext("EC50:", side=1, line=7, at=0.01, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(str(ec1[1]), side=1, line=7, at=0.1, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext("EC90:", side=1, line=8, at=0.01, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(str(ec1[2]), side=1, line=8, at=0.1, adj=0, outer=True, cex=0.75, col='brown')

    # Red Curve
    s = parseSummary(summary2)
    R.mtext(ag1+"  Red Curve, "+s['num_params']+"-param fit",  side=1, line=1, at=0.5, adj=0, outer=True, cex=0.75)
    R.mtext("Estimate", side=1, line=3, at=0.6, adj=0, outer=True, cex=0.75)
    R.mtext("p-value", side=1, line=3, at=0.68, adj=0, outer=True, cex=0.75)
    R.mtext("Lower limit", side=1, line=4, at=0.5, adj=0, outer=True, cex=0.75)
    R.mtext(s['Lower limit'][0], side=1, line=4, at=0.6, adj=0, outer=True, cex=0.75)
    R.mtext(s['Lower limit'][1], side=1, line=4, at=0.68, adj=0, outer=True, cex=0.75)
    R.mtext("Upper limit", side=1, line=5, at=0.5, adj=0, outer=True, cex=0.75)
    R.mtext(s['Upper limit'][0], side=1, line=5, at=0.6, adj=0, outer=True, cex=0.75)
    R.mtext(s['Upper limit'][1], side=1, line=5, at=0.68, adj=0, outer=True, cex=0.75)
    R.mtext("EC10:", side=1, line=6, at=0.5, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(str(ec2[0]), side=1, line=6, at=0.6, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext("EC50:", side=1, line=7, at=0.5, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(str(ec2[1]), side=1, line=7, at=.6, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext("EC90:", side=1, line=8, at=0.5, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(str(ec2[2]), side=1, line=8, at=0.6, adj=0, outer=True, cex=0.75, col='brown')
    R.mtext(R.c(strftime("%d %b %Y %H:%M",localtime())), cex=0.75, line=10, side=1, adj=1, outer=True)
    return


#####################################################################################
# drawRegCurves                                                                     #
# If the DRC model fails to converge, it drops to this function which simply graphs #
# the points and draws a simple line connecting them, along with an error msg       #
# Note: plot params with base graphics are different than drc.plot()                #
#####################################################################################
def drawRegCurves(condition,ag1,ag2,maxvalue,clone,units,x,y1,y2):
    R = ro.r
    plot = R.plot
    graphics = importr('graphics')
    # Set graph margin parameters
    R.par(mar=R.c(5.1, 6, 4.1, 12), oma=R.c(14,0,0,0), xpd=True, bty='l')

    # Draw graph then overlay second graph
    plot(x,y2, type='o',lwd=2,ylim = R.c(0,maxvalue), col='blue', main = condition+" "+ag2+" vs. "+ag1+"  "+clone, ylab = "MFI", xlab = "Conc ("+units+")", log='x')
    graphics.points(x,y1, type='o',lwd=2,ylim = R.c(0,maxvalue), col='red', log='x')

    # Draw the legend
    R.legend("topright", R.c(ag2,ag1), inset = R.c(-0.34,0.15), pch = R.c(1,2), col = R.c('blue','red'))
    R.box("figure", col = 'black')

    # Print Error msg
    R.mtext("Error: Model failed to converge, so just plotting points", side = 1, line = 1, at = 0.01, adj = 0, outer = True, cex = 0.75, col="darkred")
    return


##########################################################################
# doDRC                                                                  #
# Switch to R to calculate dose response curves. This method coordinates #
# running the model and ploting the graphs                               #
##########################################################################
def doDRC(input_file,conditions,clonelist,comparisons,units):
    # R functions
    ############################################
    grdevices = importr('grDevices')
    R = ro.r
    drc = importr('drc')

    # Import data from tmp file ito R dataframe object
    clone_data = R['read.csv'](input_file)
    
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
    
    # Perform curve-fit analysis and generate graphs for each parental/transfected antigen pair by clone and condition
    #path = MEDIA_ROOT+"/facstool/"
    path = os.getcwd()+"/"
    graphfile = path+"Rplots.pdf"
    jpg_names =[]
    for ftype in ('pdf','jpg'):
        if ftype =='pdf': 
            grdevices.pdf(onefile=True,file=graphfile,width=10,height=10)  # writes multiple graphs per single pdf file
            subprocess.call(['chmod','777',graphfile])
        for condition in subsets['Conditions']:
            for clone in subsets['Conditions'][condition]['clones']:
                for ag1,ag2 in comparisons:
                    
                    if ftype == 'jpg': # Construct filename, create file and change permissions to allow R to write to it
                        jname = path+condition+"_"+clone+'_'+ag1+'_'+ag2+'_graph.jpg'
                        grdevices.jpeg(file=jname,width=10,height=10,units='in',res=300)
                        jpg_names.append(jname)
                        subprocess.call(['chmod','777',jname])
                    
                    # set max value of yaxis from the data, add 500 so it graphs nicely
                    maxvalue = 0  
                    if max(subsets['Conditions'][condition]['clones'][clone].rx2(ag1)) >= max(subsets['Conditions'][condition]['clones'][clone].rx2(ag2)): 
                        maxvalue = max(subsets['Conditions'][condition]['clones'][clone].rx2(ag1))+100
                    else: maxvalue = max(subsets['Conditions'][condition]['clones'][clone].rx2(ag2))+100
                    
                    # Define Y~X, which columns to graph for X and Y by name, since R uses column names
                    formula1 = Formula(str(ag1)+' ~ Conc')
                    formula2 = Formula(str(ag2)+' ~ Conc')
                    
                    # DRC will take a dataframe, so don't need to do this for DRC, but the std plot() is easier to use this way with RPy2
                    # Alternatively: p, summary, ec = calculateModelFit(formula2,subsets['Conditions'][condition]['clones'][clone], L5.L5eval)
                    x = subsets['Conditions'][condition]['clones'][clone].rx2('Conc')
                    y1 = subsets['Conditions'][condition]['clones'][clone].rx2(ag1)
                    y2 = subsets['Conditions'][condition]['clones'][clone].rx2(ag2)
                   
                    # Calculate dose response model
                    p1, summary1, ec1 = calculateModelFit(formula2, ro.DataFrame({'Conc':x, ag2:y2}), L5.L5eval)
                    p2, summary2, ec2 = calculateModelFit(formula1, ro.DataFrame({'Conc':x, ag1:y1}), L5.L5eval)

                    # Draw graphs
                    cond_name = condition  # Must have a condition for parsing data, but not everything is KMD, so if only 1 cond, hide it
                    if len(subsets['Conditions']) == 1: cond_name=""
                    if p1 != 'Error' and p2 != 'Error': drawDRCCurves(ag1, ag2, p1, p2, summary1, summary2, ec1, ec2, maxvalue, units, clone, cond_name)
                    else: drawRegCurves(cond_name, ag1, ag2, maxvalue, clone, units, x, y1, y2)

                    # Close individual graph file
                    if ftype == 'jpg': grdevices.dev_off()
                    
        if ftype == 'pdf': grdevices.dev_off()    # close multi-graph file
        
    print "\n\n>>> Analysis Complete\n\n"
    return #graphfile, jpg_names


####################################
#######         Main       #########
####################################
def main():
    
    # Set up constants for this script. These variables will eventually be derived from a user input form
    # but for now are fixed to help development
    #################################################################################################################################################
    layout = '8'                                                                                                                                      #
    multiple = 0                                                                                                                                    #
    maxvalue = 0                                                                                                                                    #
    datasets = {}                                                                                                                                   #
    conc2 =[]                                                                                                                                       #
    units = 'ug'                                                                                                                                    #
    #========================== 96 well test ============                                                                                           #
    #plate_type = 96                                                                                                                                #
    #plate_layout = 'quad'                                                                                                                          #
    #multiple = 12  # if layout = 8, multiple = 12 and if layout = 12, multiple=8                                                                   #
    #conc = '25,8.33,2.77,0.92,0.30,0.10,0.03,0.01'  # Values are in ug/ml                                                                          #
    #conc = '100,25,6.25,1.56,0.39,0.10,0.02,0.01'  # Values are in nM                                                                              #
    #conditions = {'tranc_130529_40IgGs_Lep_P1.csv':"Lep+",}                                                                                        #
    #clonelist = ['Clone 1','Clone 2','Clone 3','Clone 4','Clone 5','Clone 6','Clone 7','Clone 8','Clone 9','Clone 10','Clone 11','Clone 12']       #
    #comparisons = [('Median__ss_293F','Median__ss_293F_H8'),('Median__ss_CHOK1','Median__ss_CHO_M8')]                                              #
    #========================== 384 well test ===========                                                                                           #
    plate_type = 384                                                                                                                                #
    #plate_layout = 'quad'                                                                                                                           #
    plate_layout = 'linear'                                                                                                                           #
    #multiple = 12  # if layout = 8, multiple = 12 and if layout = 12, multiple=8                                                                   #
    #multiple = 12 # for 8-pt 384 well                                                                                                               #
    #conc = '25,8.33,2.77,0.92,0.30,0.10,0.03,0.01'  # Values are in ug/ml                                                                          #
    conc = '100,25,6.25,1.56,0.39,0.10,0.02,0.01'  # Values are in nM                                                                               #
    conditions = {'384_real_test_data.csv':"Cond1",}                                                                                                #
    new_plates = {'384_real_test_data.csv':'',}                                                                                                     #
    clonelist = ['Clone 1','Clone 2','Clone 3','Clone 4','Clone 5','Clone 6','Clone 7','Clone 8','Clone 9','Clone 10','Clone 11',                   #
                 'Clone 12','Clone 13','Clone 14','Clone 15','Clone 16','Clone 17','Clone 18','Clone 19','Clone 20','Clone 21',                     #
                 'Clone 22','Clone 23','Clone 24','Clone 25','Clone 26','Clone 27','Clone 28','Clone 29','Clone 30','Clone 31',                     #
                 'Clone 32','Clone 33','Clone 34','Clone 35','Clone 36','Clone 37','Clone 38','Clone 39','Clone 40','Clone 41',                     #
                 'Clone 42','Clone 43','Clone 44','Clone 45','Clone 46','Clone 47','Clone 48']                                                      #
#    clonelist = ['Clone 1','Clone 1','Clone 2','Clone 2','Clone 3','Clone 3','Clone 4','Clone 4','Clone 5','Clone 5','Clone 6',                     #
#                 'Clone 6','Clone 7','Clone 7','Clone 8','Clone 8','Clone 9','Clone 9','Clone 10','Clone 10','Clone 11','Clone 11',                 #
#                 'Clone 12','Clone 12','Clone 13','Clone 13','Clone 14','Clone 14','Clone 15','Clone 15','Clone 16','Clone 16',                     #
#                 'Clone 17','Clone 17','Clone 18','Clone 18','Clone 19','Clone 19','Clone 20','Clone 20','Clone 21','Clone 21',                     #
#                 'Clone 22','Clone 22','Clone 23','Clone 23','Clone 24','Clone 24']                                                                 #
    comparisons = [('HEK_293F','HEK_293F_H8'),('CHO_H2','CHO_M8')]                                                                                  #
    ################################################################################################################################################
    
    files = conditions.keys()
    
    if plate_type == 96:
        multiple = getPlatePattern(layout)
        conc_range = createConcRange(conc,multiple)
        antigens = getAgList(files,conditions)
        datasets, input_file = createDataset(files,conditions,antigens,clonelist,multiple,conc_range)
        doDRC(input_file,conditions,clonelist,comparisons,units)
    
    elif plate_type == 384:
        new_conditions = {}
        new_plates = {}
        
        for fn in files:
            new_plates[fn] = merge2List(fn,split384To96(fn,plate_layout))

        for fn in files:
            new_conditions[new_plates[fn]] = conditions[fn]

        multiple = getPlatePattern(layout)
        conc_range = createConcRange(conc,multiple)
        antigens = getAgList(new_conditions.keys(),new_conditions)
        datasets, input_file = createDataset(new_conditions.keys(),new_conditions,antigens,clonelist,multiple,conc_range)
        doDRC(input_file,new_conditions,clonelist,comparisons,units)
            


if __name__ == '__main__':
    main()