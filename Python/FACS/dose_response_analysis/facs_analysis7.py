# facs_analysis5.py
# By Mark Evans
# Created 11.01.2013
# Revised 11.21.2013



import os, os.path, sys, subprocess
import math, datetime
from time import localtime, strftime
from decimal import *
import pdb

import rpy2
from rpy2.robjects.packages import importr
from rpy2 import robjects as ro
from rpy2.robjects import Formula

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# Example of log concentration values used for their dose-response curves.
# Conc = 25,8.33,2.77,0.92,0.30,0.10,0.03,0.01

CLONE_LAYOUT_384 = {'linear_layout':{'8-point':{'C1':[1, 25, 49, 73, 97, 121, 145, 169],
                                                'C2':[2, 26, 50, 74, 98, 122, 146, 170],
                                                'C3':[3, 27, 51, 75, 99, 123, 147, 171],
                                                'C4':[4, 28, 52, 76, 100, 124, 148, 172],
                                                'C5':[5, 29, 53, 77, 101, 125, 149, 173],
                                                'C6':[6, 30, 54, 78, 102, 126, 150, 174],
                                                'C7':[7, 31, 55, 79, 103, 127, 151, 175],
                                                'C8':[8, 32, 56, 80, 104, 128, 152, 176],
                                                'C9':[9, 33, 57, 81, 105, 129, 153, 177],
                                                'C10':[10, 34, 58, 82, 106, 130, 154, 178],
                                                'C11':[11, 35, 59, 83, 107, 131, 155, 179],
                                                'C12':[12, 36, 60, 84, 108, 132, 156, 180],
                                                'C13':[13, 37, 61, 85, 109, 133, 157, 181],
                                                'C14':[14, 38, 62, 86, 110, 134, 158, 182],
                                                'C15':[15, 39, 63, 87, 111, 135, 159, 183],
                                                'C16':[16, 40, 64, 88, 112, 136, 160, 184],
                                                'C17':[17, 41, 65, 89, 113, 137, 161, 185],
                                                'C18':[18, 42, 66, 90, 114, 138, 162, 186],
                                                'C19':[19, 43, 67, 91, 115, 139, 163, 187],
                                                'C20':[20, 44, 68, 92, 116, 140, 164, 188],
                                                'C21':[21, 45, 69, 93, 117, 141, 165, 189],
                                                'C22':[22, 46, 70, 94, 118, 142, 166, 190],
                                                'C23':[23, 47, 71, 95, 119, 143, 167, 191],
                                                'C24':[24, 48, 72, 96, 120, 144, 168, 192],
                                                'C25':[193, 217, 241, 265, 289, 313, 337, 361],
                                                'C26':[194, 218, 242, 266, 290, 314, 338, 362],
                                                'C27':[195, 219, 243, 267, 291, 315, 339, 363],
                                                'C28':[196, 220, 244, 268, 292, 316, 340, 364],
                                                'C29':[197, 221, 245, 269, 293, 317, 341, 365],
                                                'C30':[198, 222, 246, 270, 294, 318, 342, 366],
                                                'C31':[199, 223, 247, 271, 295, 319, 343, 367],
                                                'C32':[200, 224, 248, 272, 296, 320, 344, 368],
                                                'C33':[201, 225, 249, 273, 297, 321, 345, 369],
                                                'C34':[202, 226, 250, 274, 298, 322, 346, 370],
                                                'C35':[203, 227, 251, 275, 299, 323, 347, 371],
                                                'C36':[204, 228, 252, 276, 300, 324, 348, 372],
                                                'C37':[205, 229, 253, 277, 301, 325, 349, 373],
                                                'C38':[206, 230, 254, 278, 302, 326, 350, 374],
                                                'C39':[207, 231, 255, 279, 303, 327, 351, 375],
                                                'C40':[208, 232, 256, 280, 304, 328, 352, 376],
                                                'C41':[209, 233, 257, 281, 305, 329, 353, 377],
                                                'C42':[210, 234, 258, 282, 306, 330, 354, 378],
                                                'C43':[211, 235, 259, 283, 307, 331, 355, 379],
                                                'C44':[212, 236, 260, 284, 308, 332, 356, 380],
                                                'C45':[213, 237, 261, 285, 309, 333, 357, 381],
                                                'C46':[214, 238, 262, 286, 310, 334, 358, 382],
                                                'C47':[215, 239, 263, 287, 311, 335, 359, 383],
                                                'C48':[216, 240, 264, 288, 312, 336, 360, 384]},
                                     '12-point':{'C1':[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                                 'C2':[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
                                                 'C3':[49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60],
                                                 'C4':[73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84],
                                                 'C5':[97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108],
                                                 'C6':[121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132],
                                                 'C7':[145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156],
                                                 'C8':[169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180],
                                                 'C9':[13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
                                                 'C10':[37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48],
                                                 'C11':[61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72],
                                                 'C12':[85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96],
                                                 'C13':[109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120],
                                                 'C14':[133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,],
                                                 'C15':[157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168],
                                                 'C16':[181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192],
                                                 'C17':[193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204],
                                                 'C18':[217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228],
                                                 'C19':[241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252],
                                                 'C20':[265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276],
                                                 'C21':[289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300],
                                                 'C22':[313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324],
                                                 'C23':[337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348],
                                                 'C24':[361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372],
                                                 'C25':[205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216],
                                                 'C26':[229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240],
                                                 'C27':[253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264],
                                                 'C28':[277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288],
                                                 'C29':[301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312],
                                                 'C30':[325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336],
                                                 'C31':[349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360],
                                                 'C32':[373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384]},
                                     '12-point-vert':{'C1':[1, 25, 49, 73, 97, 121, 145, 169, 193, 217, 241, 265],
                                                'C2':[2, 26, 50, 74, 98, 122, 146, 170, 194, 218, 242, 266],
                                                'C3':[3, 27, 51, 75, 99, 123, 147, 171, 195, 219, 243, 267],
                                                'C4':[4, 28, 52, 76, 100, 124, 148, 172, 196, 220, 244, 268],
                                                'C5':[5, 29, 53, 77, 101, 125, 149, 173, 197, 221, 245, 269],
                                                'C6':[6, 30, 54, 78, 102, 126, 150, 174, 198, 222, 246, 270],
                                                'C7':[7, 31, 55, 79, 103, 127, 151, 175, 199, 223, 247, 271],
                                                'C8':[8, 32, 56, 80, 104, 128, 152, 176, 200, 224, 248, 272],
                                                'C9':[9, 33, 57, 81, 105, 129, 153, 177, 201, 225, 249, 273],
                                                'C10':[10, 34, 58, 82, 106, 130, 154, 178, 202, 226, 250, 274],
                                                'C11':[11, 35, 59, 83, 107, 131, 155, 179, 203, 227, 251, 275],
                                                'C12':[12, 36, 60, 84, 108, 132, 156, 180, 204, 228, 252, 276],
                                                'C13':[13, 37, 61, 85, 109, 133, 157, 181, 205, 229, 253, 277],
                                                'C14':[14, 38, 62, 86, 110, 134, 158, 182, 206, 230, 254, 278],
                                                'C15':[15, 39, 63, 87, 111, 135, 159, 183, 207, 231, 255, 279],
                                                'C16':[16, 40, 64, 88, 112, 136, 160, 184, 208, 232, 256, 280],
                                                'C17':[17, 41, 65, 89, 113, 137, 161, 185, 209, 233, 257, 281],
                                                'C18':[18, 42, 66, 90, 114, 138, 162, 186, 210, 234, 258, 282],
                                                'C19':[19, 43, 67, 91, 115, 139, 163, 187, 211, 235, 259, 283],
                                                'C20':[20, 44, 68, 92, 116, 140, 164, 188, 212, 236, 260, 284],
                                                'C21':[21, 45, 69, 93, 117, 141, 165, 189, 213, 237, 261, 285],
                                                'C22':[22, 46, 70, 94, 118, 142, 166, 190, 214, 238, 262, 286],
                                                'C23':[23, 47, 71, 95, 119, 143, 167, 191, 215, 239, 263, 287],
                                                'C24':[24, 48, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288]}},
                    'well_quads':{'8-point':{'C1':[1, 49, 97, 145, 193, 241, 289, 337],
                                                'C2':[2, 50, 98, 146, 194, 242, 290, 338],
                                                'C3':[3, 51, 99, 147, 195, 243, 291, 339],
                                                'C4':[4, 52, 100, 148, 196, 244, 292, 340],
                                                'C5':[5, 53, 101, 149, 197, 245, 293, 341],
                                                'C6':[6, 54, 102, 150, 198, 246, 294, 342],
                                                'C7':[7, 55, 103, 151, 199, 247, 295, 343],
                                                'C8':[8, 56, 104, 152, 200, 248, 296, 344],
                                                'C9':[9, 57, 105, 153, 201, 249, 297, 345],
                                                'C10':[10, 58, 106, 154, 202, 250, 298, 346],
                                                'C11':[11, 59, 107, 155, 203, 251, 299, 347],
                                                'C12':[12, 60, 108, 156, 204, 252, 300, 348],
                                                'C13':[13, 61, 109, 157, 205, 253, 301, 349],
                                                'C14':[14, 62, 110, 158, 206, 254, 302, 350],
                                                'C15':[15, 63, 111, 159, 207, 255, 303, 351],
                                                'C16':[16, 64, 112, 160, 208, 256, 304, 352],
                                                'C17':[17, 65, 113, 161, 209, 257, 305, 353],
                                                'C18':[18, 66, 114, 162, 210, 258, 306, 354],
                                                'C19':[19, 67, 115, 163, 211, 259, 307, 355],
                                                'C20':[20, 68, 116, 164, 212, 260, 308, 356],
                                                'C21':[21, 69, 117, 165, 213, 261, 309, 357],
                                                'C22':[22, 70, 118, 166, 214, 262, 310, 358],
                                                'C23':[23, 71, 119, 167, 215, 263, 311, 359],
                                                'C24':[24, 72, 120, 168, 216, 264, 312, 360],
                                                'C25':[25, 73, 121, 169, 217, 265, 313, 361],
                                                'C26':[26, 74, 122, 170, 218, 266, 314, 362],
                                                'C27':[27, 75, 123, 171, 219, 267, 315, 363],
                                                'C28':[28, 76, 124, 172, 220, 268, 316, 364],
                                                'C29':[29, 77, 125, 173, 221, 269, 317, 365],
                                                'C30':[30, 78, 126, 174, 222, 270, 318, 366],
                                                'C31':[31, 79, 127, 175, 223, 271, 319, 367],
                                                'C32':[32, 80, 128, 176, 224, 272, 320, 368],
                                                'C33':[33, 81, 129, 177, 225, 273, 321, 369],
                                                'C34':[34, 82, 130, 178, 226, 274, 322, 370],
                                                'C35':[35, 83, 131, 179, 227, 275, 323, 371],
                                                'C36':[36, 84, 132, 180, 228, 276, 324, 372],
                                                'C37':[37, 85, 133, 181, 229, 277, 325, 373],
                                                'C38':[38, 86, 134, 182, 230, 278, 326, 374],
                                                'C39':[39, 87, 135, 183, 231, 279, 327, 375],
                                                'C40':[40, 88, 136, 184, 232, 280, 328, 376],
                                                'C41':[41, 89, 137, 185, 233, 281, 329, 377],
                                                'C42':[42, 90, 138, 186, 234, 282, 330, 378],
                                                'C43':[43, 91, 139, 187, 235, 283, 331, 379],
                                                'C44':[44, 92, 140, 188, 236, 284, 332, 380],
                                                'C45':[45, 93, 141, 189, 237, 285, 333, 381],
                                                'C46':[46, 94, 142, 190, 238, 286, 334, 382],
                                                'C47':[47, 95, 143, 191, 239, 287, 335, 383],
                                                'C48':[48, 96, 144, 192, 240, 288, 336, 384]},
                                  '12-point':{'C1':[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23],
                                                'C2':[2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24],
                                                'C3':[25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47],
                                                'C4':[26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48],
                                                'C5':[49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71],
                                                'C6':[50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72],
                                                'C7':[73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95],
                                                'C8':[74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96],
                                                'C9':[97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119],
                                                'C10':[98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120],
                                                'C11':[121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143],
                                                'C12':[122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144],
                                                'C13':[145, 147, 149, 151, 153, 155, 157, 159, 161, 163, 165, 167],
                                                'C14':[146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168],
                                                'C15':[169, 171, 173, 175, 177, 179, 181, 183, 185, 187, 189, 191],
                                                'C16':[170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192],
                                                'C17':[193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215],
                                                'C18':[194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216],
                                                'C19':[217, 219, 221, 223, 225, 227, 229, 231, 233, 235, 237, 239],
                                                'C20':[218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240],
                                                'C21':[241, 243, 245, 247, 249, 251, 253, 255, 257, 259, 261, 263],
                                                'C22':[242, 244, 246, 248, 250, 252, 254, 256, 258, 260, 262, 264],
                                                'C23':[265, 267, 269, 271, 273, 275, 277, 279, 281, 283, 285, 287],
                                                'C24':[266, 268, 270, 272, 274, 276, 278, 280, 282, 284, 286, 288],
                                                'C25':[289, 291, 293, 295, 297, 299, 301, 303, 305, 307, 309, 311],
                                                'C26':[290, 292, 294, 296, 298, 300, 302, 304, 306, 308, 310, 312],
                                                'C27':[313, 315, 317, 319, 321, 323, 325, 327, 329, 331, 333, 335],
                                                'C28':[314, 316, 318, 320, 322, 324, 326, 328, 330, 332, 334, 336],
                                                'C29':[337, 339, 341, 343, 345, 347, 349, 351, 353, 355, 357, 359],
                                                'C30':[338, 340, 342, 344, 346, 348, 350, 352, 354, 356, 358, 360],
                                                'C31':[361, 363, 365, 367, 369, 371, 373, 375, 377, 379, 381, 383],
                                                'C32':[362, 364, 366, 368, 370, 372, 374, 376, 378, 380, 382, 384]}}}

CLONE_LAYOUT_96 = {'8-point':{'C1':[1,13,25,37,49,61,73,85],
                              'C2':[2,14,26,38,50,62,74,86],
                              'C3':[3,15,27,39,51,63,75,87],
                              'C4':[4,16,28,40,52,64,76,88],
                              'C5':[5,17,29,41,53,65,77,89],
                              'C6':[6,18,30,42,54,66,78,90],
                              'C7':[7,19,31,43,55,67,79,91],
                              'C8':[8,20,32,44,56,68,80,92],
                              'C9':[9,21,33,45,57,69,81,93],
                              'C10':[10,22,34,46,58,70,82,94],
                              'C11':[11,23,35,47,59,71,83,95],
                              'C12':[12,24,36,48,60,72,84,96] },
                   '12-point':{'C1':[1,2,3,4,5,6,7,8,9,10,11,12],
                               'C2':[13,14,15,16,17,18,19,20,21,22,23,24],
                               'C3':[25,26,27,28,29,30,31,32,33,34,35,36],
                               'C4':[37,38,39,40,41,42,43,44,45,46,47,48],
                               'C5':[49,50,51,52,53,54,55,56,57,58,59,60],
                               'C6':[61,62,63,64,65,66,67,68,69,70,71,72],
                               'C7':[73,74,75,76,77,78,79,80,81,82,83,84],
                               'C8':[85,86,87,88,89,90,91,92,93,94,95,96]} }



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
def createDataset2(files,conditions,antigens,clonelist,conc,plate_layout,plate_size,titration):
    # Parse data and write to a tmp file in modified format
    tmp = open('tmp.txt','w')
    rawdata = {}
    conds = []

    # Read raw data from input files into dictionary
    for fn in files:
        f = open(fn,'r')
        c = 1
        conds.append(conditions[fn])
        for line in f:
            vals = line.rstrip().split(',')
            
            if vals[0]=='Sample' and vals[-1:len(vals)][0].find('Ratio') != 0: c ='fulllength'  # No ratio column, data cols only
            if vals[0]!='Sample' and vals[0] !='Mean' and vals[0] != 'StdDev':
                idx = int(vals[0].split(': ')[0])   # abs well pos
                
                if rawdata.has_key(conditions[fn]):
                    if c != 'fulllength': rawdata[conditions[fn]][idx] = vals[1:-2]
                    else: rawdata[conditions[fn]][idx] = vals[1:]
                else:
                    if c != 'fulllength': rawdata[conditions[fn]] = {idx:vals[1:-2]}
                    else: rawdata[conditions[fn]] = {idx:vals[1:]}

    
    
    # Create new dataset 
    tmp.write('Conditions,Clones,Conc,'+','.join(antigens[antigens.keys()[0]])+'\n')
    for condition in conds:
        if plate_size == 384:
            if plate_layout == 'linear_layout':
                if titration == '12-point':
                    for clone_idx in range(1,33):
                        for C, pos in zip(conc,CLONE_LAYOUT_384[plate_layout][titration]['C'+str(clone_idx)]):
                            tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
                elif titration == '8-point':
                    for clone_idx in range(1,49):
                        for C, pos in zip(conc, CLONE_LAYOUT_384[plate_layout][titration]['C'+str(clone_idx)]):
                            tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
                elif titration == '12-point-vert':
                    for clone_idx in range(1,25):
                        for C, pos in zip(conc,CLONE_LAYOUT_384[plate_layout][titration]['C'+str(clone_idx)]):
                            tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
            elif plate_layout == 'well_quads': 
                if titration == '12-point':
                    for clone_idx in range(1,33):
                        for C, pos in zip(conc, CLONE_LAYOUT_384[plate_layout][titration]['C'+str(clone_idx)]):
                            tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
                elif titration == '8-point':
                    for clone_idx in range(1,49):
                        for C, pos in zip(conc, CLONE_LAYOUT_384[plate_layout][titration]['C'+str(clone_idx)]):
                            tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
        elif plate_size == 96:
            if titration == '12-point':
                for clone_idx in range(1,9):
                    for C, pos in zip(conc, CLONE_LAYOUT_96[titration]['C'+str(clone_idx)]):
                        tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
            elif titration == '8-point':
                for clone_idx in range(1,13):
                    for C, pos in zip(conc, CLONE_LAYOUT_96[titration]['C'+str(clone_idx)]):
                        tmp.write(condition+','+clonelist[clone_idx-1]+','+C+','+','.join(rawdata[condition][pos])+'\n')
    tmp.close()
    return tmp.name



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
            if s[23] != '0': var1 = str(float(s[23]))[:str(float(s[23])).index('.')+3]
            else: var1 = '0'
        if s[26] != 'NA': 
            if s[26] != '0': var2 = str(float(s[26]))[:str(float(s[26])).index('.')+5]
            else: var2 = '0'
        if s[29] != 'NA': 
            print "s[29] = ",s[29]," str Decimal= ",str(float(s[29]))
            if s[29] != '0': var3 = str(float(s[29]))[:str(float(s[29])).index('.')+3]
            else: var3 = '0'
        if s[32] != 'NA': 
            if s[32] != '0': var4 = str(float(s[32]))[:str(float(s[32])).index('.')+5]
            else: var4 = '0'

        so['num_params'] = '5'
        so['Lower limit'] = (var1, var2)
        so['Upper limit'] = (var3, var4)

    elif len(s) == 44:
        if s[22] != 'NA': 
            if s[22] != '0': var1 = str(float(s[22]))[:str(float(s[22])).index('.')+3]
            else: var1 = '0'
        if s[25] != 'NA': 
            if s[25] != '0': var2 = str(float(s[25]))[:str(float(s[25])).index('.')+5]
            else: var2 = '0'
        if s[27] != 'NA': 
            if s[27] != '0': var3 = str(float(s[27]))[:str(float(s[27])).index('.')+3]
            else: var3 = '0'
        if s[30] != 'NA': 
            if s[30] != '0': var4 = str(float(s[30]))[:str(float(s[30])).index('.')+5]
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
def drawDRCCurves(condition, clone, maxvalue, units, ag1, p1, summary1, ec1, ag2='', p2='', summary2='', ec2=''):
    R = ro.r
    plot = R.plot
    
    # Set graph margin parameters
    R.par(mar = R.c(5.1, 6, 4.1, 12), oma = R.c(11,0,0,0), xpd = True, bty = 'l')
    
    # Set param to force numbers to decimal rather than scientific notation
    R.options(scipen='1000', digits='6')

    # Draw the graphs
    if ag2 != '':
        plot(p1, type = 'all', ylim = R.c(0,maxvalue), main = condition+" "+ag2+" vs. "+ag1+"  "+clone, lwd = 2, col = "blue", ylab = "MFI", xlab = "Conc ("+units+")", log = 'x')
        plot(p2, type = 'all', ylim = R.c(0,maxvalue), add = True, lwd = 2, col = "red", pch = 2, log = 'x')
    else:
        plot(p1, type = 'all', ylim = R.c(0,maxvalue), main = condition+" "+ag1+"  "+clone, lwd = 2, col = "blue", ylab = "MFI", xlab = "Conc ("+units+")", log = 'x')
    
    # Reset param for the rest of the graphs
    R.options(scipen='0', digits='6')

    # Draw the legend
    if ag2 != '':  R.legend("topright", R.c(ag2,ag1), inset = R.c(-0.34,0.15), pch = R.c(1,2), col = R.c('blue','red'))
    else:  R.legend("topright", R.c(ag1), inset = R.c(-0.34,0.15), pch = R.c(1), col = R.c('blue'))
    R.box("figure", col='black')

    # Add Summary statistic and timestamp
    # Blue Curve
    s = parseSummary(summary1)
    R.mtext(ag1+"  Blue Curve, "+s['num_params']+"-param fit", side=1, line=1, at=0.01, adj=0, outer=True, cex=0.75)
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
    if ag2 !='':
        s = parseSummary(summary2)
        R.mtext(ag2+"  Red Curve, "+s['num_params']+"-param fit",  side=1, line=1, at=0.5, adj=0, outer=True, cex=0.75)
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
def drawRegCurves(condition,clone,maxvalue,units,ag1,x,y1,ag2='',y2=''):
    R = ro.r
    plot = R.plot
    graphics = importr('graphics')
    
    # Set graph margin parameters
    R.par(mar=R.c(5.1, 6, 4.1, 12), oma=R.c(14,0,0,0), xpd=True, bty='l')
    
    # Set param to force numbers to decimal rather than scientific notation
    R.options(scipen='1000', digits='6')
    
    # Draw graph then overlay second graph
    if ag2 !='':
        plot(x,y2, type='o',lwd=2,ylim = R.c(0,maxvalue), col='blue', main = condition+" "+ag2+" vs. "+ag1+"  "+clone, ylab = "MFI", xlab = "Conc ("+units+")", log='x')
        graphics.points(x,y1, type='o',lwd=2,ylim = R.c(0,maxvalue), col='red', log='x')
    else:
        plot(x,y1, type='o',lwd=2,ylim = R.c(0,maxvalue), col='blue', main = condition+" "+ag1+"  "+clone, ylab = "MFI", xlab = "Conc ("+units+")", log='x')

    # Reset param for the rest of the graphs
    R.options(scipen='0', digits='6')

    # Draw the legend
    if ag2 !='': R.legend("topright", R.c(ag2,ag1), inset = R.c(-0.34,0.15), pch = R.c(1,2), col = R.c('blue','red'))
    else: R.legend("topright", R.c(ag1), inset = R.c(-0.34,0.15), pch = R.c(1), col = R.c('blue'))
    R.box("figure", col = 'black')

    # Print Error msg
    R.mtext("Error: Model failed to converge, so just plotting points", side = 1, line = 1, at = 0.01, adj = 0, outer = True, cex = 0.75, col="darkred")
    return


##########################################################################
# doDRC                                                                  #
# Switch to R to calculate dose response curves. This method coordinates #
# running the model and ploting the graphs                               #
##########################################################################
def doDRC(input_file,conditions,clonelist,comparisons,units,antigens):
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
                # If there are no comparisons (single Ag, etc)
                if len(comparisons) < 2:
                    for ag1 in antigens[condition]:
                        if ftype == 'jpg': # Construct filename, create file and change permissions to allow R to write to it
                            jname = path+condition+"_"+clone+'_'+ag1+'_graph.jpg'
                            grdevices.jpeg(file=jname,width=10,height=10,units='in',res=300)
                            jpg_names.append(jname)
                            subprocess.call(['chmod','777',jname])
                        
                        # set max value of yaxis from the data, add 500 so it graphs nicely
                        maxvalue = max(subsets['Conditions'][condition]['clones'][clone].rx2(ag1)) +100
                    
                        # Define Y~X, which columns to graph for X and Y by name, since R uses column names
                        formula1 = Formula(str(ag1)+' ~ Conc')
                        
                        # DRC will take a dataframe, so don't need to do this for DRC, but the std plot() is easier to use this way with RPy2
                        # Alternatively: p, summary, ec = calculateModelFit(formula2,subsets['Conditions'][condition]['clones'][clone], L5.L5eval)
                        x = subsets['Conditions'][condition]['clones'][clone].rx2('Conc')
                        y1 = subsets['Conditions'][condition]['clones'][clone].rx2(ag1)
                   
                        # Calculate dose response model
                        p1, summary1, ec1 = calculateModelFit(formula1, ro.DataFrame({'Conc':x, ag1:y1}), L5.L5eval)

                        # Draw graphs
                        cond_name = condition  # Must have a condition for parsing data, but not everything is KMD, so if only 1 cond, hide it
                        if len(subsets['Conditions']) == 1: cond_name=""
                        if p1 != 'Error': drawDRCCurves(cond_name, clone, maxvalue, units, ag1, p1, summary1, ec1)
                        else: drawRegCurves(cond_name, clone, maxvalue, units, ag1, x, y1)

                        # Close individual graph file
                        if ftype == 'jpg': grdevices.dev_off()

                # Have at least 2 Ag to compare
                else:
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
                        if p1 != 'Error' and p2 != 'Error': drawDRCCurves(cond_name, clone, maxvalue, units, ag1, p1, summary1, ec1, ag2, p2, summary2, ec2)
                        else: drawRegCurves(cond_name, clone, maxvalue, units, ag1, x, y1, ag2, y2)

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
   # layout = '12'                                                                                                                                      #
   # multiple = 0                                                                                                                                    #
   # maxvalue = 0                                                                                                                                    #
   # datasets = {}                                                                                                                                   #
   # conc2 =[]                                                                                                                                       #
    units = 'ug/ml'                                                                                                                                    #
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
    plate_layout = 'linear_layout'            
    titration = '8-point'                                                                                                               #
    #multiple = 12  # if layout = 8, multiple = 12 and if layout = 12, multiple=8                                                                   #
    #multiple = 12 # for 8-pt 384 well                                                                                                               #
    #conc = '25,8.33,2.77,0.92,0.30,0.10,0.03,0.01'  # Values are in ug/ml                                                                          #
    conc = '100,25,6.25,1.56,0.39,0.10,0.02,0.01'.split(',')  # Values are in nM     
    #conc = '100,25,6.25,1.56,0.39,0.097,0.02,0.006,0.0015,0.00038,0.000095,0.000023'.split(',')                                                                          #
    conditions = {'384_real_test_data.csv':"Cond1",}                                                                                                #
    new_plates = {'384_real_test_data.csv':'',}                                                                                                     #
    clonelist = ['Clone 1','Clone 2','Clone 3','Clone 4','Clone 5','Clone 6','Clone 7','Clone 8','Clone 9','Clone 10','Clone 11',                   #
                 'Clone 12','Clone 13','Clone 14','Clone 15','Clone 16','Clone 17','Clone 18','Clone 19','Clone 20','Clone 21',                     #
                 'Clone 22','Clone 23','Clone 24','Clone 25','Clone 26','Clone 27','Clone 28','Clone 29','Clone 30','Clone 31',                     #
                 'Clone 32','Clone 33','Clone 34','Clone 35','Clone 36','Clone 37','Clone 38','Clone 39','Clone 40','Clone 41',                     #
                 'Clone 42','Clone 43','Clone 44','Clone 45','Clone 46','Clone 47','Clone 48']                                                      #
 #   clonelist = ['Clone 1','Clone 2','Clone 3','Clone 4','Clone 5','Clone 6','Clone 7','Clone 8','Clone 9','Clone 10','Clone 11',                   #
 #                'Clone 12','Clone 13','Clone 14','Clone 15','Clone 16','Clone 17','Clone 18','Clone 19','Clone 20','Clone 21',                     #
 #                'Clone 22','Clone 23','Clone 24','Clone 25','Clone 26','Clone 27','Clone 28','Clone 29','Clone 30','Clone 31',                     #
 #                'Clone 32']                                                                                                                        #
#    clonelist = ['Clone 1','Clone 1','Clone 2','Clone 2','Clone 3','Clone 3','Clone 4','Clone 4','Clone 5','Clone 5','Clone 6',                     #
#                 'Clone 6','Clone 7','Clone 7','Clone 8','Clone 8','Clone 9','Clone 9','Clone 10','Clone 10','Clone 11','Clone 11',                 #
#                 'Clone 12','Clone 12','Clone 13','Clone 13','Clone 14','Clone 14','Clone 15','Clone 15','Clone 16','Clone 16',                     #
#                 'Clone 17','Clone 17','Clone 18','Clone 18','Clone 19','Clone 19','Clone 20','Clone 20','Clone 21','Clone 21',                     #
#                 'Clone 22','Clone 22','Clone 23','Clone 23','Clone 24','Clone 24']                                                                 #
 #   comparisons = [('HEK_293F','HEK_293F_H8'),('CHO_H2','CHO_M8')]                                                                                  #
    comparisons = []                                                                                                                             #
    ################################################################################################################################################
    
    files = conditions.keys()
    antigens = getAgList(files, conditions)
    input_file = createDataset2(files, conditions, antigens, clonelist, conc, plate_layout, plate_type, titration)
    doDRC(input_file, conditions, clonelist, comparisons, units, antigens)
            


if __name__ == '__main__':
    main()