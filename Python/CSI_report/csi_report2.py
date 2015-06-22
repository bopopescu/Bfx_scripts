from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
import sys,os, glob, subprocess
from time import sleep

#from reportlab.lib.units import pica

# Imagemagik cmd
# convert -crop 764x664 53-C10_close_up_strip.jpg tiles_%d.jpg

idmap = {'53-C10':'XMetD247-1.19',
         '55-F8':'XMetD247-1.20',
         '56-B11':'XMetD247-1.21',
         '59-B5':'XMetD247-1.22',
         '60-F8':'XMetD247-1.23',
         '61-A9':'XMetD247-1.24',
         '61-G7':'XMetD247-1.25',
         '61-H5':'XMetD247-1.26',
         '62-G11':'XMetD247-1.27',
         '66-E11':'XMetD247-1.28',
         '72-H4':'XMetD247-1.29',
         '75-H9':'XMetD247-1.30',
         '77-F3':'XMetD247-1.31',
         '78-G7':'XMetD247-1.32',
         '79-F3':'XMetD247-1.33',
         '80-C2':'XMetD247-1.34',
         '82-C7':'XMetD247-1.35',
         '92-H4':'XMetD247-1.36',
         '95-C4':'XMetD247-1.37',
         '98-D5':'XMetD247-1.38'}
filelist = {}
filetype = {'filmstrip':'_close_up_strip.jpg',
            'green':'_green_strip.jpg',
            'plate':'_thumbnail.tiff',
            'graph':'_growth.tiff'}
assay_date = "02-03-2015"
barcode = "XOMA358"


def drawPDF(working_canvas,cid):
    stub = filelist[cid][0]
    vals = cid.split('-')
    working_canvas.drawImage(cid+"_tiles_0.jpg",5,185,width=200 , mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(cid+"_tiles_1.jpg",5,15,width=200 , mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(cid+"_tiles_2.jpg",5,-155,width=200 , mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(cid+"_tiles_3.jpg",205,185,width=200 , mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(cid+"_tiles_4.jpg",205,15,width=200 , mask=None, preserveAspectRatio=True)
    if os.path.exists(cid+"_tiles_5.jpg"):
        working_canvas.drawImage(cid+"_tiles_5.jpg",205,-155,width=200, mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(stub+filetype['plate'],437,20,width=350 , mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(stub+filetype['graph'],415,-225,width=370 , mask=None, preserveAspectRatio=True)
    working_canvas.drawImage(stub+filetype['green'],415,85, width=370, mask=None, preserveAspectRatio=True)
    working_canvas.setFont("Helvetica", 18)
    working_canvas.drawString(437,580, "Day 13 plate view")
    working_canvas.drawString(415,285,"Day 01 - Day 13 Colony growth")
    working_canvas.drawString(5,70,"Final clone ID: "+idmap[cid])
    working_canvas.drawString(5,53,"96 well plate ID: "+barcode+'-'+vals[0])
    working_canvas.drawString(5,37,"Well ID: "+vals[1])
    working_canvas.drawString(5,22,"Day 13: "+assay_date)
    working_canvas.showPage()
    return
    




def main():
    
    for infile in glob.glob (os.path.join("","*")):
        if infile.find(".py") == -1 and infile.find(".txt") == -1 and infile.find(".pdf") == -1 and os.path.isdir(infile)==False:
            fn = infile.replace(" ","_")
            os.rename(infile,fn)
            tag = fn.split('_')[0].split(barcode+'-')[1]
            if filelist.has_key(tag): filelist[tag].append(fn)
            else: filelist[tag] = [barcode+'-'+tag]

    clones = filelist.keys()
    clones.sort()
    
    c = canvas.Canvas(barcode+".pdf")
    c.setPageSize(landscape(letter))
    c.setFont("Helvetica", 18)

    for cid in clones:
        fn = filelist[cid][0]+filetype['filmstrip']
        cmd = "convert -crop 764x664 "+ fn+" tiles_%d.jpg"
        print fn
        subprocess.call(["convert","-crop","764x644",fn,cid+"_tiles_%d.jpg"])
        #sleep(3)
        drawPDF(c,cid)
        for jpg in glob.glob (os.path.join("","*")):
            if jpg.find(".py") == -1 and jpg.find(".txt") == -1 and jpg.find(".pdf") == -1 and os.path.isdir(jpg)==False:
                if jpg.find("tiles") != -1:
                    os.remove(jpg)

    c.save()



if __name__ == '__main__':
    main()