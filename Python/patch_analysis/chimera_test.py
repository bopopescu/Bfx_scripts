import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

# change to folder with data files
#os.chdir("test")

# gather the names of .pdb files in the folder
file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]

# loop through the files, opening, processing, and closing each in turn
for fn in file_names:
	replyobj.status("Processing " + fn) # show what file we're working on
	rc("open " + fn)
#	rc("align ligand ~ligand") # put ligand in front of remainder of molecule
	rc("turn 0,0,1 15")
	rc("focus") # center/zoom ligand
	rc("split #0")
	rc("surf") # surface receptor
	rc("color light gray")
#	rc("preset apply publication 3") # make everything look nice
	rc("surftransp 70") # make the surface a little bit see-through
	rc("focus")
	
	rc("color blue #0.1 :155,160,165,166,171,172,173,174,177,179,180,185,200,201,319,412")
	rc("color blue #0.2 :155,160,165,166,171,172,173,174,177,179,180,185,200,201,319,412")
	rc("color blue #0.3 :155,160,165,166,171,172,173,174,177,179,180,185,200,201,319,412")
	rc("color yellow #0.1 :157,158,162,167,168,169,170,175,176,178")
	rc("color yellow #0.2 :157,158,162,167,168,169,170,175,176,178")
	rc("color yellow #0.3 :157,158,162,167,168,169,170,175,176,178")
	png_name = fn[:-3] + ".trimer-top.png"
	rc("copy file " + png_name + " supersample 3")
	rc("turn 0,1,0 90")
	rc("turn 0,0,1 90")
	rc("turn 0,1,0 30")
	rc("~modeldisp #0.3")
	rc("focus")
	rc("color light blue #0.1")
	rc("color blue #0.1 :155,160,165,166,171,172,173,174,177,179,180,185,200,201,319,412")
	rc("color yellow #0.1 :157,158,162,167,168,169,170,175,176,178")
	# save image to a file that ends in .png rather than .pdb
	png_name = fn[:-3] + ".foo.png"
	rc("windowsize 1050 700") # 2100 1400
	rc("copy file " + png_name + " supersample 3")
	rc("close all")
# uncommenting the line below will cause Chimera to exit when the script is done
#rc("stop now")
# note that indentation is significant in Python; the fact that
# the above command is exdented means that it is executed after
# the loop completes, whereas the indented commands that 
# preceded it are executed as part of the loop.