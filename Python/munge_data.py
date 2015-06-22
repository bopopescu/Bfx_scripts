'''
Clincial data file munger for Sandra
by Mark Evans
Created 2013.05.24

This program takes as input the filename or filepath to an index file in .csv format
The input file contents should follow this format :

Subject,TP (Day),Test,Barcode,Sample,Well,IL-1b file,IL-6 File,IL-8 File,TNF-a File,INF-g File
1001,0,A,311868,U001,C01,IL-1b Unknowns 100912.csv,IL-6 Unknowns 100912.csv,IL-8 Unknowns 100912.csv,TNF-a Unknowns 100912.csv,INF-g Unknowns 100912.csv
,,,,U001,C02,IL-1b Unknowns 100912.csv,IL-6 Unknowns 100912.csv,IL-8 Unknowns 100912.csv,TNF-a Unknowns 100912.csv,INF-g Unknowns 100912.csv

The number of molecules tested (e.g. IL-6, etc) can float
The filenames in a column do not all have to be the same file
The index file and the data files DO have to be in the same place as this script

Data files should be raw csv exported from the MSD in the format:

Unknown: IL-6 (Human),,,,,,,,,,,
Sample,Well,Signal,Mean,Std. Deviation,CV,Calc. Concentration,Calc. Conc. Mean,Calc. Conc. Std. Deviation,Calc. Conc. CV,Detection Range,
U001,C01,612,682,98.28784258,14.42228064,97.0833733,107.4022335,14.59307211,13.58730785,In Detection Range,

Final output file will be called clinical_summary.csv and will be written in the same folder as this script.

'''

import os, sys

#####################################################################################
# parseDataFile                                                                     #
# This function will process a raw data file. It skips blank lines and header lines #
# Data is returned in a dictionary of lists with sample|well as the key             #
#####################################################################################
def parseDataFile(datafilename):
	datafile = open(datafilename,'r')
	data = {}

	for line in datafile:
		a = line.rstrip().split(',')
		if a[0] != '' and a[0][0:3] != 'Unk' and a[0][0:3]!='Sam':
			k = a[0]+'|'+a[1]	# create key
			data[k] = a[7:10]

			# The Detection range column is a text value, but in the summary it is 
			# a number that indicates how many samples fell outside the range
			# so here we insert a number instead of text so later we can just add them
			if a[10] == 'In Detection Range': data[k].append(0)
			elif a[10] == 'Below Detection Range': data[k].append(1)
			elif a[10] == 'Below Fit Curve Range': data[k].append(1)
			elif a[10] == 'Above Fit Curve Range': data[k].append(1)
	datafile.close()
	return data


#########################################################################################
# parseIndexFile                                                                        # 
# Processes the sample index file, handles the first 4 blank fields on every other line #
# Returns a dictionary of datafile names with test molecules as indicies                #
# Also returns a dictionary of all the samples with sample|well as the key              #
#########################################################################################
def parseIndexFile(indexfilename):
	f1 = open(indexfilename,"r")
	c = 1    # line counter just used to see the first line
	files = {}
	samples = {}
	file_keys = []
	prev_key =''

	for line in f1:
		vals = line.rstrip().split(',')
		if c == 1:
			# iterate through all test molecules and assigned them as keys. 
			# number of molec can be variable
			for j in range(6,len(vals)):
				files[vals[j].split()[0]]=[]
				file_keys.append(vals[j].split()[0])
		else:
			k = '|'.join(vals[4:6])  # create key

			# if line is not missing first four values
			if vals[0] != '':
				samples[k] = vals[0:4]

			# if it is missing first four values
			else:
				data = samples[prev_key]
				samples[k] = data
			# Loop through the test molec filenames and add them	
			i = 6
			for x in file_keys:
				if vals[i] not in files[x]:
					files[x].append(vals[i])
				i += 1 
			prev_key = k
		c += 1
	f1.close()

	return (files,samples)	


######################
# Main               #
######################
def main():
	print "Please make sure that all of the raw data files "
	print "and the index file are in the same location.\n"
	idx_filename = raw_input("Enter the name or path to the index file: ")

	r1 = open('clinical_summary.csv',"w")
	dataset = {}

	# Parse the index file
	files,samples = parseIndexFile(idx_filename)

	# Loop through all raw data files and parse them
	for test_molec in files.keys():
		if len(files[test_molec]) == 1:
			dataset[test_molec] = parseDataFile(files[test_molec][0])
		else:
			for fn in files[test_molec]:
				dataset[test_molec] = dict(dataset[test_molec].items() + parseDataFile(fn).items())

	# Begin writing summary file
	sum_header = "Barcode,Subject,TP (Day),Test,Unknown,"
	header_generic = "%s Mean,%s Std. Dev.,%s CV,%s # Below,"
	test_molecs = files.keys()
	test_molecs.sort()

	# Loop through all the test molecules and add them to the header line
	for tm in test_molecs:
		sum_header += header_generic % (tm,tm,tm,tm)

	# write header line
	r1.write(sum_header.rstrip(',')+'\n')

	# Begin writing summary data
	subject_keys = samples.keys()
	subject_keys.sort()
	dataline = ""
	prev_sample = ""
	prev_well = ""

	# Loop through all of the sample keys from the index file
	for sk in subject_keys:
		sample,well = sk.split('|')
		subject,tp,test,barcode = samples[sk][0:4]
		if sample == prev_sample:
			dataline = ','.join([barcode,subject,tp,test,sample])+','
			for tm in test_molecs:
				k = sample+'|'+well
				pk = prev_sample+'|'+prev_well
				dr = str(dataset[tm][k][3]+dataset[tm][pk][3])   # Add the values for detection range field
				dataline = dataline + ','.join([dataset[tm][k][0],dataset[tm][k][1],dataset[tm][k][2],dr])+','
			r1.write(dataline.rstrip(',')+'\n')
		else:
			prev_sample = sample
			prev_well = well
			
	r1.close()

	print "Data munging process is complete\n\n"


# Calls main to kick off process
##################################
if __name__ == '__main__':
	main()