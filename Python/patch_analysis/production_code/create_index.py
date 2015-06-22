#!/usr/bin/env python
# encoding: utf-8

import sys, os

files = os.listdir(os.getcwd())

os.mkdir('gp120')
os.mkdir('gp41')
OUTgp41 = open('gp41/result_index','w')
OUTgp120 = open('gp120/result_index','w')

for f in files:
	if f[-3:] == 'csv':
		c = 0
		region =''
		result_string=''
		INFILE = open(f,'r')
		for line in INFILE:
			if c == 2:
				a = line.split(' ')
				result_string = f+'\t'+a[1]+'\n'
				region = a[1].split('/')[-1:][0].split('_')[0]
			c += 1
		INFILE.close()
	        if region =='gp41':   
			print "moving file ",f," to gp41......" 
			os.rename(f,"gp41/"+f)
			OUTgp41.write(result_string)
		elif region =='gp120': 
			print "moving file ",f," to gp120.."
			os.rename(f,"gp120/"+f)
			OUTgp120.write(result_string)
OUTgp41.close()
OUTgp120.close()
