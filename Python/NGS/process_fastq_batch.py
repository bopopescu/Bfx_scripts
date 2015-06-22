#!/usr/bin/env python
# encoding: utf-8

import os,sys
import subprocess

fn = ["IonXpress_003_SN1-19-Gilead.fastq",
	"IonXpress_003_SN1-20-Gilead.fastq",
	"IonXpress_004_SN1-19-Gilead.fastq",
	"IonXpress_005_SN1-19-Gilead.fastq",
	"IonXpress_006_SN1-19-Gilead.fastq",
	"IonXpress_007_SN1-19-Gilead.fastq",
	"IonXpress_008_SN1-19-Gilead.fastq",
	"IonXpress_009_SN1-19-Gilead.fastq",
	"IonXpress_010_SN1-20-Gilead.fastq",
	"IonXpress_011_SN1-20_Gilead.fastq",
	"IonXpress_012_SN1-19-Gilead.fastq",
	"IonXpress_012_SN1-20-Gilead.fastq",
	"IonXpress_013_SN1-19-Gilead.fastq",
	"IonXpress_014_SN1-20-Gilead.fastq",
	"IonXpress_015_SN1-20-Gilead.fastq",
	"IonXpress_016_SN1-20-Gilead.fastq"]

path = os.getcwd()+"/"
for f in fn:
	print "Beginning processing of ",f
	subprocess.call(["python", "ns5b_fastq_trim.py", "-s",f])
	print "moving on ..."

print "All files completed\n\n"