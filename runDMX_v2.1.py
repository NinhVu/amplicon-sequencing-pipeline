#!/usr/bin/python3
# runDMX_v1.0.py, 170203, ninh - a program to demultiplex dual indexed NGS fastq file from NextSeq.
# runDMX_v2.0.py, 170511 - remove barcode error correction code
# runDMX_v2.0.1.py, 170807 - move DMX to salmonBioInfo module
# runDMX_v2.1.py, 170212 - chunking a list sample into smaller subsets
# in terminal $ runDMX_v1.0.py <barcode file> <multiplex fastq file>

import sys, re, os, math, time, multiprocessing
from multiprocessing import Process
from itertools import islice
from salmonBioInfo import DMX

def Main():
	start = time.time()
	# read & create barcode metadata
	barcodeInfo, sampleName, iCombo = [],[],[]
	with open(sys.argv[1]) as bc: barcodeFile = bc.read().splitlines()

	for i in barcodeFile[1:]:
		dataRow = re.split('\t|,|"', i)	# split either by tab, comma, " or space
		dataRow = list(filter(None, dataRow)) 	# delete empty rows in barcode list
		if len(dataRow) == 8:
			sampleName.append(dataRow[1] + dataRow[0] + '.fastq')
			iCombo.append(dataRow[5]+ '+' + dataRow[7])	# unique to bcl2fastq and NextSeq
		else: pass	# prevent error if row have less than 8 items

	# demultiplex subsets of samples = optimal is 14 processes
	if len(sampleName) <= 14: samPerProc = 1
	else: samPerProc = math.ceil(len(sampleName)/14)
	fqPerProc = [ sampleName[i:i + samPerProc] for i in range(0, len(sampleName), samPerProc)]	# divide fastq list into chunks for each processor
	indexPerProc = [ iCombo[i:i + samPerProc] for i in range(0, len(iCombo), samPerProc)]	# divide indices list into chunks for each processor

	perP = {}
	for i in range(len(fqPerProc)): perP[i] = str(i + 1)	# create dynamic variable names for each Process
	print('\nThere are {0} processes available. Each is demultiplexing ~{1:,} samples of {2:,} total samples.\n'.format(len(perP), len(fqPerProc[0]), len(sampleName)))
	for i in range(len(fqPerProc)):
		perP[i] = Process(target=DMX, args=(perP[i], fqPerProc[i], indexPerProc[i]))
		perP[i].start()
	for i in range(len(perP)): perP[i].join()	# join processes when all is done to determine run time
	
	hours, rem = divmod(time.time()-start, 3600)
	minutes, seconds = divmod(rem, 60)
	print("Total run time {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
	
Main()