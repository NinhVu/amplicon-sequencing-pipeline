#!/usr/bin/python3
# oddMark_v1.py scores indel and/or present/absent markers
# 20170725, ninh
# It requires similar input file as ProbeSeq except with targetPosition instead of correction factors
# in terminal: oddMark_v1.py

import glob, sys, re, multiprocessing, math
from multiprocessing import Process
from itertools import islice	# for reading multiple lines at once

def subHyGeno(subSetSam, markerDict):
	for marker in sorted(markerDict):
		locInfo = markerDict.get(marker)
		for file in subSetSam:
			fq = open(file, 'r')
			fwd, a1, a2 = 0,0,0
			title = fq.readline()
			fastA = fq.readline()
			
			# read and count reads
			while fastA:
				if locInfo[2] in fastA: fwd +=1
				if locInfo[2] in fastA and locInfo[0] in fastA[int(locInfo[3]):(len(locInfo[0])+int(locInfo[3])+1)]: a1 +=1
				if locInfo[2] in fastA and locInfo[1] in fastA[int(locInfo[3]):(len(locInfo[1])+int(locInfo[3])+1)]: a2 +=1
				data = list(islice(fq, 3))	# read the next three lines
				fastA = fq.readline()
		
			ratA1A2 = 0	# ratio of allele 1 & 2 is defaulted to zero
			geno = '00'	# genotype is defaulted to '00' or noCall
			
			if a1 + a2 >= 10 and locInfo[4] != 'X':	# score threshold of hybrid marker
				try: ratA1A2 = abs(a1/a2)
				except ZeroDivisionError: ratA1A2 = a1
				if ratA1A2 >= 15: geno = locInfo[4] + locInfo[4]
				elif ratA1A2 <= 10 and ratA1A2 >= 0.1: geno = locInfo[4] + locInfo[5]
				elif ratA1A2 < 0.1: geno = locInfo[5] + locInfo[5]
				
			elif fwd >= 3 and locInfo[4] == 'X' and 'NTC' not in file:	# score threshold of sex marker
				if a2 >= 20: geno = locInfo[4] + locInfo[5]	# call M
				elif a2 < 5 : geno = locInfo[4] + locInfo[4]	# call F
				
			# assign and appending calls to genos file
			genoCall = '00'	# genotype is defaulted to '00' unless 1 of 2 condition above is meet
			if geno[0] == locInfo[4] and geno[1] == locInfo[4]: genoCall = 'A1HOM'
			elif geno[0] == locInfo[4] and  geno[1] == locInfo[5]: genoCall = 'HET'
			elif geno[0] == locInfo[5] and geno[1] == locInfo[5]: genoCall = 'A2HOM'
			try: onTargetRat = ((a1+a2)/fwd)*100
			except ZeroDivisionError: onTargetRat = 0
			gName = str(file)
			gName = gName.replace('fastq', 'genos')
			genos = open(gName, 'a')
			a1T = locInfo[4] + '=' + str(a1)
			a2T = locInfo[5] + '=' + str(a2)
			genos.write(','.join((marker, a1T, a2T, str(round(ratA1A2,3)), geno, genoCall, '0,0', str(fwd), str(round(onTargetRat,1)), '0.5,Pass,0\n')))
			fq.close(), genos.close()

		
def Main():
	# read Probseq file and create lists
	with open('/home/efglserv/software/GTseq/probeSeq/CoCut_oddMarkers_ProbeSeqs.csv') as probeSeq:
		panelLociInfo = probeSeq.read().splitlines()
	while '' in panelLociInfo: panelLociInfo.remove('')	# remove empty items in list
	markerDict = {}
	for i in panelLociInfo:
		rowDat = i.split('\t')
		markerDict[rowDat[0]] = rowDat[1:]

	# divide total number of fastq files into processes
	fqFileL = []
	for file in sorted(glob.glob('*.fastq')): fqFileL.append(file)
	# multiprocess filtering subsets of fastqs
	if len(fqFileL) <= multiprocessing.cpu_count(): samPerProc = 1
	else: samPerProc = math.ceil(len(fqFileL)/multiprocessing.cpu_count())
	fqPerProc = [ fqFileL[i:i + samPerProc] for i in range(0, len(fqFileL), samPerProc)]	# divide fastq list into chunks for each processor	

	perP = {}
	for i in range(len(fqPerProc)): perP[i] = 's' + str(i)	# create dynamic variable names for each Process
	print('There are {0} processes available for scoring hybrid markers. Each is scoring ~{1:,} samples of {2:,} total samples.'.format(len(perP), len(fqPerProc[0]), len(fqFileL)))
	for i in range(len(fqPerProc)):
		perP[i] = Process(target=subHyGeno, args=(fqPerProc[i], markerDict))
		perP[i].start()
	for i in range(len(perP)): perP[i].join()	# join processes when all is done

Main()
