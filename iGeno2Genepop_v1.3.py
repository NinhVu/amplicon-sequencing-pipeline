#!/usr/bin/python3
# iGeno2Genepop_v1.3.py - Program takes iGeno file and create genepop file with Windows EOL format
# argv[1] - iGeno file w/o metadata. It must have a header row with locus names.
# agrv[2] - optional - a file with subset of panel, usually all markers except for sex or hygrid markers. One marker per line 

import sys

def Main():
	# get smaller panel loci list
	userInputChoices = ['ots', 'omy', 'sfo', 'yct', 'ots298', 'omy268', 'sfo240', 'yct161']
	args = sys.argv[2:]
	if args:
		with open(sys.argv[2]) as locus: panelLociList = locus.read().splitlines()
	else:
		print('\nWhat panel would you like to use to convert iGeno data to genepop format?')
		while True:
			try:
				userPanel = input('Your choices are omy268, ots298, sfo240 or yct161: ')
				if userPanel in userInputChoices: break
			except ValueError: continue
		if 'ots' in userPanel:
			with open('/home/efglserv/software/GTseq/ots298.txt') as locus: panelLociList = locus.read().splitlines()
		elif 'omy' in userPanel: 
			with open('/home/efglserv/software/GTseq/omy268.txt') as locus: panelLociList = locus.read().splitlines()
		elif 'sfo' in userPanel: 
			with open('/home/efglserv/software/GTseq/sfo240.txt') as locus: panelLociList = locus.read().splitlines()
		elif 'yct' in userPanel: 
			with open('/home/efglserv/software/GTseq/yct161.txt') as locus: panelLociList = locus.read().splitlines()
	panelLociList = list(filter(None, panelLociList))	# remove empty items or lines
	
	# create index of desire loci and read genotype from iGeno
	panelIndex, sampleData, pedList, sampleList=[],{},[],[]
	with open(sys.argv[1]) as geno: genoData = geno.read().splitlines()
	lociRow = genoData[0].split(',')
	lociRow = lociRow[1:]
	try:
		for i in panelLociList: panelIndex.append(lociRow.index(i))	# index desire loci
	except ValueError:
		print('\nYour panel (input) and iGeno data do not match. Try running program again with matching panel and iGeno data.\n')
		sys.exit()
	
	# put data into a dictionary but first ask user if they want to include failed samples
	print('\nWant to include failed samples? If not leave input blank.')
	while True:
		percentMiss = input('Enter percent missing data allow e.g. 0.10 for 10%: ')
		try:
			float(percentMiss)
			break
		except ValueError:
			if len(percentMiss) ==0:
				percentMiss = 0
				break
			else: continue
	count, numOfSam, noExlusion, NTCs = 0,0,0,0
	if float(percentMiss) >=1: percentMiss = 0
	for i in genoData[1:]:
		rowItems = i.split(',')
		numOfFails = rowItems[1:].count('00')
		if numOfFails/(len(panelLociList)+1) <= float(percentMiss) and 'NTC' not in rowItems[0]:
			sampleData[rowItems[0]] = rowItems[1:]	# create dictionary with failure rate samples only
			pedList.append(rowItems[0][:10])	# IDFG pedigree name convention e.g. OtsLGRA16C
			sampleList.append(rowItems[0])
			numOfSam+=1
		elif percentMiss == 0 and 'NTC' not in rowItems[0]:
			sampleData[rowItems[0]] = rowItems[1:]	# append to dictionary samples withOUT failure rate
			pedList.append(rowItems[0][:10])	# IDFG pedigree name convention e.g. OtsLGRA16C
			sampleList.append(rowItems[0])
			noExlusion+=1
		elif 'NTC' in rowItems[0]: NTCs+=1
		count+=1
	if percentMiss == 0: numOfSam = numOfSam + noExlusion
	print('\nThere are {0} NTCs in xxx_iGeno.csv file. They are excluded from genepop file.'.format(NTCs))
	if percentMiss == 0: print('There are {0} samples in genepop file.\n'.format(numOfSam))
	if percentMiss != 0:
		print('With more than {0}% missing data, {1} samples are excluded from genepop file.'.format(round(100*float(percentMiss), 1), count - NTCs - numOfSam))
		print('The remaining {0} samples in genepop file have less than or equal to {1}% missing data.\n'.format(numOfSam, round(100*float(percentMiss), 1)))
		
	# write genepop file
	pedList = list(set(pedList))	# remove duplicate pedigree
	pedList.sort()	# sort pedigree list
	if len(pedList) > 1: pedText = ' pedigrees'
	else: pedText = ' pedigree'
	genepopName = sys.argv[1].replace('iGeno.csv', 'genepop.txt')
	genePopData = open(genepopName, 'w')
	genePopData.write("".join((str(len(sampleList)),' samples and ', str(len(pedList)), pedText,'\r\n')))	# write title
	for i in panelLociList: genePopData.write("".join((i,"\r\n")))	# write list of markers
	
	for pedigree in pedList:	# write body og genepop file
		genePopData.write('pop\r\n')
		for i in sampleList:
			if pedigree in i:
				genePopData.write(''.join((i,', ')))
				sampleGeno = sampleData.get(i)
				numOfGeno=0
				for dex in panelIndex:
					genoC = sampleGeno[dex]
					geno = list(genoC)
					for allele in geno:
						if allele =='A': genePopData.write('01')
						elif allele =='C': genePopData.write('02')
						elif allele =='G': genePopData.write('03')
						elif allele =='T': genePopData.write('04')
						elif allele ==':' or allele == '-': genePopData.write('05')
						elif allele =='0': genePopData.write('00')
						else: genePopData.write('00')
					
					numOfGeno+=1
					if numOfGeno == len(panelLociList): genePopData.write('\r\n')	# do not leave a space after last allele
					else: genePopData.write(' ')	# else leave a space between genotypes

	print('Your output file {0} is ready.\n'.format(genepopName))
	genePopData.close()
	
Main()