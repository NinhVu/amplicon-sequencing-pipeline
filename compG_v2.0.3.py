#!/usr/bin/python3
# compG_v1.0.py, 170216, ninh - compiles genotype data and generates result files.
# compG requires a barcode file and all genos file with meta data in barcode file. barcode file name must begin with 'barcode', 'bc' or 'BC'.
# compG_v2.0.py adds ability to export only successful fillin genotypes, create qcGeno or fGeno files if necessary and compatibility with other labs
# compG_v2.0.1.py adds rerun to list of genotype files (170731)
# compG_v2.0.2.py move all functions into salmonBioInfo module
# in terminal $ compG_v2.0.2.py<enter>

import os, glob, re
from salmonBioInfo import *

def Main():
	# process barcode file and return meta data
	bcDict = barcode()
	
	# read geno files and store gCount, genotypes and loci list
	gCountDict, genotypeDict, lociList = readGeno()
	snpLociL = []
	for i in lociList:
		snpLociL.append(i + '-A1')
		snpLociL.append(i + '-A2')

	# calculate success rate and determine if fillin or QC samples are present
	sucRateDict = iGeno(genotypeDict, lociList)	# pass to iGeno to calculate success rate
	dataHasInitial, dataHasRerun, dataHasFillin, dataHasQC = None, None, None, None
	
	for fullName in bcDict:	# find out if data has initial
		if fullName.startswith("In") or fullName.startswith("in"):
			dataHasInitial = 'yes'
			break
	for fullName in bcDict:	# find out if data has rerun
		if fullName.startswith("R") or fullName.startswith("r"):
			dataHasRerun = 'yes'
			break
	for fullName in bcDict:	# find out if data has fillins
		if fullName.startswith("f") or fullName.startswith("F"):
			dataHasFillin = input('\nEnter percent success threshold for fillin samples you consider to be successfully genotyped e.g. 90 or 85: ')
			break
	for fullName in bcDict:	# find out if data has QC
		if fullName.startswith("qc") or fullName.startswith("QC"):
			dataHasQC = 'yes'
			break
			
	# create output file names
	print('\n\nThis program uses personalize prefix as part of output file names.')
	prefixName = input('Enter a discriptive prefix e.g. T1-T9: ')
	prefixName = re.sub(' .!@#$%^&*()', '_', prefixName)
	
	dirPath = os.getcwd().split('/')
	if dirPath[-1].startswith(("L")): dirName = dirPath[-1]
	else: dirName = dirPath[-2] + "_" + dirPath[-1]
	
	reportName = prefixName + "_" + dirName + "_report.csv"
	gCountName = prefixName + "_" + dirName + "_gCount.csv"
	dGapsName = prefixName + "_" + dirName + "_dGaps.csv"
	dGapsMapName = prefixName + "_" + dirName + "_dGapsMap.csv"
	iGenoName = prefixName + "_" + dirName + "_iGeno.csv"
	rrGenoName = prefixName + "_" + dirName + "_rrGeno.csv"
	fGenoName = prefixName + "_" + dirName + "_ALL_fGeno.csv"
	fSucGenoName = prefixName + "_" + dirName + "_Successful_fGeno.csv"
	qcGenoName = prefixName + "_" + dirName + "_qcGeno.csv"
	dFuncName = prefixName + "_" + dirName + "_Dfunc.csv"

	# write genotype fileS
	reportFile = open(reportName, 'w')
	gCountFile = open(gCountName, 'w')
	dGapsFile = open(dGapsName, 'w')
	dGapsMapFile = open(dGapsMapName, 'w')
	dFuncFile = open(dFuncName, 'w')
	genosHeader= ['sampleName'] + lociList + ['\r\n']	# header row data
	
	if dataHasInitial:
		iGenoFile = open(iGenoName, 'w')
		iGenoFile.write(','.join(genosHeader))
	if dataHasRerun:
		rrGenoFile = open(rrGenoName, 'w')
		rrGenoFile.write(','.join(genosHeader))
	if dataHasFillin:
		fGenoFile = open(fGenoName, 'w')
		fGenoFile.write(','.join(genosHeader))
		fSucGenoFile = open(fSucGenoName, 'w')
		fSucGenoFile.write(','.join(genosHeader))
	if dataHasQC:
		qcGenoFile = open(qcGenoName, 'w')
		qcGenoFile.write(','.join(genosHeader))
	
	for fullName, geno in sorted(genotypeDict.items()):
		samType = bcDict.get(fullName, None)
		progenyName = fullName[len(samType[0]):]
		rowData = [progenyName] + geno
		if dataHasFillin and (fullName.startswith(('f')) or fullName.startswith(('F'))):
			fGenoFile.write(','.join((rowData)))
			fGenoFile.write('\r\n') 
			numFailed = geno.count('0') + geno.count('00')
			if (1 - numFailed/len(geno)) >= int(dataHasFillin)/100:
				fSucGenoFile.write(','.join((rowData)))
				fSucGenoFile.write('\r\n')
		if dataHasQC and (fullName.startswith(('QC')) or fullName.startswith(('qc'))):
			qcGenoFile.write(','.join((rowData)))
			qcGenoFile.write('\r\n')
		if dataHasInitial and (fullName.startswith(('In')) or fullName.startswith(('in'))):
			iGenoFile.write(','.join((rowData)))
			iGenoFile.write('\r\n')
		if dataHasRerun and (fullName.startswith(('R')) or fullName.startswith(('r'))):
			rrGenoFile.write(','.join((rowData)))
			rrGenoFile.write('\r\n')
	
	# write ONE Dfunc file with no header
	for fullName, barcodeD in sorted(bcDict.items()):
		barcodeD = list(map(str, barcodeD))
		indiGeno = genotypeDict.get(fullName)
		rowData = [fullName] + barcodeD + indiGeno
		dFuncFile.write(','.join((rowData)))
		dFuncFile.write('\r\n')

	# write gCount
	gCountHeader= ['FullName', 'SampleType', 'SampleStatus', 'PlateID',	'i7_name', 'i7_sequence', 'i5_name', 'i5_sequence', '%Successful', '# of failed ALL loci', '%subPanel Successful', '# of failed subPanel loci', 'Raw reads', 'On-Target reads', '% On-Target', 'IFI score'] + snpLociL
	gCountFile.write(','.join(gCountHeader))
	gCountFile.write('\r\n')
	for fullName, bcMetData in sorted(bcDict.items()):
		indiSucData = sucRateDict.get(fullName)
		indiGcount = gCountDict.get(fullName)
		rowData = [fullName] + list(map(str, bcMetData)) + list(map(str, indiSucData)) + list(map(str, indiGcount))
		gCountFile.write(','.join((rowData)))
		gCountFile.write('\r\n')

	# write data gaps
	dGapsDict ={}	# for use data gaps map as well
	dGapsHeader = ['Sample name', 'Sample type', 'Sample status', '%Successful', '#failed ALL loci', '%subpanel Success', '#failed subpanel loci', 'PlateID', 'Well\r\n']
	dGapsFile.write(','.join(dGapsHeader))
	for fullName, sucRate in sorted(sucRateDict.items()):
		wellPossi = [] # for use data gaps map as well
		bcMetData = bcDict.get(fullName)
		progenyName = fullName[len(bcMetData[0]):]
		try:
			if 'NTC' not in fullName and (sucRate[0] < 90 or sucRate[2] < 90):
				rowData = [progenyName] + bcMetData[:2] + list(map(str, sucRate)) + [str(bcMetData[2])] + [str(bcMetData[5])]
				dGapsFile.write(','.join((rowData)))
				dGapsFile.write('\r\n')
				wellPossi.extend([bcMetData[2], bcMetData[5]])
				dGapsDict[fullName] = wellPossi
		except TypeError:
			if 'NTC' not in fullName and (sucRate[0] < 90 and sucRate[2] == 'na'):
				rowData = [progenyName] + bcMetData[:2] + list(map(str, sucRate)) + [str(bcMetData[2])] + [str(bcMetData[5])]
				dGapsFile.write(','.join((rowData)))
				dGapsFile.write('\r\n')
				wellPossi.extend([bcMetData[2], bcMetData[5]])
				dGapsDict[fullName] = wellPossi

	# write data gaps map
	rowName = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	trayL = []	# also for use in report file
	for fullName, loco in sorted(bcDict.items()): trayL.append(loco[2])	# grab list of trays from bcDict dictionary
	trayL = sorted(set(trayL)) # sort and remove duplicate from tray list

		# cycle through all trays for composite look
	dGapsMapFile.write(''.join((dGapsMapName,'\r\n','Composite/concatenated table of failed samples for each well\r\n')))
	compositeHeader = [dirName] + list(map(str, list(range(1,13)))) + ['\r\n']
	dGapsMapFile.write(','.join((compositeHeader)))
	for rowLetter in rowName: # cycle through rows
		dGapsMapFile.write(','.join((rowLetter)))
		dGapsMapFile.write(',')
		for wellNum in range(1,13): # cycle through body of row
			stackWellCount = 0
			for fullName, loco in dGapsDict.items():	# inner most loop - cycle through each well
				if rowLetter == loco[1][:1] and wellNum == int(loco[1][1:]):
					stackWellCount +=1
			if stackWellCount == 0: dGapsMapFile.write(',')
			else:
				dGapsMapFile.write(str(stackWellCount))
				dGapsMapFile.write(',')
		dGapsMapFile.write('\r\n')

		# cycle through each tray
	dGapsMapFile.write('\r\nTray tables with each failed sample in it\'s well position according to barcode file.\r\n')
	count=0
	for tray in trayL:
		trayHeader = [tray] + list(map(str, list(range(1,13))))
		dGapsMapFile.write(','.join((trayHeader)))
		dGapsMapFile.write('\r\n')
		for rowLetter in rowName:	# cycle through rows
			dGapsMapFile.write(','.join((rowLetter)))
			dGapsMapFile.write(',')
			for w in range(1,13):	# cycle through body of row
				aMatch = 0
				for fullName, loco in dGapsDict.items():	# inner most loop - cycle through each well
					successRate = sucRateDict.get(fullName)
					if loco[0] == tray and int(loco[1][1:]) == w and loco[1][:1] == rowLetter:
						dGapsMapFile.write(fullName)
						dGapsMapFile.write(',')
						aMatch = 1
				if aMatch ==0: dGapsMapFile.write(',')
			dGapsMapFile.write('\r\n')	# return to next row

	# write report summary
	reportFile.write('Summary of trays excluding NTCs.\r\n')
	reportFile.write('PlateID,raw reads,On-Target reads,ave. % On-Target,N,#failed samples-ALLloci,%Success-AllLoci,#failed samples-subpanel,%Success-subpanel\r\n')
	rawReads, onTreads, aveOnT, failSamAllLoci, failSamPBT = [],[],[],[],[]
	totalSample, totalOnTarget, totalLociFailed, totalPBTfailed = 0, 0, 0, 0
	for tray in trayL:	# cycle through each tray
		rawReadsTray, onTreadsTray, aveOnTTray, failSamAllLociTray, failSamPBTTray, samCount = 0,0,0,0,0,0
		for fullName, metData in bcDict.items():
			if metData[2] == tray and 'NTC' not in fullName:
				gCountData = gCountDict.get(fullName)
				sucRate = sucRateDict.get(fullName)
				rawReadsTray = rawReadsTray + gCountData[0]
				onTreadsTray = onTreadsTray + gCountData[1]
				aveOnTTray = aveOnTTray + gCountData[2]
				samCount+=1
				if sucRate[0] <90: failSamAllLociTray = failSamAllLociTray + 1
				perSucAllLoci = round((1-(failSamAllLociTray/samCount))*100, 1)
				try:
					if sucRate[2] <90: failSamPBTTray = failSamPBTTray + 1
					perSucPBT96 = round((1-(failSamPBTTray/samCount))*100, 1)
				except TypeError: perSucPBT96 = 'na'
		reportFile.write(','.join((str(tray), str(rawReadsTray), str(onTreadsTray), str(round(aveOnTTray/samCount,1)), str(samCount), str(failSamAllLociTray), str(perSucAllLoci), str(failSamPBTTray), str(perSucPBT96))))
		reportFile.write('\r\n')

		# total and average for project in library
		totalSample = totalSample + samCount
		totalOnTarget = totalOnTarget + (aveOnTTray/samCount)
		totalLociFailed = totalLociFailed + failSamAllLociTray
		totalPBTfailed = totalPBTfailed + failSamPBTTray
	reportFile.write(','.join((str(len(trayL)), 'TOTAL-AVERAGE,', str(round(totalOnTarget/len(trayL),1)), str(totalSample), str(totalLociFailed), str(round(100*(1-(totalLociFailed/totalSample)), 1)), str(totalPBTfailed), str(round(100*(1-(totalPBTfailed/totalSample)),1)))))
	reportFile.write('\n')	
		
	# concatenate contamination table to write to report.csv
	reportFile.write('\r\n\r\nContamination summary of samples: IFI score greater than 2.5 are likely contaminated.\r\n')
	reportFile.write('NTCs are included regardless of IFI scores.\r\n')
	reportFile.write('Failed samples have tendency to appear on this list as well.\r\n\r\n')
	reportFile.write('Sample name,Sample type,Sample status,%Successful,#failed ALL loci,%subpanel Success,#failed subpanel loci,PlateID,Well,IFI score\r\n')
	for fullName, gCountData in sorted(gCountDict.items()):
		bcMetData = bcDict.get(fullName)
		sucMetData = sucRateDict.get(fullName)
		progenyName = fullName[len(bcMetData[0]):]
		rowData = [progenyName] + bcMetData[:2] + list(map(str, sucMetData)) + [str(bcMetData[2])] + [str(bcMetData[5])] + [str(gCountData[3])]
		if 'NTC' in fullName or gCountData[3] > 2.50:
			reportFile.write(','.join((rowData)))
			reportFile.write('\r\n')

	# close all output files
	for file in glob.glob("*.csv"):
		with open(file) as x: x.close()
	
Main()
