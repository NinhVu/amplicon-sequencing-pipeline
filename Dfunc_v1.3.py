#!/usr/bin/python3
# Dfunc_v1.0, 170216, ninh - searches for duplicates using Dfunc input files. Dfunc assumes input files contain six columns of meta data.
# Dfunc_v1.0.4 move dSearch function into salmonBioInfo module
# Dfunc_v1.2 clean up multiprocessing and uses ALL processes if applicable
# Dfunc_v1.3 adds error rate and mismatch loci
# in terminal $ Dfunc_v1.0.4.py<enter>


import os, glob, sys, re, time, math, multiprocessing, colorama
from colorama import Fore
from multiprocessing import Process
from salmonBioInfo import *


def errorRate(allDupResults, alikeSum, unlikeSum, hetHomoSum, homoHomoSum):
	if alikeSum == 0:
		print('no duplicates found')
		pass
	
	else:
		perAlleleRate = 100 * ((hetHomoSum + (homoHomoSum*2))/(2*(alikeSum + unlikeSum)))
		if len(allDupResults) ==1:
			print('There is one duplicate pair.')
			print("The genotyping error rate for this pair of duplicates is {:5}%".format(round(100*(unlikeSum/(alikeSum + unlikeSum)),4)))
			print("Per allele error rate is {:5}%".format(round(perAlleleRate, 4)))
			
		else:
			numLength = len(str(len(allDupResults)))
			print(Fore.GREEN + "There are {num:{width},} duplicate pairs.".format(num=len(allDupResults), width=numLength + 1))
			print("For all duplicate pairs, there are a total of {num:{width},} comparisons.".format(num=alikeSum + unlikeSum, width=numLength+2))
			print("The genotyping error rate is {:5}%".format(round(100*(unlikeSum/(alikeSum + unlikeSum)),4)))
			print("Per allele error rate is {:5}%".format(round(perAlleleRate, 4)))
			
	return None


def Main():

	# get directory
	dirPath = os.getcwd().split('/')
	dirName = dirPath[-1]
	
	
	# get all Dfunc fileS into one dictionary and create sample list
	all_data, dFuncSamDict, allSam = [], {}, []
	for file in glob.glob("*Dfunc.csv"):
		with open(file) as f: oneFile = f.read().splitlines()
		all_data = all_data + oneFile

	successRate = input(Fore.GREEN + '\nEnter genotype success rate a sample must have for it to be in duplicate search e.g. 90 for 90%: ' + Fore.RED)
	if successRate == '': successRate = 0
	if float(successRate) >= 100: successRate = 100
	for oneSam in all_data:
		rDat = oneSam.split(',')
		numOfFailed = 0
		for geno in rDat[8:]:
			try:
				if int(geno) == 0: numOfFailed+=1
			except: pass
		if 1 - (numOfFailed/len(rDat[8:])) >= float(successRate)*0.01:
			dFuncSamDict[rDat[0]] = rDat[1:]
			allSam.append(rDat[0])
		else: pass
		
	totalComp = int(((len(allSam)-1)*len(allSam))/2)
	
	# multiprocess duplicate search of subsets of fastqs
	perSimilar = eval(input(Fore.GREEN + '\nEnter percent similarity threshold for a pair of samples to be consider duplicates e.g. 90 for 90%. My recommendation is 95: ' + Fore.RED))
	startTime = time.time()
	perSimilar = perSimilar/100
	if len(allSam) <= multiprocessing.cpu_count(): samPerProc = 1
	else: samPerProc = math.ceil(len(allSam)/multiprocessing.cpu_count())
	subSetPerProc = [ allSam[i:i + samPerProc] for i in range(0, len(allSam), samPerProc)]	# divide fastq list into chunks for each processor
	if len(subSetPerProc[-1]) == 1: subSetPerProc = subSetPerProc[:-1]	# if the last chuck has one item then erase - no need to compare to itself
	
	perP, smallerSamLs = {}, []
	for L in subSetPerProc:	# create smaller allSam list of lists
		if len(L) == 1: smallerSamLs.append(tuple(allSam[allSam.index(L[0]) + 1:]))	# account for list smaller than 56 samples
		else: smallerSamLs.append(tuple(allSam[allSam.index(L[1]):]))	# account for list larger than 56 samples

	for i in range(len(subSetPerProc)): perP[i] = 'pDup' + str(i + 1)	# create dynamic variable names for each Process
	print(Fore.GREEN + 'There are {0} processes available. Each is searching duplicates for ~{1:,} samples of {2:,} total samples.'.format(len(perP), len(subSetPerProc[0]), len(allSam)))
	print('There are a total of {0:,} comparisons.\n'.format(totalComp))
	for i in range(len(subSetPerProc)):
		perP[i] = Process(target=dSearch, args=(perP[i], subSetPerProc[i], smallerSamLs[i], dFuncSamDict, perSimilar))
		perP[i].start()
	for i in range(len(perP)): perP[i].join()	# join processes when all is done
		
	# calculate runtime right after processes join
	hours, rem = divmod(time.time()-startTime, 3600)
	minutes, seconds = divmod(rem, 60)


	# concatenate output subfiles into one list and delete all subfiles
	allDupResults = []
	for file in sorted(glob.glob("pDup*")):
		with open(file) as f: eachFileDup = f.read().splitlines()
		allDupResults.extend(eachFileDup)
	for file in glob.glob("pDup*"): os.remove(file)


	# account for nonmatching QC samples and get metrics for error rate calculation
	allDupSam, noMatchQC = [],[]
	alikeSum, unlikeSum, hetHomoSum, homoHomoSum = 0,0,0,0
	for i in allDupResults:
		rDat = i.split(',')
		allDupSam.append(rDat[0])
		allDupSam.append(rDat[1])
		
		alikeSum = alikeSum + int(rDat[3])
		unlikeSum = unlikeSum + int(rDat[4])
		hetHomoSum = hetHomoSum + int(rDat[6])
		homoHomoSum = homoHomoSum + int(rDat[7])
	for fullName, type in sorted(dFuncSamDict.items()):
		if (type[0] == 'qc' or type[0] == 'QC') and fullName not in allDupSam and 'NTC' not in fullName: noMatchQC.append(fullName)
	
	
	# display error rates
	errorRate(allDupResults, alikeSum, unlikeSum, hetHomoSum, homoHomoSum)

	
	# write one duplicate results file
	prefixName = input('\nEnter a descriptive prefix for result file name e.g. T1-12: ')
	if prefixName is None: dRName = dirName + "_dupSearch_results.csv"
	else: dRName = prefixName + '_' + dirName + "_dupSearch_results.csv"
	dupOut = open(dRName, 'w')
	dupOut.write('Sample:A,%Success:A,SampleType:A,SampleStatus:A,PlateID:A,Well:A,sample:B,%Success:B,SampleType:B,SampleStatus:B,PlateID:B,Well:B,%alike,#alike_loci,#unlike_loci,#non-comparisons,concordance,samNumApart,HetHomoD,HomoHomoD\n')
	for i in allDupResults:
		rDat = i.split(',')
		samA = dFuncSamDict.get(rDat[0])
		samB = dFuncSamDict.get(rDat[1])
		pLen = len(samA)-7
		sucRateA = round(((pLen - (samA[7:].count('00') + samA[7:].count('0')))/pLen)*100,1)
		sucRateB = round(((pLen - (samB[7:].count('00') + samB[7:].count('0')))/pLen)*100,1)
		proGNameA = rDat[0][len(samA[0]):]
		proGNameB = rDat[1][len(samB[0]):]
		if proGNameA == proGNameB: concordance, samNumApart = 'yes', ''
		else:
			numA, numB = proGNameA.split('_'), proGNameB.split('_')
			try: concordance, samNumApart = 'no', abs(int(numA[-1]) - int(numB[-1]))
			except ValueError: concordance, samNumApart = 'no', 0
		
		hetHOMO, homoHOMO = str(rDat[6]), str(rDat[7])
		if int(rDat[6]) == 0: hetHOMO = ''
		if int(rDat[7]) == 0: homoHOMO = ''
		rDat = [proGNameA] + [str(sucRateA)] + list(map(str, samA[:3])) + [str(samA[5])] + [proGNameB] + [str(sucRateB)] + list(map(str, samB[:3])) + [str(samB[5])] + list(map(str, rDat[2:-2])) + [concordance] + [str(samNumApart)] + [hetHOMO] + [homoHOMO]	# rowData is a list
		dupOut.write(','.join((rDat)))
		dupOut.write('\r\n')

	# write list of qc samples with no matching initial/fillin samples
	if noMatchQC:	
		dupOut.write('\r\nThese QC samples do not have matching initial or fillin samples.\r\n')
		dupOut.write('qcSample,SampleType,SampleStatus,PlateID,Well,%Success\r\n')
		for fullName in noMatchQC:
			metData = dFuncSamDict.get(fullName)
			progenyName = fullName[len(metData[0]):]
			pLen = len(metData)-7
			sucRate = round(((pLen - (metData[7:].count('00') + metData[7:].count('0')))/pLen)*100,1)
			rDat = [progenyName] + list(map(str, metData[:3])) + [metData[5]] + [str(sucRate)]
			dupOut.write(','.join(rDat))
			dupOut.write('\r\n')
	else: print('\nAll QC samples have matches or there are no QC samples in your Dfunc files.\n')
	

	# print run time
	print(Fore.BLUE + "\nTotal compute time {:0>2}:{:0>2}:{:05.2f}\n".format(int(hours),int(minutes),seconds) + Fore.WHITE)	

Main()



