#!/usr/bin/python3

# libSum_v1.1.py, 20171225 n - libSum searches for reads with i7-i5 combinations in one millions reads and display technician, sample count etc. statistics 
# in terminal $ libSum_v1.1.py

import multiprocessing, math, os, gzip, datetime, os.path, time, sys, colorama, collections
from collections import Counter
from colorama import Fore
from multiprocessing import Process
from subprocess import call
from itertools import islice


	
def i5search(reali7, comboB, i5List):
	alli5s, reali5 = [], []
	for c in comboB:
		if c[:6] == reali7 and c[7:] in i5List: alli5s.append(c[7:])
	
	for i in i5List:
		if alli5s.count(i) >10: reali5.append(i)
	
	samCountFile = open('samInTray.csv', 'a')
	samCountFile.write(','.join((reali7, str(len(reali5)))))
	samCountFile.write('\n')


def displayArchiveLib(libDat, now):
	rowDat = libDat.split(',')
	print("\n{:6} has already been archived, however here is the summary.".format(rowDat[0]))
	print('Data gathered from first one million reads.\n')
	print(Fore.RED + "{: <6}{: <8}{: ^8}{: ^11}{: <9}{: ^2}   {}".format('Tray', 'Tech', 'i7Name', 'i7Seq', 'reads', 'N', 'Species'))
	techName = rowDat[4:][::6]
	i7Name = rowDat[4:][1::6]
	i7Seq = rowDat[4:][2::6]
	numReads = rowDat[4:][3::6]
	numSam = rowDat[4:][4::6]
	species = rowDat[4:][5::6]
	
	# create dictionary of data table to be display in order of i7 names
	tableDict = {}
	for i in range(len(i7Name)): tableDict[i7Name[i]] = techName[i], i7Seq[i], numReads[i], numSam[i], species[i]
	trayCount = 1
	for k, v in sorted(tableDict.items()):
		print(Fore.WHITE + "{: <6}{: <8}{: ^8}{: ^11}{: <9}{: ^2}   {}".format(trayCount, v[0], k, v[1], v[2], v[3], v[4]))
		trayCount+=1
	print("...\nNot counting NTCs, there are{:6,} samples in {:<2} trays.\n".format(int(rowDat[3]), int(rowDat[2])))	
	sys.exit()


def callcheckTrayCount(tray2Keep, trayCount):
	try:
		if tray2Keep <= trayCount -1: return(tray2Keep)
		else: return Error
	except: return Error


def Results(uniq7Dict, techi7Dict, i7SpeciesDict, numSamDict, libOfInterest, libName):
	monthDict = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
	
	# create two dictionaries based on i7Name and read count
	i7idDict, readCountDict, totalSamInLib, trayCount = {}, {}, 0, 1
	for k,v in sorted(uniq7Dict.items()):
		i7Info = techi7Dict.get(v)
		if i7SpeciesDict.get(v) is None: speciesID = 'undetermine'
		else: speciesID = i7SpeciesDict.get(v)
		totalSamInLib = totalSamInLib + numSamDict.get(v)
		i7idDict[i7Info[1]] = [i7Info[0], v, k, str(numSamDict.get(v)), speciesID]
		readCountDict[k] = [i7Info[0], i7Info[1], v, str(numSamDict.get(v)), speciesID]
	
	# display FULL results ordered by i7 name or ID
	print(Fore.GREEN + '\nData gathered from first one million reads.', end="")
	print(' Trays with less 1,000 reads are likely not real.\n')
	print(Fore.RED + "{: <6}{: <8}{: ^8}{: ^11}{: <9}{: ^2}   {}".format('Tray', 'Tech', 'i7Name', 'i7Seq', 'reads', 'N', 'Species'))
	trayCount = 1
	for k,v in sorted(i7idDict.items()):
		print(Fore.WHITE + "{: <6}{: <8}{: ^8}{: ^11}{: <9}{: ^2}   {}".format(trayCount, v[0], k, v[1], v[2], v[3], v[4]))
		trayCount+=1
	print("...\nNot counting NTCs, there are{:6,} possible samples in {:<2} possible trays in {:>5}.\n".format(totalSamInLib, trayCount - 1, libOfInterest))

	# ask user how many trays he/she wants to keep based on read count
	attempt =0
	while attempt <=3:
		try:
			tray2Keep = eval(input('Based on READ COUNT, how many trays above are real and should be archive? '))
			userTrayInput = callcheckTrayCount(tray2Keep, trayCount)
			break
		except: 
			if attempt ==2: print('last chance')
			if attempt <3: print('enter a number no greater than the largest tray count number from table above')
			attempt+=1
			continue
	try: userTrayInput
	except:
		print('You ran out of chances. The data is not archive. Try running program again.\n')
		sys.exit()

	# create archive dictionary and new display dictionary based on what trays to keep
	tray2KeepCount, totalSam2Keep, archiveUniq7Dict, newDisplayDict = 1, 0, {}, {}
	for k,v in sorted(readCountDict.items(), reverse=True):
		archiveUniq7Dict[k] = [v[0], v[1], v[2], v[3], v[4]]
		newDisplayDict[v[1]] = [v[0], v[1], k, v[3], v[4]]
		totalSam2Keep = totalSam2Keep + int(v[3])
		if tray2KeepCount >= userTrayInput: break
		tray2KeepCount+=1

	# print new results
	trayCount = 1
	print(Fore.RED + "\n{: <10}{: <8}{: ^8}{: ^11}{: <9}{: ^2}   {}".format('Tray', 'Tech', 'i7Name', 'i7Seq', 'reads', 'N', 'Species') + Fore.WHITE)
	for k,v in sorted(newDisplayDict.items()):
		print("{: <6}{: <8}{: ^8}{: ^11}{: <9}{: ^2}   {}".format(trayCount, v[0], k, v[1], v[2], v[3], v[4]))
		trayCount+=1
	print("...\nNot counting NTCs, there are{:6,} real samples in {:<2} real trays in {:>5}.\n\n".format(totalSam2Keep, tray2KeepCount , libOfInterest))
	
	# archive results
	archiveLibData = [libOfInterest, int(os.path.getctime(libName)), userTrayInput, totalSam2Keep]
	for k, v in sorted(archiveUniq7Dict.items(), reverse=True):
		tempL = [v[0], v[1], v[2], k, v[3], v[4]]
		archiveLibData = archiveLibData + tempL
	archiveDat = open('/home/efglserv/software/GTseq/techTracker/libStatArchive.csv', 'a')
	archiveDat.write(','.join(str(x) for x in archiveLibData))
	archiveDat.write('\n')

	
def Main():

	# get supporting files from techTracker directory
	with open('/home/efglserv/software/GTseq/techTracker/i7.csv') as f: i7Dat = f.read().splitlines()
	with open('/home/efglserv/software/GTseq/techTracker/i5.csv') as f: i5Dat = f.read().splitlines()
	
	techi7Dict, orderi5Dict, i5List = {},{},[]
	for x in i7Dat:
		rowDat = x.split(',')
		techi7Dict[rowDat[0]] = rowDat[1:]
	for x in i5Dat:	# note that H12 or 96 is not included
		rowDat = x.split(',')
		orderi5Dict[rowDat[1]] = rowDat[0] + rowDat[2]
		i5List.append(rowDat[0])
	
	# get library to analyze or display archived library
	with open('/home/efglserv/software/GTseq/techTracker/libStatArchive.csv') as f: archiveDatCheck = f.read().splitlines()	# check to see if library is already in archive	
	now = datetime.datetime.now()
	count = 0
	while count <3:
		libOfInterest = input('\nEnter GTseq library you like to collect data e.g. L0150: ')
		for i in archiveDatCheck:
			if libOfInterest in i: displayArchiveLib(i, now)	# display existing library data then exit
		try:
			libName = '/media/efglserv/radspace/Library/' + str(now.year) + '/' + libOfInterest + '.fq.gz'
			libFile = gzip.open(libName, mode="rt")
			break
		except FileNotFoundError:
			print(Fore.GREEN + '\nYour library may not be a GTseq library or it may not have been archived. Try again or enter \'Ctrl+C\' to exit' + Fore.WHITE)
			count+=1
	
	# get data from new library
	timeOfCreation= "create: %s" % time.ctime(os.path.getctime(libName))	# get when the libary was created
	timeOfCreation = timeOfCreation.split(' ')

	libDat = list(islice(libFile, 16000))	 # ignore first 4,000 reads
	libDat = list(islice(libFile, 4000000))
	libFile.close()
	titleL = libDat[::4]
	fastaL = libDat[1::4]

	comboB, all7s = [], []
	for i in titleL:
		if i[-14:-8] in techi7Dict:
			all7s.append(i[-14:-8])
			comboB.append(i[-14:-1])
	
	# get only i7s with 500 reads or more
	allUniq7s, trueUniq7s, uniq7Dict = list(set(all7s)), [], {}	
	for i in allUniq7s:
		occurrence = all7s.count(i)
		if occurrence >=500:
			uniq7Dict[occurrence] = i
			trueUniq7s.append(i)
		
	# tie species to unique i7
	speciesDict = {'CGCAATGAGCCAACCCCT': 'Chinook', 'CAGGTCACAGACACACAG': 'Steelhead', 'CATCACCGGGGGCAAACT': 'Sockeye', 'CGATGGAGGATGGGAAGG': 'Cutthroat', 'AGGACTTCCCCAGGAAGA': 'Brook Trout', 'AGGTTGTTGGGAAAGCTA': 'Burbot'}
	i7SpeciesDict ={}
	for ui7 in trueUniq7s:
		tempFWDseq18bp = []
		for i in range(300000):
			if titleL[i][-14:-8] == ui7:
				tempFWDseq18bp.append(fastaL[i][:18])
		tempFWDseq18bp = sorted(tempFWDseq18bp, key=Counter(tempFWDseq18bp).get, reverse=True)	# sort list by most occurrence
		correctSpecies = ''
		for i in tempFWDseq18bp:
			try: correctSpecies = speciesDict.get(i)
			except: pass
			if correctSpecies:
				i7SpeciesDict[ui7] = speciesDict.get(i)
				break

	# get number of samples per tray
	numProc = multiprocessing.cpu_count()	
	if len(trueUniq7s) < numProc: numProc = len(trueUniq7s)
	else: numProc = 56
	perP = {}
	for i in range(numProc): perP[i] = 's' + str(i)
	for i in range(numProc):
		perP[i] = Process(target=i5search, args=(trueUniq7s[i], comboB, i5List))
		perP[i].start()
	for i in range(numProc): perP[i].join()	# join processes when all is done

	numSamDict = {}	# collect data from temp file then erase temp file, part of counting sample from above
	with open('samInTray.csv') as f: samTrayDat = f.read().splitlines()
	for i in samTrayDat:
		rowDat = i.split(',')
		numSamDict[rowDat[0]] = int(rowDat[1])
	call('rm samInTray.csv', shell=True)
	
	# display and archive results
	Results(uniq7Dict, techi7Dict, i7SpeciesDict, numSamDict, libOfInterest, libName)
	
	
Main()
