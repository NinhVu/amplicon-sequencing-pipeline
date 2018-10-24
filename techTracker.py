#!/usr/bin/python3

# techTracker.py, 20171226, n - presents GTseq library statistics for a given time period
# in terminal $ techTracker.py


import operator, time, sys, colorama
from datetime import datetime
from colorama import Fore

def checkDateFormat(date):
	try:
		date = date.split(',') 
		date = list(map(int, date))
		if date[0] > 0 and date[0] < 13 and date[1] >= 2017: return(date)
		else: return Error
	except: return Error


def dateOfInterest():
	print('\ntechTracker tracts number of samples genotyped for a given period of time')
	print('Enter a time period you like to see\n')
	
	attempt =0
	while attempt <=3:
		try:
			startTime = input('starting date e.g. 9,2017   ')
			startDate = checkDateFormat(startTime)
			break
		except: 
			if attempt ==2: print('last chance')
			if attempt <3: print('please enter correct date format and range')
			attempt+=1
			continue
	try: startDate
	except:
		print('You ran out of chances. Try rerunning program again.')
		sys.exit()

	attempt = 0 
	while attempt <=3:
		try:
			endTime = input('ending date e.g. 12,2017    ')
			endDate = checkDateFormat(endTime)
			break
		except: 
			if attempt ==2: print('last chance')
			if attempt <3: print('please enter correct date format and range')
			attempt+=1
			continue
	try: endDate
	except:
		print('You ran out of chances. Try rerunning program again.')
		sys.exit()

	if endDate[1] - startDate[1] >=0:
		return([startDate[0], startDate[1], endDate[0], endDate[1]])
	else:
		print('\nyou enter a nonlinear time range. Try running program again.\n')
		sys.exit()


def cumulativeLibStats(sTime, eTime, dateRange, archiveDat, monthDict):
	totalLibs, totalTray = 0, 0
	libList, totalTechSamTup, speciesTup = [], [], []
	for i in archiveDat:
		libDat = i.split(',')
		rightDat = []
		for z in libDat:
			try: rightDat.append(int(z))
			except: rightDat.append(z)
		if rightDat[1] >= int(sTime) and rightDat[1] < int(eTime):
			libList.append(rightDat[0])
			techName = rightDat[4:][::6]
			numSam = rightDat[4:][4::6]
			species = rightDat[4:][5::6]
			techSamTup = list(zip(techName, numSam))
			speciesCount = list(zip(species, numSam))
			totalTray = totalTray + rightDat[2]
			totalTechSamTup = totalTechSamTup + techSamTup	# create a COMPLETE list of techs-sampleS
			speciesTup = speciesTup + speciesCount

	# sum up total of samples for each tech
	libList = list(sorted(libList))
	allTechs = [x[0] for x in totalTechSamTup]	# get all tech names - with duplicates, and put in list
	uniqTechs = list(sorted(set(allTechs)))
	techSamDict= {}
	for i in uniqTechs:
		tempNsum = 0
		for n in totalTechSamTup:
			if i in n[0]: tempNsum = tempNsum + n[1]
		techSamDict[i] = [tempNsum, allTechs.count(i)]
	samArchive = sorted(techSamDict.items(), key=operator.itemgetter(1))	# create tuple in ascending order based on sample count
	
	# sum of total samples for each species and sort from least to most genotyped
	speciesSum = []
	species = [x[0] for x in speciesTup]
	species = list(set(species))
	for i in species:
		nCount=0
		for s in speciesTup:
			if s[0] == i: nCount = nCount + int(s[1])
		speciesSum.append(nCount)
	speciesCountTup = list(zip(speciesSum, species))
	speciesCountTup = list(sorted(speciesCountTup))

	# display results
	print(Fore.BLUE + '\n\nFrom the start of ' + monthDict.get(dateRange[0]) + ' ' + str(dateRange[1]) + ' to the end of ' + monthDict.get(dateRange[2]) + ' ' + str(dateRange[3]), end="")
	print(', there were ' + str(len(libList)) + ' GTseq libraries produced.')
	print(Fore.WHITE + '%s' % ', '.join(map(str, libList)))
	print(Fore.RED + "\n{:<8}{:>8}{:>11}".format('Tech', '#tray', '#sample'))
	totalAllTray, totalAllSam = 0, 0
	for k,v in samArchive:
		totalAllTray = totalAllTray + v[1]
		totalAllSam = totalAllSam + v[0]
		print(Fore.WHITE + "{:<8}{:>8}{:>11}".format(k, str(v[1]), str(v[0])))
	print('---------------------------')
	print("{:<8}{:>8}{:>11}\n\n".format('Total', str(totalAllTray), str(totalAllSam)))
	
	
	print(Fore.RED + "{:11}   {:>7}".format('species', 'N'))
	for i in speciesCountTup: print(Fore.WHITE + "{:11}   {:>7}".format(i[1], i[0]))
	print('\n')

	
def Main():
	monthDict = {1: 'January', 2: 'Febuary', 3: 'March', 4: 'April', 5: 'May', 6: 'June', 7: 'July', 8: 'August', 9: 'September', 10: 'October', 11: 'November', 12: 'December'}
	dateRange = dateOfInterest()
	formatTime = []
	if dateRange[0] < 10: formatTime.append(str(dateRange[1]) + '-0' +  str(dateRange[0]) + '-01')
	else: formatTime.append(str(dateRange[1]) + '-' +  str(dateRange[0]) + '-01')
	
	# roundup endTime
	if dateRange[2] == 12: formatTime.append(str((dateRange[3] +1)) + '-01-01')
	elif dateRange[2] < 9: formatTime.append(str(dateRange[3]) + '-0' +  str((dateRange[2] +1)) + '-01')
	else: formatTime.append(str(dateRange[3]) + '-' + str((dateRange[2] +1)) + '-01')

	startDT = datetime.strptime(formatTime[0], '%Y-%m-%d')
	sTime = time.mktime(startDT.timetuple())
	endDT = datetime.strptime(formatTime[1], '%Y-%m-%d')
	eTime = time.mktime(endDT.timetuple())

	# get archive data and display results
	with open('/home/efglserv/software/GTseq/techTracker/libStatArchive.csv') as f: archiveDat = f.read().splitlines()
	archiveDat = list(filter(None, archiveDat))
	cumulativeLibStats(sTime, eTime, dateRange, archiveDat, monthDict)

	
Main()
