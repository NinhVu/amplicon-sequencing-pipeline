#!/usr/bin/python3
# callinGeno_v1.0.py, 170216, ninh - is a wrapper program running GTseq_Genotyper_v3.pl in multiprocess mode and scoring sex and/or hybrid markers
# callinGeno_v2.0.py, 170511 - simplify and make robust panel selection and compatibility with other labs
# callinGeno_v2.1.py, 170809 - move typing function to salmonBioInfo and reorganize probeSeq files
# callinGeno_v2.2.py, 180212 - streamline multiprocessing
# in terminal $ callinGeno_v2.1.py<enter>

import sys, os, glob, math, time, multiprocessing, colorama
from colorama import Fore, Back, Style
from multiprocessing import Process
from subprocess import *
from salmonBioInfo import typing


def Main():
	start = time.time()
	# determine panel to use for scoring
	probeSeqDirPath = '/home/efglserv/software/GTseq/probeSeq/'
	probeSeqDict = {1: ['BCT-p149', 'Bear River Cutthroat', 'BCT_p149_ProbeSeq.csv'], 2: ['Llo-p97', 'Burbot', 'Llo_ProbeSeqs_v2.2.csv'], 3: ['Sfo-p242', 'Brook Trout', 'Sfo_GTseq242_ProbeSeqs_v2.1.csv'], 4: ['YCT-p156', 'Yellowstone Cutthroat', 'YCT_p156_ProbeSeq_v1.2.csv'], 5: ['One-p382', 'Sockeye', 'One_GTseq382t75_ProbeSeqs.csv'], 6: ['Omy-p192', 'Steelhead', 'Omy192_ProbeSeqs_FixedProbes_v2.csv'], 7: ['Omy-p269', 'Steelhead', 'Omy_GTseq269_ProbeSeqs_FixedProbes_v2.csv'], 8: ['Omy-p379', 'Steelhead', 'Omy_GTseq379_ProbeSeqs_75b_IDFG.csv'], 9: ['Ots-p298_75bp', 'Chinook', 'Ots_GTseq298_ProbeSeqs_v2_fixed_75b.csv']}
	print('\nPanels available for scoring GTseq data.\n')
	for p, shortName in sorted(probeSeqDict.items()):
		if p == 5: print(Fore.RED + str(p) + ' ' + shortName[0] + '  ' + shortName[1])
		elif p == 6 or p == 7 or p == 8: print(Fore.BLUE + str(p) + ' ' + shortName[0] + '  ' + shortName[1])
		elif p == 9: print(Fore.GREEN + str(p) + ' ' + shortName[0] + '  ' + shortName[1])
		elif p == 2: print(Fore.WHITE + str(p) + ' ' + shortName[0] + '   ' + shortName[1])
		else: print(Fore.WHITE + str(p) + ' ' + shortName[0] + '  ' + shortName[1])
	print(Fore.WHITE + '')
	pSeq = input('Choose a number corresponding to the panel you want to score: ')
	try:
		probeSeqFile = probeSeqDict.get(int(pSeq))
		probeFile = probeSeqDirPath + probeSeqFile[2]
	except:
		print('\nYou must enter a number. Try running callinG again.\n\n')
		sys.exit()
	
	# divide total number of fastq files into processes
	fqFileL, fqPerProc = [], []
	for file in sorted(glob.glob('*.fastq')): fqFileL.append(file)
	# multiprocess calling genotype of subsets of fastqs - optimize for 26 threads
	if len(fqFileL) <= 26: samPerProc = 1
	else: samPerProc = math.ceil(len(fqFileL)/26)
	fqPerProc = [ fqFileL[i:i + samPerProc] for i in range(0, len(fqFileL), samPerProc)]	# divide fastq list into chunks for each processor

	perP = {}
	for i in range(len(fqPerProc)): perP[i] = str(i + 1)	# create dynamic variable names for each Process
	print('There are {0} processes available. Each is scoring ~{1:,} samples of {2:,} total samples.'.format(len(perP), len(fqPerProc[0]), len(fqFileL)))
	for i in range(len(fqPerProc)):
		perP[i] = Process(target=typing, args=(perP[i], fqPerProc[i], probeFile))
		perP[i].start()
	for i in range(len(perP)): perP[i].join()	# join processes when all is done
	
	
	# sex or hybrid scoring
	if 'Omy' in probeFile or 'Ots' in probeFile: print('\r\nSex typing. Almost done!\n')
	if 'Omy' in probeFile: call('OmySEX_test_v3_IDFG.pl', shell=True)
	if 'Ots' in probeFile: call('OtsSEX_test_v2.pl', shell=True)
	if 'Sfo' in probeFile: call('moreSFOsex_v1.py', shell=True)
	if int(pSeq) == 1 or int(pSeq) == 4:
		print('\nNow scoring hybrid and sex marker. Almost done!\n')
		call('oddMark_v1.1.py', shell=True)

	hours, rem = divmod(time.time()-start, 3600)
	minutes, seconds = divmod(rem, 60)
	print("\nTotal run time {:0>2}:{:0>2}:{:05.2f}\n\n".format(int(hours),int(minutes),seconds))

Main()
