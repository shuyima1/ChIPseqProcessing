#!/usr/bin/env python
import csv
import math
import re

file1 = open('usamples.txt','r') # NonControlSamples.txt contains the names of control experiments
sbatches = file1.read()
sbatches = sbatches.split('\r')

file1.close()

#AllPeakFreq = [0]*4411532

fout = open('MTBrawChIPbindingnetwork.txt','w')
#fgout = open('MTBChIPpeakoverlaySingle.txt','w')

for dir in sbatches:
	tmp = []
	ft = open('./%s/%s.MTb.ChIPpeaks.Excel.csv' %(dir,dir),'r')
	reader = csv.reader(ft,delimiter=',',quotechar='"')
	
	x = dir.split('_')
	sTF = x[0]
	
	for row in reader:
		tmp.append(row) # tmp is a list containing all of the peak information for each experiment
	# tmp2 = tmp[1:]	# want to get rid of the first row containing row headers
	ft.close()
	
	for peak in tmp[1:]:
		try:
			Pval = float(peak[3])
		except:
			break
		try:
			VPM = float(peak[2])
		except:
			break

		Target = peak[1]
		Score = float(peak[4])
		
		Fstart = peak[32]
		#if re.search('NA',Fstart):
		#	continue
		#else:
		Fstart = int(math.floor(float(Fstart)))
		#Fstart = int(math.floor(float(peak([32])))
		
		Rstop = peak[39]
		## Convert Rstop to an integer if it's not "NA"
		#if re.search('NA',Rstop):
		#	continue
		#else:
		Rstop = int(math.ceil(float(Rstop)))
		#Rstop = int(math.ceil(float(peak[39])))
		
		# for each peak in each experiment, find the span of the peak (from Fstart to Rstop)
		# Convert Fstart to an integer if it's not "NA"
		#Find all of the bases spanned by the peak
		prange = range(Fstart,Rstop+1) 
		
		#Tally this peak base range in the variable AllPeakFreq
		#for base in prange:
		#	AllPeakFreq[base-1] = AllPeakFreq[base-1] + 1
		
		fout.write('%s\t%s\t%f\t%f\t%d\t%d\t%d\n' %(sTF,Target,Pval,Score,VPM,Fstart,Rstop))
		
		#Clear the variables used
		del Fstart
		del Rstop
		del Target
		del VPM
		del Score
		del prange
		
		del Pval
	
	del reader
	del x
	del tmp
	del sTF
	

#for basefreq in AllPeakFreq:
#	fgout.write('%d\n' %(basefreq))
#fgout.close()