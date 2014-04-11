#!/usr/bin/env python
import csv
import math
import re
import itertools

#file3 = open('samplefiles11.txt','r') # NonControlSamples.txt contains the names of control experiments
#apf = open('MTBChIPpeakoverlay6.txt','r')

#fout = open('MTB_ChIPbindingnetworkinfo11.txt','w')
#fgout = open('MTBChIPpeakoverlay11.txt','w')

file3 = open('samplefiles5.txt','r')
apf = open('MTBChIPpeakoverlay3.txt','r')

fout = open('MTB_ChIPbindingnetworkinfo5.txt','w')
fgout = open('MTBChIPpeakoverlay5.txt','w')


AllPeakFreq = apf.read()
AllPeakFreq = AllPeakFreq.split('\n')
AllPeakFreq = map(int,AllPeakFreq[:-1]) # map repeats the function in the first argument over every element of the array in the second argument

tb = file3.read()
tb = tb.split('\r')
tbatches = []
for item in tb:
	x = item.split('\t')
	tbatches.append(x)
file3.close()

#### Functions used for the main part of the code #####

# Reads the information from each .csv file
def readExcelChIP(samplename):
	ft1 = open('./%s/%s.MTb.ChIPpeaks.Excel.csv' %(samplename,samplename),'r')
	reader1 = csv.reader(ft1,delimiter=',',quotechar='"')
	chiparray = []
	for row in reader1:
		chiparray.append(row) # chiparray is a list containing all of the peak information for each experiment
	
	expname = dir[0].split('_')
	dTF = expname[0]
	ft1.close()
	return dTF, chiparray

# Extracts the relevant peak properties from each ChIP experiment	
def extractPeakParts(chiparray):
	Pval = []
	Target = []
	VPM = []
	Score = []
	Fstart = []
	Rstop = []
	Fcenter = []
	Ccenter = []
	Rcenter = []
	DNAseq = []

	for peak in chiparray[1:]:
		try:
			VPM.append(float(peak[2]))
		except:
			break
		Pval.append(float(peak[3]))
		Target.append(peak[1])
		Score.append(float(peak[4]))
		Ccenter.append(float(peak[6]))
		Fcenter.append(float(peak[14]))
		Rcenter.append(float(peak[15]))
		DNAseq.append(peak[50].split('\n')[-1])
			
		Fst = peak[32] # the start of the forward strand peak (start of the footprint)
		Fst = int(math.floor(float(Fst)))
		Fstart.append(Fst)	
		del Fst
			
		Rop = peak[39] # the end of the reverse strand peak (end of the footprint)
		Rop = int(math.ceil(float(Rop)))
		Rstop.append(Rop)
		del Rop
	return Pval,Target,VPM,Score,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq

# Outputs the peak properties to text file
def outputGoodPeak(TF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,AllPeakFreq,count1,count2):
	fout.write('%s\t%s\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%s\t%d\t%d\n' %(TF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,count1,count2))	
	prange = range(Fstart,Rstop+1) 
	#Tally this peak base range in the variable AllPeakFreq
	for base in prange:
		AllPeakFreq[base-1] = AllPeakFreq[base-1] + 1
	return AllPeakFreq

# For peaks that show up in only one experiment, outputs the peak if p-value < 0.05
def SinglePeakProcessing(Tonly,TF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,AllPeakFreq):
	for tar in Tonly:
		#ix = Target.index(tar)
		ix = [i for i,x in enumerate(Target) if x == tar] # extracts all peaks with the corresponding target gene
		for item in ix:
			if Pval[item] < 0.05:
				AllPeakFreq = outputGoodPeak(TF,Target[item],Pval[item],Score[item],VPM[item],Fstart[item],Rstop[item],Fcenter[item],Ccenter[item],Rcenter[item],DNAseq[item],AllPeakFreq,1,1)
		del ix

	return AllPeakFreq

# For peaks that show up in multiple experiments, outputs the peak with the lowest p-value if P < 0.05
def DuplicatePeakProcessing(commonTars,dir,dTF1,Target1,Pval1,Score1,VPM1,Fstart1,Rstop1,Fcenter1,Ccenter1,Rcenter1,DNAseq1,AllPeakFreq):
	
	TFdict = {}

	# iterate over all of the common target genes
	for ctar in commonTars:

		ix = [0]*len(dir)
		tardict = {}

		# iterate over each of the experiments
		for sample in range(len(dir)):

			# find all instances (peaks) of the target gene in each experiment
			ix[sample] = [i for i,x in enumerate(Target1[sample]) if x ==ctar]

			# iterate over each peak of the target gene
			for idx in ix[sample]:
				if Pval1[sample][idx] < 0.05:

					#compile a dictionary of the relevant parameters for each peak
					tardict['-'.join([str(sample),str(idx),str(Fstart1[sample][idx])])] = [ctar,sample,idx,Fstart1[sample][idx],Rstop1[sample][idx],Pval1[sample][idx],2,1]
		
		# compare each pair of peaks 			
		for pair in itertools.combinations(tardict.keys(),2):
			
			#extract the binding footprint of each pair of peaks under consideration
			try:
				bindrange0 = range(tardict[pair[0]][3],tardict[pair[0]][4]+1)
			except:
				bindrange0 = [0]
			try:
				bindrange1 = range(tardict[pair[1]][3],tardict[pair[1]][4]+1)
			except:
				bindrange1 = [-1]

			# calculate the overlap
			overlap = len(set(bindrange0).intersection(set(bindrange1)))
			
			#del bindrange0
			del bindrange1
			
			# if the overlap is greater than threshold, then delete the peak with the greater p-value
			if overlap > 0.5*len(bindrange0):
				if tardict[pair[0]][5] < tardict[pair[1]][5]:
					del tardict[pair[1]]
					tardict[pair[0]][7] = 2
				else:
					del tardict[pair[0]]
					tardict[pair[1]][7] = 2
			
			del bindrange0
		
		for key in tardict:
			AllPeakFreq = outputGoodPeak(dTF1[tardict[key][1]],Target1[tardict[key][1]][tardict[key][2]],Pval1[tardict[key][1]][tardict[key][2]],Score1[tardict[key][1]][tardict[key][2]],VPM1[tardict[key][1]][tardict[key][2]],Fstart1[tardict[key][1]][tardict[key][2]],Rstop1[tardict[key][1]][tardict[key][2]],Fcenter1[tardict[key][1]][tardict[key][2]],Ccenter1[tardict[key][1]][tardict[key][2]],Rcenter1[tardict[key][1]][tardict[key][2]],DNAseq1[tardict[key][1]][tardict[key][2]],AllPeakFreq,tardict[key][6],tardict[key][7])

	return AllPeakFreq

######################################

# iterate over each TF
for dir in tbatches: #dir is separate for each experiment
	tmp1 = []
	dTF1 = []
	Pval1 = []
	Target1 = []
	VPM1 = []
	Score1 = []
	Fstart1 = []
	Rstop1 = []
	Fcenter1 = []
	Ccenter1 = []
	Rcenter1 = []
	DNAseq1 = []
	
	#iterate over each experiment of each TF
	for fnum in range(len(dir)):
		
		##print(dir[fnum])
		tmp1.append([])
		dTF1.append([])
		Pval1.append([])
		Target1.append([])
		VPM1.append([])
		Score1.append([])
		Fstart1.append([])
		Rstop1.append([])
		Fcenter1.append([])
		Ccenter1.append([])
		Rcenter1.append([])
		DNAseq1.append([])
			
		# extract the relevant parameters for every peak
		dTF1[fnum],tmp1[fnum] = readExcelChIP(dir[fnum])

		Pval1[fnum],Target1[fnum],VPM1[fnum],Score1[fnum],Fstart1[fnum],Rstop1[fnum],Fcenter1[fnum],Ccenter1[fnum],Rcenter1[fnum],DNAseq1[fnum]=extractPeakParts(tmp1[fnum])
		
	# want to find the target genes that only show up in one experiment vs. the genes that show up in multiple experiments
	
	# bigIntersect contains the union of target genes in all but one of the experiments
	intersects = list(itertools.combinations(range(len(dir)),r=len(dir)-1)) #combinations of n choose (n-1)
	bigIntersect = [[] for item in range(len(dir))]
	for item,itemr in zip(range(len(dir)),reversed(range(len(dir)))):
		for element in intersects[item]:
			bigIntersect[itemr] = list(set(bigIntersect[itemr]) | set(Target1[element]))
		
	# uniqueTargets contains target genes found in only one of the experiments
	uniqueTargets = [[] for item in range(len(dir))]
	for item in range(len(dir)):		
		uniqueTargets[item] = list(set(Target1[item]).difference(bigIntersect[item]))
		AllPeakFreq = SinglePeakProcessing(uniqueTargets[item],dTF1[item],Target1[item],Pval1[item],Score1[item],VPM1[item],Fstart1[item],Rstop1[item],Fcenter1[item],Ccenter1[item],Rcenter1[item],DNAseq1[item],AllPeakFreq)
		
	del bigIntersect

	# commonTargets contains target genes found in more than one experiments
	commonTargets = []
	for item in range(len(dir)):
		commonTargets = list(set(commonTargets) | set(Target1[item]).difference(set(uniqueTargets[item])))
	AllPeakFreq = DuplicatePeakProcessing(commonTargets,dir,dTF1,Target1,Pval1,Score1,VPM1,Fstart1,Rstop1,Fcenter1,Ccenter1,Rcenter1,DNAseq1,AllPeakFreq)
	
	del commonTargets			
	del Fstart1
	del Rstop1
	del Target1
	del VPM1
	del Score1
	del Pval1
	del Fcenter1
	del Ccenter1
	del Rcenter1
	del DNAseq1
	del uniqueTargets	

for basefreq in AllPeakFreq:
	fgout.write('%d\n' %(basefreq))
fgout.close()