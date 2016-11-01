#!/usr/bin/env python

import os, re, sys, string

fcbs='./fcbs'
rs='Rs'

class Run:
	def __init__(self,data,d,m,D2=False):
		self.data = data
		self.d = d
		self.D2 = D2
		self.m = m
		self.args = '-v 2 -b 1000'
		if self.D2: self.args = self.args + ' -D2'
		self.prog=fcbs

	def execline(self):
		return '%s -r ../../%s -d %s -m %i %s' \
			%(self.prog,self.data,self.d,self.m,self.args)

	def rundir(self):
		rname = self.data
		rundir = '%s.%s' %(rname,self.d)
		if(self.D2): rundir = rundir + '.D2'
		rundir = '%s.%i' %(rundir,self.m)
		if not os.path.exists(rundir): os.mkdir(rundir)
		return rundir

	def script(self):
		self.fname = 'sge.%s.csh' %self.rundir()
		f = open(self.fname,'w')
		f.write('''
#!/bin/bash
#$ -cwd
#$ -l h_rt=120:00:00
#$ -M justin.ashworth@uts.edu.au
#$ -m n
#$ -N fcbs
#$ -o log
#$ -e err
#$ -S /bin/bash\n
''')

		f.write('''
r=`head -n $SGE_TASK_ID %s | tail -n 1`

cd %s

# multiscale bootstrap levels (run in separate dirs)
mkdir $r
cd $r

%s -s $r >& log

''' %(rs,self.rundir(), self.execline()))

runs = [
#	Run('data.tab','euclidean',1),
#	Run('data.tab','euclidean',4),
	Run('data.tab','euclidean',4,True),
#	Run('data.tab','pearson1',1),
#	Run('data.tab','pearson1',4),
#	Run('data.tab','pearson1',4,True),
#	Run('data.tab','pearson2',1),
#	Run('data.tab','pearson2',4),
#	Run('data.tab','pearson2',4,True),
#	Run('data.tab','spearman',1),
#	Run('data.tab','spearman',4),
#	Run('data.tab','spearman',4,True),
]

for run in runs:
	run.script()
