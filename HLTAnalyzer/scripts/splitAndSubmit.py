#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import time

usage = "usage: python splitAndSubmit.py -t HLTAnalyzer_template_cfg.py -p ./ -q 1nh -i fileList.txt -o myTree --split 10"

parser = optparse.OptionParser(usage)
parser.add_option("-t", "--template", dest="template",
    help="name of the template config to be used",
    )

parser.add_option("-p", "--templatePath", dest="templatePath",
    help="path where to find the template",
    )

parser.add_option('-q', '--queue',       action='store',     dest='queue',       
    help='run in batch in queue specified as option (default -q cmslong)', 
    default='cmsan',
    metavar="QUEUE")

parser.add_option("-i", "--input", dest="input",
    help="path and name of the fileList",
    )

parser.add_option("-o", "--output", dest="output",
    help="the root file outName",
    metavar="OUTDIR")

parser.add_option("--split", dest="filesperjob", type=int,
    help="files to analyze per job ",
    default=10)

parser.add_option('-I', '--interactive',      
    action='store_true',
    dest='interactive',      
    help='run the jobs interactively, 2 jobs at a time',
    default=False)

(opt, args) = parser.parse_args()
################################################


###
pwd = os.environ['PWD']
current_time = datetime.datetime.now()
simpletimeMarker = "_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
timeMarker = "submit_%04d%02d%02d_%02d%02d%02d" % (current_time.year,current_time.month,current_time.day,current_time.hour,current_time.minute,current_time.second) 
workingDir = pwd+"/batch/"+timeMarker

os.system("mkdir -p "+workingDir)
os.system("cp "+opt.template+" "+workingDir)

template = workingDir+"/"+opt.template

inputlist = []
njobs_list = []

num_lines = sum(1 for line in open(opt.input, "r"))

ins = open(opt.input, "r")
##loop over lists (one for datasets) to create splitted lists
count = 0
jobCount = 0
for line in  ins:
    count = count+1
    line = "root://eoscms.cern.ch//eos/cms"+line.rstrip('\n')
    inputlist.append(line)
    if count%opt.filesperjob == 0 or count==num_lines:
        jobCount = jobCount+1
        os.system("mkdir "+workingDir+"/"+str(jobCount))
        
        with open(template) as fi:
            contents = fi.read()
            replaced_contents = contents.replace('INPUTLIST', str(inputlist))
            replaced_contents = replaced_contents.replace('OUTPUTFILE', "\""+opt.output+"_"+str(jobCount)+".root\"")
        with open(workingDir+"/"+str(jobCount)+"/config.py", "w") as fo:
            fo.write(replaced_contents)

        inputlist = []
        
        os.system("echo cd "+pwd+" > launch.sh")
        os.system("echo 'eval `scramv1 runtime -sh`\n' >> launch.sh")
        os.system("echo cd - >> launch.sh")
        os.system("echo cmsRun "+workingDir+"/"+str(jobCount)+"/config.py >> launch.sh")
        os.system("echo mv "+opt.output+"_"+str(jobCount)+".root "+workingDir+" >> launch.sh")
        os.system("chmod 755 launch.sh")
        os.system("mv launch.sh "+workingDir+"/"+str(jobCount))
        njobs_list.append("bsub -q"+opt.queue+" -o "+workingDir+"/"+str(jobCount)+"/log.out -e "+workingDir+"/"+str(jobCount)+"/log.err "+workingDir+"/"+str(jobCount)+"/launch.sh")
        
for job in njobs_list:
    print job






