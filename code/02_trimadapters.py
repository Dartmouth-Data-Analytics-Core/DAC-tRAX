#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
import re
import subprocess
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool

# -*- coding: de_DE.UTF-8 -*- 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
Function declarations below
"""
#----- Function to update the manifest file
def updatemanifest(indexfilename, runname, runfile):
    filelocs = dict()
    try:
        indexfile = open(indexfilename)
        for currline in indexfile:
            fields = currline.split()
            if len(fields) > 1:
                filelocs[fields[0]] = fields[1]
        indexfile.close()
    except IOError as e:
        pass
    filelocs = filelocs
    filelocs[runname] = runfile
    outfile = open(indexfilename, "w")
    for currname in filelocs.iterkeys():
           print >>outfile, currname+"\t"+filelocs[currname]

#----- Function to call an RScript
def runrscript(*script):
    print >>sys.stderr, "Rscript "+" ".join(script)

    retcode = subprocess.call("Rscript "+" ".join(script), shell=True,  stderr = subprocess.STDOUT)

    if retcode > 0:
        print >>sys.stderr, "R script "+script[0]+" failed"

        # Uncommenting this line will kill the script if the R code fails
        #sys.exit()
    return retcode

#----- Function to extract sample name and subprocess arguments from argnames    
def subprocesspool(argnames):
	# First argument is treated as sample name
    samplename = argnames[0]
    # Second argument is expected to be a list with subprocess arguments
    args = argnames[1]
    # Launch a subprocess using the provided arguments
    process = subprocess.Popen(*args[0], **args[1])
    # Wait for the subprocess to complete before moving on
    process.wait()
    # Return the sample name, the executed command (joined into a string), and the process object
    return samplename," ".join(args[0]),  process

"""
Function to take any number of position (*args) and keyword args (*kwargs)
and returns a tuple containing the positional arguments as a list and the keyword
argument as a dictionary.
"""
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])

"""
Function to process the output of SeqPrep command line tool. It extracts specific information
from the tool's standard error output (stderr) which contains details about the total
number of reads, the number of merged reads, and discarded reads. It returns a dictionary of 
these counts along with the error output for reference.
"""
def readseqprep(processoutput):
	# Initialize a dictionary to store read counts
    seqprepcounts = dict()

    # Retrieve the output (stdout and stderr) from the process execution
    output = processoutput.communicate()
    errinfo = output[1] # This contains stderr information

    # Check if the process failed (return code is non-zero)
    if processoutput.returncode != 0:
        print >>sys.stderr, "seqprep failed"
        print >>sys.stderr, errinfo

    # Parse each line of stderr output
    for line in errinfo.split("\n"):
    	# Match patterns to extract total reads, merged reads, and discarded reads
        totalmatch = rereadtotal.match(line)
        mergematch = rereadmerge.match(line)
        discardmatch = rereaddiscard.match(line)
        # Extract reads based on pattern matching and remove commas
        if totalmatch:
            totalreads = int(totalmatch.group(1).replace(",",""))
        elif mergematch: 
            merged = int(mergematch.group(1).replace(",",""))
        elif discardmatch:
            discard = int(discardmatch.group(1).replace(",",""))

    # Calculate the counts for merged, unmerged, and discarded reads
    seqprepcounts["merged"] = merged
    seqprepcounts["unmerged"] = totalreads - (merged + discard)
    seqprepcounts["discarded"] = discard

    # Return the dictionary of counts and the stderr output for reference
    return seqprepcounts, errinfo

"""
This function is similar to the above, but processess the output of the command line tool
cutadapt, which is used for trimming adapter sequences from reads. The stderr file is 
parsed to extract information about total reads, trimmed reads, discarded reads, and
written reads.
"""
def readcutadapt(processoutput):
	# Initialize a dictionary to store cutadapt counts
    cutadaptcounts = dict()

    # Retrieve the output (stdout and stderr) from the process execution
    output = processoutput.communicate()
    errinfo = output[1] # This contains stderr information

    # Check if the process failed (return code is non-zero)
    if processoutput.returncode != 0:
        print >>sys.stderr, "cutadapt failed"
        print >>sys.stderr, errinfo

    # Parse each line of stderr output
    for line in errinfo.split("\n"):

    	# Match patterns to extract total, trimmed, discarded, and written reads
        totalmatch = recutadapttotal.match(line)
        trimmatch = recutadapttrimmed.match(line)
        discardmatch = recutadaptshort.match(line)
        writtenmatch = recutadaptwritten.match(line)
        #trimmed = None

    	# Extract based on pattern matching and remove commas.
        if totalmatch:
            totalreads = int(totalmatch.group(1).replace(",",""))
        elif trimmatch:
            trimmed = int(trimmatch.group(1).replace(",",""))
        elif discardmatch:
            discard = int(discardmatch.group(1).replace(",",""))
        elif writtenmatch:
            written = int(writtenmatch.group(1).replace(",",""))

    # Calculate the counts for trimmed, untrimmed and discarded reads
    cutadaptcounts["trimmed"] = trimmed - discard
    cutadaptcounts["untrimmed"] = totalreads - (trimmed)
    cutadaptcounts["discarded"] = discard

    # Return the dictionary of counts and the stderr output for reference
    return cutadaptcounts, errinfo

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define the directory of the script
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"

"""
Below, we will set up the argument parsers which will hold arguments passed to the script
via the command-line.
"""
# Initialize an argument parser for the command-line input
parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')

# Add required argument for run name
parser.add_argument('--runname',required=True,
                   help='run name to be used')

# Add required argument for the input file (run file)
parser.add_argument('--runfile',required=True,
                   help='run file')

# Add an optional argument for the first adapter sequence
# (default adapter for trimming)
parser.add_argument('--firadapter',default ='AGATCGGAAGAGCACACGTC' ,
                   help='subset file')

# Add an optional argument for the second adapter sequence
# (default adapter for trimming)
parser.add_argument('--secadapter',default ='GATCGTCGGACTGTAGAACTC' ,
                   help='subset file')

# Add an optional argument for minimum sequence length (after trimming)
parser.add_argument('--minlength',default ='15' ,
                   help='minimum length of sequence')

# Add a flag (boolean option) for single-end mode
# (default is paired-end) (FALSE)
parser.add_argument('--singleend', action="store_true", default=False,
                   help='single-end mode (uses cutadapt)')

# Add an optional argument to specify the length of the UMI (Unique Molecular Identifier)
# (default = 0)
parser.add_argument('--umilength',default ='0',
                   help='length of UMI (uses umi_tools to extract if present)')

# Add a flag to indicate if theUMI is at the 3' end of the sequence
# (default = FALSE)
parser.add_argument('--umithreeprime', action="store_true", default=False,
                   help='umi is at the three prime end')

# Add an optional argument to specify the number of processor cores to use
parser.add_argument('--cores',
                   help='number of processors to use')

# Parse the command-line arguments
args = parser.parse_args()

# Assign the parsed arguments to variables
runname = args.runname
seqprepfile = args.runfile
firadapter = args.firadapter
secadapter =   args.secadapter
minlength = args.minlength
singleendmode = args.singleend
threeprimeumi = args.umithreeprime
umilength = int(args.umilength)
cores = cpu_count()

# If the cores argument is passed, override the default cpu_count
if args.cores is not None:
    cores = args.cores

#firadapter = 'AGATCGGAAGAGCACACGTC'
#secadapter = 'GATCGTCGGACTGTAGAACTC'

"""
 Use the python re module to compile regular expressions (regex) that match specific patterns 
 in the output of tools such as SeqPrep and CutAdapt. These compiled regex patterns are used 
 later to search for, extract, and process numerical data (read counts) from the output of
 these tools.
 """

 # For SeqPrep
rereadtotal = re.compile(r'Pairs Processed:\s+(\d+)')
rereadmerge = re.compile(r'Pairs Merged:\s+(\d+)' )
rereadapter = re.compile(r'Pairs With Adapters:\s+(\d+)')
rereaddiscard = re.compile(r'Pairs Discarded:\s+(\d+)' )

# For CutAdapt
recutadapttotal = re.compile(r'Total reads processed:\s+([\d\,]+)')
recutadapttrimmed = re.compile(r'Reads with adapters:\s+([\d\,]+)' )
recutadaptshort = re.compile(r'Reads that were too short:\s+([\d\,]+)')
recutadaptwritten = re.compile(r'Reads written (passing filters):\s+([\d\,]+)' )

# Initialize dictionaries and lists
samplefiles = dict() # Dict to store sample names and their associated file paths
sampleorder = list() # List to maintain the order of sample names
cutadaptorder = list() # List to track sample names when using single-end mode

# Open the input file 'seqprepfile' and process each line
for currline in open(seqprepfile):

    # Split the current line by whitespace to extract fields
    fields = currline.split()

    # Check for paired-end mode (when not in SE Mode and there are at least 3 fields)
    if not singleendmode and len(fields) > 2:

    	# Store the sample name and its paired file path as a tuple
        samplefiles[fields[0]] = tuple([fields[1],fields[2]])
        # Append the same name to maintain its order
        sampleorder.append(fields[0])

    # Check that you are in SE Mode and there are at least 2 fields
    elif singleendmode and len(fields) > 1:

    	# Store the sample name and its paired file path as a tuple
        samplefiles[fields[0]] = tuple([fields[1]])
        cutadaptorder.append(fields[0])
        sampleorder.append(fields[0])

# Check if both sampleorder and cutadapt order are empty (i.e., no samples processed)
if len(sampleorder) < 1 and len(cutadaptorder) < 1:
    print >>sys.stderr, "Failed to read "+seqprepfile
    print >>sys.stderr, "Perhaps you failed to specify --singleend mode?"
    sys.exit(1)

# Initialize variables to store read counts and statistics
totalreads = None # Total number of reads processed
merged = None # Total number of merged reads (for Paired-End)
discard = None # Total number of discarded reads
samplenum = 1 # Counter to track the samples being processed
seqprepcounts = dict() # Dict to store the counts for each sample processed (merged, discarded)
cutadaptcounts = dict() # Dict to store the counts (trimmed, untrimmed)
allsamples = set() # Set to store the unique sample names that have been processed

# Initialize default parameters and list to hold the output file names
minsize = 15
outputfiles = list()

"""
Check the existence of files for each sample listed in samplefiles. This code
could be used for error logging as well.
"""
# Initialize an empty string to store preparation output
prepout = ""

# Iterate through each sample in the samplefiles directory
for currsample in samplefiles.iterkeys():

	# For each sample, iterate through its associated files
    for currfile in samplefiles[currsample]:

    	# Check if the current file does not exist
        if not os.path.isfile(currfile):
            print >>sys.stderr, currfile +" does not exist"
            sys.exit(1)

# Initialize empty dictionaries to store SeqPrep and CutAdapt runs for each sample
seqprepruns = dict() # Used to track the processess or outcomes of SeqPrep per sample
cutadaptruns = dict() # Used to track the processess or outcomes of CutAdapt per sample


#print >>sys.stderr, cores

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
Below is the code that executes the trimming of the sample
"""

# Create a pool of worker processes for parallel execution
trimpool = Pool(processes=int(cores))

# Initialize UMI args
umiargs = ""

# Check if 3' UMI extraction is requested
if threeprimeumi:

	# If TRUE, add the flag for 3' UMI to the arguments
    umiargs += " --3prime "

# Iterate through the samples in the specified order
for currsample in sampleorder:

    # Check if the mode is not single-end (paired-end)
    if not singleendmode:

    	# Initialize the command for sequence preparation
        seqprepcommmand = ""

        # If a UMI length is specified, prepare the SeqPrep command
        if umilength  > 0:
            seqprepcommmand = 'SeqPrep -L '+str(minsize)+ ' -A '+firadapter+' -B '+secadapter +' -f '+samplefiles[currsample][0]+'  -r '+samplefiles[currsample][1]+' -1 '+currsample+'_left.fastq.gz     -2 '+currsample+'_right.fastq.gz   -s '+currsample+'_m.fastq.gz; '

            # Command to extract UMI from merged fastq file
            seqprepcommmand += "umi_tools extract --stdin="+currsample+"_m.fastq.gz "+umiargs+"--bc-pattern="+("N"*int(umilength))+" --stdout="+currsample+'_merge.fastq.gz 1>&2'
        else:

            # Prepare a seqprep command without UMI extraction
            seqprepcommmand = 'SeqPrep -L '+str(minsize)+ ' -A '+firadapter+' -B '+secadapter +' -f '+samplefiles[currsample][0]+'  -r '+samplefiles[currsample][1]+' -1 '+currsample+'_left.fastq.gz     -2 '+currsample+'_right.fastq.gz   -s '+currsample+'_merge.fastq.gz'

		# Append the merged output file to the list
        outputfiles.append(currsample+'_merge.fastq.gz')

        # Initialize the sequence preparation run command
        seqprepruns[currsample] = None

        # Execute the command and capture outputs and errors
        seqprepruns[currsample] = compressargs(seqprepcommmand, shell = True, stderr = subprocess.PIPE)

	# Run in single-end mode if not paired-end
    else:

    	# Initialize a cutadapt command
        cutadaptcommand = ""

        # If a UMI length is specified, prepare the cutadapt command and command to extract UMI
        if umilength > 0:
            cutadaptcommand = 'cutadapt -m '+str(minsize)+ ' --adapter='+firadapter+' '+samplefiles[currsample][0]  +' | gzip -c >'+ currsample+'_t.fastq.gz;'
            cutadaptcommand += "umi_tools extract --stdin="+currsample+"_t.fastq.gz "+umiargs+"--bc-pattern="+("N"*int(umilength))+" --stdout="+currsample+'_trimmed.fastq.gz 1>&2'

        # Prepare cutadapt command without UMI extraction
        else:
            cutadaptcommand = 'cutadapt -m '+str(minsize)+ ' --adapter='+firadapter+' '+samplefiles[currsample][0]  +' | gzip -c >'+ currsample+'_trimmed.fastq.gz'

	# Append the output dile to the list
        outputfiles.append(currsample+'_trimmed.fastq.gz')

        # Initialize the cutadapt run command
        cutadaptruns[currsample] = None

        # Execute the cutadapt command and capture the outputs/errors
        cutadaptruns[currsample] = compressargs(cutadaptcommand, shell = True, stderr = subprocess.PIPE)
# check if the mode is paired-end
if not singleendmode:

    # Use the process pool to apply the seqprep commands in parallel
    results = trimpool.imap_unordered(subprocesspool, list(tuple([currsample, seqprepruns[currsample]]) for currsample in sampleorder))

    # Iterate over the results returned from the process pool
    for samplename, command, spoutput in results:
        print >>sys.stderr, samplename +" merged"
        #print >>sys.stderr, spoutput

        # Read the output from the SeqPrep command and capture the count of reads and any error information
        seqprepcounts[samplename], errinfo = readseqprep(spoutput)

        # Append the sample name, command used, and any error information to the prepout string
        prepout += samplename +"\n"
        prepout += command+"\n"
        prepout += errinfo+"\n"

# For single-end mode, do the same as above
else:
    results = trimpool.imap_unordered(subprocesspool, list(tuple([currsample, cutadaptruns[currsample]]) for currsample in sampleorder))
    for samplename,  command, caoutput in results:

        print >>sys.stderr, samplename +" trimmed"
        cutadaptcounts[samplename], errinfo = readcutadapt(caoutput)
        prepout += samplename+"\n"
        prepout += command+"\n"
        prepout += errinfo+"\n"

#print >>sys.stderr,  cutadaptcounts.keys()
#sys.exit()
''' 
for currsample in sampleorder:
    if not singleendmode:
        output = seqprepruns[currsample].communicate()
        errinfo = output[1]
        if seqprepruns[currsample].returncode != 0:
            print >>sys.stderr, "seqprep failed"
            print >>sys.stderr, errinfo
        
        prepout += errinfo
        for line in errinfo.split("\n"):
        
            
            totalmatch = rereadtotal.match(line)
            mergematch = rereadmerge.match(line)
            discardmatch = rereaddiscard.match(line)
            if totalmatch:
                totalreads = int(totalmatch.group(1).replace(",",""))
            elif mergematch: 
                merged = int(mergematch.group(1).replace(",",""))
            elif discardmatch:
                discard = int(discardmatch.group(1).replace(",",""))
        
        seqprepcounts[currsample]["merged"] = merged
        seqprepcounts[currsample]["unmerged"] = totalreads - (merged + discard)
        seqprepcounts[currsample]["discarded"] = discard
    else:
        output = cutadaptruns[currsample].communicate()
        errinfo = output[1]
        if cutadaptruns[currsample].returncode != 0:
            print >>sys.stderr, "seqprep failed"
            print >>sys.stderr, errinfo
        print >>sys.stderr, errinfo
        prepout += errinfo
        for line in errinfo.split("\n"):
            totalmatch = recutadapttotal.match(line)
            trimmatch = recutadapttrimmed.match(line) 
            discardmatch = recutadaptshort.match(line) 
            writtenmatch = recutadaptwritten.match(line) 
            #trimmed = None

            if totalmatch:
                totalreads = int(totalmatch.group(1).replace(",",""))
            elif trimmatch: 
                trimmed = int(trimmatch.group(1).replace(",",""))
            elif discardmatch:
                discard = int(discardmatch.group(1).replace(",",""))
            elif writtenmatch:
                written = int(writtenmatch.group(1).replace(",",""))
               
        cutadaptcounts[currsample]["trimmed"] = trimmed
        cutadaptcounts[currsample]["untrimmed"] = totalreads - (trimmed + discard)
        cutadaptcounts[currsample]["discarded"] = discard
'''

#print >>sys.stderr, cutadaptcounts

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Check if the mode is paired-end
if not singleendmode:

	# Open files for output: sample file, log file, and manifest file for replicates
    samplefile = open(runname+"_sp.txt", "w")
    logfile = open(runname+"_log.txt", "w")
    replicatefile = open(runname+"_manifest.txt", "w")

    # Write sample names and corresponding output files to the manifest file
    for i, curr in enumerate(sampleorder):
        print >>replicatefile, curr+"\t"+outputfiles[i]
    print >>samplefile,"\t".join(sampleorder)

    # Write counts of merged, unmerged, and discarded reads for each sample
    for currtype in ["merged","unmerged","discarded"]:
        print >>samplefile,currtype+"\t"+"\t".join(str(seqprepcounts[currsample][currtype]) for currsample in sampleorder)

    # Close the sample file after writing
    samplefile.close()

    # Write various information to the log file
    print >>logfile ,"samplefile: "+seqprepfile
    print >>logfile ,"first adapter: "+ firadapter
    print >>logfile ,"second adapter: "+ secadapter
    print >>logfile ,"output files: "+ ",".join(outputfiles)
    print >>logfile ," ".join(sys.argv)
    print >>logfile ,"************************"
    print >>logfile ,prepout

	# Close the log file after writing
    logfile.close()

    # Run an Rscript for feature types using the sample file
    runrscript(scriptdir+"/02b_featuretypesreal.R",runname+"_sp.txt",runname+"_sp.pdf")

# Do the same for single-end mode
else:
    samplefile = open(runname+"_ca.txt", "w")
    logfile = open(runname+"_log.txt", "w")
    replicatefile = open(runname+"_manifest.txt", "w")
    for i, curr in enumerate(cutadaptorder):
        print >>replicatefile, curr+"\t"+outputfiles[i]
    print >>samplefile,"\t".join(cutadaptorder)
    for currtype in ["trimmed","untrimmed","discarded"]:
        print >>samplefile,currtype+"\t"+"\t".join(str(cutadaptcounts[currsample][currtype]) for currsample in cutadaptorder)

    samplefile.close()
    print >>logfile ,"samplefile: "+seqprepfile
    print >>logfile ,"adapter: "+ firadapter
    print >>logfile ,"output files: "+ ",".join(outputfiles)
    print >>logfile ," ".join(sys.argv)
    print >>logfile ,"************************"
    print >>logfile ,prepout

    logfile.close()
    runrscript(scriptdir+"/02b_featuretypesreal.R",runname+"_ca.txt",runname+"_ca.pdf")

# Update the manifest
updatemanifest("trimindex.txt",runname,seqprepfile)

