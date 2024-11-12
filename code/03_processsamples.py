#!/usr/bin/env python

#----- IMPORT LIBRARIES
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import subprocess
import os
import time
import mapreads
import countreads
import getcoverage
import getends
import countreadtypes
import maketrackhub
import traxqc
from distutils.version import LooseVersion, StrictVersion
from multiprocessing import Pool, cpu_count

#expname is experiment name
#dbname is database name
#samplefile is sample file
#$4 is bed feature for other sRNAs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- ARGUMENTS
# Initialize an argument parser to pass command line arguments to the code
parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')

#----- Required arguments
# Argument for experiment name
parser.add_argument('--experimentname',required=True,
                   help='experiment name to be used')
# Argument for database name (formatted as /path/to/tRNAdb/<db_prefix>)
parser.add_argument('--databasename',required=True,
                   help='name of the tRNA database')
# Path to sample file (sample name, trimmed fastq file name, no column names)
parser.add_argument('--samplefile',required=True,
                   help='sample file')

# Optional arguments
parser.add_argument('--ensemblgtf',
                   help='The ensembl gene list for that species')
parser.add_argument('--exppairs',
                   help='List of sample pairs to compare')
parser.add_argument('--bedfile', nargs='*',
                   help='Additional bed files for feature list')
parser.add_argument('--lazyremap', action="store_true", default=False,
                   help='Skip mapping reads if bam files exit')
parser.add_argument('--nofrag', action="store_true", default=False,
                   help='Omit fragment determination (Used for TGIRT mapping)')
parser.add_argument('--olddeseq', action="store_true", default=False,
                   help='Use old DESeq1 for analysis')
parser.add_argument('--nosizefactors', action="store_true", default=False,
                   help='Don\'t use Deseq size factors in plotting')
parser.add_argument('--maxmismatches',
                   help='Maximum allowed mismatches')
parser.add_argument('--mincoverage',
                   help='Minimum read count for coverage plots')
parser.add_argument('--minnontrnasize',type=int,default=20,
                   help='Minimum read length for non-tRNAs')
parser.add_argument('--paironly', action="store_true", default=False,
                   help='Generate only pair files (for adding a pair file after initial processing)')
parser.add_argument('--makehub', action="store_true", default=False,
                   help='make a track hub')
parser.add_argument('--hubonly', action="store_true", default=False,
                   help='Only make the track hub')
parser.add_argument('--maponly', action="store_true", default=False,
                   help='Only do the mapping step')
parser.add_argument('--dumpother', action="store_true", default=False,
                   help='Dump "other" features when counting gene types')
parser.add_argument('--local', action="store_true", default=False,
                   help='use local bam mapping')
parser.add_argument('--cores',
                   help='number of cores to use')
parser.add_argument('--bamdir',
                   help='directory for placing bam files (default current working directory)')

# Additional argument options that were commented out by the code author
"""
parser.add_argument('--maketdr', action="store_true", default=False,
                   help='create tdrs')
parser.add_argument('--makeall', action="store_true", default=False,
                   help='make both track hub and tdrs')
parser.add_argument('--splittypecounts', action="store_true", default=False,
                   help='Split type counts into tRNA types')
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- HANDLING R SCRIPTS
# Initialize an Rlog file and open it as writeable
rlogname = "Rlog.txt"
rlogfile = open(rlogname, "w")

# Define a function to run R scripts
def runrscript(*script):
    print >>sys.stderr, "Rscript "+" ".join(script)
    print >>rlogfile, "*******************************************************" 
    print >>rlogfile, "Rscript "+" ".join(script)
    rlogfile.flush()
    retcode = subprocess.call("Rscript "+" ".join(script), shell=True, stdout = rlogfile, stderr = subprocess.STDOUT)

    if retcode > 0:
        print >>rlogfile, script[0]+" failed"
        print >>sys.stderr, "R script "+script[0]+" failed"
        print >>sys.stderr, "Check "+rlogname+" for details"
        
        #sys.exit()
    return retcode
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- DEFINING CLASSES FOR OOP
# Define a class for trnadatabase
"""
The trnadb class manages various files and operations realted to tRNA geonmic data
"""
class trnadatabase:
    # Intitialize the class with a database name (dbname), which is used to construct various file paths
    def __init__(self, dbname):
        self.dbname = dbname
        self.trnatable = dbname+"-trnatable.txt" # Files storing a table of tRNA information
        self.bowtiedb = dbname+"-tRNAgenome" # Bowtie2 database for aligning tRNA sequences to the genome
        self.locifile = dbname+"-trnaloci.bed" # BED files containing tRNA loci positions in the genome
        self.maturetrnas=dbname+"-maturetRNAs.bed" # BED file for mature tRNA sequences
        self.trnaalign = dbname+"-trnaalign.stk" # Stockholm alignment file for tRNA sequences
        self.locialign = dbname+"-trnaloci.stk" # Stockholm alignment file for tRNA loci
        self.trnanums = dbname+"-alignnum.txt" # File containing the number of aligned tRNA sequences
        self.locinums = dbname+"-locusnum.txt" # File containing the number of aligned tRNA loci
        self.trnafasta = dbname+"-maturetRNAs.fa" # FASTA file containing mature tRNA sequences
        self.modomics = dbname+"-modomics.txt" # File realted to post-transcriptional tRNA modifications (modomics)
        self.otherseqs = dbname+"-otherseqs.txt" # File containing other sequecnes such as non-tRNA
        self.dbinfo = dbname+"-dbinfo.txt" # Metadata or information file about the database
    
    # Method to run a Bowtiw2 alignment job on a test FASTQ file
    def test(self):
        # Test run as a subprocess
        bowtie2job = subprocess.Popen(
            ["bowtie2","-x", self.bowtiedb, "-U", scriptdir+"test.fq"],
            stdout = subprocess.PIPE,stderr = subprocess.STDOUT
            )
        
        # Retrieve and process the output of Bowtie2
        rstatsresults = rstatsjob.communicate()[0]
        
        # Check if Bowtie 2 returned a non-zero exit code
        if bowtie2job.returncode  != 0:
                print >>sys.stderr, "bowtie2 failed to run"
    
    # Method to determine the organism type (e.g, eukaryote) from the database info file
    def getorgtype(self):
        
        # Default organism type is eukaryotic ("euk")
        orgtype = "euk"
        
        # Open the database information file and look for the "orgmode" field
        for currline in open(self.dbinfo):
            fields = currline.split()
            if fields[0] == "orgmode":
                orgtype = fields[1]
        return orgtype

# Define a class for expdatabase  
"""
The expdatabase class manages various files and results related to a tRNA expression experiment
Each file is organized under a directory named after the experiment (expname)
"""         
class expdatabase:
    def __init__(self, expname):
        
        # Initialize the class with an experment name which is used to construt various file paths
        self.expname = expname
        
        # BED file containing all genomic features identified in the experiment
        self.allfeats = expname+"/"+expname+"-allfeats.bed" 
        
        # Mapping statistics and visualization files for the entire genome
        self.mapinfo = expname+"/"+expname+"-mapinfo.txt" 
        self.mapplot = expname+"/"+expname+"-mapinfo.pdf"

        # Mapping statistics and plots specifically for tRNA regions
        self.trnamapfile = expname+"/"+expname+"-trnamapinfo.txt"
        self.trnamapplot = expname+"/"+expname+"-trnamapinfo.pdf"
        
        # Log file capturing mapping statistics
        self.maplog = expname+"/"+expname+"-mapstats.txt"

        # Files containing gene type classification and read counts for the entire experiment
        self.genetypes = expname+"/"+expname+"-genetypes.txt"
        self.genecounts = expname+"/"+expname+"-readcounts.txt"
        self.trnacounts = expname+"/"+expname+"-trnacounts.txt"
        
        # File containing normalized read counts after size factor correction
        self.normalizedcounts = expname+"/"+expname+"-normalizedreadcounts.txt"
        self.sizefactors = expname+"/"+expname+"-SizeFactors.txt"

        # Gene type-specific read counts and corresponding plots
        self.genetypecounts=expname+"/"+expname+"-typecounts.txt"
        self.genetypeplot=expname+"/"+expname+"-typecounts.pdf"

        # Files containing "real" gene type counts
        self.genetyperealcounts=expname+"/"+expname+"-typerealcounts.txt"
        self.genetyperealplot=expname+"/"+expname+"-typerealcounts.pdf"
        
        # Files containing amino acid counts for tRNA and corresponding plots
        self.trnaaminofile=expname+"/"+expname+"-aminocounts.txt"
        self.trnaaminoplot=expname+"/"+expname+"-aminocounts.pdf"
        self.trnaaminorealplot=expname+"/"+expname+"-aminorealcounts.pdf"
        
        # Files containing anticodon-specific counts for tRNAs
        self.trnaanticodonfile=expname+"/"+expname+"-anticodoncounts.txt"
        
        # Files tracking the length of reads mapped to tRNA regions and corresponding plots
        self.trnalengthfile=expname+"/"+expname+"-readlengths.txt"
        self.trnalengthplot=expname+"/"+expname+"-readlengths.pdf"
        
        # Files tracking mismatch counts and corresponding plots
        self.mismatchcountfile=expname+"/"+expname+"-mismatches.txt"
        self.mismatchcountplot=expname+"/"+expname+"-mismatches.pdf"
        
        # Files tracking tRNA coverage (how much of the tRNA region is covered by reads) and corresponding plots
        self.trnacoveragefile=expname+"/"+expname+"-coverage.txt"
        self.trnacoverageplot=expname+"/"+expname+"-coverage.pdf"
        self.trnacombinecoverageplot=expname+"/"+expname+"-combinecoverage.pdf"
        
        # File containing unique coverage information for tRNA regions
        self.trnauniqcoveragefile=expname+"/"+expname+"-uniqcoverage.txt"

        # Files tracking coverage of pre-tRNA loci and corresponding plots
        self.locicoveragefile=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.txt"
        self.locicoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.pdf"
        self.locicombinecoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcombinecoverage.pdf"
        
        # Files related to mismatches in tRNA sequences and deletions, along with plots and reports
        self.trnamismatchfile = expname+"/mismatch/"+expname+"-mismatchcoverage.txt"
        self.trnamismatchplot = expname+"/mismatch/"+expname+"-mismatchcoverage.pdf"
        
        # Files realted to tRNA deletion events and corresponding plots
        self.trnadeletefile = expname+"/mismatch/"+expname+"-deletecoverage.txt"
        self.trnadeleteplot = expname+"/mismatch/"+expname+"-deletecoverage.pdf"
        
        # Report summarizing mismatch statistics and unique counts and end positions
        self.trnamismatchreport = expname+"/mismatch/"+expname+"-mismatchreport.txt"
        self.trnauniquefile=expname+"/unique/"+expname+"-trnauniquecounts.txt"
        self.trnaendfile=expname+"/"+expname+"-trnaendcounts.txt"
        
        # Plots for prinicipal component analysis
        self.pcaplot = expname+"/"+expname+"-pca.pdf"
        self.pcatrnaplot = expname+"/"+expname+"-pcatrna.pdf"
        
        # Quality assessment report output as html file
        self.qaoutputname = expname+"/"+expname+"-qa.html"
        
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- FUNCTION DEFINITIONS
# Function to make a bed file
"""
The makefeaturebed function consolidates genomic features from various sources.
(e.g., tRNA information, GTF files, BED files) into a single output BED file.
This output BED file will contain all features related to the experiment
"""
def makefeaturebed(trnainfo, expinfo, ensgtf, bedfiles):
    
    # Open the output file for writing
    allfeatfile = open(expinfo.allfeats, "w")
    
    # Process the mature tRNAs, tRNA loci, gene features, and additional features from the trnainfo objects
    for currfeature in readbed(trnainfo.maturetrnas):
        print >>allfeatfile, currfeature.bedstring()
    for currfeature in readbed(trnainfo.locifile):
        print >>allfeatfile, currfeature.bedstring()
    for currfeature in readgtf(ensgtf):
        print >>allfeatfile, currfeature.bedstring(name = currfeature.data["genename"])
    for currbed in bedfiles:
        for currfeature in readbed(currbed):
            print >>allfeatfile, currfeature.bedstring()
    allfeatfile.close()

# Function to map sequencing reads to tRNA and other genomic sequences
"""
The mapsamples function calls the testmain method from the mapreads module to run the mapping process.
The method maps the reads in the sample file to the tRNA database and other non-tRNA sequences, while 
logging the results and creating output files.
"""
def mapsamples(samplefile, trnainfo,expinfo, lazyremap, bamdir = "./",  cores = 8, minnontrnasize = 20, local = False):
    
    # Calls the testmain method from the mapreads module to run the mapping process.
    mapreads.testmain(
        samplefile = samplefile, # Input file containing sample sequencing reads
        trnafile = trnainfo.trnatable, # Path to the tRNA table file, which contains tRNA annotations
        bowtiedb = trnainfo.bowtiedb,  # Path to the Bowtie2 index datavase for tRNA mapping
        bamdir = bamdir, # Directory where BAM files (mapped reads) will be stored
        otherseqs = trnainfo.otherseqs, # File containing non-tRNA sequences for additional mapping
        logfile = expinfo.maplog, # Path to the log file where mapping statistics are recorded
        mapfile = expinfo.mapinfo, # Path to the output file where mapping info (e.g. statistics) will be stores
        trnamapfile = expinfo.trnamapfile, # Path to the output file for tRNA-soecific mapping statistics
        lazy = lazyremap, # Flag to enable/disable lazy remapping (skipping redundant steps.)
        cores = cores, # Number of CPU cores to use for the mapping process (default = 8)
        minnontrnasize = minnontrnasize, # Minimum read size for mapping to non-tRNA sequences (default = 20 bp)
        local = local # Bolean flag to enable local alignment (if set to True)
        ) 
    
# Function to  count reads
"""
countfeatures function counts reads aligned to various features (tRNA, genes, etc.)
in sequencing data and outputs counts for further analysis, using parallel processing
and customizable mismatch thresholds.
Calls the testmain method from the countreads module to count reads for different genomic features.
This includes tRNA features, gene types (from GTF), and other regions specified in BED files
"""    
def countfeatures(samplefile, trnainfo,expinfo, ensgtf, bedfiles,  bamdir = "./", cores = 8, maxmismatches = None):
    countreads.testmain(
        samplefile = samplefile, # Input file with sample reads to be counted
        ensemblgtf = ensgtf, # Ensembl GTF file containing gene annotations
        maturetrnas = [trnainfo.maturetrnas], # List containing path to BED file with mature tRNA annotations
        bamdir = bamdir, # Directory containing BAM files of mapped reads
        otherseqs = trnainfo.otherseqs, # File containing non-tRNA sequences for additional mapping
        trnaloci = [trnainfo.locifile], # List containing ED file with tRNA loci (pre-tRNA regions)
        removepseudo = True, # Flag to remove pseudogenes from the counts (default = True)
        genetypefile = expinfo.genetypes, # Path to file where gene type annotations will be stored
        trnatable = trnainfo.trnatable, # File containing tRNA table for counting tRNA-related reads
        countfile = expinfo.genecounts, # output file for storing gene-level read counts
        bedfile = bedfiles, # Additional BED files containing regions of intrest to count
        trnacounts = expinfo.trnacounts, # Output file for tRNA-specific read counts
        trnaends = expinfo.trnaendfile, # Output file for read counts at tRNA 3' and 5' ends
        trnauniquecounts = expinfo.trnauniquefile, # Output file for unique tRNA read counts
        nofrag = nofrag, # Option to skip fragment counting, a potential feature in countreads
        cores = cores, # Number of CPU cores to use for parallel processing (default = 8)
        maxmismatches = maxmismatches # maximum allowed mismathes per read (if specified)
        )
    
    # Optional R scripts to run PCA
    #runrscript(scriptdir+"pcareadcounts.R",expinfo.normalizedcounts,samplefile,expinfo.pcaplot)
    #runrscript(scriptdir+"pcareadcounts.R",expinfo.trnacounts,samplefile,expinfo.pcatrnaplot)

# Function to count types
"""
The counttypes function handles counting and analyzing different types of genomic features,
particularly focusing on gene types, tRNA amino acids, and mismatch rates using mapped reads
"""
def counttypes(samplefile, trnainfo, expinfo, ensgtf, bedfiles, bamdir = "./",  ignoresizefactors = False, countfrags = False, bamnofeature = False, cores = 8):
    
    # If size factors should be considered in normalization
    if not ignoresizefactors:
        
        # Perform the counting of different read types
        countreadtypes.main(
            sizefactors = expinfo.sizefactors, # Use pre-computed size factors to normalize counts
            combinereps = True,  # Combine replicates of the samples for counting
            bamdir = bamdir, # Directory containing BAM files (aligned reads)
            otherseqs = trnainfo.otherseqs, # File with non-tRNA sequences
            samplefile = samplefile, # Input file containing the sample reads
            maturetrnas = [trnainfo.maturetrnas], # BED file of mature tRNA annotations
            trnatable = trnainfo.trnatable, # Table with detailed tRNA annotations
            trnaaminofile = expinfo.trnaaminofile, # File to store counts of reads aligned to tRNA amino acids
            trnaanticodonfile = expinfo.trnaanticodonfile, # File to sotre counts of reads aligned to tRNA anticodons
            ensemblgtf = ensgtf, # Ensembl GTF file with gene annotations
            trnaloci = [trnainfo.locifile], #BED file with tRNA loci (pre-tRNA regions)
            countfile = expinfo.genetypecounts, # Output file for counts by gene type
            realcountfile = expinfo.genetyperealcounts, # output file for read counts (non-normalized) by gene type
            mismatchfile = expinfo.mismatchcountfile, # File for storing counts of mismatches in the reads
            bedfile = bedfiles, # Additional BED files for regions of intrest
            readlengthfile = expinfo.trnalengthfile, # File for read length distributions
            countfrags = countfrags, # Boolean to count fragments instead of reads (optional)
            bamnofeature = bamnofeature, # Boolean for handling BAM files with no feature tag
            cores = cores # NUmber of CPU cores for parallel processing (defaul = 8)
            )
        
        #Optional R script to plot reads by gene type and tRNAs by amino acid
        #runrscript(scriptdir+"genefeatures.R",expinfo.genetypecounts,expinfo.genetypeplot)
        #runrscript(scriptdir+"featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot, "all")
        #runrscript(scriptdir+"featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminorealplot, "all")
        #
        #runrscript(scriptdir+"featuretypesreal.R",expinfo.genetyperealcounts,expinfo.genetyperealplot)
        #
        #runrscript(scriptdir+"readlengthhistogram.R",expinfo.trnalengthfile,samplefile,expinfo.trnalengthplot)
        #runrscript(scriptdir+"plotreadmismatches.R",expinfo.mismatchcountfile,expinfo.mismatchcountplot)
        
    else:
        countreadtypes.testmain(combinereps = True,
                                samplefile = samplefile,
                                maturetrnas = [trnainfo.maturetrnas],
                                otherseqs = expinfo.otherseqs, 
                                bamdir = bamdir, 
                                trnatable = trnainfo.trnatable,
                                trnaaminofile = expinfo.trnaaminofile,
                                trnaanticodonfile = expinfo.trnaanticodonfile,
                                ensemblgtf = ensgtf,
                                trnaloci = [trnainfo.locifile],
                                countfile = expinfo.genetypecounts,
                                realcountfile = expinfo.genetyperealcounts,
                                bedfile = bedfiles,
                                readlengthfile = expinfo.trnalengthfile,
                                countfrags = countfrags, 
                                cores = cores
                                )
        
    # Plot reads by gene type and tRNAs by amino acid
    runrscript(scriptdir+"03b_genefeatures.R",expinfo.genetypecounts,expinfo.genetypeplot)
    runrscript(scriptdir+"03c_featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot, "all")
    runrscript(scriptdir+"02b_featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminorealplot, "all")
    runrscript(scriptdir+"02b_featuretypesreal.R",expinfo.genetyperealcounts,expinfo.genetyperealplot)
    runrscript(scriptdir+"03d_readlengthhistogram.R",expinfo.trnalengthfile,samplefile,expinfo.trnalengthplot)
    runrscript(scriptdir+"03e_plotreadmismatches.R",expinfo.mismatchcountfile,expinfo.mismatchcountplot)

# Function to get tRNA coverage
"""
the gettrnacoverage function calculates coverage across tRNA ergions from BAM files,
normalizes the data (if size factors are provided), and generates coverage and mismatch
plots using R scripts
"""        
def gettrnacoverage(samplefile, trnainfo,expinfo, bamdir = "./",  orgtype = "euk",ignoresizefactors = False, cores = 8, mincoverage = None):
    #print >>sys.stderr, orgtype
    
    # IF size factors are to be considered (normalization applied)
    if not ignoresizefactors:
        getcoverage.testmain(
            samplefile = samplefile, # Input sample file with aligned reads
            bedfile = [trnainfo.maturetrnas], # BED file with mature tRNA regions
            locibed = [trnainfo.locifile], # BED file with tRNA loci (pre-tRNA regions)
            locistk = trnainfo.locialign, # Pre-aligned file for the tRNA loci
            bamdir = bamdir, # Directory containing BAM files
            lociedgemargin = 30, # Margin for the tRNA loci region
            sizefactors = expinfo.sizefactors, # Size factors for normalization
            orgtype = orgtype, # Organism type
            locicoverage = expinfo.locicoveragefile, # Output file for tRNA loci coverage
            stkfile = trnainfo.trnaalign, # Pre-aligned file for mature tRNAs
            numfile = trnainfo.trnanums, # File with tRNA numbering information
            locinums = trnainfo.locinums, # File with loci numbering information
            allcoverage = expinfo.trnacoveragefile, # Output file for overall tRNA coverage
            trnafasta = trnainfo.trnafasta, #FASTA file with tRNA sequences
            cores = cores, # NUmber of CPU cores for parallel processing (default = 8)
            uniqcoverage = expinfo.trnauniqcoveragefile,  # Output file for unique tRNA coverage
            mincoverage = mincoverage # Minimum coverage threshold for reporting
            )
        
        # Run Rscripts to plot coverage and mismatches
        runrscript(
            scriptdir+"03f_newcoverageplots.R",
            "--cov="+expinfo.trnacoveragefile,
            "--locicov="+expinfo.locicoveragefile,
            "--trna="+trnainfo.trnatable,
            "--samples="+samplefile,
            "--allcov="+expinfo.trnacoverageplot,
            "--runname="+expname,
            "--modomics="+trnainfo.modomics,
            "--combinecov="+expinfo.trnacombinecoverageplot,
            "--directory="+expname
            )
        runrscript(
            scriptdir+"03g_boxplotmismatches.R",
            "--runname="+expinfo.expname,
            "--mismatch="+expinfo.trnacoveragefile,
            "--trna="+trnainfo.trnatable,
            "--samples="+samplefile,
            "--directory="+expname+"/mismatch/"
            )
    else:
        getcoverage.testmain(
            samplefile = samplefile,
            bedfile = [trnainfo.maturetrnas],
            stkfile = trnainfo.trnaalign,
            uniquename = expname+"/"+expname,
            orgtype = orgtype, 
            bamdir = bamdir, 
            allcoverage = expinfo.trnacoveragefile,
            trnafasta = trnainfo.trnafasta, 
            cores = cores, 
            uniqcoverage = expinfo.trnauniqcoveragefile, 
            mincoverage = mincoverage
            )
        
        # Run Rscripts to plot coverage and mismatches
        runrscript(
            scriptdir+"03f_newcoverageplots.R",
            "--cov="+expinfo.trnacoveragefile,
            "--locicov="+expinfo.locicoveragefile,
            "--trna="+trnainfo.trnatable,
            "--samples="+samplefile,
            "--allcov="+expinfo.trnacoverageplot,
            "--runname="+expname,
            "--modomics="+trnainfo.modomics,
            "--combinecov="+expinfo.trnacombinecoverageplot,
            "--directory="+expname
            )
        runrscript(
            scriptdir+"03g_boxplotmismatches.R",
            "--runname="+expinfo.expname,
            "--mismatch="+expinfo.trnacoveragefile,
            "--trna="+trnainfo.trnatable,
            "--samples="+samplefile,
            "--directory="+expname+"/mismatch/"
            )
'''
def getendscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getends.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],sizefactors=expinfo.sizefactors,stkfile=trnainfo.trnaalign,uniquename=expname+"/mismatch/"+expname, allmismatch=expinfo.trnamismatchfile,trnafasta = trnainfo.trnafasta,mismatchfile=expinfo.trnamismatchfile,mismatchreport=expinfo.trnamismatchreport, indelfile=expinfo.trnadeletefile)
        runrscript(scriptdir+"endplots.R","--cov="+expinfo.trnamismatchfile,"--mismatchcov="+expinfo.trnamismatchfile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnamismatchplot,"--uniquename="+expname+"/mismatch/"+expname,"--modomics="+trnainfo.modomics,"--directory="+expname+"/mismatch/")
        runrscript(scriptdir+"boxplotmismatches.R","--mismatch="+expinfo.trnamismatchreport,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
    else:
        getends.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/mismatch/"+expname, allmismatch=expinfo.trnamismatchfile,trnafasta = trnainfo.trnafasta,mismatchfile=expinfo.trnamismatchfile,mismatchreport=expinfo.trnamismatchreport )
        runrscript(scriptdir+"endplots.R","--cov="+expinfo.trnamismatchfile,"--mismatchcov="+expinfo.trnamismatchfile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnamismatchplot,"--uniquename="+expname+"/mismatch/mismatch/"+expname,"--modomics="+trnainfo.modomics,"--directory="+expname+"/mismatch/")
        runrscript(scriptdir+"boxplotmismatches.R","--mismatch="+expinfo.trnamismatchreport,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
def getlocuscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],sizefactors=expinfo.sizefactors,stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile) #removed minextend = 5
        runrscript(scriptdir+"locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
    else:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile)
        runrscript(scriptdir+"locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
'''

# Function to extract tRNA-derived fragment information
"""
The gettdrinfo function runs a shell script to extract tDR (tRNA-derived fragments)
information from the given sample file using a specific tDR database and experiment name
"""
def gettdrinfo(samplefile, dbname,expname):
    
    # Contruct bash command
    tdrcommand = " ".join([
        "bash",
        scriptdir+""+"tdrtrax.bash", 
        samplefile, 
        dbname,
        expname+"/"+expname+"tdrs"
        ])
    
    # Print the generated command to the standard error for debugging purposes
    print >>sys.stderr, tdrcommand

    # Execute the bash script using the constructed command
    tdrjob = subprocess.Popen(
        tdrcommand,
        stdout = subprocess.PIPE,
        stderr = subprocess.STDOUT, 
        shell = True
        )
    print >>sys.stderr, tdrjob.communicate()[0]

# Function to make trackhub plot    
def createtrackhub(samplefile, trnainfo,expinfo):
    maketrackhub.main(
        genomedatabase = trnainfo, 
        samplefile = samplefile,
        expname = expinfo.expname
        )

# Function to get traxqc
def gettraxqc(samplefile, trnainfo,expinfo, tgirtmode = False):
    traxqc.main(
        samplefile = samplefile,
        databasename = trnainfo.dbname,
        experimentname = expinfo.expname,
        tgirt = tgirtmode, 
        output = expinfo.qaoutputname
        )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- PARSE ARGUMENTS
args = parser.parse_args()
dbname = args.databasename
expname = args.experimentname
pairfile =  args.exppairs
ensgtf = args.ensemblgtf
samplefilename = args.samplefile
lazyremap = args.lazyremap
bedfiles= args.bedfile
nofrag= args.nofrag
nosizefactors = args.nosizefactors
olddeseq = args.olddeseq
bamdir = args.bamdir
maponly = args.maponly
local = args.local
maxmismatches = args.maxmismatches
mincoverage = args.mincoverage
mismatch = False
paironly= args.paironly
splittypecounts = False
bamnofeature = args.dumpother
minnontrnasize = args.minnontrnasize
if args.cores is None:
    cores = min(8,cpu_count())
else:
    cores = int(args.cores)
hubonly = args.hubonly
makehubs = args.makehub 
maketdrs=  False # args.maketdr

'''
if args.makeall:
    makehubs = True
    maketdrs = True
'''

# Path to script dir
#scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- SOFTWARE TESTING
# Function to test see if samtools exists
def testsamtools(): #Version: 1.6
    samversionre = re.compile(r"Version\:\s*([\.\d]+)")
    samtoolsloc = get_location("samtools")
    if samtoolsloc is None:
            print >>sys.stderr, "Cannot find samtools in path"
            print >>sys.stderr, "Make sure samtools is installed"
    samtoolsjob = subprocess.Popen([samtoolsloc,"--help"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    samtoolsresults = samtoolsjob.communicate()[0]
    if samtoolsjob.returncode  != 0:
            print >>sys.stderr, "Samtools failed to run"
            print >>sys.stderr, "Make sure samtools is functioning" 
    samtoolsres = samversionre.search(samtoolsresults)
    if samtoolsres:
        if LooseVersion(samtoolsres.group(1)) < LooseVersion("1.0.0"):
            print >>sys.stderr, "Old samtools version "+samtoolsres.group(1)+" found"
            print >>sys.stderr, "Upgrade to latest version"
            sys.exit(1)
    else:
        print >>sys.stderr, "Could not find samtools version number"

# Function to see if R exists       
def testrstats():
    rstatsversionre = re.compile(r"R\s+version\s+((\d+)\.(\d+)\.(\d+))")
    rstatsloc = get_location("R")
    if rstatsloc is None:
            print >>sys.stderr, "Cannot find R in path"
            print >>sys.stderr, "Make sure R is installed"
            sys.exit(1)
    rstatsjob = subprocess.Popen([rstatsloc, "--version"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    rstatsresults = rstatsjob.communicate()[0]
    if rstatsjob.returncode  != 0:
            print >>sys.stderr, "R failed to run"
            print >>sys.stderr, "Make sure R is functioning" 
    rstatsres = rstatsversionre.search(rstatsresults)
    if rstatsres:
        if LooseVersion(rstatsres.group(1)) < LooseVersion("3.1.2"):
            print >>sys.stderr, "Old R version "+rstatsres.group(1)+" found"
            print >>sys.stderr, "Upgrade to latest version"
            sys.exit(1)
    else:
        print >>sys.stderr, "Could not find R version number"

# Test if R exists       
testrstats()
get_location("Rscript")

# Test if samtools and Bowtie2 exist
testsamtools()
get_location("bowtie2")
gitversion, gitversionhash = getgithash(scriptdir)

#trnainfo.test(trnainfo)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- PARSING INPUT FILES
# Read in the sampledata file
sampledata = samplefile(samplefilename)
samples = sampledata.getsamples()

# Check for duplicate names in the sample file
if len(samples) < len(set(samples)):
    print >>sys.stderr, "duplicate sample names in sample file"
    sys.exit(1)

# Check for hypens and leading numerical characters in sample names
for currsample in samples:
    if '-' in currsample:
        print >>sys.stderr, "Sample names containing '-' character are not allowed"
        sys.exit(1)
    if currsample[0].isdigit():
        print >>sys.stderr, "Sample names starting with digits are not allowed"
        sys.exit(1)

# Same logic as above for replicates column
replicates = sampledata.allreplicates()
for currsample in replicates:
    if '-' in currsample:
        print >>sys.stderr, "Sample names containing '-' character are not allowed"
        sys.exit(1)
    if currsample[0].isdigit():
        print >>sys.stderr, "Sample names starting with digits are not allowed"
        sys.exit(1)

# Set replicates       
replicates = set(replicates)        
       
# If pair file exists, check for missing pairs
if pairfile is not None:
    missingnames = set()
    for fir, sec in getpairfile(pairfile):
        #print >>sys.stderr, "**"
        if fir not in replicates:
            missingnames.add(fir)
        if sec not in replicates:
            missingnames.add(sec)
    if len(missingnames) > 0:
        print >>sys.stderr, "Pair names "+",".join(missingnames)+" are not present in sample file"
        sys.exit(1)
        
#sys.exit(0)

# Set DESeq version
deseqversion = "DESeq2"
if olddeseq:
    deseqversion = "DESeq"
if runrscript(scriptdir+"03h_checkRmodules.R",deseqversion) > 0:
    print >>sys.stderr, "Not all R modules needed are installed"
    print >>sys.stderr, "check README for needed R modules"
    sys.exit(1)
    

# Ensure that directories and paths exist and make directories if necessary
#mkdir -p expname
if not os.path.exists(expname):
    os.makedirs(expname)
if not os.path.exists(expname+"/mismatch"):
    os.makedirs(expname+"/mismatch")
if not os.path.exists(expname+"/pretRNAs"):
    os.makedirs(expname+"/pretRNAs")
if not os.path.exists(expname+"/unique"):
    os.makedirs(expname+"/unique")

'''
if not os.path.exists(expname+"/indiv"):
    os.makedirs(expname+"/indiv")
'''

# File handling   
if bedfiles is None:   
    bedfiles = list()
dbname = os.path.expanduser(dbname)
if ensgtf is not None:
    ensgtf = os.path.expanduser(ensgtf)
if bedfiles is not None:    
    bedfiles = list(os.path.expanduser(curr) for curr in bedfiles)
trnainfo = trnadatabase(dbname)
orgtype = trnainfo.getorgtype()
#print >>sys.stderr, orgtype
expinfo = expdatabase(expname)
getsamples = samplefile(samplefilename)
if len(getsamples.getsamples()) == 1:
    nosizefactors = True
    

# Check if a pairing file is provided and if the analysis should focus on pairs
if pairfile and paironly:
    
    # If using the old DESeq method
    if olddeseq:

        # Run the DESeq analysis script with provided parameters and capture the return value
        deseqret = runrscript(
            scriptdir+"03i_deseq1.R",
            expname,
            expinfo.genecounts,
            samplefilename
            )
        
        # Check if the DESeq analysis was unsuccessful
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        print >>sys.stderr, scriptdir+"03j_analyzecounts.R", expname, expinfo.genecounts, samplefilename

        # Run the alternatuve counts analysis script and capture the return value
        deseqret = runrscript(
            scriptdir+"03j_analyzecounts.R",
            expname,
            expinfo.genecounts,
            samplefilename, 
            pairfile)

        # Check if the analysis failed
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)
    
    # Generate a scatter plot based on the normalized counts and other parameters
    runrscript(
        scriptdir+"03k_makescatter.R",
        expname,
        expinfo.normalizedcounts,
        trnainfo.trnatable,
        expinfo.genetypes,
        samplefilename,
        pairfile)
    sys.exit(0)

# If the analysis mode is set to pair-only but no pairing file is provided
elif paironly:
    print >>sys.stderr, "pair only mode used but no --pairfile used"
    sys.exit(1)

# If the analysis is set to create only a track hub
if hubonly:
    print >>sys.stderr, "Creating trackhub"      

    # Create the track hub
    createtrackhub(samplefilename, dbname,expinfo)
    sys.exit(0)

"""
# Commented out code from the original code author
getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)
gettdrinfo(samplefilename, dbname,expname)
        tdrtrax.bash samplefile.txt traxdb outputname
        coverage plot of tRNAs
"""

      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- MAP THE READS
# Get the current runtime in seconds since the epoch
runtime = time.time()

# Convert the runtime into a local time structure
loctime = time.localtime(runtime)
print >>sys.stderr, "Mapping Reads"

# need to check here for names with dashes
# Call the map samples function
mapsamples(
    samplefilename, # The name of the file containing the sample reads to be mapped
    trnainfo, # Information about the tRNA database and related parameters
    expinfo, # Experimental information including settings and paths
    lazyremap, # A flag indicating whether to lazily remap reads
    bamdir = bamdir, # Directory for BAM files; set to the variable bamdir
    cores = cores, # NUmber of cores to use for parallel processing; set to the variable cores
    minnontrnasize = minnontrnasize,  # Minimum size for non-tRNA reads; set to the variable minnontrnasize
    local = local # A flag inficating if the mapping should be done locally; set to the variable local
    )

# Construct the filename for the runtime information, combining the experiment name with a suffic
runinfoname = expname+"/"+expname+"-runinfo.txt"
dbinfo = None

# Check if lazy remapping is not enabled
if not lazyremap:
    
    # If lazy remapping is not enabled, open the run info file for writing
    dbinfo = open(runinfoname,"w")
    print >>dbinfo, "Starting"

else:
    
    # If lazy remapping is enabled, open the file for appending
    dbinfo = open(runinfoname,"a")
    print >>dbinfo, "---------------------------------------------------------"
    print >>dbinfo, "redoing"
    
print >>dbinfo, "expname\t"+expname
print >>dbinfo, "time\t"+str(runtime)+" ("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")"
print >>dbinfo, "samplefile\t"+os.path.realpath(samplefilename)
print >>dbinfo, "dbname\t"+os.path.realpath(dbname)
print >>dbinfo, "git version\t"+gitversion
print >>dbinfo, "git version hash\t"+gitversionhash
print >>dbinfo, "command\t"+" ".join(sys.argv)
dbinfo.close()

# Run Rscripts
runrscript(
    scriptdir+"03c_featuretypes.R",
    expinfo.mapinfo,
    expinfo.mapplot
    )
runrscript(
    scriptdir+"03c_featuretypes.R",
    expinfo.trnamapfile,
    expinfo.trnamapplot)

#print >>sys.stderr, "Counting Read Types"
#counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = nosizefactors)

if maponly:
    sys.exit()
    
# Run the makefeaturebed function    
makefeaturebed(
    trnainfo,expinfo, 
    ensgtf, 
    bedfiles
    )  
  
# Count the reads for DEseq2 and scatter plots
print >>sys.stderr, "Counting Reads"

# Run the countfeatures function
countfeatures(
    samplefilename, 
    trnainfo,expinfo, 
    ensgtf, 
    bedfiles, 
    bamdir = bamdir,
    cores = cores, 
    maxmismatches = maxmismatches
    )

# Create a plot of mapped reads                                
print >>sys.stderr, "Analyzing counts"



#Analyze counts and create scatter plots if pair file is provided
if pairfile:
    if olddeseq:
        deseqret = runrscript(
            scriptdir+"03i_deseq1.R",
            expname,
            expinfo.genecounts,
            samplefilename
            )
        if deseqret >  0:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(
            scriptdir+"03j_analyzecounts.R",
            expname,
            expinfo.genecounts,
            samplefilename, 
            pairfile
            )
        print >>sys.stderr, scriptdir+"03j_analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile

        if deseqret > 0:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            print >>sys.stderr, "Likely that a sample did not contain enough reads"
            sys.exit(1)

    # Run plotting scripts
    runrscript(
        scriptdir+"03l_pcareadcounts.R",
        expinfo.normalizedcounts,
        samplefilename,
        expinfo.pcaplot
        )
    runrscript(
        scriptdir+"03l_pcareadcounts.R",
        expinfo.trnacounts,
        samplefilename,
        expinfo.pcatrnaplot
        )
    runrscript(
        scriptdir+"03k_makescatter.R",
        expname,
        expinfo.normalizedcounts,
        trnainfo.trnatable,
        expinfo.genetypes,
        samplefilename,
        pairfile
        )
    runrscript(
        scriptdir+"03m_ccaendplot.R",
        "--ends="+expinfo.trnaendfile,
        "--trna="+trnainfo.trnatable,
        "--samples="+samplefilename,
        "--directory="+expname+"/",
        "--runname="+expname
        )
    
    

elif not nosizefactors:
    if olddeseq:
        deseqret = runrscript(
            scriptdir+"03i_deseq1.R",
            expname,expinfo.genecounts,
            samplefilename
            )
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(
            scriptdir+"03j_analyzecounts.R",
            expname,expinfo.genecounts,
            samplefilename
            ) 
        print >>sys.stderr, scriptdir+"03j_analyzecounts.R",expname,expinfo.genecounts,samplefilename
        
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)

    # Run plotting scripts
    runrscript(
        scriptdir+"03l_pcareadcounts.R",
        expinfo.normalizedcounts,
        samplefilename,
        expinfo.pcaplot
        )
    runrscript(
        scriptdir+"03l_pcareadcounts.R",
        expinfo.trnacounts,
        samplefilename,
        expinfo.pcatrnaplot
        )
    runrscript(
        scriptdir+"03m_ccaendplot.R",
        "--ends="+expinfo.trnaendfile,
        "--trna="+trnainfo.trnatable,
        "--samples="+samplefilename,
        "--directory="+expname+"/",
        "--runname="+expname
        )

# Count the reads by gene type
print >>sys.stderr, "Counting Read Types"

# Run the counttypes function
counttypes(
    samplefilename, 
    trnainfo,expinfo, 
    ensgtf, 
    bedfiles,
    bamdir = bamdir, 
    ignoresizefactors = nosizefactors,
    countfrags = splittypecounts, 
    bamnofeature = bamnofeature, 
    cores = cores
    )


#coverage plot of tRNAs
print >>sys.stderr, "Generating Read Coverage plots"      
#print >>sys.stderr, orgtype

# Run the gettrnacoverage function
gettrnacoverage(
    samplefilename, 
    trnainfo,
    expinfo, 
    bamdir = bamdir, 
    orgtype = orgtype, 
    ignoresizefactors = nosizefactors, 
    cores = cores, 
    mincoverage = mincoverage
    )

"""
# Code commented out by the original code authors
#coverage plot of pre-tRNAs
getlocuscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

#coverage plot of mismatches
getmismatchcoverage(samplefilename, trnainfo,expinfo, nosizefactors)
print >>sys.stderr, "Counting mismatches"      

getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)
"""

runrscript(
    scriptdir+"03n_summed_mismatch_heatmap.R",
    "--runname="+expinfo.expname,
    "--mismatch="+expname+"/"+expname+"-coverage.txt",
    "--pairfile="+samplefilename,
    "--directory="+expname+"/mismatch/"
)
sys.exit(0)


# Run the gettraxqc function
gettraxqc(
    samplefilename, 
    trnainfo, 
    expinfo, 
    tgirtmode = nofrag
    )

if makehubs:
    print >>sys.stderr, "Creating trackhub"      

    # Run the trackhub function
    createtrackhub(samplefilename, dbname,expinfo)


if (os.path.isfile(scriptdir+""+"tdrtrax.bash") and maketdrs):
    print >>sys.stderr, "Creating tdrs"      

    gettdrinfo(samplefilename, dbname,expname)
        #tdrtrax.bash samplefile.txt traxdb outputname

    