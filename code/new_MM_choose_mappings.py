#!/usr/bin/env python

import re
import sys
import os.path
import itertools
import pysam
import subprocess
import argparse
from collections import defaultdict
from trnasequtils import *

# Constants
defminnontrnasize = 20
maxmaps = 50
BATCH_SIZE = 10000

def getgithash(scriptdir):
    try:
        gitversion = subprocess.check_output(['git', 'describe', '--always'], cwd=scriptdir).strip()
        githash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=scriptdir).strip()
        return gitversion, githash
    except:
        return "unknown", "unknown"

def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)

def process_batch(batch, trnadata, trnatranscripts, bamfile, outfile, stats, minnontrnasize):
    for pairedname, allmaps in batch:
        # Skip unmapped reads
        mappings = [m for m in allmaps if m.tid != -1]
        if not mappings:
            continue

        stats['total_reads'] += 1

        # Find best scoring mappings
        best_score = None
        best_mappings = set()
        for mapping in mappings:
            score = mapping.get_tag("AS")
            if best_score is None or score > best_score:
                best_score = score
                best_mappings = {mapping}
            elif score == best_score:
                best_mappings.add(mapping)
        
        if len(best_mappings) > 1:
            stats['multimaps'] += 1
        
        # Process tRNA mappings
        trna_mappings = [m for m in best_mappings if bamfile.getrname(m.tid) in trnatranscripts]
        if trna_mappings:
            stats['trna_reads'] += 1
            
            # Get unique anticodons and amino acids
            anticodons = frozenset(trnadata.getanticodon(bamfile.getrname(m.tid)) for m in trna_mappings)
            aminos = frozenset(trnadata.getamino(bamfile.getrname(m.tid)) for m in trna_mappings)
            
            # Get all loci
            loci = set()
            for mapping in trna_mappings:
                ref_name = bamfile.getrname(mapping.tid)
                loci.update(trnadata.transcriptdict.get(ref_name, []))
            
            # Track statistics for ambiguous mappings
            if len(anticodons - {'NNN'}) > 1:
                stats['amb_anticodon'] += 1
            if len(aminos - {'Und'}) > 1:
                stats['amb_amino'] += 1
            if len(trna_mappings) > 1:
                stats['amb_trna'] += 1
            
            # Add tRNA-specific tags
            for m in trna_mappings:
                m.tags.extend([
                    ("YA", len(anticodons)),
                    ("YM", len(aminos)),
                    ("YR", len(trna_mappings)),
                    ("YL", len(loci))
                ])
                outfile.write(m)
        else:
            # Process non-tRNA mappings
            if len(mappings[0].seq) < minnontrnasize:
                continue
                
            if len(best_mappings) > maxmaps:
                stats['dupe_remove'] += 1
                continue
                
            if len(best_mappings) > 1:
                stats['nonunique_nontrnas'] += 1
            else:
                stats['unique_nontrnas'] += 1
                
            for m in best_mappings:
                outfile.write(m)

def getbesttrnamappings(trnafile, bamout=True, logfile=sys.stderr, progname=None, fqname=None, libname=None, setcountfile=None, extraseqfilename=None, minnontrnasize=defminnontrnasize):
    # Initialize data structures
    trnadata = transcriptfile(trnafile)
    trnatranscripts = set(trnadata.gettranscripts())
    stats = defaultdict(int)
    
    # Setup input/output
    bamfile = pysam.Samfile("-", "r")
    newheader = bamfile.header.to_dict()
    
    # Setup header
    newheader["RG"] = [{"ID": fqname}] if fqname else []
    if progname:
        scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
        gitversion, githash = getgithash(scriptdir)
        newheader["PG"] = [{"PN": progname, "ID": progname, "VN": gitversion}]
    
    outfile = pysam.Samfile("-", "wb", header=newheader)
    
    # Process reads in batches
    batch = []
    for pairedname, allmaps in itertools.groupby(bamfile, lambda x: x.qname):
        batch.append((pairedname, allmaps))
        if len(batch) >= BATCH_SIZE:
            process_batch(batch, trnadata, trnatranscripts, bamfile, outfile, stats, minnontrnasize)
            batch = []
    
    # Process remaining reads
    if batch:
        process_batch(batch, trnadata, trnatranscripts, bamfile, outfile, stats, minnontrnasize)
    
    # Write statistics
    print >>logfile, "tRNA Reads with multiple transcripts: %d" % stats['amb_trna']
    print >>logfile, "tRNA Reads with multiple anticodons: %d" % stats['amb_anticodon']
    print >>logfile, "tRNA Reads with multiple aminos: %d" % stats['amb_amino']
    print >>logfile, "Total tRNA Reads: %d" % stats['trna_reads']
    print >>logfile, "Single mapped non-tRNAs: %d" % stats['unique_nontrnas']
    print >>logfile, "Multiply mapped non-tRNAs: %d" % stats['nonunique_nontrnas']
    print >>logfile, "Duplicate reads removed: %d" % stats['dupe_remove']
    
    outfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process tRNA mappings efficiently')
    parser.add_argument('trnaname',
                    help='name of tRNA database')
    parser.add_argument('--progname',
                       help='program name')
    parser.add_argument('--fqname',
                       help='fastq file name')
    parser.add_argument('--expname',
                       help='library name')
    parser.add_argument('--trnasetcounts',
                       help='Counts for all sets of tRNAs')
    parser.add_argument('--minnontrnasize', type=int, default=20,
                       help='Minimum read length for non-tRNAs')
    
    args = parser.parse_args()
    getbesttrnamappings(args.trnaname, 
                       progname=args.progname, 
                       fqname=args.fqname, 
                       libname=args.expname, 
                       setcountfile=args.trnasetcounts, 
                       minnontrnasize=args.minnontrnasize)