#!/usr/bin/env python3
"""
ronaldo calculate false positive rates for different assays

ronaldo calculate false positive rates for different assays

### CHANGE LOG ### 
2020-09-03 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build - split from dirty scripts
"""
import collections
import logging
import csv 
import datetime
import matplotlib.pyplot as plt
from os import path, mkdir, listdir
import numpy as np
import argparse
import meta
import sys
import time
from sam_util import get_genome_recovery, get_genome_coverage
import csv

epi = "Licence: " + meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig()
log = logging.getLogger()

def check_blanks(blank_list, ref_length=29000):
    max_coverage = 0
    max_recovery_10 = 0
    max_recovery_20 = 0
    max_reads = 0
    # TODO: Logic to handle Illumina or Nanopore
    for bam_file in blank_list:
        current_recovery_10, current_recovery_20,  = get_genome_recovery(bam_file)
        current_coverage, current_reads = get_genome_coverage(bam_file)
        if current_recovery_10 > max_recovery_10:
            max_recovery_10 = current_recovery_10
            max_recovery_20 = current_recovery_20
        if current_coverage > max_coverage:
            max_coverage = current_coverage
        if current_reads > max_reads:
            max_recovery = current_reads                        
    return max_coverage, max_recovery_10, max_recovery_20, max_reads

def get_sample_metrics(bam_file):
    coverage = 0
    # TODO: Logic to handle Illumina or Nanopore
    recovery_10 = 0
    recovery_20 = 0
    reads = 0
    recovery_10, recovery_20 = get_genome_recovery(bam_file)
    coverage, reads = get_genome_coverage(bam_file)
    return coverage, recovery_10, recovery_20, reads

def calculate_metrics(args): 
    # Create path for blanks
    blank_paths = [path.join(args.bamfolder, blank_file) for blank_file in args.blankbam]
    input_good = True
    for blank_path in blank_paths:
        if not path.exists(blank_path):
            log.error(f'Filepath to Blank does not exist: {blank_path}')
            input_good = False
    if input_good:
        #Create output dir 
        if not path.exists(args.db):
            mkdir(args.db)
        # Check blanks
        blank_coverage, blank_recovery_10, blank_recovery_20, blank_reads = check_blanks(blank_paths)
        log.info(f'Max blank genome coverage: {blank_coverage}')
        log.info(f'Max blank genome recovery >10X: {blank_recovery_10}')
        log.info(f'Max blank genome recovery >20X: {blank_recovery_20}')
        log.info(f'Max blank number of mapped reads: {blank_reads}')
        cov_cut_not_ok = False
        if blank_recovery_20 > args.blank_recovery_cutoff and args.ont:
            cov_cut_not_ok = True
        elif blank_recovery_10 > args.blank_recovery_cutoff and not args.ont:
            cov_cut_not_ok = True
        if cov_cut_not_ok or blank_reads > args.blank_read_cutoff:
            log.info('RUN SKIPPED: Specified Blanks from this run have too much SARSCOV2 content')
        else:
            log.info('BLANKS OK!')
            all_samples = {}
            for bam_file in [path.join(args.bamfolder, bam_file) for bam_file in listdir(args.bamfolder) if bam_file not in args.blankbam]: 
                coverage, recovery_10, recovery_20, reads = get_sample_metrics(bam_file)
                # TODO: Include mapping/lookup table
                all_samples[path.basename(bam_file)] = dict(sample_name=path.basename(bam_file), mean_cov=coverage, pc_pos_gte_20=recovery_20, pc_pos_gte_10=recovery_10, no_reads=reads)
            db_path = path.join(args.db, f'ronaldo.db.{path.basename(args.bamfolder)}.csv')
            db_out = csv.DictWriter(open(db_path, 'w'), fieldnames=list(all_samples.values())[0].keys())
            db_out.writeheader()
            db_out.writerows(all_samples.values())


def assess_run(args):
    print(args)


def is_valid_dir(parser, arg):
    if not path.exists(arg):
        parser.error("The directory %s does not exist!" % arg)
    else:
        return arg


if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[1].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    subparsers = parser.add_subparsers(help='commands')

    # Metrics parser
    metrics_parser = subparsers.add_parser('calculate', help='Calculate genome cov metrics')
    metrics_parser.add_argument('-d','--db', action='store',help='DB directory', default='ronaldo_db')
    metrics_parser.add_argument('--ctdata', action='store',  help='Path to table with assay information')
    metrics_parser.add_argument('--blank_read_cutoff', action='store',  default=500,  help='Run skipped if blanks have number of mapped reads')
    metrics_parser.add_argument('--blank_recovery_cutoff', action='store', default=4.0, help='Run skipped if blanks have higher perc. genome recovery')
    metrics_parser.add_argument('--ont', action='store_true', default=False, help='Data is OXFORD NANOPORE')
    metrics_parser.add_argument('bamfolder', action='store', help='Folder of SARSCOV2 BAM files', type=lambda x: is_valid_dir(metrics_parser, x))
    metrics_parser.add_argument('blankbam', metavar='N', nargs='+', help='Negative control BAM file')
    metrics_parser.set_defaults(func=calculate_metrics)

    # Filter parser
    filter_parser = subparsers.add_parser('filter', help='Update Lineage values')
    filter_parser.add_argument('-d','--db', action='store',help='DB directory', default='ronaldo_db')
    filter_parser.add_argument('-o','--output',action='store',help='output directory', default='ronaldo_out')
    filter_parser.add_argument('-c','--coverage',action='store',help='Minimum coverage for fully mapped reads', default=1)
    filter_parser.add_argument('-n','--readcount',action='store',help='Minimum number of fully mapped reads', default=10)

    args = parser.parse_args()
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())    
    if hasattr(args, 'func'):
        args.func(args)
    else: 
        parser.print_help()
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))
    sys.exit(0)